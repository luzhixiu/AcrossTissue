



library(AnaCoDa)
library(profmem)
library(argparse)
rm(list=ls())


parser <- ArgumentParser()
parser$add_argument("-i","--input",type="character",default="./")
parser$add_argument("-o","--output",type="character",default="./")
parser$add_argument("-d","--div",type="integer",default=0)
parser$add_argument("-s","--samp",type="integer",default=1000)
parser$add_argument("-a","--adapt",type="integer",default=100)
parser$add_argument("-t","--thin",type="integer",default=20)
parser$add_argument("-n","--threads",type="integer",default=1)


args <- parser$parse_args()
div <- args$div
input <- args$input
directory <- args$output
thin <- args$thin
adapt <- args$adapt
samp <- args$samp
num_threads <- args$threads

#locally defined parameters

init_run_number <- 2
max_run_number <- 4

# Helper functions
createParameterOutput <- function(parameter,numMixtures,samples,mixture.labels,samples.percent.keep=1,relative.to.optimal.codon=F,report.original.ref=T)
{
  for (i in 1:numMixtures)
  {
    getCSPEstimates(parameter,paste(dir_name,"Parameter_est",mixture.labels[i],sep="/"),i,samples*samples.percent.keep,relative.to.optimal.codon=relative.to.optimal.codon,report.original.ref = report.original.ref)
  }
}

createTracePlots <- function(trace, model,genome,numMixtures,samples,mixture.labels,samples.percent.keep=1)
{
  for (i in 1:numMixtures)
  {
    plot(trace, what = "Mutation", mixture = i)
    plot(trace, what = "Selection", mixture = i)
    plot(model, genome, samples = samples*samples.percent.keep, mixture = i,main = mixture.labels[i])
  }
}

fasta.folders <- input #, "../data/cds/sampled/",  "../data/cds/sampled/", "../data/cds/filtered/")
fasta.files <- list.files(path=fasta.folders,pattern="*.fasta",full.names = F)
print(fasta.files)
mixture.labels <- unlist(strsplit(fasta.files,split=".fasta"))
fasta.paths <- paste0(fasta.folders, fasta.files)
numMixtures <- length(fasta.files)
mixture.sizes <- rep(0, numMixtures)

## Note: writing a for loop to deal with all mixtures (1 - n.mixtures) is tricky.
## Part of the issue is the appending of the object defined in the command and the assignment of the output
mixture.index <- 1;

genome <- initializeGenomeObject(file=fasta.paths[mixture.index],match.expression.by.id = FALSE,append = FALSE)

mixture.sizes[mixture.index] <- length(genome)
if(numMixtures > 1){
  for(mixture.index in 2:numMixtures)
  {
    tmp.length <- length(genome)
    genome <- initializeGenomeObject(file=fasta.paths[mixture.index],genome=genome,match.expression.by.id = FALSE,append = TRUE,positional = T)
    mixture.sizes[mixture.index] <- length(genome) - tmp.length
  }
}

if(length(genome) != sum(mixture.sizes)){
  stop("length(genomeObj) != sum(mixture.sizes), but it should.")
}else{
  print("FASTA successfully files loaded:");
  print(fasta.files[1:numMixtures])
}

cat("Genome loaded\n")
#initialize parameter object

sphi_init <- rep(1,numMixtures)
with.phi <- F
mixDef <- "allUnique"
percent.to.keep <- 1
size <- length(genome)
cat(size,"\n")
index <- c(1:size)
#geneAssignment <- c(rep(1,size.tmp),rep(2,size.tmp.2-size.tmp),rep(3,size-size.tmp.2))
geneAssignment <- rep(1:numMixtures, mixture.sizes)

# init_phi <- c()
# for (i in phi.path)
# {
#   segment_exp <- read.table(file=i,sep=",",header=TRUE)
#   init_phi <- c(init_phi,segment_exp[,2])
# }

# if(length(genome) != length(init_phi)){
#   stop("length(genomeObj) != length(init_phi), but it should.")
# }else{
#   print("Initial Phi values successfully files loaded:");
# }


parameter <- initializeParameterObject(genome,model="ROC",sphi_init,numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
# Do we need to reissue this command for each restart?
parameter$initMutationCategories("/data1/compbio/zlu21/AcrossTissue/RunResults/Crei_Mutation_NoRef.csv",1,TRUE)
# parameter$initSelectionCategories(c(sel),1)



#initialize MCMC object
samples <-samp
thinning <- thin
adaptiveWidth <-adapt
mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                             est.expression=T, est.csp=TRUE, est.hyper=T,est.mix = FALSE)

# get model object
model <- initializeModelObject(parameter, "ROC", with.phi)

run_number <- init_run_number

if(run_number == 1){ 
dir.create(directory)
dir_name <- paste0(directory,"/run_",run_number)
dir.create(dir_name)
dir.create(paste(dir_name,"Graphs",sep="/"))
dir.create(paste(dir_name,"Restart_files",sep="/"))
dir.create(paste(dir_name,"Parameter_est",sep="/"))
dir.create(paste(dir_name,"R_objects",sep="/"))
setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, F)
#run mcmc on genome with parameter using model
sys.runtime<-system.time(
  runMCMC(mcmc, genome, model, num_threads,divergence.iteration = div)
)


sys.runtime <- data.frame(Value=names(sys.runtime),Time=as.vector(sys.runtime))
write.table(sys.runtime,file=paste(dir_name,"mcmc_runtime.csv",sep="/"),sep=",",col.names = T,row.names = T,quote=F)

print(paste0("Done run_number ", run_number, "run time info: ", sys.runtime))

createParameterOutput(parameter = parameter,numMixtures = numMixtures,mixture.labels = mixture.labels,samples = samples,samples.percent.keep = percent.to.keep,relative.to.optimal.codon = F,report.original.ref = T)
expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*percent.to.keep)
write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)

pdf(paste(dir_name,"Graphs/Parameter_comparisons.pdf",sep="/"), width = 11, height = 12)
plot(parameter,what="Mutation",samples=samples*percent.to.keep,mixture.name=mixture.labels)
plot(parameter,what="Selection",samples=samples*percent.to.keep,mixture.name=mixture.labels)
dev.off()

trace <- parameter$getTraceObject()
pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
createTracePlots(trace=trace,model=model,genome=genome,numMixtures=numMixtures,samples=samples,samples.percent.keep = 1,mixture.labels = mixture.labels)
dev.off()

#plots different aspects of trace

pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
plot(mcmc,what = "LogPosterior")
plot(trace, what = "ExpectedPhi")
aa <- aminoAcids()
done.adapt <- TRUE
for(a in aa)
{
  if (a=="M"||a=="X"||a=="W") next
  accept.trace <- trace$getCodonSpecificAcceptanceRateTraceForAA(a)
  len <- length(accept.trace)
  mean.acceptance <- mean(accept.trace[(len-len*0.5):len])
  if (mean.acceptance < 0.1 || mean.acceptance > 0.44) done.adapt <- FALSE
  plot(accept.trace,main=paste0("Acceptace Rate for ",a),xlab="Samples",ylab="Acceptance Rate",type="l")
}
acfCSP(parameter,csp="Selection",numMixtures = numMixtures,samples=samples*percent.to.keep)
acfCSP(parameter,csp="Mutation",numMixtures = numMixtures,samples=samples*percent.to.keep)
dev.off()
Step

for (i in 1:numMixtures)
{
  param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Selection",mixture=i,frac1=0.1)
  z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
  if (length(z.scores) > 0)
  {
    param.conv <- FALSE
  }
  write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_eta_",i,".txt"),ncolumns = 1)
}


for (i in 1:numMixtures)
{
  param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Mutation",mixture=i,frac1=0.1)
  z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
  if (length(z.scores) > 0)
  {
    param.conv <- FALSE
  }
  write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_m_",i,".txt"),ncolumns = 1)
}


writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))


diag <- convergence.test(mcmc,samples = samples*percent.to.keep,thin=thinning,frac1=0.2)
z<-abs(diag$z)
done <- (z > 1.96) && param.conv
rm(parameter)
rm(trace)
rm(model)

run_number <- run_number + 1

} # end run_number = 1


while((!done) && (run_number < max_run_number))
{
  parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"),model="ROC")
  # Restarts crash when \DeltaM is fixed in initial run
  # Insert fix \DeltaM command here?
  # parameter$initMutationCategories("/data1/compbio/zlu21/AcrossTissue/RunResults/Crei_Mutation_NoRef.csv",1,TRUE)

  dir_name <- paste0(directory,"/run_",run_number)
  dir.create(dir_name)
  dir.create(paste(dir_name,"Graphs",sep="/"))
  dir.create(paste(dir_name,"Restart_files",sep="/"))
  dir.create(paste(dir_name,"Parameter_est",sep="/"))
  dir.create(paste(dir_name,"R_objects",sep="/"))
  
  mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                               est.expression=T, est.csp=TRUE, est.hyper=T,est.mix=FALSE)
  
  model <- initializeModelObject(parameter, "ROC", with.phi)
  setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, F)
  sys.runtime <- system.time(
    runMCMC(mcmc, genome, model, num_threads,div=0)
  )
  sys.runtime <- data.frame(Value=names(sys.runtime),Time=as.vector(sys.runtime))
  write.table(sys.runtime,file=paste(dir_name,"mcmc_runtime.csv",sep="/"),sep=",",col.names = T,row.names = T,quote=F)
  
  print(paste0("Done run_number ", run_number, "run time info: ", sys.runtime))

  createParameterOutput(parameter = parameter,numMixtures = numMixtures,samples = samples,mixture.labels = mixture.labels,samples.percent.keep = percent.to.keep,relative.to.optimal.codon = F,report.original.ref = T)
  
  expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*percent.to.keep)
  write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)
  
  #
  #plots different aspects of trace
  trace <- parameter$getTraceObject()
  pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
  plot(mcmc,what = "LogPosterior")
  plot(trace, what = "ExpectedPhi")
  aa <- aminoAcids()
  done.adapt <- TRUE
  for(a in aa)
  {
    if (a=="M"||a=="X"||a=="W") next
    accept.trace <- trace$getCodonSpecificAcceptanceRateTraceForAA(a)
    len <- length(accept.trace)
    mean.acceptance <- mean(accept.trace[(len-len*0.5):len])
    if (mean.acceptance < 0.1 || mean.acceptance > 0.44) done.adapt <- FALSE
    plot(accept.trace,main=paste0("Acceptace Rate for ",a),xlab="Samples",ylab="Acceptance Rate",type="l")
  }
  acfCSP(parameter,csp="Selection",numMixtures = numMixtures,samples=samples*percent.to.keep)
  acfCSP(parameter,csp="Mutation",numMixtures = numMixtures,samples=samples*percent.to.keep)
  dev.off()
  
  
  for (i in 1:numMixtures)
  {
    param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Selection",mixture=i,frac1=0.1)
    z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
    if (length(z.scores) > 0)
    {
      param.conv <- FALSE
    }
    write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_eta_",i,".txt"),ncolumns = 1)
  }
  
  
  for (i in 1:numMixtures)
  {
    param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Mutation",mixture=i,frac1=0.1)
    z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
    if (length(z.scores) > 0)
    {
      param.conv <- FALSE
    }
    write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_m_",i,".txt"),ncolumns = 1)
  }
  pdf(paste(dir_name,"Graphs/Parameter_comparisons.pdf",sep="/"), width = 11, height = 12)
  plot(parameter,what="Mutation",samples=samples*percent.to.keep,mixture.name=mixture.labels)
  plot(parameter,what="Selection",samples=samples*percent.to.keep,mixture.name=mixture.labels)
  dev.off()
  
  
  pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
  createTracePlots(trace=trace,model=model,genome=genome,numMixtures=numMixtures,samples=samples,samples.percent.keep = percent.to.keep,mixture.labels = mixture.labels)
  dev.off()
  writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
  writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))
  
  diag <- convergence.test(mcmc,samples = samples*percent.to.keep,thin=thinning,frac1=0.1)
  z<-abs(diag$z)
  done <- (z > 1.96) && param.conv
  rm(parameter)
  rm(trace)
  rm(model)

  run_number <- run_number + 1
} #end run_number > 1 & < max_run_number loop

print(paste0("Done printing run_number ", run_number))


samples <- 10000
thinning <- thin 
parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"),model="ROC")
# Restarts crash when \DeltaM is fixed in initial run
# Insert fix \DeltaM command here?
# parameter$initMutationCategories("/data1/compbio/zlu21/AcrossTissue/RunResults/Crei_Mutation_NoRef.csv",1,TRUE)
run_number <- run_number + 1
dir_name <- paste0(directory,"/final_run")
dir.create(dir_name)
dir.create(paste(dir_name,"Graphs",sep="/"))
dir.create(paste(dir_name,"Restart_files",sep="/"))
dir.create(paste(dir_name,"Parameter_est",sep="/"))
dir.create(paste(dir_name,"R_objects",sep="/"))


mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE,est.mix=FALSE)


#this part set steps 
adptiveRatio=0.5
adaptiveSamples=samples*thinning*adaptiveRatio
mcmc$setStepsToAdapt(0)

model <- initializeModelObject(parameter, "ROC", with.phi)
setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, F)
#run mcmc on genome with parameter using model
#p<-profmem({
sys.runtime <- system.time(
  runMCMC(mcmc, genome, model, num_threads)
)
sys.runtime <- data.frame(Value=names(sys.runtime),Time=as.vector(sys.runtime))
write.table(sys.runtime,file=paste(dir_name,"mcmc_runtime.csv",sep="/"),sep=",",col.names = T,row.names = T,quote=F)

print(paste0("Done run_number ", run_number, "run time info: ", sys.runtime))


createParameterOutput(parameter = parameter,numMixtures = numMixtures,samples = samples,mixture.labels = mixture.labels,samples.percent.keep = 1,relative.to.optimal.codon = F,report.original.ref = T)

# mixtureAssignment <- getMixtureAssignmentEstimate(parameter,c(1:size),samples*0.5)
expressionValues <- getExpressionEstimates(parameter,c(1:size),samples)
write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)


#plots different aspects of trace
trace <- parameter$getTraceObject()
pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
plot(mcmc,what = "LogPosterior")
plot(trace, what = "ExpectedPhi")
acfCSP(parameter,csp="Selection",numMixtures = numMixtures,samples=samples)
acfCSP(parameter,csp="Mutation",numMixtures = numMixtures,samples=samples)
dev.off()

pdf(paste(dir_name,"Graphs/Parameter_comparisons.pdf",sep="/"), width = 11, height = 12)
plot(parameter,what="Mutation",samples=samples,mixture.name=mixture.labels)
plot(parameter,what="Selection",samples=samples,mixture.name=mixture.labels)
dev.off()


pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
createTracePlots(trace=trace,model=model,genome=genome,numMixtures=numMixtures,samples=samples,samples.percent.keep = 1,mixture.labels = mixture.labels)
dev.off()

for (i in 1:numMixtures)
{
  param.diag<-convergence.test(trace,samples=samples,thin = thinning,what="Selection",mixture=i,frac1=0.1)
  z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
  if (length(z.scores) > 0)
  {
    param.conv <- FALSE
  }
  write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_eta_",i,".txt"),ncolumns = 1)
}

for (i in 1:numMixtures)
{
  param.diag<-convergence.test(trace,samples=samples,thin = thinning,what="Mutation",mixture=i,frac1=0.1)
  z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
  if (length(z.scores) > 0)
  {
    param.conv <- FALSE
  }
  write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_m_",i,".txt"),ncolumns = 1)
}
rm(trace)
rm(model)
writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))


