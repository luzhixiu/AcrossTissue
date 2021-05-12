## Description:
## WARNING: There are mulitple copies of this file.
## WARNING: Master version is in 
## Code orignally form Alex Cope, modified for work on C elegans
## More info here!
## File called from script run.sh



## Load libraries
print("Loading Libraries")
library(AnaCoDa)
library(profmem)
library(argparse)

## Clear previous definitions
print("Clearing local object definitions prior to using script arguments")
rm(list=ls())


## Global parameters

## Create list indicating which sets of parameters to estimate 
## Values could be "Selection", "Mutation", or both of them.

nRuns = 2; # number of model fit runs 
cspToEstimate = c("Selection")
# Values to use for mutation
file.mutation.values = "/data1/compbio/zlu21/AcrossTissue/RunResults/Crei_Mutation_NoRef.csv"

with.phi <- FALSE
mixDef <- "allUnique"
proportion.to.keep <- 1 ## note function uses percent.to.keep, but it's really proportion
adptive.proportion <- 0.5 # proportion of run that uses the adaptive MCMC
initial.phi <- "random" # either "random", a function name, or files 
initial.run.number <- 1; # uses restart files for run > 1
run.convergence.tests <- FALSE
## Process arguments passed from script file

## NOTE: We should have the following command print the script name
print("Processing arguments passed from script")

parser <- ArgumentParser()
parser$add_argument("-i","--input",type="character",default="./")
parser$add_argument("-o","--output",type="character",default="./")
parser$add_argument("-d","--div",type="integer",default=0)
parser$add_argument("-s","--samp",type="integer",default=5000)
parser$add_argument("-a","--adapt",type="integer",default=100)
parser$add_argument("-t","--thin",type="integer",default=20)
parser$add_argument("-n","--threads",type="integer",default=1)

args <- parser$parse_args()
div <- args$div
input <- args$input
directory <- args$output
thinning <- args$thin
adaptiveWidth <- args$adapt
samples <- args$samp
num_threads <- args$threads


## Local Functions

print("Defining local functions.")
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
    for(what in cspToEstimate){
      plot(trace, what = what, mixture = i)
    }
    plot(model, genome, samples = samples*samples.percent.keep, mixture = i,main = mixture.labels[i])
  }
}

## Process arguments
print("Processing arguments")

fasta.folders <- input #, "../data/cds/sampled/",  "../data/cds/sampled/", "../data/cds/filtered/")
fasta.files <- list.files(path=fasta.folders,pattern="*.fasta",full.names = F)
print(fasta.files)
mixture.labels <- unlist(strsplit(fasta.files,split=".fasta"))
fasta.paths <- paste0(fasta.folders, fasta.files)
numMixtures <- length(fasta.files)
mixture.sizes <- rep(0, numMixtures)
mixture.indexes <- 1:numMixtures; ## This should be defined based on previous arguments, such as numMixtures

## Fit Model
##
## Note: writing a for loop to deal with all mixtures (1 - n.mixtures) is tricky.
## Part of the issue is the appending of the object defined in the command and the assignment of the output



## Loop over mixtures

for(mixture.index in mixture.indexes){
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
}

# Check final genome object and make sure it's as expected
if(length(genome) != sum(mixture.sizes)){
  stop("length(genomeObj) != sum(mixture.sizes), but it should.")
}else{
  print("FASTA successfully files loaded:");
  print(fasta.files[1:numMixtures])
}
print("Genome loaded\n")


## Initialize parameter objects

sphi_init <- rep(1,numMixtures)
size <- length(genome)
cat(size,"\n")
index <- c(1:size)
#geneAssignment <- c(rep(1,size.tmp),rep(2,size.tmp.2-size.tmp),rep(3,size-size.tmp.2))
geneAssignment <- rep(1:numMixtures, mixture.sizes)



## Should have flag here for 

if(initialPhi !="random"){
  print("Don't know how to handle initial phi values. Exiting")
  break;
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
}


done=FALSE; #Flag used in while() for model fits
run.number <- initial.run.number

while((!done) && (run.number <= 3))
{

  ## Set up directories
  dir_name <- paste0(directory,"/run_",run.number)
  old_dir_name <- paste0(directory,"/run_",run.number-1) ## TODO: Test to make sure file exists 

  dir.create(directory)
  dir.create(dir_name)
  dir.create(paste(dir_name,"Graphs",sep="/"))
  dir.create(paste(dir_name,"Restart_files",sep="/"))
  dir.create(paste(dir_name,"Parameter_est",sep="/"))
  dir.create(paste(dir_name,"R_objects",sep="/"))

  if(run.number==1){
    ## initial run, not using restart files
    parameter <- initializeParameterObject(genome,model="ROC",sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
  }else{
    ## use a restart file
    parameter<-initializeParameterObject(init.with.restart.file = paste(old_dir_name,"Restart_files/rstartFile.rst_final",sep="/"),model="ROC")

  }

  ## Fix mutation based on values in file
  ## May need to only evaluate for run.number = 1
  if(!("Mutation" %in% cspToEstimate) ) parameter$initMutationCategories(file.mutation.values,1,TRUE)

  ## Initialize MCMC object
  mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                               est.expression=T, est.csp=TRUE, est.hyper=T,est.mix = FALSE)

  ## Get model object
  model <- initializeModelObject(parameter, "ROC", with.phi)

  setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, F)
  #run mcmc on genome with parameter using model
  sys.runtime<-system.time(
    runMCMC(mcmc, genome, model, num_threads,divergence.iteration = div)
  )

  print(paste0("Run ", run.number, " completed."))
  ## Writing objects

  print("Writing objects.")

  writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
  writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))



  sys.runtime <- data.frame(Value=names(sys.runtime),Time=as.vector(sys.runtime))
  write.table(sys.runtime,file=paste(dir_name,"mcmc_runtime.csv",sep="/"),sep=",",col.names = T,row.names = T,quote=F)

  createParameterOutput(parameter = parameter,numMixtures = numMixtures,mixture.labels = mixture.labels,samples = samples,samples.percent.keep = proportion.to.keep,relative.to.optimal.codon = F,report.original.ref = T)
  expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*proportion.to.keep)
  write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)

  ## Create Plots

  ## Compare parameters between mixtures
  if(numMixtures > 1)
  {
    print("Printing comparisons of estimated CSP parameters between mixtures")
    pdf(paste(dir_name,"Graphs/Parameter_comparisons.pdf",sep="/"), width = 11, height = 12)
    for(what in cspToEstimate){
      plot(parameter,
           what=what,
           samples=samples*proportion.to.keep,
           mixture.name=mixture.labels)
    }
    dev.off()
  }
  
  ## Print Traces
  print("Printing Traces")
  trace <- parameter$getTraceObject()
  pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
  createTracePlots(trace=trace,model=model,genome=genome,numMixtures=numMixtures,samples=samples,samples.percent.keep = 1,mixture.labels = mixture.labels)
  dev.off()

  #plots different aspects of trace
  pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
  plot(mcmc, what = "LogPosterior")
  plot(trace, what = "ExpectedPhi")
  aa <- aminoAcids()

  for(a in aa)
  {
    if (a=="M"||a=="X"||a=="W") next
    accept.trace <- trace$getCodonSpecificAcceptanceRateTraceForAA(a)
    len <- length(accept.trace)
    mean.acceptance <- mean(accept.trace[(len-len*proportion.to.keep):len])
    plot(accept.trace,main=paste0("Acceptace Rate for ",a),xlab="Samples",ylab="Acceptance Rate",type="l")
  }
  acfCSP(parameter,csp="Selection",numMixtures = numMixtures,samples=samples*proportion.to.keep)
  acfCSP(parameter,csp="Mutation",numMixtures = numMixtures,samples=samples*proportion.to.keep)
  dev.off()

  if(run.convergence.tests){
    ## Run convergence tests
    for (i in 1:numMixtures)
    {
      param.diag<-convergence.test(trace,samples=samples*proportion.to.keep,thin = thinning,what="Selection",mixture=i,frac1=0.1)
      z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
      if (length(z.scores) > 0)
      {
        param.conv <- FALSE
      }
      write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_eta_",i,".txt"),ncolumns = 1)
    }


    for (i in 1:numMixtures)
    {
      param.diag<-convergence.test(trace,samples=samples*proportion.to.keep,thin = thinning,what="Mutation",mixture=i,frac1=0.1)
      z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
      if (length(z.scores) > 0)
      {
        param.conv <- FALSE
      }
      write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_m_",i,".txt"),ncolumns = 1)
    }



    diag <- convergence.test(mcmc,samples = samples*proportion.to.keep,thin=thinning,frac1=0.2)
    z<-abs(diag$z)
    done <- (z > 1.96) && param.conv
  }
  
  # clear out objects
  rm(parameter)
  rm(trace)
  rm(model)

  print(paste0("Completed running and plotting run ", run.number)) 
  # Update run index
  run.number = run.number+1
}




  
while((!done) && (run.number <= 3))
{
  parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"),model="ROC")
  run.number <- run.number + 1
  dir_name <- paste0(directory,"/run_",run.number)
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
  
  createParameterOutput(parameter = parameter,numMixtures = numMixtures,samples = samples,mixture.labels = mixture.labels,samples.percent.keep = proportion.to.keep,relative.to.optimal.codon = F,report.original.ref = T)
  
  expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*proportion.to.keep)
  write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)
  
  #
  #plots different aspects of trace
  trace <- parameter$getTraceObject()
  pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
  plot(mcmc,what = "LogPosterior")
  plot(trace, what = "ExpectedPhi")
  aa <- aminoAcids()
  for(a in aa)
  {
    if (a=="M"||a=="X"||a=="W") next
    accept.trace <- trace$getCodonSpecificAcceptanceRateTraceForAA(a)
    len <- length(accept.trace)
    mean.acceptance <- mean(accept.trace[(len-len*proportion.to.keep):len])
    plot(accept.trace,main=paste0("Acceptace Rate for ",a),xlab="Samples",ylab="Acceptance Rate",type="l")
  }
  acfCSP(parameter,csp="Selection",numMixtures = numMixtures,samples=samples*proportion.to.keep)
  acfCSP(parameter,csp="Mutation",numMixtures = numMixtures,samples=samples*proportion.to.keep)
  dev.off()
  
  
  for (i in 1:numMixtures)
  {
    param.diag<-convergence.test(trace,samples=samples*proportion.to.keep,thin = thinning,what="Selection",mixture=i,frac1=0.1)
    z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
    if (length(z.scores) > 0)
    {
      param.conv <- FALSE
    }
    write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_eta_",i,".txt"),ncolumns = 1)
  }
  
  
  for (i in 1:numMixtures)
  {
    param.diag<-convergence.test(trace,samples=samples*proportion.to.keep,thin = thinning,what="Mutation",mixture=i,frac1=0.1)
    z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
    if (length(z.scores) > 0)
    {
      param.conv <- FALSE
    }
    write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_m_",i,".txt"),ncolumns = 1)
  }
  pdf(paste(dir_name,"Graphs/Parameter_comparisons.pdf",sep="/"), width = 11, height = 12)
  plot(parameter,what="Mutation",samples=samples*proportion.to.keep,mixture.name=mixture.labels)
  plot(parameter,what="Selection",samples=samples*proportion.to.keep,mixture.name=mixture.labels)
  dev.off()
  
  
  pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
  createTracePlots(trace=trace,model=model,genome=genome,numMixtures=numMixtures,samples=samples,samples.percent.keep = proportion.to.keep,mixture.labels = mixture.labels)
  dev.off()
  writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
  writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))
  
  diag <- convergence.test(mcmc,samples = samples*proportion.to.keep,thin=thinning,frac1=0.1)
  z<-abs(diag$z)
  done <- (z > 1.96) && param.conv
  rm(parameter)
  rm(trace)
  rm(model)
}

samples <- 10000
thinning <- 5
parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"),model="ROC")
run.number <- run.number + 1
dir_name <- paste0(directory,"/final_run")
dir.create(dir_name)
dir.create(paste(dir_name,"Graphs",sep="/"))
dir.create(paste(dir_name,"Restart_files",sep="/"))
dir.create(paste(dir_name,"Parameter_est",sep="/"))
dir.create(paste(dir_name,"R_objects",sep="/"))

mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE,est.mix=FALSE)


#this part set steps 
adaptiveSamples=samples*thinning*adaptive.proportion
mcmc$setStepsToAdapt(adaptiveSamples)

model <- initializeModelObject(parameter, "ROC", with.phi)
setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, F)
#run mcmc on genome with parameter using model
#p<-profmem({
sys.runtime <- system.time(
  runMCMC(mcmc, genome, model, num_threads)
)
sys.runtime <- data.frame(Value=names(sys.runtime),Time=as.vector(sys.runtime))
write.table(sys.runtime,file=paste(dir_name,"mcmc_runtime.csv",sep="/"),sep=",",col.names = T,row.names = T,quote=F)


createParameterOutput(parameter = parameter,numMixtures = numMixtures,samples = samples,mixture.labels = mixture.labels,samples.percent.keep = 1,relative.to.optimal.codon = F,report.original.ref = T)

# mixtureAssignment <- getMixtureAssignmentEstimate(parameter,c(1:size),samples*proportion.to.keep)
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


