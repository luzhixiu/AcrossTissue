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
}

samples <- 10000
thinning <- 5
parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"),model="ROC")
run_number <- run_number + 1
dir_name <- paste0(directory,"/final_run")
dir.create(dir_name)
dir.create(paste(dir_name,"Graphs",sep="/"))
dir.create(paste(dir_name,"Restart_files",sep="/"))
dir.create(paste(dir_name,"Parameter_est",sep="/"))
dir.create(paste(dir_name,"R_objects",sep="/"))

mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE,est.mix=FALSE)

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

