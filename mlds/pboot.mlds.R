pboot.mlds <- function(x, nsim, no.warn = TRUE, workers="localhost", master="localhost", username='guille') {
  
  
  # working with snow
  cl <- makeSOCKcluster(workers, user=username, master= master)
  #
  #sfInit(parallel= parallel, cpus= cpus)  # initialize cluster and 
  
  #sfLibrary(MLDS)   # loads MLDS package on each of the workers
  clusterEvalQ(cl, library(MLDS))
  clusterEvalQ(cl, library(psyphy))
  
  if (no.warn){
    old.opt <- options(warn = -1)
    on.exit(old.opt)
  }
  d <- if (x$method == "glm") 
    as.matrix(x$obj$data[, -1]) else
      as.matrix(make.ix.mat(x$data)[, -1])
  #				as.mlds.df(x$obj$data) else
  #				x$data
  p <- fitted(x)
  rsim <- matrix(rbinom(length(p) * nsim, 1, p), 
                 nrow = length(p), ncol = nsim)
  
  #sfExport("x", "rsim", "d") # export needed variables to each worker 
  clusterExport(cl, c("x", "rsim", "d"), envir= environment() )

  bts.samp <- parApply(cl, rsim, 2, function(y, dd) {
    #		dd$resp <- x
    #		psct <- mlds(dd, ...)$pscale
    
    psct <- glm.fit(dd, y, family = binomial(x$link))$coefficients
    names(psct) <- x$stimulus[-1]
    c(psct, sigma = 1)/psct[length(psct)]
  },  dd = d)
  
  #sfStop() # after all done, stops cluster
  stopCluster(cl)
  
  res <- list(boot.samp = bts.samp,
              bt.mean = apply(bts.samp, 1, mean),
              bt.sd = apply(bts.samp, 1, sd),
              N = nsim
  )
  class(res) <- c("mlds.bt", "list")
  res
}
