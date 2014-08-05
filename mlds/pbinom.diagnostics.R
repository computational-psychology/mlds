pbinom.diagnostics <- function(obj, nsim=200, type = "deviance", no.warn = TRUE, workers="localhost", master="localhost", username='guille') {
  
  
  # working with snow
  cl <- makeSOCKcluster(workers, user=username, master= master)
  # initialize cluster and 
  
  # loads MLDS package on each of the workers
  clusterEvalQ(cl, library(MLDS))
  clusterEvalQ(cl, library(psyphy))
  
  if (no.warn){
    old.opt <- options(warn = -1)
    on.exit(options(old.opt))
  }
  n <- length(fitted(obj))
  d <- as.mlds.df(obj$obj$data)  
  
  # export needed variables to each worker 
  clusterExport(cl, c("obj", "n", "d", "type"), envir= environment() )

  # apply run in cluster
  res <- parSapply(cl, seq_len(nsim), function(x, obj){
    ys <- rbinom(n, 1, fitted(obj))
    d$resp <- ys
    br <- mlds(d, lnk= obj$link)
    rs <- residuals(br$obj, type = type)
    rsd <- sort(rs)
    fv.sort <- sort(fitted(br), index.return = TRUE)
    rs <- rs[fv.sort$ix]
    rs <- rs > 0
    runs <- sum(rs[1:(n-1)] != rs[2:n])
    list(resid = as.vector(rsd), NumRuns = runs)
  }, obj = obj
  )
  
  # after all done, stops cluster
  stopCluster(cl)
  
  fres <- list(NumRuns = sapply(seq(2, length(res), 2), 
                                function(x) res[[x]]))
  fres$resid <- t(do.call("cbind", 
                          list(sapply(seq(1, length(res), 2), 
                                      function(x) res[[x]]))))
  fres$resid <- apply(fres$resid, 2, sort)
  fres$Obs.resid <- residuals(obj$obj, type = type)
  rs <- residuals(obj$obj, type = type)
  fv.sort <- sort(fitted(obj), index.return = TRUE)
  rs <- rs[fv.sort$ix]
  rs <- rs > 0
  obs.runs <- sum(rs[1:(n-1)] != rs[2:n])
  nr <- sum(fres$NumRuns > obs.runs) 
  fres$ObsRuns <- obs.runs
  fres$p <- 1 - nr/nsim
  class(fres) <- c("mlds.diag", "list")
  fres
  
  
}
