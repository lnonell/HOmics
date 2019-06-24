#' filters results using p-value thresholds
#'
#' @param res resultslist
#' @param param parameter of the table to filter with (vars_n, aic, Cor2, p or adj.p). Default is "adj.p"
#' @param comp comparison param comp threshold. Default is "below"
#' @param threshold numerical, the threshold to select associated variables related to the specified param. Default is 0.05
#' @param as.data.frame collapse as data.frame. Default = TRUE

#' @export getsignif


getsignif <- function (res, param = "p.pos", threshold = 0.05, as.data.frame = T)
{
  
  if (!param %in% c("p.pos","p.neg"))
    stop("param has to be one of the following parameters: p.pos or p.neg")
  
  if (!(is.null(threshold))) {
   # res.f <- sapply(res, function(x) {x[x[,param]<threshold,]})
    res.f <- map(res, function(x) filter(x,!!as.name(param)<threshold))
    res.f <- keep(res.f,function(x) nrow(x)>0)
    if (as.data.frame)  res.f <- bind_rows(res.f,.id="gene")
    return(res.f)
  } 
}

