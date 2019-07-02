#' plots results in terms of 95CI
#'
#' @param res resultslist
#' @param element number or name of group of features to plot. Default = 1

#' @import dplyr
#' @import ggplot2
#' @export plot.res



plot.res <- function (res, element = 1)
{
  if (class(res)!="HOmics") stop("res must be an HOmics class object")
  
   resi <- res[[element]]
   p <- ggplot(data=resi) +
     geom_segment(aes(x=`2.5%`,y=feature,xend=`97.5%`,yend=feature),
                  arrow=arrow(length=unit(0.1,"cm"),
                              ends='both'),size=0.5) +
      geom_vline(linetype   ='dashed',  xintercept = 0) +
      xlab ("95 % probability") + 
      theme(panel.border = element_blank(),
            axis.ticks.y = element_blank())
   p

}

