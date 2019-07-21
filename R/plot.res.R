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
  if (!is.numeric(element) & !(element %in% names(res))) stop("element must be a numeric value or one of the element names of res")

   resi <- res[[element]]
   resi$feature <- factor(resi$feature, levels=sort(unique(resi$feature),decreasing = TRUE)) 
   p <- ggplot(data=resi) +
     geom_segment(aes(x=`2.5%`,y=feature,xend=`97.5%`,yend=feature),
                  arrow=arrow(length=unit(0.15,"cm"),
                              ends='both')) +
      geom_vline(linetype   ='dashed',  xintercept = 0) +
      xlab ("95 % probability") + 
      theme(panel.border = element_blank(),
            axis.ticks.y = element_blank())
   p

}

