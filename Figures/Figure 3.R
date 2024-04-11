

load("~/Documents/Projects_since_2021/dSOMNiBUS-add-application/Rcode/Summary/RES_all_methods.RData")
#RES_Biseq_GlobalTest
#RES_smooth_smsc
load("~/Documents/Projects_since_2021/dSOMNiBUS-add-application/Rcode/Summary/dSOMNiBUS_old_testStats.RData")
#RES_old
load("~/Documents/Projects_since_2021/dSOMNiBUS-add-application/Rcode/Summary/dSOMNiBUS_old_testStats_with_error.RData")
#RES_error
load("~/Documents/Projects_since_2021/dSOMNiBUS-add-application/Rcode/Summary/dSOMNiBUS_old_testStats_with_error_use_trc_median.RData")
#RES_error_trc_mean

qqPlot <- function(pval, truncate = FALSE, ylim=NULL, thinThreshold=NULL, ci=TRUE, ...) 
{
  # pvalue is a vector of p-values
  # truncate is T/F whether to truncate the y-axis (observed) to the same limits as the x-axis (expected)
  # thin is whether to thin insignificant p-values
  
  pval <- -log10(sort(pval)) # sort() removes NAs
  n <- length(pval)
  a <- 1:n
  b <- a/n
  x <- -log10(b)
  
  if (!is.null(thinThreshold)){
    breaks <- seq(from=0, to=thinThreshold, length.out=11)
    pval.cut <- cut(pval, breaks=breaks, right=FALSE)
    #quant <- quantile(pval[pval < thinThreshold], probs=1:10/10)
    #pval.cut <- cut(pval, breaks=c(0, quant), right=FALSE)
    ind.list <- list()
    for (level in levels(pval.cut)){
      sel <- which(pval.cut == level)
      ind.list[[level]] <- sel[sample.int(length(sel), min(1000, length(sel)))]
    }
    ind.list[["max"]] <- which(pval >= thinThreshold)
    ind <- unlist(ind.list, use.names=FALSE)
    ind <- sort(ind) # sorting necessary for polygon, below
  } else {
    ind <- 1:n
  }
  
  
  char <- rep(19,n)
  if(!is.logical(truncate) | truncate){
    if (is.logical(truncate)){
      maxx <- max(x)+2  
    } else {
      maxx <- min(truncate, max(pval))
    }
    
    ylm <- c(0,maxx)
    ylb <- expression(paste(-log[10], "(observed P) - truncated"))
    nx <- length(which(pval > maxx))
    if(nx > 0){
      pval[1:nx] <- maxx
      char[1:nx] <- 2
    }
    
  } else {
    ylm <- ylim
    ylb <- expression(Observed~~-log[10](italic(p)))
  }
  plot(x[ind], pval[ind], type = "n", ylim = ylm, ylab = ylb,pch = 19,
       xlab = expression(Expected~~-log[10](italic(p))), ...)
  # upper and lower have already been subset
  if (ci){
    upper <- qbeta(0.025, a[ind], rev(a)[ind])
    lower <- qbeta(0.975, a[ind], rev(a)[ind])
    
    polygon(-log10(c(b[ind], rev(b[ind]))), -log10(c(upper, rev(lower))), density=NA, col="gray")
  }
  points(x[ind], pval[ind], pch = char, ...)  
  abline(0,1,col="red")  
}
library(QCEWAS)


load("~/Documents/Projects_since_2021/dSOMNiBUS-add-application/Rcode/Summary/Tab_all/Res_13004_genes_6_methods_trc_median_autosomes_only.RData")
setwd("~/Documents/Projects_2024/dSOMNiBUS_simu/Figures-code")
#setwd("~/Documents/Projects_since_2021/dSOMNiBUS-add-application/Rcode/Summary/Fig2-qqplot-all-methods")
pdf(file ="Allqqplots_with_lams_autosomes_only.pdf", width = 10, height = 5)

postscript(file = "Figure 3.eps", width = 10, height = 5, horizontal = FALSE, onefile = FALSE, paper ="special",
           pointsize = 12)

par(mfrow = c(2, 3),     # 2x2 layout
    oma = c(0, 0, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 2, 0.5), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0) # axis label at 2 rows distance, tick labels at 1 row
) 


pval1 = RES_old$`p-value`[which(RES_old$gene %in% anno1$gene)] 
qqPlot(pval1, main = "dSOMNiBUS (zero error)")
lam = round(P_lambda(pval1), 2 )
text(0.5, max(-log10(pval1), na.rm= T), bquote(lambda["GC"] == .(lam)))
pval2=RES_error_trc_mean$`p-value`[which(RES_error_trc_mean$gene%in% anno1$gene)] 
qqPlot(pval2,  main = "dSOMNiBUS (non-zero error)")
lam = round(P_lambda(pval2), 2 )
text(0.5, max(-log10(pval2), na.rm= T), bquote(lambda["GC"] == .(lam)))
pval3 = RES_Biseq_GlobalTest$BiSeq[which(RES_Biseq_GlobalTest$gene %in% anno1$gene)] 
qqPlot(pval3,  main = "BiSeq")
lam = round(P_lambda(pval3), 2 )
text(0.5, max(-log10(pval3), na.rm= T), bquote(lambda["GC"] == .(lam)))
pval4=RES_Biseq_GlobalTest$GlobalTest[which(RES_Biseq_GlobalTest$gene %in% anno1$gene)] 
qqPlot(pval4,  main = "GlobalTest")
lam = round(P_lambda(pval4), 2 )
text(0.5, max(-log10(pval4), na.rm= T), bquote(lambda["GC"] == .(lam)))
pval5=RES_smooth_smsc$BSmooth[which(RES_smooth_smsc$gene%in% anno1$gene)] 
qqPlot(pval5, main = "BSmooth")
lam = round(P_lambda(pval5), 2 )
text(0.5, max(-log10(pval5[pval5>0]), na.rm= T), bquote(lambda["GC"] == .(lam)))
pval6 = RES_smooth_smsc$SMSC[which(RES_smooth_smsc$gene%in% anno1$gene)] 
qqPlot(pval6, main = "SMSC")
lam = round(P_lambda(pval6), 2 )
text(0.5, max(-log10(pval6[pval6>0]), na.rm= T), bquote(lambda["GC"] == .(lam)))

dev.off()

