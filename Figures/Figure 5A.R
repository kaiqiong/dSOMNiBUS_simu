









myqq <- function(pval.1){
  pval.1_full <- pval.1[!is.na(pval.1)]
  # p.exp <-  -log10((1:length(pval.1_full ))/(length(pval.1_full )+1 ))
  p.obs <- -log10(pval.1_full )[(order(pval.1_full ))]
  
  p.obs <- c(p.obs, rep(NA, length(pval.1) -  length(pval.1_full)))
  return((p.obs))
}

myqq_exp <- function(pval.1){
  pval.1_full <- pval.1[!is.na(pval.1)]
  p.exp <-  -log10((1:length(pval.1_full ))/(length(pval.1_full )+1 ))
  #p.obs <- -log10(pval.1_full )[(order(pval.1_full ))]
  
  p.obs <- c(p.exp, rep(NA, length(pval.1) -  length(pval.1_full)))
  return((p.obs))
}



load("~/scratch/kaiqiong.zhao/Projects/Dispersion_simu_results/WithError/qqplot/FinalPval_mlog10_all_6_methods_fix_transform_correct.RData")

RES_ALL <- RES_ALL[RES_ALL$Method!= "dmrseq-minp",]
#------------------
# Try use QQ plots
#------------------

#library(data.table)
library(ggplot2)
#library(gridExtra)
#library(pubr)


#postscript(file = "QQplot.eps", width = 7, height = 5, horizontal = FALSE, onefile = FALSE,
#           paper ="special", pointsize = 12)
#trop = c('#1b9e77','hotpink')
#trop = c("#e6ab02", "#e7298a","#66a61e")
#trop = c('#1b9e77', '#e66101','#fdb863','#b2abd2','#5e3c99', "red", "blue")
#trop = c('#1b9e77','#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd')
#trop = c('#1b9e77','#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462')
#trop=c( '#1b9e77', '#d7191c','#fdae61','#ffffbf','#a6d96a','#1a9641','#fdb462')


cols= c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99')



cols = c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')
#cols = c(rev(sequential_hcl(5, "Red-Purple")[1:3]), '#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
#cols <- rev(cols)


cols = c('#1b9e77', cols[c(2,4,6,8,10, 12)])
library(grDevices)
library(RColorBrewer)
library(colorspace)
#pal <- choose_palette()
#cols = diverge_hcl(14)

RES_ALL$obs_p_ori <- 10^(-RES_ALL$obs_p)
RES_ALL$exp_p_ori <- 10^(-RES_ALL$exp_p)
RES_ALL$Method <- factor(RES_ALL$Method, levels = c(
  "SOMNiBUS_Fletcher", "GlobalTest", "dmrseq-fisherp",  "BSmooth", "SMSC", "BiSeq"
))
#trop = c('#1b9e77',"darkorange", "dodgerblue", "hotpink","limegreen", "#984ea3","red")
trop = c('#1b9e77',"darkorange", "red", "dodgerblue","hotpink", "#984ea3")

postscript(file = "Figure 5A.eps", width = 8, height = 5, horizontal = FALSE, onefile = FALSE, paper ="special",
           pointsize = 12)

ggplot(data = RES_ALL, 
       aes(x = exp_p, y = obs_p,  color = Method)) +
  geom_point(size = 0.2, alpha =1) + 
  geom_point(data = subset(RES_ALL, Method == 'SOMNiBUS_Fletcher'), size = 0.3, alpha =1)+

  xlab(NULL) +
  facet_grid(phi~sigma2, labeller = label_bquote(cols = sigma[0]^2==.(sigma2),
                                                 rows = phi==.(phi)))+
  scale_y_continuous(name=expression(Observed~~-log[10](italic(p))), limits = c(0,4))+
  scale_x_continuous(name=expression(Expected~~-log[10](italic(p))), limits = c(0,3))+
  #  scale_y_continuous(limits = c(0,4))+
  scale_color_manual(labels = c("dSOMNiBUS", "GlobalTest", "dmrseq", "BSmooth", "SMSC", "BiSeq"), 
                     values = trop, name = "")+
  geom_abline(linetype = 4)+
  theme(legend.position="bottom")+
  guides(col = guide_legend(nrow = 1))

dev.off()
#setwd("~/myscratch/kaiqiong.zhao/Projects/Dispersion_simu_results/WithError/qqplot")
ggsave("3-qqplot-6-methods-fixedEM-transform-correct-dmrseq-fisherp.pdf", width = 8, height = 5)


