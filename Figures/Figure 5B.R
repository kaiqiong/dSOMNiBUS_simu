load("~/scratch/kaiqiong.zhao/Projects/Dispersion_simu_results/WithError/power/Power_fixedEM_transform_6methods.RData")
#------------------
# Try use QQ plots
#------------------

#library(data.table)
library(ggplot2)
#library(gridExtra)
#library(pubr)


trop = c('#1b9e77',"darkorange", "red", "dodgerblue","hotpink", "#984ea3")

#postscript(file = "QQplot.eps", width = 7, height = 5, horizontal = FALSE, onefile = FALSE,
#           paper ="special", pointsize = 12)

#trop = c("darkorange", "dodgerblue","hotpink")
#trop = c('#d73027','#fc8d59','#fee08b','#d9ef8b','#91cf60','#1a9850')
#trop = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')[rev(c(1:6))]
#trop[c(5,6,2)] <-  c('#1b9e77','#d95f02','#7570b3')
#trop[4] <- '#66a61e'

#            Biseq,              "GlobalTest" "SMSC"
#trop <- c('#1b9e77', "#e6ab02", "#7570b3", "#e7298a", "#66a61e" ,"#1b9e77", "#d95f02")
#trop = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33')[rev(c(1:6))]

#ggsave("Power_6methods_fix_transforme_dmrseq_fisher.pdf", width = 10, height = 5)
postscript(file = "Figure 5B.eps", width = 10, height = 5, horizontal = FALSE, onefile = FALSE, paper ="special",
           pointsize = 12)

ggplot(data = Power_long, 
       aes(x = dif, y = power,  color = Method)) +
  geom_line(size = 0.5)+
  geom_line(data = subset(Power_long, Method == 'SOMNiBUS_Fletcher'), size = 0.7)+
  xlab(NULL) +
  facet_grid(phi~sigma2, labeller = label_bquote(cols = sigma[0]^2==.(sigma2),
                                                 rows = phi==.(phi)))+
  scale_y_continuous(name="Power")+
  scale_x_continuous(name="Methylation Differences")+
  #  scale_y_continuous(limits = c(0,4))+
  scale_color_manual(values = trop, labels = c('dSOMNiBUS','GlobalTest', "dmrseq",
                                               'BSmooth', 'SMSC', 'BiSeq'), name="" )+
  theme(legend.position="bottom")+
  guides(col = guide_legend(nrow = 1))
dev.off()
#setwd("~/myscratch/kaiqiong.zhao/Projects/Dispersion_simu_results/WithError/power")
ggsave("Power_6methods_fix_transforme_dmrseq_fisher.pdf", width = 10, height = 5)


# The code ends here

