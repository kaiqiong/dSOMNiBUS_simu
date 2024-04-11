

#setwd("~/scratch/GREENWOOD_SCRATCH/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/Summary_parametric_bootstrap")

# Sasha's 2 result
load("~/scratch/kaiqiong.zhao/Projects/Sasha_new/Sasha_new/Analysis1_RE/Results/Summary/Ordered_res_12983.RData")

RES2 <- RES_order


# Sasha's 1 result

load("~/scratch/kaiqiong.zhao/Projects/Robust_SOMNiBUS/Sasha_data/Analysis1_keep_RE/Results/Summary/Sasha_1_10759.RData")

RES1 <- RES_order



res1 <-RES1[,c("sigma00", "phi")]
res2 <- RES2[,c("sigma00", "phi")]


res_all <- data.frame(cbind(
  rbind(res1, res2), Data = 
    c(rep("Data 1", nrow(res1)), rep("Data 2", nrow(res2)))
))



res_all$logsig <- log(res_all$sigma00)
res_all$logphi <- log(res_all$phi)

#####---------


library(patchwork)



bin2d = ggplot(data=res_all[res_all$Data=="Data 2",], aes(x=logphi, y=logsig) ) +
  #geom_bin2d() +
  scale_y_continuous(name=expression(log(hat(sigma[0]^2))),limits = c(-20,15))+
  scale_x_continuous(name=expression(log(hat(phi))), limits = c(-0.5,3))+
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  geom_hline(yintercept=log(0.05), linetype = 4)+
  geom_vline(xintercept=0, linetype = 4)+
  theme_bw()+
  theme(legend.position=c(0.85, 0.78), legend.text = element_text(size = 7), 
        legend.title=element_text(size=7),
        legend.key.size = unit(0.5, 'cm'))
bin2d


hist_phi<- ggplot(res_all[res_all$Data=="Data 2",], aes(x=logphi)) + 
  geom_histogram(color = "black", fill = "white") + 
  scale_x_continuous(limits = c(-0.5,3))+
  theme_void()


hist_sig <- ggplot(res_all[res_all$Data=="Data 2",], aes(x=logsig)) + 
  geom_histogram(color = "black", fill = "white") + 
  scale_x_continuous(limits = c(-20,15))+
  theme_void()+  coord_flip()
hist_phi + plot_spacer() + bin2d + hist_sig + 
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(4, 1),
    heights = c(1, 4)
  ) 
ggsave("data2_new.pdf", width = 5, height = 5)
postscript(file = "Figure 2A.eps", width = 5, height = 4, horizontal = FALSE, onefile = FALSE, paper ="special",
           pointsize = 12)
hist_phi + plot_spacer() + bin2d + hist_sig + 
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(4, 1),
    heights = c(1, 4)
  ) 
dev.off()

