

#setwd("~/myscratch/kaiqiong.zhao/Projects/Dispersion_simu_results/WithError/qqplot")
library(ggplot2)



load("~/scratch/kaiqiong.zhao/Projects/Dispersion_simu_results/WithError/qqplot/Pval_mlog10_long_quasi_only_fixed_transform_correct.RData")

res_q = Pval_mlog10_long

load("~/scratch/kaiqiong.zhao/Projects/Dispersion_simu_results/WithError/qqplot/Pval_mlog10_long_RE_only.RData")

res_re = Pval_mlog10_long

load("~/scratch/kaiqiong.zhao/Projects/Dispersion_simu_results/WithError/qqplot/Pval_mlog10_long_both_fix_transform_correct.RData")

res_both = Pval_mlog10_long

#load("~/scratch/kaiqiong.zhao/Projects/Dispersion_simu_results/WithError/qqplot/Pval_mlog10_long_neither.RData")

#res_neither = Pval_mlog10_long



RES = data.frame(cbind(rbind(res_q[res_q$Method=="SOMNiBUS_Fletcher",],
                             res_re[res_re$Method=="SOMNiBUS_Fletcher",],
                             res_both[res_both$Method=="SOMNiBUS_Fletcher",]),
                       Model = rep(c("Quasi only", "RE only", "Both"), each = nrow(res_q[res_q$Method=="SOMNiBUS_Fletcher",]))
))

RES$Model = factor(RES$Model, levels = c("Both",  "Quasi only" , "RE only"))
trop = c('#1b9e77','#d95f02','#7570b3', "dodgerblue")


postscript(file = "Figure 4B.eps", width = 8, height = 5, horizontal = FALSE, onefile = FALSE, paper ="special",
           pointsize = 12)

ggplot(data = RES[RES$Model!="Both",], 
       aes(x = exp_p, y = obs_p,  color = Model)) +
  geom_point( size = 0.7, alpha =1) + 
  geom_point(data = subset(RES, Model == 'Both'), size = 0.5, alpha = 1)+
  xlab(NULL) +
  facet_grid(phi~sigma2, labeller = label_bquote(cols = sigma[0]^2==.(sigma2),
                                                 rows = phi==.(phi)))+
  scale_y_continuous(name=expression(Observed~~-log[10](italic(p))), limits = c(0,3))+
  scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))+
  geom_abline(linetype = 4)+
  scale_color_manual(values = trop, labels = c('dSOMNiBUS','multiplicative-dispersion-only model',
                                               'additive-dispersion-only model', "SOMNiBUS"), name="" )+
  theme(legend.position="bottom")

dev.off()
ggsave("1-impact-of-dispersion_fixedScale_transform_correct.pdf", width = 8, height = 5)


# change alpha = 0.5 to alpha = 1 for make sure correct postscript eps export




