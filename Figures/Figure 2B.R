myqq <- function(pval.1){
  p.exp <-  -log10((1:length(pval.1))/(length(pval.1)+1 ))
  p.obs <- -log10(pval.1)[(order(pval.1))]
  return((p.obs))
}



load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/Summary_parametric_bootstrap/Sasha1_all_res.RData")



reg_both<- data.frame(obs_p = c(Region_p_all[,1], Region_p_all[,2], Region_p_all[,3]),
                     obs_p_log10 = c(myqq(Region_p_all[,1]), myqq(Region_p_all[,2]), myqq(Region_p_all[,3]) ),
                     exp_p_log10 = c( rep(-log10((1:nrow(Region_p_all))/(nrow(Region_p_all)+1 )) , 3)),
                     pboot_p = c(Region_emp_p_all[,1], Region_emp_p_all[,2], Region_emp_p_all[,3]),
                     method = rep(c("SOMNiBUS_Fletcher", "SOMNIBUS_REML", "mgcv_gam"), each = nrow(Region_p_all)))
#ordered(sizes, levels = c("small", "medium", "large"))


load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/Summary_parametric_bootstrap/Sasha1_Quasionly_all_res.RData")


reg_quasi<- data.frame(obs_p = c(Region_p_all[,1], Region_p_all[,2], Region_p_all[,3]),
                      obs_p_log10 = c(myqq(Region_p_all[,1]), myqq(Region_p_all[,2]), myqq(Region_p_all[,3]) ),
                      exp_p_log10 = c( rep(-log10((1:nrow(Region_p_all))/(nrow(Region_p_all)+1 )) , 3)),
                      pboot_p = c(Region_emp_p_all[,1], Region_emp_p_all[,2], Region_emp_p_all[,3]),
                      method = rep(c("SOMNiBUS_Fletcher", "SOMNIBUS_REML", "mgcv_gam"), each = nrow(Region_p_all)))



load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/Summary_parametric_bootstrap/Sasha1_REonly_all_res.RData")

reg_RE<- data.frame(obs_p = c(Region_p_all[,1], Region_p_all[,2], Region_p_all[,3]),
                       obs_p_log10 = c(myqq(Region_p_all[,1]), myqq(Region_p_all[,2]), myqq(Region_p_all[,3]) ),
                       exp_p_log10 = c( rep(-log10((1:nrow(Region_p_all))/(nrow(Region_p_all)+1 )) , 3)),
                       pboot_p = c(Region_emp_p_all[,1], Region_emp_p_all[,2], Region_emp_p_all[,3]),
                       method = rep(c("SOMNiBUS_Fletcher", "SOMNIBUS_REML", "mgcv_gam"), each = nrow(Region_p_all)))




RES_1_ALL <- data.frame(rbind(reg_quasi, reg_RE, reg_both), variation = c(rep("Dispersion only", nrow(reg_quasi)),
                                                                          rep("RE only", nrow(reg_RE)),
                                                                          rep("Dispersion + RE", nrow(reg_both))
                                                                          ))

RES_1_ALL$method <- as.factor(RES_1_ALL$method)
RES_1_ALL$method <-relevel(RES_1_ALL$method, ref = "SOMNiBUS_Fletcher")
trop = c('#1b9e77','#d95f02','#7570b3')
ggplot(data = RES_1_ALL[RES_1_ALL$method=="SOMNiBUS_Fletcher",], 
       aes(x = exp_p_log10, y = obs_p_log10,  color = variation)) +
  geom_point() + 
  xlab(NULL) +
#  facet_grid(cols= vars(variation), labeller = label_both)+
  scale_y_continuous(name=expression(Observed~~-log[10](italic(p))), limits = c(0,10))+
  scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))+
  #  scale_y_continuous(limits = c(0,4))+
  geom_abline(linetype = 4)+
  scale_color_manual(values = trop)+
  theme(legend.position="bottom")



#----

# Data 2

load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/Summary_parametric_bootstrap/Sasha2_all_res.RData")



reg_both<- data.frame(obs_p = c(Region_p_all[,1], Region_p_all[,2], Region_p_all[,3]),
                      obs_p_log10 = c(myqq(Region_p_all[,1]), myqq(Region_p_all[,2]), myqq(Region_p_all[,3]) ),
                      exp_p_log10 = c( rep(-log10((1:nrow(Region_p_all))/(nrow(Region_p_all)+1 )) , 3)),
                      pboot_p = c(Region_emp_p_all[,1], Region_emp_p_all[,2], Region_emp_p_all[,3]),
                      method = rep(c("SOMNiBUS_Fletcher", "SOMNIBUS_REML", "mgcv_gam"), each = nrow(Region_p_all)))
#ordered(sizes, levels = c("small", "medium", "large"))


load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/Summary_parametric_bootstrap/Sasha2_Quasionly_all_res.RData")


reg_quasi<- data.frame(obs_p = c(Region_p_all[,1], Region_p_all[,2], Region_p_all[,3]),
                       obs_p_log10 = c(myqq(Region_p_all[,1]), myqq(Region_p_all[,2]), myqq(Region_p_all[,3]) ),
                       exp_p_log10 = c( rep(-log10((1:nrow(Region_p_all))/(nrow(Region_p_all)+1 )) , 3)),
                       pboot_p = c(Region_emp_p_all[,1], Region_emp_p_all[,2], Region_emp_p_all[,3]),
                       method = rep(c("SOMNiBUS_Fletcher", "SOMNIBUS_REML", "mgcv_gam"), each = nrow(Region_p_all)))



load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/Summary_parametric_bootstrap/Sasha2_REonly_all_res.RData")

reg_RE<- data.frame(obs_p = c(Region_p_all[,1], Region_p_all[,2], Region_p_all[,3]),
                    obs_p_log10 = c(myqq(Region_p_all[,1]), myqq(Region_p_all[,2]), myqq(Region_p_all[,3]) ),
                    exp_p_log10 = c( rep(-log10((1:nrow(Region_p_all))/(nrow(Region_p_all)+1 )) , 3)),
                    pboot_p = c(Region_emp_p_all[,1], Region_emp_p_all[,2], Region_emp_p_all[,3]),
                    method = rep(c("SOMNiBUS_Fletcher", "SOMNIBUS_REML", "mgcv_gam"), each = nrow(Region_p_all)))




RES_2_ALL <- data.frame(rbind(reg_quasi, reg_RE, reg_both), variation = c(rep("Dispersion only", nrow(reg_quasi)),
                                                                          rep("RE only", nrow(reg_RE)),
                                                                          rep("Dispersion + RE", nrow(reg_both))
))

RES_2_ALL$method <- as.factor(RES_2_ALL$method)
RES_2_ALL$method <-relevel(RES_2_ALL$method, ref = "SOMNiBUS_Fletcher")



RES_12_ALL <- data.frame(rbind(RES_1_ALL, RES_2_ALL), data = c(rep("Data 1", nrow(RES_1_ALL)),
                                                               rep("Data 2", nrow(RES_2_ALL))))





#---
# No RE no Quasi


load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/Summary_parametric_bootstrap/Sasha1_noquai_nore_all_res.RData")


reg_neither1<- data.frame(obs_p = c(Region_p_all[,1], Region_p_all[,2], Region_p_all[,3]),
                      obs_p_log10 = c(myqq(Region_p_all[,1]), myqq(Region_p_all[,2]), myqq(Region_p_all[,3]) ),
                      exp_p_log10 = c( rep(-log10((1:nrow(Region_p_all))/(nrow(Region_p_all)+1 )) , 3)),
                      pboot_p = rep(NA, nrow(Region_p_all)*3),
                      method = rep(c("SOMNiBUS_Fletcher", "SOMNIBUS_REML", "mgcv_gam"), each = nrow(Region_p_all)))

load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/Summary_parametric_bootstrap/Sasha2_noquai_nore_all_res.RData")

reg_neither2<- data.frame(obs_p = c(Region_p_all[,1], Region_p_all[,2], Region_p_all[,3]),
                          obs_p_log10 = c(myqq(Region_p_all[,1]), myqq(Region_p_all[,2]), myqq(Region_p_all[,3]) ),
                          exp_p_log10 = c( rep(-log10((1:nrow(Region_p_all))/(nrow(Region_p_all)+1 )) , 3)),
                          pboot_p = rep(NA, nrow(Region_p_all)*3),
                          method = rep(c("SOMNiBUS_Fletcher", "SOMNIBUS_REML", "mgcv_gam"), each = nrow(Region_p_all)))

RES_12_ALL_vars <- data.frame(
rbind(RES_12_ALL, cbind(reg_neither1, 
                        variation = rep("No dispersion; No RE", nrow(reg_neither1)),
                        data = rep("Data 1", nrow(reg_neither1))
                        ),
      cbind(reg_neither2, 
            variation = rep("No dispersion; No RE", nrow(reg_neither2)),
            data = rep("Data 2", nrow(reg_neither2))
      )
))
#ordered(sizes, levels = c("small", "medium", "large"))


table(RES_12_ALL_vars$variation)


trop = c('#1b9e77','#d95f02','#7570b3', "dodgerblue")

RES_2only <- RES_12_ALL_vars[RES_12_ALL_vars$data=="Data 2",]

RES_2only$variation = factor(RES_2only$variation, levels = c("Dispersion + RE"   , "Dispersion only", "RE only", "No dispersion; No RE"))
trop = c('#1b9e77','#d95f02','#7570b3', "dodgerblue")


library(magrittr)
library(QCEWAS)

lam_d = round(P_lambda(RES_2only$obs_p[RES_2only$variation=="Dispersion + RE"]) ,2)
lam_mdo = round( P_lambda(RES_2only$obs_p[RES_2only$variation=="Dispersion only"]) ,2)
lam_ado = round( P_lambda(RES_2only$obs_p[RES_2only$variation=="RE only"]) ,2)
lam_no = round( P_lambda(RES_2only$obs_p[RES_2only$variation=="No dispersion; No RE"]) ,2)

postscript(file = "Figure 2B.eps", width = 5, height = 5, horizontal = FALSE, onefile = FALSE, paper ="special",
           pointsize = 12)

ggplot(data = RES_2only[RES_2only$method=="SOMNiBUS_Fletcher",], 
       aes(x = exp_p_log10, y = obs_p_log10,  color = variation)) +
  geom_point() + 
  xlab(NULL) +
  scale_y_continuous(name=expression(Observed~~-log[10](italic(p))), limits = c(0,10))+
  scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))+
  #  scale_y_continuous(limits = c(0,4))+
  geom_abline(linetype = 4)+
scale_color_manual(values = trop, labels = c('dSOMNiBUS','MDO model',
                                             'ADO model', "SOMNiBUS"), name="" )+
  theme_minimal()+
  theme(legend.position = "top")+
  annotate("text", x= 1.1, y = 10, label = bquote(lambda["GC"] == .(lam_no)), col=trop[4])+
  annotate("text", x= 1.5, y = 8.8, label = bquote(lambda["GC"] == "4.20"), col=trop[3])+
annotate("text", x= 1.7, y = 7.2, label = bquote(lambda["GC"] == .(lam_mdo)), col=trop[2])+
annotate("text", x= 2.4, y = 6.5, label = bquote(lambda["GC"] == .(lam_d)), col=trop[1])
  
dev.off()



ggsave("Excess_heterogeneity_dat2_only.pdf", width = 5, height = 5)
#ggsave("Excess_heterogeneity.pdf", width = 8, height = 5)



trop = c('#1b9e77','#d95f02','#7570b3', "dodgerblue")
ggplot(data = RES_12_ALL_vars, 
       aes(x = exp_p_log10, y = obs_p_log10,  color = variation)) +
  geom_point() + 
  xlab(NULL) +
  facet_grid(cols= vars(data), rows = vars(method))+
  scale_y_continuous(name=expression(Observed~~-log[10](italic(p))), limits = c(0,10))+
  scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))+
  #  scale_y_continuous(limits = c(0,4))+
  geom_abline(linetype = 4)+
  scale_color_manual(values = trop)+
  theme(legend.position="bottom")
