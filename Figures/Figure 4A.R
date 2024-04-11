# Load the output from the simulation
load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/fixedScaleEM_transform_correct_piS/Exp_phi_3_re_3_quasi_RE_newb0_em_correct_F/Res_all_Samp_100_Simu1000.RData")
load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BANK1betas.RData")
beta.0 <- beta.0/4 + 2

beta_all <- cbind(beta.0, beta.1, beta.2, rep(0, length(beta.0)))

beta_all[,1] <- beta_all[,1]-3.5
beta.0 <- beta.0 - 3.5
# Plot type is line


#setwd("~/myscratch/kaiqiong.zhao/Projects/Dispersion_simu_results/WithError")


#postscript(file = "Finalshadow-plot_phi_3_re_3_quasi_RE.eps", width = 7, height = 5, horizontal = FALSE, onefile = FALSE, paper ="special",
#           pointsize = 12)

postscript(file = "Figure 4A.eps", width = 7, height = 5, horizontal = FALSE, onefile = FALSE, paper ="special",
           pointsize = 12)

#pdf(file="Finalshadow-plot_phi_3_re_3_quasi_RE.pdf", width = 7, height = 5)
par(mfrow = c(2, 2),     # 2x2 layout
    oma = c(0, 0, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 2, 0.5), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0) # axis label at 2 rows distance, tick labels at 1 row
) 
plot(pos[order(pos)], beta.0[order(pos)], col="red", xaxt ="n",  type="l", ylim=c(min(Beta.est[,,1], na.rm = T), max(Beta.est[,,1], na.rm = T)),
     xlab="Genomic Position", ylab=" ", main = expression(paste(beta[0], "(t)")))
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)
for( i in 1:M){
  lines(  pos.index[order(pos.index[,i]),i], Beta.est[order(pos.index[,i]),i , 1],pch = as.character(i),col="gray") 
}
lines(pos[order(pos)], beta.0[order(pos)], col="red")

plot(pos[order(pos)], beta.1[order(pos)],col=2, xaxt ="n", ylim=c(min(Beta.est[,,2], na.rm = T), max(Beta.est[,,2], na.rm = T)), type="l",
     xlab="Genomic Position", ylab=" ", main = expression(paste(beta[1], "(t)")))
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)
for( i in 1:M){
  lines(pos.index[order(pos.index[,i]),i], Beta.est[order(pos.index[,i]),i,2 ], pch = as.character(i), col="gray") 
}
lines(pos[order(pos)], beta.1[order(pos)], col="red")
abline(h=0, lty=4)



plot(pos[order(pos)], beta.2[order(pos)], col=2, xaxt ="n",  type="l",
     xlab="Genomic Position", ylab=" ",  main = expression(paste(beta[2], "(t)")) , ylim=c(min(Beta.est[,,3], na.rm = T), max(Beta.est[,,3], na.rm = T)))  #ylim = c(min(beta.2.est), max(beta.2.est))
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)

axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)

for( i in 1:M){
  lines(pos.index[order(pos.index[,i]),i], Beta.est[order(pos.index[,i]),i,3], pch = as.character(i), col="gray") 
}
lines(pos[order(pos)], beta.2[order(pos)], col="red")
abline(h=0, lty=4)


beta.3 = rep(0, length(pos))
plot(pos[order(pos)], beta.3[order(pos)], col=2, xaxt ="n", ylim=c(min(Beta.est[,,4], na.rm = T), max(Beta.est[,,4], na.rm = T)), type="l",
     xlab="Genomic Position", ylab=" ", main = expression(paste(beta[3], "(t)")))
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = seq(round(min(pos)), round(max(pos)),length.out = 10 ) , tck= -0.02)

for( i in 1:M){
  lines(pos.index[order(pos.index[,i]),i], Beta.est[order(pos.index[,i]),i, 4], pch = as.character(i), col="gray") 
}
lines(pos[order(pos)], beta.3[order(pos)], col="red")
abline(h=0, lty=4)

dev.off()
