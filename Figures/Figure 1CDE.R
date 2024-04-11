load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/manu_plots/Scripts/Final_plots/Uni_dispersion_beta1_actual.RData")
#load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/manu_plots/Scripts/Final_plots/ERICH1_uni_result.RData")
#setwd("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/manu_plots/Scripts")
pdf("Fig1-lower_panels_dispersion_real_data_simulation.pdf", width = 12/4*3, height = 3)
postscript(file = "Figure 1CED.eps", width = 12/4*3, height = 3, horizontal = FALSE, onefile = FALSE, paper ="special",
           pointsize = 12)
par(mfrow = c(1, 3),     # 2x2 layout
    oma = c(0, 0, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 3.5, 0.5), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0) # axis label at 2 rows distance, tick labels at 1 row
) 

load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/manu_plots/Scripts/Final_plots/ERICH1_uni_result.RData")
plot(pos, Res_uni[,3,2], pch = 19, 
     xlab = "Genomic Positions",  xaxt = "n", ylab ="Multiplicative Dispersion Estimates")
#grid()
axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = round(seq((min(pos)), (max(pos)),length.out = 10)) , tck= -0.02)
abline(h=1, lty = 4)
mtext(side=3, line = 0.3 , at = 637080, adj =0, cex = 1.2, "Data application")


load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/manu_plots/Scripts/Final_plots/Uni_dispersion_beta1_actual.RData")
plot(pos, Res_uni[,3,2,2], pch = 19, xlab = "Genomic Positions", 
     ylim = c(0,20), ylab = "Multiplicative Dispersion Estimates", xaxt = "n")

axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = round(seq((min(pos)), (max(pos)),length.out = 10)) , tck= -0.02)

abline(h = 1, lty = 4)
mtext(side=3, line = 0.3 , at = 637080, adj =0, cex = 1.2, 
      "A simulated dataset \nwithout subject-level RE")

plot(pos, Res_uni[,3,2,8], pch = 19, xlab = "Genomic Positions", 
     ylim = c(0,20), ylab = "Multiplicative Dispersion Estimates", xaxt = "n")

axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = round(seq((min(pos)), (max(pos)),length.out = 10)) , tck= -0.02)

abline(h = 1, lty = 4)
mtext(side=3, line = 0.3 , at = 637080, adj =0, cex = 1.2, "A simulated dataset \nwith subject-level RE")


dev.off()

