
# Run Fig1-methyprop-variation-6samps-or-allsamps.R first

load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/manu_plots/Scripts/Final_plots/Uni_199_samp_permutations.RData")
load("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/manu_plots/Scripts/Final_plots/ERICH1_uni_result.RData")
#pdf("Fig1-2-univariate_obs_emp_pvals.pdf", width = 4, height = 4)
postscript(file = "Figure 1B.eps", width = 4, height = 4, horizontal = FALSE, onefile = FALSE, paper ="special",
           pointsize = 12)
par(mfrow = c(1, 1),     # 2x2 layout
    oma = c(0, 0, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 2, 0.5), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0) # axis label at 2 rows distance, tick labels at 1 row
) 

plot(Res_uni[,2,1], emp_p[,1], col = "darkorange", xlab ="Observed P-value", 
     ylab= "Empirical P-value", pch = 19, main = "Single-site P-values")
points( Res_uni[,2,2], emp_p[,2], col = "dodgerblue", pch = 19)
abline(0,1, lty = 4)
mtext(side=3, line = 0.3 , at = 637080, adj =0, cex = 1.2, 
      "Estimated multiplicative dispersion \nfrom the data")

legend("bottomright", fill = c("darkorange", "dodgerblue") ,
       legend= c("Binomial", "Quasi-binomial"), bty ='n')

dev.off()

