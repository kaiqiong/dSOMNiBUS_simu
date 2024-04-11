# use some plots to illustrate 
library(magrittr)
source("~/scratch/kaiqiong.zhao/Projects/Robust_SOMNiBUS/functions/Functions_better_approx_Quasi_error_Rand_effect_Refine_CI.R")


source("~/scratch/kaiqiong.zhao/Projects/SOMNiBUS_RE_Simu/functions/BSMethEM_April_update_scale_estp0p1_export_Q_known_error_with_phi_Mstep_quasi_correct_Fletcher.R")

# how about the tope 1 region?

# Sasha's 2 result
load("~/scratch/kaiqiong.zhao/Projects/Sasha_new/Sasha_new/Analysis1_RE/Results/Summary/ALL_combine.RData")
colnames(res_re)


head(res_re.o)



which.min(res_re.o$NoRE)

#plot(res_re.o$p.value, res_re.o$NoRE)
see = res_re.o
# The top region has phi as 2.85 and sigma00: 2.048


# focus on the ERICH region


#id = 104 # the top region identified by model without RE
id = 1  
see[id,]


load("~/scratch/kaiqiong.zhao/Projects/Sasha_new/Sasha_new/Cov_dat/RA102.RData")

#--------


#-----------------------------
# Load the Methylation data
#------------------------------
chr <-strsplit(see$seqnames[id]%>% as.character(), "chr")[[1]][2] %>% as.numeric()

int_num = strsplit(see$fileName[id], "int")[[1]][2] %>% as.numeric()

DATA_PATH <- "/project/greenwood/share/Sasha_Bernatsky/Methylation_TCBS/04-03-2019/NOVOALIGN_INTERVAL_DATA/"

load(file.path(paste0(DATA_PATH,"CHR", chr, "/", paste0("SB.TCBS.Methylation_NovoMethyl_stranded_dedup_HG19.", see$fileName[id], ".RData"))))



methMat<-datMeth[[1]][, colnames(datMeth[[1]])%in% Covs$PROJECT_CODE ]#rows are CpG sites; and columns are samples
totalMat<-datMeth[[2]][, colnames(datMeth[[2]])%in% Covs$PROJECT_CODE  ]#rows are CpG sites; and columns are samples



cpgR <- rownames(methMat)
tmp1 <-lapply(strsplit(cpgR,"-"),function(x)x[1]) %>% unlist
position <- lapply(strsplit(tmp1,":"),function(x)x[2]) %>% unlist %>% as.numeric




## Descriptive Analysis


print(all(colnames(methMat)==Covs$PROJECT_CODE))
print(all(colnames(totalMat)==Covs$PROJECT_CODE))

pos = position
propMat <- methMat/totalMat

setwd("~/dSOMNiBUS_revision/Figures")
pdf("methyprop-variation-allsamps.pdf", width = 8, height = 4)


postscript(file = "Figure 1A.eps", width = 8, height = 4, horizontal = FALSE, onefile = FALSE, paper ="special",
          pointsize = 12)
par(mfrow = c(1, 1),     # 2x2 layout
    oma = c(0, 0, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 2, 0.5), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0) # axis label at 2 rows distance, tick labels at 1 row
) 

blue_cols = c('#edf8fb','#b2e2e2','#66c2a4','#238b45')

orange_cols = c('#ffffb2','#fecc5c','#fd8d3c','#e31a1c')

plot(1, pch = 20, xlim = c(min(pos), max(pos)),
     ylab = "Methylation proportion", ylim = c(0,1), xaxt = "n", xlab = "Genomic Position", col = "dodgerblue")

axis(side = 1, at = pos[order(pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
axis(side=1, at = round(seq((min(pos)), (max(pos)),length.out = 10)) , tck= -0.02)
for( i in 1:sum(Covs$ACPA4==0)){
  points(position,propMat[,which(Covs$ACPA4==0)[i]], pch = 19, col = blue_cols[2] )
}


for( i in 1:sum(Covs$ACPA4==1)){
  points(position,propMat[,which(Covs$ACPA4==1)[i]], pch = 20, col = orange_cols[2] )
}

#abline(v= 637562, lty = 4)
#abline(v= 637734, lty = 4)
legend("bottomleft", c("ACPA Negative", "ACPA Positive"),
       fill = c(blue_cols[2], orange_cols[2]), bty="n", cex = 0.8)

dev.off()

