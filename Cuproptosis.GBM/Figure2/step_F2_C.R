# CNV

# 设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)

# 读入数据
rt=read.table("Input_data/CNV/cnvMatrix.txt", header=T, sep="\t", check.names=F, row.names=1)   
GAIN=rowSums(rt> 0)      
LOSS=rowSums(rt< 0)      
GAIN=GAIN/ncol(rt)*100    
LOSS=LOSS/ncol(rt)*100  
data=cbind(GAIN, LOSS)
data=data[order(data[,"GAIN"],decreasing = T),]

data.max = apply(data, 1, max)
pdf(file="Figure2/out_polt/F2_C.pdf", width=9, height=6)
cex=1.3
par(cex.lab=cex, cex.axis=cex, font.axis=2, las=1, xpd=T)
bar=barplot(data.max, col="grey80", border=NA,
            xlab="", ylab="CNV.frequency(%)", space=1.5,
            xaxt="n", ylim=c(0,1.2*max(data.max)))
points(bar,data[,"GAIN"], pch=20, col=2, cex=3)
points(bar,data[,"LOSS"], pch=20, col=3, cex=3)
legend("top", legend=c('GAIN','LOSS'), col=2:3, pch=20, bty="n", cex=2, ncol=2)
par(srt=45)
text(bar, par('usr')[3]-0.2, rownames(data), adj=1)
dev.off()

