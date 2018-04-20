setwd('/Users/boxia/PycharmProjects/Yeast_Greg/insertion_genome_location_result/')

pos_df = read.csv('/Users/boxia/PycharmProjects/Yeast_Greg/insertion_genome_location_result/insertions.xls', sep='\t')
neg_df = read.csv('/Users/boxia/PycharmProjects/Yeast_Greg/insertion_genome_location_result/controls.xls', sep='\t')

pos_df = log10(pos_df)
neg_df = log10(neg_df)

pdf('ARS_insertions_distance_box_plot.pdf')
par(mfrow = c(2, 2))
boxplot(pos_df[,'All'],neg_df[,'All'], outline=FALSE,names=c('',''),col=c("red","cyan"), ylim=c(0,3.2))
boxplot(pos_df[,'Confirmed'],neg_df[,'Confirmed'], outline=FALSE,names=c('',''),col=c("red","cyan"), ylim=c(0,3.2))
boxplot(pos_df[,'Dubious'],neg_df[,'Dubious'], outline=FALSE,names=c('',''),col=c("red","cyan"), ylim=c(0,3.2))
boxplot(pos_df[,'Likely'],neg_df[,'Likely'], outline=FALSE,names=c('',''),col=c("red","cyan"), ylim=c(0,3.2))
dev.off()

pdf('other_insertions_distance_box_plot.pdf')
par(mfrow = c(2, 2))
boxplot(pos_df[,'rH2A'],neg_df[,'rH2A'], outline=FALSE,names=c('',''),col=c("red","cyan"), ylim=c(0,3.2))
boxplot(pos_df[,'Rloops'],neg_df[,'Rloops'], outline=FALSE,names=c('',''),col=c("red","cyan"), ylim=c(0,3.2))
dev.off()

wilcox.test(pos_df$All, neg_df$All, alternative="less")$p.value
wilcox.test(pos_df$Confirmed, neg_df$Confirmed, alternative="less")$p.value
wilcox.test(pos_df$Dubious, neg_df$Dubious, alternative="less")$p.value
wilcox.test(pos_df$Likely, neg_df$Likely, alternative="less")$p.value
wilcox.test(pos_df$rH2A, neg_df$rH2A, alternative="less")$p.value
wilcox.test(pos_df$Rloops, neg_df$Rloops, alternative="less")$p.value
