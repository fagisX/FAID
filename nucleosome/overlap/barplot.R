setwd('/Users/boxia/PycharmProjects/Yeast_Greg/')

df = read.csv('/Users/boxia/PycharmProjects/Yeast_Greg/insertion_distance_to_TSSs.xls', sep='\t')

pdf('TSS_insertions_distance_box_plot.pdf')
par(mfrow = c(2, 2))
boxplot(df[,'insertions'],df[,'controls'], outline=FALSE,names=c('',''),col=c("red","cyan"))
dev.off()


wilcox.test(df$insertions, df$controls, alternative="less")$p.value

setwd('/Users/boxia/PycharmProjects/Yeast_Greg/nucleosome/overlap/')

pdf('nucleosome_occupacy.pdf')
par(mfrow = c(2, 2))
df = read.csv('/Users/boxia/PycharmProjects/Yeast_Greg/nucleosome/overlap/Dextrose_nucleosome_occupacy_detail.xls', sep='\t')
boxplot(df[,'insertions'],df[,'controls'], outline=FALSE,names=c('',''),col=c("red","cyan"))
wilcox.test(df$insertions, df$controls, alternative="less")$p.value
df = read.csv('/Users/boxia/PycharmProjects/Yeast_Greg/nucleosome/overlap/Dextrose_nucleosome_occupacy.xls', sep='\t')
boxplot(df[,'insertions'],df[,'controls'], outline=FALSE,names=c('',''),col=c("red","cyan"))
wilcox.test(df$insertions, df$controls, alternative="less")$p.value
df = read.csv('/Users/boxia/PycharmProjects/Yeast_Greg/nucleosome/overlap/Galactose_nucleosome_occupacy_detail.xls', sep='\t')
boxplot(df[,'insertions'],df[,'controls'], outline=FALSE,names=c('',''),col=c("red","cyan"))
wilcox.test(df$insertions, df$controls, alternative="less")$p.value
df = read.csv('/Users/boxia/PycharmProjects/Yeast_Greg/nucleosome/overlap/Galactose_nucleosome_occupacy.xls', sep='\t')
boxplot(df[,'insertions'],df[,'controls'], outline=FALSE,names=c('',''),col=c("red","cyan"))
wilcox.test(df$insertions, df$controls, alternative="less")$p.value
dev.off()


