yeast = read.table('/Users/guangyu/Work/Yeast/Data/donor_anno.bed',header = F,sep = '\t')
H = read.table('/Users/guangyu/Work/Yeast/Data/Hi-C_rep1.bed',header = F,sep = '\t')
yeast = na.omit(yeast)
random = sample(H[,4],length(yeast[,5]))
P = data.frame(1,2)
for(i in 1:1000){
  r = sample(H[,4],length(yeast[,5]))
  random = data.frame(random,r)
  P[i,1] = wilcox.test(yeast[,5],r)$p.value
  P[i,2] = i
}

medium=P[order(P[,1])==500,2]
P[order(P[,1])==500,1]
random=random[,medium+1]

pdf('/Users/guangyu/Work/Yeast/Data/boxplot1.pdf')
boxplot(yeast[,5],random,outline = F,ylim=c(0,100))
dev.off()



yeast = read.table('/Users/guangyu/Work/Yeast/Data/donor_anno2.bed',header = F,sep = '\t')
H = read.table('/Users/guangyu/Work/Yeast/Data/Hi-C_rep1.bed',header = F,sep = '\t')
yeast = na.omit(yeast)
random = sample(H[,4],length(yeast[,5]))
P = data.frame(1,2)
for(i in 1:1000){
  r = sample(H[,4],length(yeast[,5]))
  random = data.frame(random,r)
  P[i,1] = wilcox.test(yeast[,5],r)$p.value
  P[i,2] = i
}

medium=P[order(P[,1])==500,2]
P[order(P[,1])==500,1]
random=random[,medium+1]

pdf('/Users/guangyu/Work/Yeast/Data/boxplot2.pdf')
boxplot(yeast[,5],random,outline = F,ylim=c(0,100))
dev.off()



yeast = read.table('/Users/guangyu/Work/Yeast/Data/donor_anno_MATa2_1.bed',header = F,sep = '\t')
H = read.table('/Users/guangyu/Work/Yeast/Data/Hi-C_rep1.bed',header = F,sep = '\t')
yeast = na.omit(yeast)
random = sample(H[,4],length(yeast[,5]))
P = data.frame(1,2)
for(i in 1:1000){
  r = sample(H[,4],length(yeast[,5]))
  random = data.frame(random,r)
  P[i,1] = wilcox.test(yeast[,5],r)$p.value
  P[i,2] = i
}

medium=P[order(P[,1])==500,2]
P[order(P[,1])==500,1]
random=random[,medium+1]

pdf('/Users/guangyu/Work/Yeast/Data/boxplot3.pdf')
boxplot(yeast[,5],random,outline = F,ylim=c(0,100))
dev.off()

yeast = read.table('/Users/guangyu/Work/Yeast/Data/donor_anno_MATa2_2.bed',header = F,sep = '\t')
H = read.table('/Users/guangyu/Work/Yeast/Data/Hi-C_rep1.bed',header = F,sep = '\t')
yeast = na.omit(yeast)
random = sample(H[,4],length(yeast[,5]))
P = data.frame(1,2)
for(i in 1:1000){
  r = sample(H[,4],length(yeast[,5]))
  random = data.frame(random,r)
  P[i,1] = wilcox.test(yeast[,5],r)$p.value
  P[i,2] = i
}

medium=P[order(P[,1])==500,2]
P[order(P[,1])==500,1]
random=random[,medium+1]

pdf('/Users/guangyu/Work/Yeast/Data/boxplot4.pdf')
boxplot(yeast[,5],random,outline = F,ylim=c(0,100))
dev.off()



