f = open('/Users/guangyu/Work/Yeast/Data/Hi-C_rep2.txt')
nf =open('/Users/guangyu/Work/Yeast/Data/Hi-C_rep2.bed','w')

f.readline()


for line in f:
	replacements = ('|', ':', '-', '\t')
	for r in replacements:
		line = line.replace(r, ' ')
	l = line.split()
	
	nf.write(l[3]+'\t'+l[4]+'\t'+l[5]+'\t'+l[6]+'\n')

f.close()
nf.close()

f = open('/Users/guangyu/Work/Yeast/Data/donor_posi_MATa2.txt')
nf =open('/Users/guangyu/Work/Yeast/Data/donor_posi_MATa2.bed','w')

f.readline()

for line in f:
	l = line.strip()
	nf.write('chr'+l + '\t' +'1\n')

f.close()
nf.close()
