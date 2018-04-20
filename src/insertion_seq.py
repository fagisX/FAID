from itertools import product
from collections import defaultdict
import pandas as pd
from time import time

def reverse_complement(seq):
    # return seq
    new_seq = ''
    map = {'A':'T', 'T':"A", 'G':'C', 'C':'G'}
    for s in seq:
        new_seq = map[s] + new_seq
    return new_seq

if __name__ == "__main__":
    if False:
        f = open('insertion_seq.txt', 'r')
        info = f.readlines()
        f.close()

        f = open('insertion_seq.txt', 'w')
        for i in range(len(info)):
            f.write('>'+str(i)+'\n')
            f.write(info[i])
        f.close()


    if False:
        f = open('insertion_seq.txt', 'r')
        info = f.readlines()
        f.close()
        seq = []
        for line in info:
            line = line.strip()
            if line.startswith('>'):
                continue
            else:
                seq.append(line)
        motifs =[''.join(a) for a in product('ATGC', repeat=7)] + [''.join(a) for a in product('ATGC', repeat=6)] + [''.join(a) for a in product('ATGC', repeat=5)] + [''.join(a) for a in product('ATGC', repeat=4)] + [''.join(a) for a in product('ATGC', repeat=3)] + [''.join(b) for b in product('ATGC', repeat=2)] + [''.join(c) for c in product('ATGC', repeat=1)]
        # [''.join(a) for a in product('ATGC', repeat=7)] + [''.join(a) for a in product('ATGC', repeat=6)] + [''.join(a) for a in product('ATGC', repeat=5)] + [''.join(a) for a in product('ATGC', repeat=4)] +
        n = 10
        motifs = motifs[n*2000:(n+1)*2000]
        # if 'CATCATC' in motifs:
        #
        #     lala
        # else:
        #     print 'wrong'
        #     lala
        out = 'motif_count'+str(n)+'.xls'
        reverse_comp_motif = {}
        for m in motifs:
            reverse_comp_motif[m] = reverse_complement(m)

        results = {}
        results['insertion'] = defaultdict(int)
        results['genome'] = defaultdict(int)
        insertion_length = 0
        print len(seq)
        for s in seq:
            for m in motifs:
                cur_count = s.count(m)
                results['insertion'][m] += cur_count
                results['insertion'][reverse_comp_motif[m]] += cur_count
            insertion_length += len(s)
        print 'first part done!'
    # get yeast genome at gc composition
        genome_size = 0
        seq = ['chrmt.fsa']
        for i in range(1, 17):
            if i < 10:
                seq.append('chr0'+str(i)+'.fsa')
            else:
                seq.append('chr'+str(i)+'.fsa')
        for f in seq:
            f = open('ref_data/'+f, 'r')
            info = f.readlines()[1:]
            info = ''.join([line.strip() for line in info])
            for m in motifs:
                cur_count = info.count(m)
                results['genome'][m] += cur_count
                results['genome'][reverse_comp_motif[m]] += cur_count
            genome_size += len(info)
        df = pd.DataFrame(results)
        print df
        df['insertion_percentage'] = df['insertion']/(insertion_length*1.)
        df['genome_percentage'] = df['genome']/(genome_size*1.)
        df.to_csv(out, sep='\t')

    if False:
        df = None
        results = {}
        results['insertion'] = defaultdict(int)
        results['genome'] = defaultdict(int)
        for i in range(11):
            cur_df = pd.read_csv('motif_count'+str(i)+'.xls', sep='\t', index_col=0)
            for i in cur_df.index:
                results['insertion'][i] += cur_df.ix[i, 'insertion']
                results['genome'][i] += cur_df.ix[i, 'genome']
        df = pd.DataFrame(results)
        df.to_csv('motif_count'+'.xls', sep='\t')



