import os, pandas as pd

if False:
    peaks = ['Galactose.bed', 'Dextrose.bed']
    insertions = ['insertion_total.bed'] + ['insertion_total_control'+str(i)+'.bed' for i in range(0, 10000, 100)]

    for p in peaks:
        for i in insertions:
            os.system('/Users/boxia/tools/bedtools2/bin/bedtools intersect -wa -a '+p +' -b '+i+' > '+p[:-4]+'_'+i)

if True:
    peaks = ['Galactose', 'Dextrose']
    insertions = ['insertion_total.bed'] + ['insertion_total_control'+str(i)+'.bed' for i in range(0, 10000, 100)]

    for p in peaks:
        results = {}
        results['insertions'] = []
        results['controls'] = []
        for i in insertions:
            f = p + '_' + i
            try:
                df = pd.read_csv(f, sep='\t', header=None)
            except:
                results['controls'].append(0)
                continue
            if i == 'insertion_total.bed':
                results['insertions'].append(df.iloc[:, 4].sum())
            else:
                if df.shape[0] == 0:
                    results['controls'].append(0)
                else:
                    results['controls'].append(df.iloc[:, 4].sum())
        print len(results['controls']), len(results['insertions'])
        results['insertions'] = results['insertions']+ ['']*(len(results['controls'])-len(results['insertions']))
        final_df = pd.DataFrame(results)
        final_df.to_csv(p+'_nucleosome_occupacy.xls', sep='\t', index=None)

if False:
    peaks = ['Galactose', 'Dextrose']
    insertions = ['insertion_total.bed'] + ['insertion_total_control'+str(i)+'.bed' for i in range(0, 10000, 100)]
    for p in peaks:
        results = {}
        results['insertions'] = []
        results['controls'] = []
        for i in insertions:
            f = p + '_' + i
            try:
                df = pd.read_csv(f, sep='\t', header=None)
            except:
                results['controls'].append(0)
                continue
            if i == 'insertion_total.bed':
                results['insertions']+= list(df.iloc[:, 4])
            else:
                if df.shape[0] == 0:
                    results['controls'].append(0)
                else:
                    results['controls'] += list(df.iloc[:, 4])
        print len(results['controls']), len(results['insertions'])
        results['insertions'] = results['insertions']+ ['']*(len(results['controls'])-len(results['insertions']))
        final_df = pd.DataFrame(results)
        final_df.to_csv(p+'_nucleosome_occupacy_detail.xls', sep='\t', index=None)

