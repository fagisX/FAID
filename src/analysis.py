"""
Main port for analyzing the insertion DNAs and yeast genome features
"""

import pandas as pd, os, numpy as np
from collections import defaultdict
import random

def int_to_roman(input):
   """
   Convert an integer to Roman numerals.

   Examples:
   >>> int_to_roman(0)
   Traceback (most recent call last):
   ValueError: Argument must be between 1 and 3999

   >>> int_to_roman(-1)
   Traceback (most recent call last):
   ValueError: Argument must be between 1 and 3999

   >>> int_to_roman(1.5)
   Traceback (most recent call last):
   TypeError: expected integer, got <type 'float'>

   >>> for i in range(1, 21): print int_to_roman(i)
   ...
   I
   II
   III
   IV
   V
   VI
   VII
   VIII
   IX
   X
   XI
   XII
   XIII
   XIV
   XV
   XVI
   XVII
   XVIII
   XIX
   XX
   >>> print int_to_roman(2000)
   MM
   >>> print int_to_roman(1999)
   MCMXCIX
   """
   if type(input) != type(1):
      raise TypeError, "expected integer, got %s" % type(input)
   if not 0 < input < 4000:
      raise ValueError, "Argument must be between 1 and 3999"
   ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
   nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
   result = ""
   for i in range(len(ints)):
      count = int(input / ints[i])
      result += nums[i] * count
      input -= ints[i] * count
   return result

def roman_to_int(input):
    """
    Convert a roman numeral to an integer.

    >>> r = range(1, 4000)
    >>> nums = [int_to_roman(i) for i in r]
    >>> ints = [roman_to_int(n) for n in nums]
    >>> print r == ints
    1

    >>> roman_to_int('VVVIV')
    Traceback (most recent call last):
     ...
    ValueError: input is not a valid roman numeral: VVVIV
    >>> roman_to_int(1)
    Traceback (most recent call last):
     ...
    TypeError: expected string, got <type 'int'>
    >>> roman_to_int('a')
    Traceback (most recent call last):
     ...
    ValueError: input is not a valid roman numeral: A
    >>> roman_to_int('IL')
    Traceback (most recent call last):
     ...
    ValueError: input is not a valid roman numeral: IL
    """
    if type(input) != type(""):
        raise TypeError, "expected string, got %s" % type(input)
    input = input.upper()
    nums = ['M', 'D', 'C', 'L', 'X', 'V', 'I']
    ints = [1000, 500, 100, 50, 10, 5, 1]
    places = []
    for c in input:
        if not c in nums:
            raise ValueError, "input is not a valid roman numeral: %s" % input
    for i in range(len(input)):
        c = input[i]
        value = ints[nums.index(c)]
        # If the next place holds a larger number, this value is negative.
        try:
            nextvalue = ints[nums.index(input[i + 1])]
            if nextvalue > value:
                value *= -1
        except IndexError:
            # there is no next place.
            pass
        places.append(value)
    sum = 0
    for n in places: sum += n
    # Easiest test for validity...
    if int_to_roman(sum) == input:
        return sum
    else:
        raise ValueError, 'input is not a valid roman numeral: %s' % input


def ori_to_index(df, bin=100, distance=0):
    results = {}
    f = open(df, 'r')
    for line in f.readlines()[1:]:
        # print line
        line = line.split()
        if len(line) <=1:
            continue
        t = (int(line[0]),
             int(line[1])-distance,
             int(line[2])+distance)
        if t[0] not in results:
            results[t[0]] = {}

        for i in range(t[1]/bin, t[2]/bin+1):
            if i in results[t[0]]:
                results[t[0]][i].add(t)
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(t)
    f.close()
    return results

def feature_to_index(df, bin=100):
    results = {}
    f = open(df, 'r')
    for line in f.readlines()[1:]:
        line = line.split()
        if len(line) <=1:
            continue
        t = (int(line[0]),
             int(line[1]),
             int(line[2]),
             )
        if t[0] not in results:
            results[t[0]] = {}

        for i in range(t[1]/bin, t[2]/bin+1):
            if i in results[t[0]]:
                results[t[0]][i].add(t)
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(t)
    f.close()
    return results

def gtf_to_index(df, bin=100):
    results = {}
    f = open(df, 'r')
    for line in f.readlines()[1:]:
        line = line.split()
        if len(line) <= 1:
            continue
        t = [line[2],
             line[3],
             int(line[4]),
             int(line[5]),
             ]
        t[0] = roman_to_int(t[0].replace('chr',''))
        if t[0] not in results:
            results[t[0]] = {}
        for i in range(t[2] / bin, t[3] / bin + 1):
            if i in results[t[0]]:
                results[t[0]][i].add(tuple(t))
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(tuple(t))
    f.close()
    return results

def gtf_to_index_genename(df, bin=100):
    results = {}
    f = open(df, 'r')
    for line in f.readlines()[1:]:
        line = line.split()
        if len(line) <= 1:
            continue
        t = [line[2],
             line[3],
             int(line[4]),
             int(line[5]),
             line[1],
             ]
        t[0] = roman_to_int(t[0].replace('chr',''))
        if t[0] not in results:
            results[t[0]] = {}
        for i in range(t[2] / bin, t[3] / bin + 1):
            if i in results[t[0]]:
                results[t[0]][i].add(tuple(t))
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(tuple(t))
    f.close()
    return results

def random_location(df_path, chromosome_size="./ref_data/sacCer_chrom_size.xlsx", random_state=0, rDNA=False):
    df = pd.read_csv(df_path, sep='\t')
    chr_sizes = pd.read_excel(chromosome_size, index_col=0)

    results = []
    random.seed(random_state)

    hotspots_map = defaultdict(set)
    rDNA_hotspots_map = defaultdict(set)

    for i in range(df.shape[0]):
        try:
            if rDNA:
                chr, start, end, hotspots, region, rDNA_label = list(df.iloc[i, :])
            else:
                chr, start, end, hotspots, region = list(df.iloc[i, :])
            chr = int(str(chr).replace('chr', ''))
            # print region, not pd.isnull(region)
            if not pd.isnull(region):
                region = region.strip()
                region = region[1:-1]
                region = [int(x) for x in region.split(',')]
                region = tuple(region)
            if hotspots and rDNA and rDNA_label=='rDNA':
                rDNA_hotspots_map[region].add((chr, start, end))
                continue
            elif hotspots:
                hotspots_map[region].add((chr, start, end))
                continue
            if chr not in chr_sizes.index:
                continue
            chr_size = chr_sizes.ix[chr, 'Length (bp)'] - (end - start)

            new_start = int(random.uniform(1, chr_size))
            new_end = new_start + (end - start)
            if rDNA and rDNA_label == 'rDNA':
                results.append([chr, new_start, new_end, 0, '', 'rDNA'])
            elif rDNA:
                results.append([chr, new_start, new_end, 0, '', ''])
            else:
                results.append([chr, new_start, new_end, 0, ''])
        except:
            continue

    for key, value in hotspots_map.items():
        value = list(value)
        target_chr = key[0]
        target_length = key[2] - key[1]
        target_number = len(value)
        chr_size = chr_sizes.ix[target_chr, 'Length (bp)'] - target_length
        new_region_start = int(random.uniform(1, chr_size))
        new_region_end = new_region_start + target_length

        for i in range(target_number):
            new_start = int(random.uniform(new_region_start, new_region_end))
            new_length = value[i][2] - value[i][1]
            new_end = new_start + new_length
            if rDNA:
                results.append([target_chr, new_start, new_end, 1, (target_chr, new_region_start, new_region_end), ''])
            else:
                results.append([target_chr, new_start, new_end, 1, (target_chr, new_region_start, new_region_end)])

    for key, value in rDNA_hotspots_map.items():
        value = list(value)
        target_chr = key[0]
        target_length = int(key[2]) - int(key[1])
        target_number = len(value)
        chr_size = chr_sizes.ix[target_chr, 'Length (bp)'] - target_length
        new_region_start = int(random.uniform(1, chr_size))
        new_region_end = new_region_start + target_length

        for i in range(target_number):
            new_start = int(random.uniform(new_region_start, new_region_end))
            new_length = value[i][2] - value[i][1]
            new_end = new_start + new_length
            results.append([target_chr, new_start, new_end, 1, (target_chr, new_region_start, new_region_end), 'rDNA'])
    final_df = pd.DataFrame(results)
    if rDNA:
        final_df.columns = ['chr', 'start', 'end', 'hotspot', 'hotspot_region', 'rDNA_label']
    else:
        final_df.columns = ['chr', 'start', 'end', 'hotspot', 'hotspot_region']
    final_df.to_csv('./controls/'+df_path[df_path.rfind('/')+1:-4]+'_control'+str(random_state)+'.xls', sep='\t', index=False)

if __name__ == "__main__":
    ## Q0s: separate OriDB to different catogory
    if False:
        ori_df = pd.read_csv('All_OriDB.xls', sep='\t')
        confirmed_df = ori_df[ori_df['status']=='Confirmed']
        likely_df = ori_df[ori_df['status']=='Likely']
        dubious_df = ori_df[ori_df['status']=='Dubious']
        confirmed_df.to_csv('Confirmed_OriDB.xls', sep='\t', index=False)
        likely_df.to_csv('Likely_OriDB.xls', sep='\t', index=False)
        dubious_df.to_csv('Dubious_OriDB.xls', sep='\t', index=False)

    ## Q0s: separate insertions table to get just the chromosome, start, end
    if False:
        xls = pd.ExcelFile('./ref_data/insertion-V3.0 -20180205-Bo.xlsx')
        if False:
            df_total = pd.read_excel(xls, 'total but no rDNA', index_col=0)
            df_total = df_total[['donor chromosome(s)', "5' template end", "3' template end", "insert sequnce"]]
            df_total.columns = ['chromosome', 'start', 'end', 'sequence']
            df_total.to_csv('insertion_total.xls', sep='\t', index=False)
            seq = list(df_total.sequence)
            f = open('insertion_seq.txt', 'w')
            i = 1
            for s in seq:
                f.write('>'+str(i) + '\n')
                f.write(s+'\n')
                i += 1
            f.close()
            df_total = df_total[['chromosome', 'start', 'end']]
            df_total = df_total.sort_values(['chromosome', 'start'], ascending=[True, True])
            df_total.to_csv('insertion_total.xls', sep='\t', index=False)

        df_total = pd.read_excel(xls, 'total with rDNA', index_col=0)
        df_total = df_total[['donor chromosome(s)', "5' template end", "3' template end", "insert sequnce", "rDNA"]]
        df_total.columns = ['chromosome', 'start', 'end', 'sequence', 'rDNA_label']

        seq = list(df_total.sequence)
        f = open('insertion_seq_rDNA.txt', 'w')
        i = 1
        for s in seq:
            f.write('>' + str(i) + '\n')
            f.write(s + '\n')
            i += 1
        f.close()
        df_total = df_total[['chromosome', 'start', 'end', 'rDNA_label']]
        df_total = df_total.sort_values(['chromosome', 'start'], ascending=[True, True])
        df_total.to_csv('insertion_total_rDNA.xls', sep='\t', index=False)

    ## Get the hotspot insertions
    if False:
        final = {}
        hotspots = []
        hot_chr, hot_start, hot_end = None, None, None
        df_total = pd.read_csv('insertion_total.xls', sep='\t')
        for i in range(df_total.shape[0]):
            chr, start, end = list(df_total.ix[i, :])
            if hot_chr is None:
                hot_chr, hot_start, hot_end = chr, start, end
                hotspots = [[chr, start, end]]
            elif hot_chr != chr:
                if len(hotspots) > 1:
                    print hotspots
                    for h in hotspots:
                        final[tuple(h)] = (hot_chr, hot_start, hot_end)
                hot_chr, hot_start, hot_end = chr, start, end
                hotspots = [[chr, start, end]]
            elif end - hot_end > 3000:
                if len(hotspots) > 1:
                    print hotspots
                    for h in hotspots:
                        final[tuple(h)] = (hot_chr, hot_start, hot_end)
                hot_chr, hot_start, hot_end = chr, start, end
                hotspots = [[chr, start, end]]
            elif end - hot_end <= 3000:
                if end > hot_end:
                    hot_end = end
                hotspots.append([chr, start, end])
        hotspots_column1 =[]
        hotspots_column2 = []
        for i in range(df_total.shape[0]):
            chr, start, end = list(df_total.ix[i, :])
            if (chr, start, end) in final.keys():
                hotspots_column1.append(1)
                hotspots_column2.append(final[(chr, start, end)])
            else:
                hotspots_column1.append(0)
                hotspots_column2.append('')
        df_total['hotspot'] = hotspots_column1
        df_total['hotspot_region'] = hotspots_column2
        df_total.to_csv('insertion_total.xls', sep='\t', index=False)

    if False:
        if False:
            final = {}
            hotspots = []
            hot_chr, hot_start, hot_end = None, None, None
            df_total = pd.read_csv('insertion_total_rDNA.xls', sep='\t')
            for i in range(df_total.shape[0]):
                chr, start, end, rDNA = list(df_total.ix[i, :])
                if hot_chr is None:
                    hot_chr, hot_start, hot_end = chr, start, end
                    hotspots = [[chr, start, end]]
                elif hot_chr != chr:
                    if len(hotspots) > 1:
                        print hotspots
                        for h in hotspots:
                            final[tuple(h)] = (hot_chr, hot_start, hot_end)
                    hot_chr, hot_start, hot_end = chr, start, end
                    hotspots = [[chr, start, end]]
                elif end - hot_end > 3000:
                    if len(hotspots) > 1:
                        print hotspots
                        for h in hotspots:
                            final[tuple(h)] = (hot_chr, hot_start, hot_end)
                    hot_chr, hot_start, hot_end = chr, start, end
                    hotspots = [[chr, start, end]]
                elif end - hot_end <= 3000:
                    if end > hot_end:
                        hot_end = end
                    hotspots.append([chr, start, end])
            hotspots_column1 = []
            hotspots_column2 = []
            for i in range(df_total.shape[0]):
                chr, start, end, rDNA = list(df_total.ix[i, :])
                if (chr, start, end) in final.keys():
                    hotspots_column1.append(1)
                    hotspots_column2.append(final[(chr, start, end)])
                else:
                    hotspots_column1.append(0)
                    hotspots_column2.append('')
            df_total['hotspot'] = hotspots_column1
            df_total['hotspot_region'] = hotspots_column2
            df_total.to_csv('insertion_total_rDNA.xls', sep='\t', index=False)

        if False:
            df_total = pd.read_csv('insertion_total_rDNA.xls', sep='\t')
            for i in range(df_total.shape[0]):
                chr, start, end, rDNA = df_total.ix[i, 0], df_total.ix[i, 1], df_total.ix[i, 2], df_total.ix[i, 5]
                # print rDNA
            df_total.to_csv('insertion_total_rDNA_1.xls', sep='\t', index=False)


    ## Get the random control without rDNA
    if False:
        for r in range(0, 10000, 10):
            random_location('./ref_data/insertion_total.xls', random_state=r, rDNA=False)

    ## Get the random control with rDNA
    if False:
        for r in range(0, 10000, 100):
            random_location('./ref_data/insertion_total_rDNA.xls', random_state=r, rDNA=True)

    ## Create ARSes, Rloops, rH2As, CENs, rDNA veresion
    if False:
        # features = ['All_OriDB.xls', 'Confirmed_OriDB.xls', 'Dubious_OriDB.xls', 'Likely_OriDB.xls', 'Rloops.xls', 'rH2A.xls', 'CEN_features.xls']
        # features = ['TEL_features.xls']
        features = ['Ty_features.xls']
        for feature in features:
            df = pd.read_csv('./ref_data/'+feature, sep='\t')
            print feature
            for i in range(df.shape[0]):
                if df.ix[i, 'chr'] != 12:
                    continue
            df.to_csv('./ref_data/' + feature[:-4]+'_rDNA.xls', sep='\t', index=False)

    # For TEL
    if False:
        TEL = pd.read_csv('./ref_data/TEL_features.xls', sep='\t')
        results = set()
        for i in range(TEL.shape[0]):
            chr, start, end = TEL.ix[i, :]
            if start >= end:
                results.add((chr, end, start))
            else:
                results.add((chr, start, end))
        results = list(results)
        final_df = pd.DataFrame(results)
        final_df.columns = ['chr', 'start', 'end']
        final_df.to_csv('TEL_features.xls', sep='\t', index=None)


    ## Q1: Is there enrichment of inserts within douious or likely ARSes, or proximity of ARSes?
    if False:
        insert_df_path = './ref_data/insertion_total.xls'
        # oris = ['Confirmed_OriDB.xls', 'Likely_OriDB.xls', 'Dubious_OriDB.xls',  'Rloops.xls', 'All_OriDB.txt','rH2A.txt','CEN_features.txt']
        # oris = ['TEL_features.xls']
        # oris = ['Ty_features.xls']
        oris = ['tRNA_feature.txt']
        final_result = []
        for ori_df_path in oris:
            ori_df = ori_to_index('./ref_data/' + ori_df_path, distance=0)
            # print ori_df
            result = 0
            distance = []
            insertion_df = pd.read_csv(insert_df_path, sep='\t')
            for i in range(insertion_df.shape[0]):
                chr, start, end, hotspot = int(insertion_df.iloc[i, 0]), \
                                           int(insertion_df.iloc[i, 1]), \
                                           int(insertion_df.iloc[i, 2]), \
                                           int(insertion_df.iloc[i, 3])
                # print name
                if chr not in ori_df.keys():
                    # print type(chr), ori_df.keys(), type(ori_df.keys()[3]), 'lala'
                    continue
                else:
                    add = False
                    min_distance = float('inf')
                    for j in range(start/100, end/100+1):
                        if j not in ori_df[chr].keys():
                            continue
                        else:
                            target_ori = ori_df[chr][j]
                            for t in ori_df[chr][j]:
                                if (t[1] <= start <= t[2] or t[1] <= end <= t[2] or (t[1] <= start <= t[2] and t[1] <= end <= t[2]) or (start <= t[1] <= end and start <= t[2] <= end)):
                                    if not add:
                                        result+=1
                                        add = True
                                    cur_distance = 0
                                else:
                                    cur_distance = min(abs(t[1]-start), abs(t[1]-end), abs(t[2]-start), abs(t[2]-end))
                                if min_distance > cur_distance:
                                    min_distance = cur_distance
                    if not add:
                        candidates = set()
                        for k in ori_df[chr].keys():
                            candidates = candidates.union(ori_df[chr][k])
                        for t in candidates:
                            cur_distance = min(abs(t[1] - start), abs(t[1] - end), abs(t[2] - start),
                                               abs(t[2] - end))
                            if cur_distance < min_distance:
                                min_distance = cur_distance
                    distance.append([chr, start, end, hotspot, min_distance])
            final_result.append([ori_df_path[:-4], result])
            distance_df = pd.DataFrame(distance)
            # print ori_df_path
            # print distance_df
            distance_df.columns=['chr', 'start', 'end', 'hotspot', 'distance']
            distance_df = distance_df.sort_values(by='distance')
            distance_df['rank'] = range(1, distance_df.shape[0]+1)
            distance_df.to_csv(ori_df_path[:-4]+'_'+insert_df_path[insert_df_path.rfind('/')+1:-4]+'_distance.xls', sep='\t', index=False)
        final_df = pd.DataFrame(final_result)
        final_df.to_csv('ORI_overlap_change.xls', sep='\t')

    if False:
        insert_df_path = './ref_data/insertion_total_rDNA.xls'
        # oris = ['Confirmed_OriDB_rDNA.xls', 'Likely_OriDB_rDNA.xls', 'Dubious_OriDB_rDNA.xls',  'Rloops_rDNA.xls', 'All_OriDB_rDNA.xls','rH2A_rDNA.xls','CEN_features_rDNA.xls']
        oris = ['TEL_features_rDNA.xls']
        # oris = ['Ty_features_rDNA.xls']
        final_result = []
        for ori_df_path in oris:
            ori_df = ori_to_index('./ref_data/' + ori_df_path, distance=0)
            # print ori_df
            result = 0
            distance = []
            insertion_df = pd.read_csv(insert_df_path, sep='\t')
            for i in range(insertion_df.shape[0]):
                chr, start, end, hotspot = int(insertion_df.iloc[i, 0]), \
                                           int(insertion_df.iloc[i, 1]), \
                                           int(insertion_df.iloc[i, 2]), \
                                           int(insertion_df.iloc[i, 3])
                # print name
                if chr not in ori_df.keys():
                    # print type(chr), ori_df.keys(), type(ori_df.keys()[3]), 'lala'
                    continue
                else:
                    add = False
                    min_distance = float('inf')
                    for j in range(start/100, end/100+1):
                        if j not in ori_df[chr].keys():
                            continue
                        else:
                            target_ori = ori_df[chr][j]
                            for t in ori_df[chr][j]:
                                if (t[1] <= start <= t[2] or t[1] <= end <= t[2] or (t[1] <= start <= t[2] and t[1] <= end <= t[2]) or (start <= t[1] <= end and start <= t[2] <= end)):
                                    if not add:
                                        result+=1
                                        add = True
                                    cur_distance = 0
                                else:
                                    cur_distance = min(abs(t[1]-start), abs(t[1]-end), abs(t[2]-start), abs(t[2]-end))
                                if min_distance > cur_distance:
                                    min_distance = cur_distance
                    if not add:
                        candidates = set()
                        for k in ori_df[chr].keys():
                            candidates = candidates.union(ori_df[chr][k])
                        for t in candidates:
                            cur_distance = min(abs(t[1] - start), abs(t[1] - end), abs(t[2] - start),
                                               abs(t[2] - end))
                            if cur_distance < min_distance:
                                min_distance = cur_distance
                    distance.append([chr, start, end, hotspot, min_distance])
            final_result.append([ori_df_path[:-4], result])
            distance_df = pd.DataFrame(distance)
            # print ori_df_path
            # print distance_df
            distance_df.columns=['chr', 'start', 'end', 'hotspot', 'distance']
            distance_df = distance_df.sort_values(by='distance')
            distance_df['rank'] = range(1, distance_df.shape[0]+1)
            distance_df.to_csv(ori_df_path[:-4]+'_'+insert_df_path[insert_df_path.rfind('/')+1:-4]+'_distance.xls', sep='\t', index=False)
        final_df = pd.DataFrame(final_result)
        final_df.to_csv('ORI_overlap_change_rDNA.xls', sep='\t')


    ## random insertions as control
    if False:
        insert_df_path = './controls/insertion_total_control'
        # oris = ['Confirmed_OriDB.xls', 'Rloops.xls']
        # oris = ['Confirmed_OriDB.xls', 'Likely_OriDB.xls', 'Dubious_OriDB.xls', 'Rloops.xls', 'All_OriDB.txt', 'rH2A.txt',
        #         'CEN_features.txt']
        # oris = ['TEL_features.xls', 'CEN_features.txt']
        # oris = ['Ty_features.xls']

        oris = ['tRNA_feature.txt']

        for ori_df_path in oris:
            ori_df = ori_to_index('./ref_data/' + ori_df_path, distance=0)
            # print ori_df
            distances = {}
            distances_hotspot = {}
            final_result = []
            for random_state in range(0, 10000, 10):
                result = 0
                distance = []
                insertion_df = pd.read_csv(insert_df_path+str(random_state)+'.xls', sep='\t')
                for i in range(insertion_df.shape[0]):
                    chr, start, end, hotspot = int(insertion_df.iloc[i, 0]), \
                                               int(insertion_df.iloc[i, 1]), \
                                               int(insertion_df.iloc[i, 2]), \
                                               int(insertion_df.iloc[i, 3])
                    # print name
                    if chr not in ori_df.keys():
                        # print type(chr), ori_df.keys(), type(ori_df.keys()[3]), 'lala'
                        continue
                    else:
                        add = False
                        min_distance = float('inf')
                        for j in range(start / 100, end / 100 + 1):
                            if j not in ori_df[chr].keys():
                                continue
                            else:
                                target_ori = ori_df[chr][j]
                                for t in ori_df[chr][j]:
                                    if (t[1] <= start <= t[2] or t[1] <= end <= t[2] or (
                                                t[1] <= start <= t[2] and t[1] <= end <= t[2]) or (
                                                start <= t[1] <= end and start <= t[2] <= end)):
                                        if not add:
                                            result += 1
                                            add = True
                                        cur_distance = 0
                                    else:
                                        cur_distance = min(abs(t[1] - start), abs(t[1] - end), abs(t[2] - start),
                                                           abs(t[2] - end))
                                    if min_distance > cur_distance:
                                        min_distance = cur_distance
                        if not add:
                            candidates = set()
                            for k in ori_df[chr].keys():
                                candidates = candidates.union(ori_df[chr][k])
                            for t in candidates:
                                cur_distance = min(abs(t[1] - start), abs(t[1] - end), abs(t[2] - start),
                                                   abs(t[2] - end))
                                if cur_distance < min_distance:
                                    min_distance = cur_distance
                        distance.append([chr, start, end, hotspot, min_distance])
                final_result.append(['random_state_'+str(random_state), result])
                distances['random_state_'+str(random_state)] = sorted([d[-1] for d in distance])
                distances_hotspot['random_state_'+str(random_state)] = sorted([d[-1] for d in distance if d[-2]])
            distance_df = pd.DataFrame(distances)
            distance_df['average'] = distance_df.mean(axis=1)
            distances_hotspot_df = pd.DataFrame(distances_hotspot)
            distances_hotspot_df['average'] = distances_hotspot_df.mean(axis=1)


            distance_df = distance_df.sort_values(by='average')
            distance_df['rank'] = range(1, distance_df.shape[0] + 1)
            distance_df.to_csv(ori_df_path[:-4] + '_' + insert_df_path[insert_df_path.rfind('/') + 1:] + '_distance.xls',
                               sep='\t', index=False)

            distances_hotspot_df = distances_hotspot_df.sort_values(by='average')
            distances_hotspot_df['rank'] = range(1, distances_hotspot_df.shape[0] + 1)
            distances_hotspot_df.to_csv(ori_df_path[:-4] + '_' + insert_df_path[insert_df_path.rfind('/') + 1:] + 'hotspot_distance.xls',
                               sep='\t', index=False)
            final_df = pd.DataFrame(final_result)
            final_df.to_csv(ori_df_path[:-4] + 'control_overlap.xls', sep='\t')

    if False:
        insert_df_path = './controls/insertion_total_rDNA_control'
        # oris = ['Confirmed_OriDB.xls']
        # oris = ['Confirmed_OriDB_rDNA.xls', 'Likely_OriDB_rDNA.xls', 'Dubious_OriDB_rDNA.xls', 'Rloops_rDNA.xls',
        #         'All_OriDB_rDNA.xls',
        #         'rH2A_rDNA.xls', 'CEN_features_rDNA.xls']
        oris = ['TEL_features_rDNA.xls']
        # oris = ['Ty_features_rDNA.xls']
        for ori_df_path in oris:
            ori_df = ori_to_index('./ref_data/' + ori_df_path, distance=0)
            # print ori_df
            distances = {}
            distances_hotspot = {}
            final_result = []
            for random_state in range(0, 10000, 100):
                result = 0
                distance = []
                insertion_df = pd.read_csv(insert_df_path + str(random_state) + '.xls', sep='\t')
                for i in range(insertion_df.shape[0]):
                    chr, start, end, hotspot = int(insertion_df.iloc[i, 0]), \
                                               int(insertion_df.iloc[i, 1]), \
                                               int(insertion_df.iloc[i, 2]), \
                                               int(insertion_df.iloc[i, 3])
                    # print name
                    if chr not in ori_df.keys():
                        # print type(chr), ori_df.keys(), type(ori_df.keys()[3]), 'lala'
                        continue
                    else:
                        add = False
                        min_distance = float('inf')
                        for j in range(start / 100, end / 100 + 1):
                            if j not in ori_df[chr].keys():
                                continue
                            else:
                                target_ori = ori_df[chr][j]
                                for t in ori_df[chr][j]:
                                    if (t[1] <= start <= t[2] or t[1] <= end <= t[2] or (
                                                        t[1] <= start <= t[2] and t[1] <= end <= t[2]) or (
                                                        start <= t[1] <= end and start <= t[2] <= end)):
                                        if not add:
                                            result += 1
                                            add = True
                                        cur_distance = 0
                                    else:
                                        cur_distance = min(abs(t[1] - start), abs(t[1] - end), abs(t[2] - start),
                                                           abs(t[2] - end))
                                    if min_distance > cur_distance:
                                        min_distance = cur_distance
                        if not add:
                            candidates = set()
                            for k in ori_df[chr].keys():
                                candidates = candidates.union(ori_df[chr][k])
                            for t in candidates:
                                cur_distance = min(abs(t[1] - start), abs(t[1] - end), abs(t[2] - start),
                                                   abs(t[2] - end))
                                if cur_distance < min_distance:
                                    min_distance = cur_distance
                        distance.append([chr, start, end, hotspot, min_distance])
                final_result.append(['random_state_' + str(random_state), result])
                distances['random_state_' + str(random_state)] = sorted([d[-1] for d in distance])
                distances_hotspot['random_state_' + str(random_state)] = sorted([d[-1] for d in distance if d[-2]])
            distance_df = pd.DataFrame(distances)
            distance_df['average'] = distance_df.mean(axis=1)
            distances_hotspot_df = pd.DataFrame(distances_hotspot)
            distances_hotspot_df['average'] = distances_hotspot_df.mean(axis=1)

            distance_df = distance_df.sort_values(by='average')
            distance_df['rank'] = range(1, distance_df.shape[0] + 1)
            distance_df.to_csv(
                ori_df_path[:-4] + '_' + insert_df_path[insert_df_path.rfind('/') + 1:] + '_distance_rDNA.xls',
                sep='\t', index=False)

            distances_hotspot_df = distances_hotspot_df.sort_values(by='average')
            distances_hotspot_df['rank'] = range(1, distances_hotspot_df.shape[0] + 1)
            distances_hotspot_df.to_csv(
                ori_df_path[:-4] + '_' + insert_df_path[insert_df_path.rfind('/') + 1:] + 'hotspot_distance_rDNA.xls',
                sep='\t', index=False)
            final_df = pd.DataFrame(final_result)
            final_df.to_csv(ori_df_path[:-4] + 'control_overlap_rDNA.xls', sep='\t')

    ## Q4: Is there distribution of insertions random or not, does the number of insertion corresponds to size of chromosome?
    if False:
        df_total = pd.read_csv('./ref_data/insertion_total.xls', sep='\t')

        count_df = df_total['chromosome'].value_counts().to_frame()
        print count_df
        size_df = pd.read_excel('./ref_data/sacCer_chrom_size.xlsx', index_col=0)
        count_df['length'] = [np.nan] * count_df.shape[0]
        for i in size_df.index:
            if i in count_df.index:
                count_df.ix[i, 'length'] = size_df.ix[i, 'Length (bp)']
        count_df.columns = ['number of insertions', 'length']
        count_df.index.name = 'chromosome'
        print count_df.corr(method='spearman')
        count_df.to_excel('chrsize_vs_insertions.xlsx')

        df_total = pd.read_csv('./ref_data/insertion_total_rDNA.xls', sep='\t')

        count_df = df_total['chromosome'].value_counts().to_frame()
        size_df = pd.read_excel('./ref_data/sacCer_chrom_size.xlsx', index_col=0)
        count_df['length'] = [np.nan] * count_df.shape[0]
        for i in size_df.index:
            if i in count_df.index:
                count_df.ix[i, 'length'] = size_df.ix[i, 'Length (bp)']
        count_df.columns = ['number of insertions', 'length']
        count_df.index.name = 'chromosome'
        print count_df.corr(method='spearman')
        count_df.to_excel('chrsize_vs_insertions_rDNA.xlsx')








