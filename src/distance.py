import pandas as pd, numpy as np

if __name__ == "__main__":
    sizes = pd.read_excel('./ref_data/sacCer_chrom_size.xlsx', index_col=0)

    results = {}

    ori_df = pd.read_csv('ALL_OriDB.xls', sep='\t')

    for chrom in sizes.index:
        cur_ori_df = ori_df[ori_df.chr==chrom]
        cur_result = []
        for i in range(cur_ori_df.shape[0]):
            if len(cur_result) == 0:
                cur_result.append(cur_ori_df.iloc[i, 1] - 1)
                if cur_ori_df.iloc[i, 1] - 1 < 0:
                    print 1
            elif i == cur_ori_df.shape[0] - 1:
                cur_result.append(sizes.ix[chrom, 'Length (bp)'] - cur_ori_df.iloc[i, 2])
                if sizes.ix[chrom, 'Length (bp)'] - cur_ori_df.iloc[i, 2] < 0:
                    print 2
            else:
                cur_result.append(cur_ori_df.iloc[i+1, 1] - cur_ori_df.iloc[i, 2])
                if cur_ori_df.iloc[i+1, 1] - cur_ori_df.iloc[i, 2] < 0:
                    print 3
        results[chrom] = cur_result
    final = []

    for k, v in results.items():
        final += v

    for f in final:
        print f
