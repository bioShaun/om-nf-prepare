import pandas as pd
from pandas import DataFrame
import fire
import os


def split_bed(input_file, outdir, chr_list=None, fai=False):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    bed_df = pd.read_table(input_file, header=None, index_col=0)
    bed_df.index = bed_df.index.map(unicode)
    if chr_list is not None:
        chr_df = pd.read_table(chr_list, header=None, index_col=0)
        chr_df.index = chr_df.index.map(unicode)
        bed_df = DataFrame(bed_df.loc[chr_df.index])
    if fai:
        bed_df = DataFrame(bed_df.loc[:, 1])
        bed_df.columns = ['end']
        bed_df.loc[:, 'start'] = 0
        bed_df = bed_df.loc[:, ['start', 'end']]
    for each_idx in bed_df.index.unique():
        each_split_df = bed_df.loc[each_idx]
        if fai:
            each_split_df = DataFrame(each_split_df).T
        each_split_file = os.path.join(outdir, '{}.bed'.format(each_idx))
        each_split_df.to_csv(each_split_file, sep='\t', header=False)


if __name__ == '__main__':
    fire.Fire(split_bed)
