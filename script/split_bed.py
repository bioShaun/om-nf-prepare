import pandas as pd
import fire
import os


def split_bed(bed_file, outdir, chr_list=None):
    bed_df = pd.read_table(bed_file, header=None, index_col=0)
    if chr_list is not None:
        chr_df = pd.read_table(chr_list, header=None, index_col=0)
        bed_df = bed_df.loc[chr_df.index]
    for each_idx in bed_df.index:
        each_split_df = bed_df.loc[each_idx]
        each_split_file = os.path.join(outdir, '{}.bed'.format(each_idx))
        each_split_df.to_csv(each_split_file, sep='\t', header=False,
                             index=False)


if __name__ == '__main__':
    fire.Fire(split_bed)