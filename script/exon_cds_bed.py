#!/usr/bin/env python

import gtfparse
import fire
from pathlib import PurePath
import envoy


def feature_bed(gtf_file, features=['cds', 'exon']):
    gtf_df = gtfparse.read_gtf(gtf_file)
    gtf_df.loc[:, 'feature_lower'] = gtf_df.feature.map(str.lower)
    gtf_df.loc[:, 'bed_start'] = gtf_df.start - 1
    gtf_file = PurePath(gtf_file)
    for feature in features:
        feature_gtf_df = gtf_df[gtf_df.feature_lower == feature]
        feature_bed_file = gtf_file.with_suffix('.{}.bed'.format(feature))
        feature_gtf_df.to_csv(str(feature_bed_file), index=False,
                              columns=['seqname', 'bed_start', 'end'],
                              sep='\t', header=False)
        sort_bed_file = gtf_file.with_suffix('.{}.sort.bed'.format(feature))
        sort_bed = 'sort -k1,1 -k2,2n {bed}'.format(
            bed=feature_bed_file
        )
        sort_bed_response = envoy.run(sort_bed)
        with open(str(sort_bed_file), 'w') as bed_inf:
            bed_inf.write(sort_bed_response.std_out)
        merge_bed_file = gtf_file.with_suffix('.{}.merged.bed'.format(feature))
        merge_bed = 'bedtools merge -i {sorted_bed}'.format(
            sorted_bed=sort_bed_file
        )
        merge_bed_response = envoy.run(merge_bed)
        with open(str(merge_bed_file), 'w') as m_bed_inf:
            m_bed_inf.write(merge_bed_response.std_out)


if __name__ == '__main__':
    fire.Fire(feature_bed)
