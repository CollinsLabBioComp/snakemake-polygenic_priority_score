#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="Calculate genetic correlation matrix for a given vcf."
    )

    parser.add_argument(
        '-in', '--in_list',
        action='store',
        dest='inlst',
        default='',
        help='''
        Comma-separated list of files to merge
        '''
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='merged_file.tsv.gz',
        help='''
        Out file
        '''
    )
    options = parser.parse_args()

    inlist = [ pd.read_csv(x, sep='\t', header=0) for x in options.inlst.split(',') ]
    final_df = pd.concat(inlist, axis=0, ignore_index=True)
  
    # Save data
    final_df.to_csv(
        options.of,
        index=False,
        header=True,
        sep='\t',
        compression='gzip'
    )
    return


if __name__ == '__main__':
    main()
