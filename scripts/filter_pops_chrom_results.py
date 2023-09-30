#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import pandas as pd


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description=""
    )

    parser.add_argument(
        '-pr', '--pops_results',
        action='store',
        dest='pr',
        required=True,
        help='''
        PoPs reults. Should contain the following column: `ENSGID`
        '''
    )

    parser.add_argument(
        '-gl', '--gene_locations',
        action='store',
        dest='gl',
        required=True,
        help='''
        Gene locations. Should contain the following columns: `ENSGID` and `CHR`
        '''
    )

    parser.add_argument(
        '-chr', '--chromosome',
        action='store',
        dest='chr',
        required=True,
        help='''
        Chromosome
        '''
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='merged_df.tsv.gz',
        help='Output file with merged dataframes.'
    )
    options = parser.parse_args()

    pops = pd.read_csv(
        options.pr,
        sep='\t',
        header=0,
        dtype=str 
    )
    genes = pd.read_csv(
        options.gl,
        sep='\t',
        header=0,
        dtype=str 
    )

    # subset
    genes = genes[genes.CHR == options.chr]
    pops = pops[[ x in genes.ENSGID.tolist() for x in pops.ENSGID]]
    pops.to_csv(
        options.of,
        index=False,
        header=True,
        sep='\t'
    )
    return


if __name__ == '__main__':
    main()
