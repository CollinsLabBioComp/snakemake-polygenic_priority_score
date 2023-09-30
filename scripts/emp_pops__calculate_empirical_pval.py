#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd

def calculate_empirical_pvalue(ix, obs_df, null_df, basis, score_col, comparison):
    obs_val = obs_df.loc[ix, score_col]

    if basis == 'grouped':
        null_scores = null_df.loc[ix, score_col].values
    elif basis == 'pooled':
        null_scores = null_df[score_col].values

    if comparison == '<':
        null_scores_bool = null_scores < obs_val
    elif comparison == '>':
        null_scores_bool = null_scores > obs_val
    elif comparison == '>=':
        null_scores_bool = null_scores >= obs_val
    elif comparison == '<=':
        null_scores_bool = null_scores <= obs_val
    else:
        raise ValueError("Unknown comparison operator")
    
    # n should be 18,000 genes x 1,000 perms = 18,000,000
    n = len(null_scores_bool)
    # r should be how many times that vector is greater than observed score across ALL genes  
    r = np.count_nonzero(null_scores_bool)
    
    pval = (r + 1) / (n + 1)
    return(pval)

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description=""
    )

    parser.add_argument(
        '-obs', '--observed',
        action='store',
        dest='obs',
        required=True,
        help='''
        Observed PoPS
        '''
    )

    parser.add_argument(
        '-null', '--null',
        action='store',
        dest='null',
        required=True,
        nargs='+',
        help='''
        Null PoPS
        '''
    )

    parser.add_argument(
        '-ic', '--index_col',
        action='store',
        dest='ic',
        required=True,
        help='''
        Index column
        '''
    )

    parser.add_argument(
        '-sc', '--score_col',
        action='store',
        dest='sc',
        required=True,
        help='''
        Score column
        '''
    )

    parser.add_argument(
        '-co', '--comparison_operator',
        action='store',
        dest='co',
        default='<',
        choices=['<', '>', '>=', '<='],
        help='''
        Comparison operator. Compares null_scores `operator` observed_score.
        If evals to true, it increases the p-value (becomes less significant)
        '''
    )

    parser.add_argument(
        '-b', '--basis',
        action='store',
        dest='basis',
        default='pooled',
        choices=['pooled', 'grouped'],
        help='''
        Select the basis for comparison:
        'pooled' uses the entire population
        'grouped' uses subsetted data on the condition
        '''
    )

    parser.add_argument(
        '-uav', '--use_absolute_value',
        action='store_true',
        dest='use_abs',
        default=False,
        help='''
        Use absolute value in empirical p-value calculation
        '''
    )

    parser.add_argument(
        '-smn', '--save_merged_null',
        action='store',
        dest='smn',
        default='merged_df.tsv.gz',
        help='Output file with merged dataframes.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='merged_df.tsv',
        help='Output file with merged dataframes.'
    )
    options = parser.parse_args()

    ix_col = options.ic
    score_col = options.sc

    observed = pd.read_csv(options.obs, sep='\t', header = 0, index_col = ix_col)
    if len(observed) != len(observed.index.unique()):
        raise Exception('Non-unique elements in the observed dataframe.')

    print("Reading null results...")
    null_results = []
    for f in options.null:
        try:
            null_results.append(pd.read_csv(f, sep='\t', header=0))
        except Exception as e:
            print('There was an error. Skipping', f)
            continue

    # Concatenate all dataframes
    print('Concatenating dataframes...')
    null_results = pd.concat(null_results, axis=0)
    if options.smn is not None:
        print("Path passed to save_merged_null, saving...")
        null_results.to_csv(options.smn, sep='\t', header=True)

    # Absolute value if requested
    if options.use_abs:
        observed[score_col] = np.abs(observed[score_col])
        null_results[score_col] = np.abs(null_results[score_col])

    # Calculate empirical p-value
    print("Calculating empirical p-values...")
    null_results_ix = null_results.set_index(ix_col)
    observed['pvalue'] = [
        calculate_empirical_pvalue(
            ix = x,
            obs_df = observed,
            null_df = null_results_ix,
            basis = options.basis,
            score_col = score_col,
            comparison = options.co
        ) for x in observed.index
    ]

    # Write updated observed
    observed.to_csv(options.of, sep='\t', index=True, header=True)
    return

if __name__ == '__main__':
    main()
