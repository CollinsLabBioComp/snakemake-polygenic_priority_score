#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import random

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description=""
    )

    parser.add_argument(
        '-fr', '--feature_rows',
        action='store',
        dest='fr',
        required=True,
        help='''
        Feature rows
        '''
    )

    parser.add_argument(
        '-fm', '--feature_matrices',
        action='store',
        dest='fm',
        required=True,
        help='''
        Comma-separated list of *.npy Feature matrices.
        '''
    )

    parser.add_argument(
        '-s', '--seed',
        action='store',
        dest='seed',
        required=True,
        help='''
        Seed
        '''
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='shuffled.npy',
        help='Output'
    )
    options = parser.parse_args()

    with open(options.fr, 'r') as f:
        rows = [ x.strip('\n') for x in f ]

    random.seed(options.seed)
    
    # get index list then shuffle
    ix = list(range(0, len(rows)))
    random.shuffle(ix)

    for each in options.fm.split(','):
        mtx = np.load(each)
        mtx = mtx[ix, ].copy()
        np.save(each, mtx)

    return

if __name__ == '__main__':
    main()
