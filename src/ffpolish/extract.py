"""Extract"""


import os
import pandas as pd
from . import deepsvr_utils as dp

import numpy as np
from loguru import logger


def convert_coordinates(row: np.array):
    """Converts 0-based coordinates to 1-based"""
    row['start'] += 1
    return row


def extract(ref, bam, outdir, prefix, labels):
    """Extract"""

    if not prefix:
        prefix = os.path.basename(bam.split('.')[0])

    # Generate matrix of true variants
    if labels or not os.path.exists(os.path.join(outdir, 'true_vars.pkl')):
        logger.info('Converting ground truth variants to 1-based coordinates')
        true_vars = pd.read_csv(labels, sep='\t', index_col=None, header=None, dtype={0: str})
        true_vars.columns = ['chr', 'start', 'end', 'ref', 'alt']
        # true_vars = true_vars.progress_apply(convert_one_based, axis=1)
        true_vars.to_pickle(os.path.join(outdir, 'true_vars.pkl'))

    logger.info('Preparing data')
    prep_data = dp.PrepareData(prefix, bam, bed_file_path, ref, outdir)
    
    df = prep_data.training_data

    true_vars = pd.read_pickle(os.path.join(DATA, 'train_data', 'true_vars.pkl'))
    df['real'] = 0

    sample = df.index[0].split('~')[0]
    true_vars_set = set(df.index.str.replace(sample + '~', ''))

    for index,row in true_vars.iterrows():
        progress(index, true_vars.shape[0])
        var = "{0}:{1}-{2}{3}>{4}".format(row.chr, row.start, row.end, row.ref, row.alt)
        if var in true_vars_set:
            df.loc[sample + '~' + var, 'real'] = 1

    # TODO: Save dataframe in hdf5 format?
    df.to_hdf(outdir / "train.hdf5")
    # vaex_df = vaex.from_pandas(df)
    # vaex_df.export(os.path.join(outdir, 'train.hdf5'))
