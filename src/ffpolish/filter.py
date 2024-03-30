"""Filter"""

from importlib import resources
from pathlib import Path
from typing import Union

import os
import pandas as pd
from . import models
from . import deepsvr_utils as dp
import joblib
import gzip

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import GridSearchCV

from loguru import logger


def vcf2bed(vcf: Union[str, Path]):
    """Convert vcf file to bed"""

    if isinstance(vcf, str):
        vcf = Path(vcf)

    with vcf.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            (chr, start, name, span) = line.strip().split("\t")
            start -= start
            span = len(span)
            end = start + span
            logger.info(f"{chr=} {start=} {name=} {span=} {end=}")
            print("\t".join([chr, start, end, name]))

def filter(ref: Union[str, Path], vcf: Union[str, Path], bam: Union[str, Path], outdir: Union[str, Path], prefix: Union[str, Path], retrain: bool = False, grid_search: bool = True, cores: int = 1, seed: int = 24159):
    """"""

    logger.info("Running FFPolish prediction")

    if not prefix:
        prefix = bam.stem
    
    if retrain:
        orig_train_df = pd.read_hdf("orig_train.hdf5")
        # OLD CODE
        # orig_train_df = vaex.open(os.path.join(BASE, 'orig_train.hdf5'))
        # train_features = orig_train_df.get_column_names(regex='tumor')
        # new_train_df = vaex.open(retrain)

        orig_train_df = orig_train_df.to_pandas_df(train_features + ['real'])
        new_train_df = new_train_df.to_pandas_df(train_features + ['real'])

        clf = LogisticRegression(penalty='l2', random_state=seed, solver='saga', max_iter=10000, 
                                 class_weight='balanced', C=0.001)
        scaler = MinMaxScaler()

        logger.info('Concatenate old and new training feature matrices')
        train_df = pd.concat([orig_train_df, new_train_df], ignore_index=True)
        X = train_df[train_features]
        y = train_df['real']

        logger.info('Scaling training feature matrix')
        X = scaler.fit_transform(X)

        logger.info('Training model')
        if grid_search:
            logger.info('Training using grid search')
            param = {'C': [1e-3, 1e-2, 1e-1, 1]}
            metrics = ['accuracy', 'precision', 'recall', 'roc_auc', 'f1']
            gs = GridSearchCV(clf, param, n_jobs=cores, scoring=metrics, cv=10, refit='f1')
            gs.fit(X, y)
            clf = gs.best_estimator_
        else:
            logger.info('Training using previous optimized parameters')
            clf.fit(X, y)
    else:
        clf = joblib.load(resources.files(models) / "trained.clf")
        scaler = joblib.load(resources.files(models) / "trained.scaler")

    if not outdir.exists():
        outdir.mkdir(parents=True)

    logger.info('Converting VCF to bed file')
    bed_file_path = outdir / prefix.with_suffix(".bed")
    os.system(f"zcat {vcf} | grep PASS | vcf2bed | cut -f 1-3,6-7 > {bed_file_path}")
    
    logger.info('Preparing data')
    prep_data = dp.PrepareData(prefix, bam, bed_file_path, ref, outdir)

    df = prep_data.training_data

    logger.info('Hard-filter variants')
    df = df[df.tumor_VAF > 0.05]
    df = df[df.tumor_depth > 10]
    df = df[df.tumor_var_num_minus_strand + df.tumor_var_num_plus_strand > 4]

    df = df.drop(['ref', 'var'], axis=1)

    logger.info('Scaling features')
    df_scaled = scaler.transform(df)

    logger.info('Obtaining predictions')
    preds = clf.predict(df_scaled)
    df['preds'] = preds
    df = df[df.preds == 1]

    logger.info('Filtering VCF')
    kept_vars = set(df.index.str.replace(prefix + '~', ''))

    vcf_file_path = os.path.join(outdir, prefix + '_filtered.vcf')
    with open(vcf_file_path, 'w') as f_out:
        with gzip.open(vcf, 'rt') as f_in:
            for line in f_in:
                if not line.startswith('#'):
                    split = line.split('\t')
                    var = '{}:{}-{}{}>{}'.format(split[0], split[1], split[1], split[3], split[4])
                    if var in kept_vars:
                        f_out.write(line)
                else:
                    f_out.write(line)
