#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
import os
import argparse
from scipy.stats import hypergeom
import statsmodels.stats.multitest as mt
from tqdm import tqdm
import time

import ld_proxy
import logger

def parse_gwas(gwas_fp, logger):
    #print('Parsing GWAS associations...')
    cols = ['SNPS', 'SNP_ID_CURRENT', 'DISEASE/TRAIT', #'MAPPED_GENE',
            'P-VALUE',  'OR or BETA', '95% CI (TEXT)']
    if gwas_fp is None:
        fp = 'https://www.ebi.ac.uk/gwas/api/search/downloads/full'
        df = pd.read_csv(fp, sep='\t', usecols=cols, low_memory=False)
    else:
        if not os.path.isfile(gwas_fp):
            return
        df = pd.read_csv(gwas_fp, sep='\t', usecols=cols, low_memory=False)
    df['SNP_ID_CURRENT'] = df['SNP_ID_CURRENT'].fillna('').apply(lambda x: clean_snp_ids(x))
    df = df.assign(SNPS=df['SNPS'].str.split(';')).explode('SNPS')
    df = df.assign(SNPS=df['SNPS'].str.split(',')).explode('SNPS')
    df['SNPS'] = df['SNPS'].str.strip()
    return df

def clean_snp_ids(snp_id):
    try:
        return f"rs{int(snp_id)}"
    except:
        if snp_id == '':
            return snp_id
        else:
            snp_id = snp_id.split('-')[0]
            snp_id = snp_id.split('_')[0]
            return f"rs{snp_id}" if not snp_id[:-1].isdigit() else f"rs{snp_id[:-1]}"

#Called within comorbid.py's pipeline() as such:
#    sig_res = find_snp_disease.find_disease(
#        gwas, output_dir, output_dir, args.ld, args.correlation_threshold,
#        args.window, args.window_size, args.window_control, args.population, args.ld_dir, logger, bootstrap=bootstrap)

#multimorbid3dm: I add a new/additional argument specification for this function -> window_size AND window_control
def find_disease(gwas, ppin_dir, out, ld, corr_thresh, window, window_size, window_control, population, ld_dir,
                 logger, disable_pg=True, bootstrap=False):
    #logger.write('Identifying GWAS traits.')
    sig_res = []
    probs_res = []
    eqtl_fps = [fp for fp in sorted(os.listdir(ppin_dir)) if fp.endswith('snp_gene.txt')]
    cols = ['level', 'trait', 'total_gwas_snps', 'trait_snps',
                  'eqtls_in_catalog', 'trait_eqtls', 'pval']
    # Total GWAS SNPs.
    M = gwas['SNPS'].nunique()
    for level_fp in eqtl_fps:
        level = f"{os.path.basename(level_fp).split('.')[0].split('_')[0]}"
        logger.write(f"{level}")          ##########TEMPORARY CHECK#############
        #logger.write(f"\t{level}")
        df = pd.read_csv(os.path.join(ppin_dir, level_fp), sep='\t')
        if df.empty:
            continue
        snps = df['snp'].drop_duplicates().tolist()
        if ld:
            #multimorbid3dm: I add "window_size" AND "window_control" here
            ld_snps = ld_proxy.ld_proxy(snps, corr_thresh, window, window_size, window_control, population, ld_dir, logger, bootstrap)
            if not ld_snps.empty:
                snps = ld_snps['rsidt'].drop_duplicates().tolist()

                logger.write('START 2nd merge') ##########TEMPORARY CHECK#############
                df = (df.merge(ld_snps, left_on='snp', right_on='rsidq')
                      .sort_values(by=['snp', 'gene', 'dprime'])
                      .drop_duplicates(subset=['snp', 'rsidt', 'gene', 'dprime']))
                logger.write('DONE 2nd merge') ##########TEMPORARY CHECK#############

                write_results(df,  f'{level}_snp_gene.txt', out)
        overlap = gwas.loc[gwas['SNPS'].isin(snps)].drop_duplicates()
        # Total level eQTLs that are in the GWAS Catalog. 
        N = overlap['SNPS'].nunique()
        probs_df = []
        snp_trait_df = []
        for trait in tqdm(overlap['DISEASE/TRAIT'].drop_duplicates(),
                          disable=disable_pg):
            # Total trait-associated SNPs in GWAS Catalog 
            n = gwas[gwas['DISEASE/TRAIT']==trait]['SNPS'].nunique()
            trait_overlap = overlap[overlap['DISEASE/TRAIT'] == trait][['SNPS', 'DISEASE/TRAIT']]
            # Total trait-associated eQTLs
            X = trait_overlap['SNPS'].nunique()
            # See https://alexlenail.medium.com/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458
            pval = hypergeom.sf(X-1, M, n, N)
            probs_df.append((level, trait, M, n, N, X, pval))
            snp_trait_df.append(trait_overlap.drop_duplicates())
        probs_df = pd.DataFrame(probs_df, columns=cols)
        probs_res.append(probs_df)
        try:
            probs_df['adj_pval'] = mt.multipletests(probs_df['pval'], method='bonferroni')[1]
        except:
            continue
        sig_df = probs_df[probs_df['adj_pval'] < 0.05]
        sig_res.append(sig_df)
        snp_trait_df = pd.concat(snp_trait_df)
        snp_trait_df = snp_trait_df.rename(columns={'DISEASE/TRAIT': 'trait',
                                                    'SNPS': 'snp'})
        snp_trait_df = snp_trait_df[snp_trait_df['trait'].isin(sig_df['trait'])]
        if not bootstrap:
            if ld:
                write_results(
                    snp_trait_df.merge(df.drop(columns=['snp']),
                                       how='inner', left_on='snp', right_on='rsidt'),
                    f'{level}_sig_interactions.txt', out)
            else:
                write_results(
                    snp_trait_df.merge(df, how='inner'),
                    f'{level}_sig_interactions.txt', out)
            write_results(probs_df,  f'{level}_enrichment.txt', out)
    if len(probs_res) > 0 and not bootstrap:
        probs_res = pd.concat(probs_res)
        write_results(probs_res, 'enrichment.txt',  out)
    if len(sig_res) == 0:
        sig_res = pd.DataFrame(columns=cols)
    else:
        sig_res = pd.concat(sig_res)
    #if not bootstrap:
    write_results(sig_res, 'significant_enrichment.txt',  out)
    return sig_res

def write_results(res, fp, out):
    #logger.write(f'\tWriting {fp}...')
    os.makedirs(out, exist_ok=True)
    res.to_csv(os.path.join(out, fp), sep='\t', index=False)

        
def parse_args():
    parser = argparse.ArgumentParser(
        description='Find disease associated with PPIN eQTLs.')
    parser.add_argument(
        '-p', '--ppin-dir', required=True,
        help='Filepath to directory containing eQTL-gene pairs for PPIN levels.')
    parser.add_argument(
        '-o', '--output-dir', required=True, help='Directory to write results.')
    parser.add_argument(
        '--gwas', default=None,
        help='''Filepath to GWAS associations. 
        Default: Associations from the GWAS Catalog 
        (https://www.ebi.ac.uk/gwas/api/search/downloads/full) ''')
    parser.add_argument(
        '--ld', action='store_true', default=False,
        help='Include LD SNPs in identifying eQTLs and GWAS traits.')
    parser.add_argument(
        '-c', '--correlation-threshold', default=0.8, type=int,
        help='The r-squared correlation threshold to use.')
    parser.add_argument(
        '-w', '--window', default=5000, type=int,
        help='The genomic window (+ or - in bases) within which proxies are searched.')
    parser.add_argument(
        '--population', default='EUR', choices=['EUR'],
        help='The ancestral population in which the LD is calculated.')
    parser.add_argument(
        '--ld-dir', default=os.path.join(os.path.dirname(__file__),'data/ld/dbs/super_pop/'),
        help='Directory containing LD database.')
    return parser.parse_args()


if __name__=='__main__':
    args = parse_args()
    start_time = time.time()
    os.makedirs(args.output_dir, exist_ok=True)
    logger = logger.Logger(logfile=os.path.join(args.output_dir, 'find_snp_disease.log'))
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {getattr(args, arg)}')
    gwas = parse_gwas(args.gwas, logger)
    find_disease(gwas, args.ppin_dir, args.output_dir,
                 args.ld, args.correlation_threshold, args.window, args.population,
                 args.ld_dir, logger)
    #write_results(res, args.output_dir)
    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')
