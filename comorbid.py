#! /usr/bin/env python

import pandas as pd
import istarmap
import argparse
import requests
import sys
import os
import stringdb_params as sdb
from io import StringIO
import time
import shutil
import multiprocessing as mp
from tqdm import tqdm
from itertools import repeat

import query_grn
import ppin
import query_string
import query_proper
import get_ppi_eqtls
import find_snp_disease
import logger


def write_results(df, output_fp, logger):
    out_dir = os.path.dirname(output_fp)
    os.makedirs(out_dir, exist_ok=True)
    logger.write('Writing output...')
    df.to_csv(output_fp, sep='\t', index=False)



def parse_input(inputs):
    '''Return a dataframe of gene input.'''
    logger.write('Parsing input...')
    df = pd.DataFrame()
    if os.path.isfile(inputs[0]): # Input is file.
        df = pd.read_csv(inputs[0], sep='\t')
        df = df[['gene']].drop_duplicates()
    else: # Input probably space-separated list of genes.
        df = pd.DataFrame({'gene': [i.upper() for i in inputs]})
        df = df[['gene']].drop_duplicates()
    return df


def parse_args():
    parser = argparse.ArgumentParser(
        description='Identify multimorbid traits based on eQTL associations and protein-protein interactions.')
    parser.add_argument(
        '-g', '--genes', nargs='+',
        help='''A space-separated list of gene symbols or filepath to a file 
        containing gene symbols in the 'gene' column.''' )
    parser.add_argument(
        '-s', '--snps', nargs='+',
        help='''A space-separated list of SNP rsIDs or filepath to a file 
        containing SNP rsids in the 'snp' column.''' )
    parser.add_argument(
        '--trait',
        help='''GWAS trait to query. Note: this flag is mutually exclusive with
        the --snps and --pmid flags''')
    parser.add_argument(
        '--pmid',
        help='''PubMed ID of the GWAS to query. Note: this flag is mutually exclusive with 
        the --snps and --trait flag''')
    parser.add_argument(
        '--grn-dir', required=True,
        help='''Directory containing tissue gene regulatory network.
        The subdirectories should contain significant_eqtls.txt for each chromosome.''')
    parser.add_argument(
        '--gwas', default=None,
        help='''Filepath to GWAS associations.
        Default: Associations from the GWAS Catalog
        (https://www.ebi.ac.uk/gwas/api/search/downloads/full) ''')
    parser.add_argument(
        '-o', '--output-dir', required=True,
        help='Directory to write results.')
    parser.add_argument(
        '-l', '--levels', default=1, type=int,
        help='Path length (i.e. number of nodes) to query. Default = 1')
    parser.add_argument(
        '-p', '--ppin', nargs='+', default=['string', 'proper'],
        choices=['string', 'proper'],
        help='''The protein-protein-interaction database(s) to use.
        Default: ['string', 'proper']''')
    parser.add_argument(
        '--string-score', default=0.7, type=float,
        help='Cut-off score for STRING interactions. Default = 0.7')
    parser.add_argument(
        '--bootstrap', default=False, action='store_true',
        help='Perform a bootstrap. Default = False')
    parser.add_argument(
        '--bootstraps', default=1000, type=int,
        help='Number of bootstrap datasets. Default: 1000')
    parser.add_argument(
        '--keep-bootstraps', action='store_true', default=False,
        help='Keep bootstrap results. Default: False ')
    parser.add_argument(
        '--non-spatial', action='store_true', default=False,
        help='Include non-spatial eQTLs. Default = False')
    parser.add_argument(
        '--non-spatial-dir', default=os.path.join(os.path.dirname(__file__), 'data/GTEx/'),
        help='Filepath to non-spatial eQTLs.')
    parser.add_argument(
        '--snp-ref-dir', default=os.path.join(os.path.dirname(__file__), 'data/snps/'),
        help='Filepath to SNP BED databases.')
    parser.add_argument(
        '--gene-ref-dir', default=os.path.join(os.path.dirname(__file__), 'data/genes/'),
        help='Filepath to gene BED.')
    parser.add_argument(
        '--ld', action='store_true', default=False,
        help='Include LD SNPs in identifying eQTLs and GWAS traits. Default = False')
    parser.add_argument(
        '-c', '--correlation-threshold', default=0.8, type=int,
        help='The r-squared correlation threshold to use.')
    parser.add_argument(
        '-w', '--window', action='store_true', default=False,
        help='Use genomic distance as filter for ld')
    parser.add_argument(
        '-ws', '--window-size', default=5000, type=int,
        help='The genomic window (+ or - in bases) within which proxies are searched. Default = 5000')
    parser.add_argument(
        '-wc', '--window-control', default=5000, type=int,
        help='The genomic window (+ or - in bases) within which proxies are searched in the problematic locus on chromosome 17, highly recommended to keep this number low. Default = 5000')
    parser.add_argument(
         '--population', default='EUR', choices=['EUR'],
        help='The ancestral population in which the LD is calculated. Default = "EUR"')
    parser.add_argument(
        '--ld-dir', default=os.path.join(os.path.dirname(__file__), 'data/ld/dbs/super_pop/'),
        help='Directory containing LD database.')
    return parser.parse_args()


#Called within the execution part as such:
#        snps, genes = parse_snps(
#            args.snps, args.trait, args.pmid, gwas, grn, args.output_dir,
#            args.non_spatial, args.non_spatial_dir, args.snp_ref_dir, args.gene_ref_dir,    #multimorbid3dm homework, implement args.window, args.window_size, args.window_control here somehow!!!
#            args.ld, args.correlation_threshold, args.window, args.window_size, args.window_control, args.population, args.ld_dir,
#            logger)

#multimorbid3dm hw realization: added window_size AND window_control here!
def parse_snps(snp_arg, trait_arg, pmid_arg, gwas, grn, output_dir,
               non_spatial, non_spatial_dir, snp_ref_dir, gene_ref_dir,
               ld, corr_thresh, window, window_size, window_control, population, ld_dir, logger):
    if (snp_arg and trait_arg) or (snp_arg and pmid_arg) or (trait_arg and pmid_arg):
        sys.exit('Only one of --snps, --trait, or --pmid is required.\nExiting.')
    snps = pd.DataFrame()
    logger.write('Parsing SNP input...')
    if snp_arg:
        if os.path.isfile(snp_arg[0]):
            df = pd.read_csv(snp_arg[0], sep='\t', header=None, names=['snp'])
            df = df[df['snp'] != 'snp']
            df['snp'] = df['snp'].str.strip()
            snps = df['snp'].drop_duplicates()
        else:
            df = pd.DataFrame({'snp':  snp_arg})
            snps = df['snp'].drop_duplicates()
    elif trait_arg:
        snps = query_grn.extract_trait_snps(trait_arg, gwas, logger)
    elif pmid_arg:
        snps = query_grn.extract_pmid_snps(pmid_arg, gwas, logger)
    #multimorbid3dm hw realization: added window_size AND window_control here!!
    eqtls = query_grn.get_eqtls(snps, grn, output_dir,
                                non_spatial, non_spatial_dir, snp_ref_dir, gene_ref_dir,
                                ld, corr_thresh, window, window_size, window_control, population, ld_dir, logger)
    return snps, eqtls

def parse_genes(genes_args, logger):
    logger.write('Parsing gene input...')
    df = pd.DataFrame()
    if os.path.isfile(genes_args[0]):
        df = pd.read_csv(genes_args[0], sep='\t')
        df = df[['gene']].drop_duplicates()
    else:
        df = pd.DataFrame({'gene': [i.upper() for i in genes_args]})
        df = df[['gene']].drop_duplicates()
    return df

def join_path(*args):
    fp = ''
    for arg in args:
        fp = os.path.join(fp, arg)
    return fp

#Called within the execution part as such:
#sig_res = pipeline(genes, gwas, args.output_dir,  args, logger)
def pipeline(genes, gwas, output_dir, args, logger, bootstrap=False):
    # PPIN
    gene_list = []
    graph = pd.DataFrame()
    if not bootstrap:
        logger.write('Identifying protein interactions...')
    gene_list = ppin.make_ppin(
        genes, args.levels, output_dir, args.ppin, args.string_score, logger, bootstrap=bootstrap)
    if sum([len(level) for level in gene_list]) == 0:
        return pd.DataFrame()
    # PPIN eQTLs
    if not bootstrap:
        logger.write('Identifying gene eQTLs...')
    query_grn.get_gene_eqtls(
        gene_list, grn, output_dir,
        args.non_spatial, args.non_spatial_dir, args.snp_ref_dir, args.gene_ref_dir,
        logger, bootstrap=bootstrap)
    # Traits
    if not bootstrap:
        logger.write('Identifying GWAS traits...')
    #multimorbid3dm: I add TWO additional argument for the input of the find_snp_disease.find_disease() function -> args.window_size AND args.window_control 
    sig_res = find_snp_disease.find_disease(
        gwas, output_dir, output_dir, args.ld, args.correlation_threshold,
        args.window, args.window_size, args.window_control, args.population, args.ld_dir, logger, bootstrap=bootstrap)
    return sig_res

def prep_bootstrap(sim, gene_num, sims_dir, res_dict, grn_genes, gwas, args):
    sim_output_dir = join_path(sims_dir, sim)
    sim_genes = pd.DataFrame(
        {'gene': grn_genes.sample(gene_num, random_state=int(sim)).tolist()})
    sim_res = pipeline(sim_genes, gwas, sim_output_dir, args, logger, bootstrap=True)
    if not sim_res is None:
        for i, row in sim_res.iterrows():
            k = row['level'] + '__' + row['trait']
            try:
                res_dict[k] += 1
            except:
                pass

            
def bootstrap_genes(sig_res, genes, gwas, num_sims, grn, args):
    logger.write('Preparing simulations...')
    gene_num = genes[genes['gene'].isin(grn['gene'])]['gene'].nunique()
    sims_dir = join_path(args.output_dir, 'bootstrap')
    manager = mp.Manager()
    res_dict = manager.dict()
    for _, row in sig_res.iterrows():
        k = row['level'] + '__' + row['trait']
        if not k in res_dict.keys():
            res_dict[k] = 0
    desc = 'Boostrapping'
    bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
    sims = [str(i) for i in range(num_sims)]
    if args.ld:
        for sim in tqdm(
                sims,
                total=len(sims), desc=desc, bar_format=bar_format,
                unit='simulations', ncols=80
        ):
            #logger.write(sim)
            prep_bootstrap(
                sim, gene_num, sims_dir, res_dict,
                grn['gene'].drop_duplicates(), gwas, args)
    else:
        with mp.Pool(4) as pool:
            for _  in tqdm(
                    pool.istarmap(
                        prep_bootstrap,
                        zip(
                            sims,
                            repeat(gene_num),
                            repeat(sims_dir),
                            repeat(res_dict),
                            repeat(grn['gene'].drop_duplicates()),
                            repeat(gwas),
                            repeat(args))
                    ),
                    total=len(sims), desc=desc, bar_format=bar_format,
                    unit='simulations', ncols=80
            ):
                pass
    sim_df = []
    for k in res_dict:
        ks = k.split('__')
        sim_df.append([ks[0], ks[1], res_dict[k]])
    if len(sim_df) == 0:
        logger.write('Warning: No output for simulations.')
        return
    sim_df = pd.DataFrame(sim_df, columns=['level', 'trait', 'sim_count'])
    # Using (count + 1) / (num_sims +1) See https://doi.org/10.1086/341527
    sim_df['sim_pval'] = round((sim_df['sim_count'] + 1) / (num_sims + 1), 3)
    res = (sig_res
           .merge(sim_df, on=['level', 'trait'], how='left')
           .fillna(0)
    )
    if not args.keep_bootstraps:
        shutil.rmtree(sims_dir)
    write_results(res,
                  join_path(args.output_dir, 'significant_enrichment_bootstrap.txt'),
                  logger)
    
        
if __name__=='__main__':
    pd.options.mode.chained_assignment = None
    args = parse_args()
    if not args.genes and not args.snps and not args.trait and not args.pmid:
        sys.exit('FATAL: One of --genes, --snps, --trait, or --pmid is required.\nExiting.')
    start_time = time.time()
    if args.genes and (args.snps or args.trait or args.pmid):
        sys.exit('Only one of --genes, --snps, --trait, or --pmid is required.\nExiting.')
    os.makedirs(args.output_dir, exist_ok=True)
    global logger
    logger = logger.Logger(logfile=os.path.join(args.output_dir, 'comorbid.log'))
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {getattr(args, arg)}')
    logger.write('\n')
    grn = query_grn.parse_grn(args.grn_dir, logger)
    gwas = query_grn.parse_gwas(args.gwas, logger)
    snps = []
    genes = pd.DataFrame()
    if args.genes:
        genes = parse_genes(args.genes, logger)
    else:
        #multimorbid3dm hw realization: added args.window_size AND args.window_control here!
        snps, genes = parse_snps(
            args.snps, args.trait, args.pmid, gwas, grn, args.output_dir,
            args.non_spatial, args.non_spatial_dir, args.snp_ref_dir, args.gene_ref_dir,
            args.ld, args.correlation_threshold, args.window, args.window_size, args.window_control, args.population, args.ld_dir,
            logger)
        if genes.empty:
            logger.write('Exiting: No gene targets for SNPs found.')
            sys.exit()
    sig_res = pipeline(genes, gwas, args.output_dir,  args, logger)
    write_results(sig_res, os.path.join(args.output_dir, 'significant_enrichment.txt'), logger)
    if sig_res.empty:
        logger.write('Exiting: No results found.')
        exit()
    # Bootstrap
    if args.bootstrap:
        bootstrap_genes(sig_res, genes, gwas, args.bootstraps, grn, args)
    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')
