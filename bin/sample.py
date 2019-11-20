#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""sample.py: some description"""

import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
import re
import ast
from multiprocessing import Pool
import os
from functools import partial


__author__ = "Qiongzi Qiu"
__copyright__ = "Copyright 2019, The xx Project"
__credits__ = ["Qiongzi Qiu"]
__license__ = "GPL"
__version__ = "3.0.1"
__maintainer__ = "Qiongzi Qiu"
__email__ = "3110102477@zju.edu.cn"
__status__ = "Production"


def main():
    args = usage()
    if args.separate == 'replace':
        separate = ['replace']
    elif args.separate != '':
        separate = [args.separate]
    else:
        separate = ''

    if not args.dataframe:
        mutation_df_outfile = args.outfile.replace('.out', 'replace.out')
        print(mutation_df_outfile)
        mutation_df = generate_mutation_df(args.mutation, args.LoF, separate)
        write_table(mutation_df, mutation_df_outfile)
    else:
        mutation_df = pd.read_table(args.dataframe, sep='\t')
        rc = re.compile(r"\[|\]|\'|,")
        mutation_df[['replace']] = mutation_df[['replace']].apply(lambda x: x.apply(lambda x: rc.sub('', x).split()))

    num_processes = max(len(os.sched_getaffinity(0)) - 10, 5)
    chunk_size = int(mutation_df.shape[0]/num_processes)
    chunks = [mutation_df.ix[mutation_df.index[i:i + chunk_size]] for i
              in range(0, mutation_df.shape[0], chunk_size)]
    pool = Pool(processes=num_processes)
    gene_pair_df = pd.concat(pool.map(partial(sample_sort, paralogous=args.paralogous), chunks))

    if not args.paralogous:
        gene_pair_df = gene_pair_df.groupby(['replace'])['replace'].apply(list).reset_index()

    write_table(gene_pair_df, args.outfile)


def usage():
    parser = argparse.ArgumentParser()
    parser.add_argument('-O', '--outfile', dest='outfile', help='')
    parser.add_argument('-DF', '--dataframe', dest='dataframe', default=False, help='')
    parser.add_argument('-M', '--mutation', dest='mutation', help='')
    parser.add_argument('-L', '--LoF', dest='LoF', default=False, help='')
    parser.add_argument('-S', '--separate', dest='separate', default='', help='')
    parser.add_argument('-P', '--paralogous', dest='paralogous', default=False, help='')
    args = parser.parse_args()
    return args


def write_table(gene_pair_df, outfile):
    gene_pair_df.to_csv(outfile, sep='\t', index=False)


def generate_mutation_df(mutation_file, LoF, separate):
    with open(mutation_file) as f:
        line = next((l for l in f if 'replace' in l), None)
        vep_field_names = line.split('replace')[-1].strip('">\n').split('|')
    with open(mutation_file) as f:
        line = next((l for l in f if 'replace' in l), None)
        cell_mutation_columns = line.strip('\n').split('\t')
    chunk_reader = pd.read_table(mutation_file, sep='\t', header=None, comment='#', iterator=True, chunksize=10000)
    cell_mutation = pd.concat([chunk for chunk in chunk_reader], ignore_index=True)
    cell_mutation.rename(columns=pd.Series(cell_mutation_columns), inplace=True)
    if LoF:
        cell_mutation['replace'] = cell_mutation['replace'].apply(
            lambda x: 'replace' if [l.split('|')[vep_field_names.index('replace')] for l in
                                    x.split("replace")[-1].split(',')].count('replace') > 0 else '')
    else:
        mutation_class = pd.read_table('replace', sep='\t', comment='#')
        cell_mutation['replace'] = mutation_class.iloc[0:mutation_class.shape[0], ]
    cell_mutation = cell_mutation[cell_mutation['replace'] == 'replace']
    cell_mutation['replace'] = cell_mutation['replace'].apply(
        lambda x: x.split("replace")[-1].split(',')[0].split('|')[vep_field_names.index('replace')])
    cell_mutation = pd.melt(cell_mutation.iloc[:, 9:], id_vars=['replace', 'replace'],
                                 value_vars=cell_mutation.iloc[:, 9:].columns[:-2].tolist(),
                                 var_name='replace', value_name='replace')
    cell_mutation = cell_mutation[~cell_mutation['replace'].str.contains(re.compile(r"(^0/0|\./\.)"))]
    cell_mutation = cell_mutation[['replace', 'replace', 'replace']].drop_duplicates()
    ko_list = ['replace', 'replace', 'replace', 'replace']
    cell_annotation = pd.read_table('replace', sep='\t')
    cell_annotation['replace'] = cell_annotation[cell_annotation['replace'] == 'replace'].groupby('replace')[
        'replace'].transform('size')
    cell_annotation['replace'] = cell_annotation['replace'].apply(lambda x: re.match(r"^(\w+-\w+)", x).group(0))
    mutation_df = cell_mutation.groupby(
            ['replace'])['replace'].apply(lambda x: x.unique()).reset_index(name='replace')
    if separate != '':
        missense_mutation_filter = cell_mutation[cell_mutation['replace'].isin(['replace'])].groupby(
            ['replace'])['replace'].apply(lambda x: x.unique()).reset_index(name='replace')
        mutation_df = mutation_df.merge(missense_mutation_filter[['replace', 'replace']],
                                        left_on='replace', right_on='replace', how='left')
        mutation_df = mutation_df.drop(['replace'], axis=1)
    if LoF:
        ko_mutation_filter = cell_mutation[cell_mutation['replace'] == 'replace'].groupby(
            ['replace'])['replace'].apply(lambda x: x.unique()).reset_index(name='ko_mutation_sample')
    else:
        if separate == '':
            ko_mutation_filter = cell_mutation[cell_mutation['replace'].isin(ko_list)].groupby(
                ['replace'])['replace'].apply(lambda x: x.unique()).reset_index(name='replace')
        else:
            ko_mutation_filter = cell_mutation[cell_mutation['replace'].isin(separate)].groupby(
                ['replace'])['replace'].apply(lambda x: x.unique()).reset_index(name='replace')
    mutation_df = mutation_df.merge(ko_mutation_filter[['replace', 'replace']],
                          left_on='replace', right_on='replace', how='left')
    mutation_df['replace'] = mutation_df.apply(lambda x: np.setdiff1d(x['replace'], x['replace']).tolist(), axis=1)
    mutation_df['replace'] = mutation_df['replace'].apply(lambda x: len(x))
    mutation_df = mutation_df.drop(['replace'], axis=1)
    return mutation_df


def sample_sort(mutation_df, paralogous):
    mutation_df = mutation_df[mutation_df['replace'] >= 3]
    if paralogous:
        gene_pair_df = pd.DataFrame(columns=['replace'])
    else:
        gene_pair_df = pd.DataFrame(columns=['replace'])
    cell_annotation = pd.read_table('file', sep='\t')
    cell_annotation['replace'] = cell_annotation[cell_annotation['replace'] == 'replace'].groupby('replace')[
        'replace'].transform('size')
    sample_type = cell_annotation['replace'].drop_duplicates()
    gene_family = pd.read_table('replace', sep='\t')
    paralogous_gene = pd.read_table('replace', sep='\t')
    for index, row in mutation_df.iterrows():
        mutation_tmp = row[0]
        for type_tmp in sample_type:
            ko_mutation_sample_list_tmp = cell_annotation[(cell_annotation['replace'] == 'replace') &
                                              (cell_annotation['replace'] == type_tmp) &
                                              (cell_annotation['replace'].isin(row['replace']))]['replace'].tolist()
            other_mutation_sample_list_tmp = cell_annotation[(cell_annotation['replace'] == 'replace') &
                                              (cell_annotation['replace'] == type_tmp) &
                                              (cell_annotation['replace'].isin(row['replace']))]['replace'].tolist()
            non_mutation_sample_list_tmp = cell_annotation[(cell_annotation['replace'] == 'replace') &
                                              (cell_annotation['replace'] == type_tmp) &
                                              (cell_annotation['replace'].isin(row['replace']))]['replace'].tolist()
            ko_mutation_sample_num_tmp = int(len(ko_mutation_sample_list_tmp))
            other_mutation_sample_num_tmp = int(len(other_mutation_sample_list_tmp))
            non_mutation_sample_num_tmp = int(len(non_mutation_sample_list_tmp))
            if paralogous:
                homologous_gene_list = paralogous_gene[paralogous_gene['replace'] == mutation_tmp][
                                           'replace'].tolist() + [mutation_tmp]
                add_tmp = [mutation_tmp, type_tmp, ko_mutation_sample_num_tmp,
                           other_mutation_sample_num_tmp, non_mutation_sample_num_tmp,
                           ','.join(ko_mutation_sample_list_tmp), ','.join(other_mutation_sample_list_tmp),
                           ','.join(non_mutation_sample_list_tmp)]
                gene_pair_df = gene_pair_df.append(pd.DataFrame([add_tmp[:1] + [x] + add_tmp[1:] for
                                                                 x in homologous_gene_list],
                                                                columns=gene_pair_df.columns), ignore_index=True)
            else:
                gene_family_list = pd.unique(
                            gene_family[gene_family['replace'] == mutation_tmp]['replace']).tolist()
                for gene_family_tmp in gene_family_list:
                    homologous_gene_list = pd.unique(gene_family[gene_family['replace'] == gene_family_tmp]['replace'])
                    add_tmp = [mutation_tmp, gene_family_tmp, type_tmp, ko_mutation_sample_num_tmp,
                               other_mutation_sample_num_tmp, non_mutation_sample_num_tmp,
                               ','.join(ko_mutation_sample_list_tmp), ','.join(other_mutation_sample_list_tmp),
                               ','.join(non_mutation_sample_list_tmp)]
                    gene_pair_df = gene_pair_df.append(pd.DataFrame([add_tmp[:1] + [x] + add_tmp[1:] for
                                                                     x in homologous_gene_list],
                                                                    columns=gene_pair_df.columns), ignore_index=True)
    return gene_pair_df


if __name__ == "__main__":
    main()
