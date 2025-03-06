import pandas as pd
import os
import pybedtools
from pybedtools import chromsizes
from itertools import chain
import numpy as np
import math
import pyBigWig
import fnmatch
import gzip
import re
from timeit import default_timer as clock
import subprocess
#from multiprocessing import Pool
import multiprocessing
import argparse
import shutil
import json
import sys
import GenerateInput.TSS_Fetcher as TSS_Fetcher
import GenerateInput.BigWig_Counter as BigWig_Counter

"""Script to write the input files (samples X region features filled with bigwig signals) for the pretrained models.
Please refer to the README at https://github.com/SchulzLab/ExpressionPredictionModels."""

# --------------------------------------------------------------------------------------------------
# Process input arguments
# --------------------------------------------------------------------------------------------------
multiprocessing.set_start_method('fork')  # Multiprocess otherwise can cause issues on MacOS.
parser = argparse.ArgumentParser()
parser.add_argument('json_file', help='Path to the JSON file with the run specifications.')
args = parser.parse_args()

input_dict = json.load(open(args.json_file))

for entry in ['bigwigs', 'gene_file', 'mode', 'out_folder', 'provided_input']:
    if entry not in input_dict:
        print("ERROR: missing entry in JSON file for " + entry)
        sys.exit()

# Get the file with the model performance per gene.
if os.path.isfile(input_dict['provided_input'] + '/' + 'ModelPerformances.txt.gz'):
    performance_file = input_dict['provided_input'] + '/' + 'ModelPerformances.txt.gz'
else:
    if os.path.isfile(input_dict['provided_input'] + '/' + 'ModelPerformances_chr21Genes.txt.gz'):
        performance_file = input_dict['provided_input'] + '/' + 'ModelPerformances_chr21Genes.txt.gz'
    else:
        print("ERROR: missing the mapping file for CREs to genes in the provided input folder. Expected" +
              input_dict['provided_input'] + '/' + 'ModelPerformances.txt.gz')
        sys.exit()

if 'cores' not in input_dict:
    input_dict['cores'] = 1

if 'correlation_cutoff' not in input_dict:
    input_dict['correlation_cutoff'] = 0

if not os.path.isdir(input_dict['out_folder']):
    os.mkdir(input_dict['out_folder'])

# Copy the used JSON file to the output directory.
shutil.copyfile(args.json_file, input_dict['out_folder'] + '/' + args.json_file.split('/')[-1])

# Fetch the bigwig file paths and how to name them in the output files.
if type(input_dict['bigwigs']) == list:
    bigwigs = {x: x.split('/')[-1] for x in input_dict['bigwigs']}
elif os.path.isdir(input_dict['bigwigs']):  # Take all bigwigs in the folder.
    bigwigs = {input_dict['bigwigs'] + '/' + x: x.split('/')[-1] for x in os.listdir(input_dict['bigwigs'])
               if x.lower().endswith('bigwig') or x.lower().endswith('.bw')}
else:  # Assume pattern matching.
    bw_folder = '/'.join(input_dict['bigwigs'].split('/')[:-1])
    bw_pattern = input_dict['bigwigs'].split('/')[-1]
    re_pattern = re.compile(bw_pattern.replace('*', '(.*)'))
    bigwigs = {bw_folder + '/' + x: re_pattern.search(x).group(1)
               for x in os.listdir(bw_folder) if fnmatch.fnmatch(x, bw_pattern)}

gene_set = set([x.strip().split('\t')[0] for x in open(input_dict['gene_file']) if not x.startswith('#')])

# Check for potential performance cutoff and for failed models.
gene_misses = {g: [] for g in gene_set}  # As list, so we can append strings per mode.
writable_genes = {m: set() for m in ['CRE-RF', 'Binned-CNN']}
performance_df = pd.read_table(performance_file, sep='\t', header=0, index_col='Ensembl ID')
performance_df = performance_df[performance_df.index.isin(gene_set)]
for gene in gene_set:
    # If a gene is not in the performance table, it didn't pass the filtering and was not considered for training.
    if gene not in performance_df.index:
        gene_misses[gene].append("Not considered for training")
    else:
        for mode in ['CRE-RF', 'Binned-CNN']:
            if input_dict['mode'] == 'all' or mode.lower() in input_dict['mode'].lower():
                if performance_df.loc[gene][mode] == 0:  # A 0 means the model failed entirely.
                    gene_misses[gene].append(mode+" model training failed")
                elif performance_df.loc[gene][mode] < input_dict['correlation_cutoff']:
                    gene_misses[gene].append(mode+" model performance below correlation cutoff")
                else:
                    writable_genes[mode].add(gene)
open(input_dict['out_folder'] + '/FailedGenes_Input.txt', 'w').write('\n'.join([g + '\t' + ', '.join(val)
                                                                          for g, val in gene_misses.items() if val]))

# --------------------------------------------------------------------------------------------------
# CRE-mode
# --------------------------------------------------------------------------------------------------
if input_dict['mode'] == 'all' or 'cre' in input_dict['mode'].lower():
    cre_out_folder = input_dict['out_folder'] + '/' + 'CRE_input/'
    if not os.path.isdir(cre_out_folder):
        os.mkdir(cre_out_folder)
    # Get the path to the gene-CRE mapping, only allow the name of the full one and the mini example.
    if os.path.isfile(input_dict['provided_input'] + '/' + 'CRE_GenePeakMap.txt.gz'):
        gene_cre_file = input_dict['provided_input'] + '/' + 'CRE_GenePeakMap.txt.gz'
    else:
        if os.path.isfile(input_dict['provided_input'] + '/' + 'CRE_GenePeakMap_Example.txt.gz'):
            gene_cre_file = input_dict['provided_input'] + '/' + 'CRE_GenePeakMap_Example.txt.gz'
        else:
            print("ERROR: missing the mapping file for CREs to genes in the provided input folder. Expected" +
                  input_dict['provided_input'] + '/' + 'CRE_GenePeakMap.txt.gz')
            sys.exit()

    startcre = clock()
    print("Writing input files for CRE mode")
    # First get the mapping of regions to genes.
    gene_peaks_map = {}
    with gzip.open(gene_cre_file, 'rt') as map_in:
        for row in map_in:
            gene = row.strip().split('\t')[0]
            if gene in writable_genes['CRE-RF']:
                gene_peaks_map[row.strip().split('\t')[0]] = row.strip().split('\t')[1].split(',')  # Fixed order.
    # Get the set of regions that are within reach of genes and get their signal in the bigwig files.
    window_regions = pybedtools.BedTool('\n'.join(['chr'+x.replace('-', '\t') for x in set(chain(*gene_peaks_map.values()))]),
                                        from_string=True)
    peaks_counts, errors = BigWig_Counter.bigwig_counts(window_regions, list(bigwigs.keys()), n_cores=input_dict['cores'])
    peaks_counts.columns = ['#chr', 'start', 'end'] + [bigwigs[x] for x in peaks_counts.columns[3:]]
    # We use hyphens for the peak name to later be able to write them to output without accidental tabs.
    peaks_counts_dict = {'-'.join([str(y).replace('chr', '') for y in val[:3]]): val[3:] for val in peaks_counts.values}

    def write_gene(this_gene):
        window_out = cre_out_folder + this_gene + '.txt'
        these_peaks = gene_peaks_map[this_gene]
        if these_peaks:
            with open(window_out, 'w') as window:
                window.write('Sample\t' + '\t'.join(these_peaks).replace('chr', '') + '\n')
                for s_idx, sample in enumerate(bigwigs.values()):
                    window.write(sample + '\t')
                    window.write('\t'.join([f'{round(peaks_counts_dict[peak][s_idx], 5):g}' for peak in these_peaks]) + '\n')
        subprocess.call("gzip " + window_out, shell=True)

    process_pool = multiprocessing.Pool(processes=input_dict['cores'])
    process_pool.map(write_gene, [g for g in writable_genes['CRE-RF'] if g in gene_peaks_map])
    process_pool.close()
    print('Input for CRE mode written', clock() - startcre)

# --------------------------------------------------------------------------------------------------
# Binned-mode
# --------------------------------------------------------------------------------------------------
if input_dict['mode'] == 'all' or 'binned' in input_dict['mode'].lower():
    binned_out_folder = input_dict['out_folder'] + '/' + 'Binned_input/'
    if not os.path.isdir(binned_out_folder):
        os.mkdir(binned_out_folder)
    # Get the path to the gtf-fle, only allow the name of the full one and the mini example.
    if os.path.isfile(input_dict['provided_input'] + '/' + 'gencode.v38.annotation.gtf.gz'):
        gtf_file = input_dict['provided_input'] + '/' + 'gencode.v38.annotation.gtf.gz'
    else:
        if os.path.isfile(input_dict['provided_input'] + '/' + 'gencode.v38.annotation_chr21Genes.gtf'):
            gtf_file = input_dict['provided_input'] + '/' + 'gencode.v38.annotation_chr21Genes.gtf'
        else:
            print("ERROR: missing the gtf-file in the provided input folder. Expected" +
                  input_dict['provided_input'] + '/' + 'gencode.v38.annotation.gtf.gz')
            sys.exit()

    start_bin = clock()
    print("Writing input files for Binned mode")
    window_size = 1000000
    bin_size = 100
    tss_locs = TSS_Fetcher.gene_window_bed(gtf_file, tss_type='5', dict_only=True, gene_set=writable_genes['Binned-CNN'])
    chrom_boundaries = chromsizes('hg38')

    def get_gene_bins(bin_gene):
        """Bin 1 MB window around a gene's 5' into bins of size 100 bp and get the signals from each bigwig file.
        Bins outside the chromosome boundaries are skipped, meaning that not all genes will have 10,000 entries."""
        locs = tss_locs[bin_gene]
        gene_out = binned_out_folder + bin_gene + ".txt"
        # locs["tss"] is a set with only the 5'.
        upper_end = min(list(locs["tss"])[0] + window_size // 2, chrom_boundaries[locs['chr']][1])
        lower_end = max(list(locs["tss"])[0] - window_size // 2, chrom_boundaries[locs['chr']][0])
        possible_bins = math.floor((upper_end - lower_end) // bin_size)
        bins = ['-'.join([locs['chr'].replace('chr', ""), str(lower_end + b * bin_size), str(lower_end + (b + 1) * bin_size)]) for b in
                range(possible_bins)]
        gene_counts = np.zeros((len(bigwigs), possible_bins), dtype="float64")
        for n_bw, bw_file in enumerate(bigwigs.keys()):
            this_bw = pyBigWig.open(bw_file)
            sample_counts = this_bw.stats(locs['chr'], lower_end, upper_end, type="mean", nBins=possible_bins)
            gene_counts[n_bw][:] = sample_counts

        # Now add the expression and overwrite the bed file.
        gene_df = pd.DataFrame(gene_counts, index=bigwigs.values(), columns=bins).round(5)
        gene_df.fillna(0, inplace=True)
        gene_df.index.name = 'Sample'
        gene_df.to_csv(gene_out, header=True, index=True, sep='\t')

        subprocess.call("gzip " + gene_out, shell=True)

    process_pool = multiprocessing.Pool(processes=input_dict['cores'])
    process_pool.map(get_gene_bins, list(writable_genes['Binned-CNN']))
    process_pool.close()
    print('Input for Binned mode written', clock() - start_bin)

