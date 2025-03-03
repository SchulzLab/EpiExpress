import gzip
import os
import pandas as pd
import numpy as np
import random
import pybedtools
import gseapy as gp
from itertools import chain
from timeit import default_timer as clock
import sys
from multiprocessing import Pool
# sys.path.append("*ISP/")  # Append the directory with the scripts if it's not the working directory.
import eQTL_Helpers

random.seed(1234)

"""From the different models get a specified number of top scored interactions in a sample and 
calculate the recall for the eQTL-gene pairs in the matching GTEx tissue. The matching of samples and 
input files for the models were prepared beforehand."""


def match_colname(c_name):
    if 'ensembl' in c_name.lower():
        return 'Ensembl ID'
    elif 'prediction' in c_name.lower():
        return 'predicted'
    elif 'ism' in c_name.lower():
        return 'ISM'
    elif c_name == 'chrom':
        return 'chr'
    else:
        return c_name


tag = 'RecallRun'
out_dir = "*/eQTLs_" + tag + '_'
input_folder = "*/ValidateInteractions/"  # Directory with a folder written for each model, with the in-silico mutagenesis simulation.

kept_genes_file = "*/kept_genes.txt"
kept_genes = set([x.strip().split('\t')[1] for x in open(kept_genes_file).readlines()[1:]])

models = ['ENCODE_RF_ism_all_reverse', 'ENCODE_MLP_ism_all_reverse', 'STITCHIT_ism_all_eqtl_backscaled',
          'binnedRF_os_ism_all_expand_len10', 'CNN_bestwarm_os_ism_all_expand_len10']

hg38_annotation = "*/gencode.v38.annotation.gtf"
hg38_tss = eQTL_Helpers.gene_window_bed(hg38_annotation, tss_type='5', dict_only=True)
ihec_mapping_file = "*/ValidateInteractions/IHEC_ValidationSamples.txt"
ihec_mapping_df = pd.read_table(ihec_mapping_file, sep='\t', header=0)
gtex_tissues = {'Thyroid', 'Whole_Blood', 'Brain_Cortex', 'Nerve_Tibial'}
tissue_order = sorted(list(gtex_tissues))

# Get a more direct map of GTEx tissues to EpiRR IDs for easier iterations.
tissue_sample_map = {}
for tissue in gtex_tissues:
    match_samples = ihec_mapping_df[ihec_mapping_df['MatchedOntology'] == tissue]
    tissue_sample_map[tissue] = set(match_samples['EpiRR'])
sample_tissue_map = {x[0]: x[1] for x in chain(*[[[s, k] for s in vals] for k, vals in tissue_sample_map.items()])}

base_folder = "*/"  # Folder with eQTL files from GTEx.
eqtl_types = {'CAVIAR': base_folder + "/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz",
              'CaVEMaN': base_folder + "/GTEx_v8_finemapping_CaVEMaN/GTEx_v8_finemapping_CaVEMaN.txt.gz",
              'DAP-G': base_folder + "/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.CS95.txt.gz"}

max_distance = 500000
eqtl_beds = {t: {} for t in gtex_tissues}
tissue_genes = {t: set() for t in gtex_tissues}  # To know the genes with eQTL in a tissue.

# Get a mapping of the minimum number of eQTL hits for each sample and eQTL type from a previous run.
minimum_hits_file = '/projects/apog/work/IHEC/ValidateInteractions/eQTLs_DFs/eQTLs_MinimumOverlap.txt'
minimum_hits = pd.read_table(minimum_hits_file, sep='\t', header=0, index_col=0).to_dict()

for tissue in gtex_tissues:
    print(tissue)
    these_tissues = tissue.split(',')  # We have cases where we combine multiple GTEx tissues.
    for eqtl in eqtl_types:
        eqtl_list = []

        if eqtl == "CAVIAR":
            with gzip.open(eqtl_types[eqtl], 'rt') as caviar:
                caviar_head = {x: i for i, x in enumerate(caviar.readline().strip().split('\t'))}
                for entry in caviar:
                    entry = entry.strip().split('\t')
                    if entry[caviar_head['TISSUE']] in these_tissues:
                        this_eqtl = 'chr' + entry[caviar_head['eQTL']]
                        this_gene = entry[caviar_head['GENE']].split('.')[0]
                        if this_gene in hg38_tss and \
                                abs(int(this_eqtl.split('_')[1]) - list(hg38_tss[this_gene]['tss'])[0]) <= max_distance:
                            eqtl_list.append(this_gene + '\t' + this_eqtl)

        if eqtl == "CaVEMaN":
            with gzip.open(eqtl_types[eqtl], 'rt') as caveman:
                cave_head = {x: i for i, x in enumerate(caveman.readline().strip().split('\t'))}
                for entry in caveman:
                    entry = entry.strip().split('\t')
                    if entry[cave_head['TISSUE']] in these_tissues:
                        this_eqtl = '_'.join(entry[cave_head['eQTL']].split('_')[:2])
                        this_gene = entry[cave_head['GENE']].split('.')[0]
                        if this_gene in hg38_tss and \
                                abs(int(this_eqtl.split('_')[1]) - list(hg38_tss[this_gene]['tss'])[0]) <= max_distance:
                            eqtl_list.append(this_gene + '\t' + this_eqtl)

        if eqtl == "DAP-G":
            with gzip.open(eqtl_types[eqtl], 'rt') as dapg:
                for entry in dapg:
                    entry = entry.strip().split('\t')
                    this_eqtl = '_'.join(entry[2].split('_')[:2])
                    for sub_entry in entry[5].split('|'):
                        if sub_entry and sub_entry.split('@')[1].split('=')[0] in these_tissues:
                            this_gene = sub_entry.split('.')[0]
                            if this_gene in hg38_tss and \
                                    abs(int(this_eqtl.split('_')[1]) - list(hg38_tss[this_gene]['tss'])[
                                        0]) <= max_distance:
                                eqtl_list.append(this_gene + '\t' + this_eqtl)

        # Set still necessary if we have multiple GTEx tissues.
        eqtl_list = [x for x in eqtl_list if x.split('\t')[0] in kept_genes]
        eqtl_bed = pybedtools.BedTool('\n'.join(
            set([x.split('\t')[0] + '\t' + str(int(x.split('_')[1]) - 1) + '\t' + x.split('_')[1] for x in eqtl_list])),
                                      from_string=True)
        eqtl_beds[tissue][eqtl] = eqtl_bed
        tissue_genes[tissue] |= set([x.fields[0] for x in eqtl_bed])
print('eQTL-beds collected')


def file_header(file_path):
    """Check number of reader lines."""
    header_rows = []
    if file_path.endswith('.gz'):
        file_opener = gzip.open(file_path, 'rt')
    else:
        file_opener = open(file_path)
    with file_opener as read_header:
        for row in read_header:
            if not row.startswith('#'):
                break
            header_rows.append(row.strip())

    if not header_rows:
        print("No headers in", model, file_path)
        return []

    return header_rows


def get_sample_recall(args):
    """Get the recall for eQTLs for a sample. The function allows parallelization with multiprocess."""
    model = args[0]
    sample = args[1]
    plot_gsea = args[2]
    print(sample)
    match_folder = [input_folder + '/' + x for x in os.listdir(input_folder) if x == model][0] + '/' + sample
    sample_scores = []
    df_list = []
    for file in [match_folder + '/' + x for x in os.listdir(match_folder) if
                 x.endswith('.bed') or x.endswith('.txt') or x.endswith('.gz') or x.endswith('.zip')]:
        this_gene = file.split('/')[-1].split('.')[0].split('_')[0]
        if this_gene not in kept_genes or this_gene not in tissue_genes[sample_tissue_map[sample]]:
            continue

        header = file_header(file)
        this_df = pd.read_table(file, header=0, sep='\t', skiprows=len(header))
        this_df.columns = [match_colname(c) for c in this_df.columns]
        this_df['chr'] = [str(x).replace('chr', '') for x in this_df['chr'].values]
        this_df['Ensembl ID'] = this_gene
        this_df = this_df.astype({'start': 'int', 'end': 'int', "ISM": 'float64', "predicted": 'float64'})
        if this_df['ISM'].min() < 0 or this_df['predicted'].min() < 0:
            this_df['ISM'] = this_df['ISM'] + abs(min([this_df['ISM'].min(), this_df['predicted'].min()]))
            this_df['predicted'] = this_df['predicted'] + abs(min([this_df['ISM'].min(), this_df['predicted'].min()]))

        this_df['score'] = np.log2((this_df['ISM'] + 1) / (this_df['predicted'] + 1)).abs()
        if this_df['score'].sum() > 0:
            this_df["score normG"] = this_df['score'] / this_df['score'].sum()
        else:
            this_df['score normG'] = 0
        df_list.append(this_df)

    if df_list:
        collected_df = pd.concat(df_list)
    else:
        print("No remaining genes", model)
        return []

    print(sample, model, len(collected_df), 'total interactions')

    for score_col in ['score', 'score normG']:
        sorted_score_df = collected_df.sort_values(by=score_col, ascending=False)
        this_score_df = sorted_score_df.iloc[:100000]
        this_score_df.index = ['\t'.join([str(x) for x in val]) for val in
                               this_score_df[['Ensembl ID', 'start', 'end']].values]
        this_score_df = this_score_df[[score_col]]
        collected_bed = pybedtools.BedTool('\n'.join(this_score_df.index), from_string=True)

        for eqtl in eqtl_types:
            start_e = clock()
            inter_hits = collected_bed.intersect(eqtl_beds[sample_tissue_map[sample]][eqtl], u=True)
            downsampled_hits = set(
                random.sample(set(['\t'.join(x.fields[:3]) for x in inter_hits]), minimum_hits[eqtl][sample]))
            re_res = gp.prerank(rnk=this_score_df,  # CARE they sort themselves again by the first column.
                                gene_sets={
                                    model + '; ' + score_col + '; ' + sample_tissue_map[sample]: downsampled_hits},
                                min_size=1,
                                max_size=len(eqtl_beds[sample_tissue_map[sample]][eqtl]) + 1,
                                permutation_num=100,  # reduce number to speed up testing
                                # If set creates a directory at that path.
                                outdir=None if not plot_gsea else out_dir + model + "_" + score_col.replace(" ",
                                                                                                            '') + "_" + sample + "_" + eqtl,
                                seed=1234,
                                threads=1,
                                weight=0,
                                no_plot=not plot_gsea,
                                verbose=False,
                                ).res2d
            print(clock() - start_e, eqtl)
            sample_scores.append({"EpiRR": sample,
                                  'eqtl': eqtl,
                                  'model': model + (' normG') if 'normG' in score_col else model,
                                  'ranks': 100000,
                                  'eQTL-gene pair hits': len(set(['\t'.join(x.fields[:3]) for x in inter_hits])),
                                  'ES': re_res.iloc[0]['ES'],  # We only have one entry.
                                  'NES': re_res.iloc[0]['NES'],
                                  'NOM p-val': re_res.iloc[0]['NOM p-val']})

    return sample_scores


# Get the interactions from the models.
collected_recalls = []
for model in models:
    print(model)
    startw = clock()
    process_pool = Pool(processes=20)
    model_recalls = process_pool.map(get_sample_recall,
                                     [[model, sample, False] for sample in list(sample_tissue_map.keys())])
    process_pool.close()
    collected_recalls += list(chain(*model_recalls))
    print(clock() - startw, model)

recall_df = pd.DataFrame.from_dict(collected_recalls, orient='columns')
recall_df['GTEx tissue'] = [sample_tissue_map[s] for s in recall_df['EpiRR'].values]
recall_df.to_csv(out_dir + "CollectedRecall.tsv", sep='\t', header=True, index=False)



