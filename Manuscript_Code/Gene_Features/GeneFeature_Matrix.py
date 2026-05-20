import pandas as pd
import numpy as np
import os
from timeit import default_timer as clock
from multiprocessing import Pool
from pybedtools import BedTool
from sklearn.preprocessing import MinMaxScaler
from collections import Counter
import sys
sys.path.append("/home/dhecker/FuFis/src/")
import GeneFeature_Helpers


"""Collect multiple features on gene level into one DataFrame."""

n_cores = 40
annotation = "*/gencode.v38.annotation.gtf"
ihec_meta_file = "*IHEC_metadata_harmonization.v1.1.csv"
experiment_meta_file = "*230723_epiatlas_metadata.csv"
expression_tpm_file = "*genes_TPM.csv"
expression_counts_file = "*genes_expected_count_DESeq2_H3K27acFormatted.tsv"
tad_file = "*TADMap_scaffold_hs_hg38.bed"  # From TADMap https://cb.csail.mit.edu/tadmap/
encode_pecres_file = '*ENCODE_peCREsV3_GRCh38.bed'
out_path = "*GeneFeatureMat.txt.gz"

partition_file = "*partition0_*.csv"
test_samples = set(pd.read_table(partition_file.replace('*', 'test'), sep=',', header=0)['epirr_id_without_version'].values)
train_samples = set(pd.read_table(partition_file.replace('*', 'train'), sep=',', header=0)['epirr_id_without_version'].values)

# Start from the kept genes.
kept_genes_file = "*kept_genes.txt"
kept_genes = set([x.strip().split('\t')[1] for x in open(kept_genes_file).readlines()[1:]])
gene_set = kept_genes

gene_df = GeneFeature_Helpers.gene_feature_table(gtf_file=annotation, gene_set=gene_set, extend=500000)

# ------------------------------------------------------------------------------------------------
# Add expression-based features
# ------------------------------------------------------------------------------------------------
sample_meta = pd.read_csv(ihec_meta_file, header=0)
experiment_df = pd.read_csv(experiment_meta_file, header=0)
meta_df = pd.merge(experiment_df[['uuid', 'epirr_id_without_version']], sample_meta, on="epirr_id_without_version", how='inner')
uuid_cell_map = {x['uuid']: x['harmonized_sample_ontology_intermediate'] for x in meta_df.to_dict(orient='records')}
expression_df = pd.read_table(expression_tpm_file, sep=',', header=0, index_col='id_col')
expression_df.columns = [uuid_cell_map[c.split('.')[-1]] for c in expression_df.columns]
expression_mean_ct = expression_df.groupby(by=expression_df.columns, axis=1).mean()
ubi_expression = {k.split('.')[0]: val/expression_mean_ct.shape[1] for k, val in (expression_mean_ct >= 0.5).sum(axis=1).items()}
gene_df['Expression ubiquitousness'] = [ubi_expression[g] for g in gene_df.index]

# Get the fraction of samples in train/test with non-zero counts.
counts_df = pd.read_table(expression_counts_file, sep='\t', header=0, index_col='id_col')

counts_df = counts_df[[c for c in counts_df.columns if c in test_samples or c in train_samples]]
counts_std = counts_df.std(axis=1).to_dict()
gene_df['Expression counts std'] = [None if g not in counts_std else counts_std[g] for g in gene_df.index]

for sample_tag, sample_set in [['train', train_samples], ['test', test_samples]]:
    sample_counts = counts_df[[c for c in counts_df.columns if c in sample_set]]
    nonzero_samples = (sample_counts > 0).sum(axis=1)
    gene_df['Fraction nonzero in '+sample_tag] = nonzero_samples / len(sample_set)
    sample_std = sample_counts.std(axis=1).to_dict()
    gene_df['Expression counts std in ' + sample_tag] = [None if g not in sample_std else sample_std[g] for g in gene_df.index]
    if sample_tag == 'train':
        scaler = MinMaxScaler()
        scaled_counts = pd.DataFrame(scaler.fit_transform(np.log2(sample_counts + 1).T), columns=sample_counts.index).T
        scaled_counts.columns = sample_counts.columns
        scaled_counts_std = scaled_counts.std(axis=1).to_dict()
        gene_df['Expression scaled counts std in ' + sample_tag] = [None if g not in scaled_counts_std else scaled_counts_std[g] for g in gene_df.index]

# ------------------------------------------------------------------------------------------------
# Distance to the closest TAD border
# ------------------------------------------------------------------------------------------------
# Get the distance to the closest TAD border. First need to rewrite the TADs from bin-bin to two 1-bp positions.
tad_borders = []
for tad in [x.strip().split('\t') for x in open(tad_file).readlines()]:
    tad_borders.append('\t'.join([tad[0], tad[1], str(int(tad[1])+1)]))
    tad_borders.append('\t'.join([tad[0], str(int(tad[2]) - 1), tad[2]]))
tad_borders_bed = BedTool('\n'.join(tad_borders), from_string=True).sort()
open(tad_file.replace('.bed', '_borders.bed'), 'w').write(str(tad_borders_bed))
tads = BedTool(tad_file)
print(np.mean([x.length for x in tads]), np.min([x.length for x in tads]), np.max([x.length for x in tads]))

tss_bed = GeneFeature_Helpers.gene_window_bed(gtf_file=annotation, extend=0, gene_set=gene_set, tss_type='5').sort()
closest_tad = tss_bed.closest(tad_borders_bed, d=True)
closest_tad_dict = {x.fields[3]: x.fields[-1] for x in closest_tad}
gene_df['Distance TAD border'] = [None if g not in closest_tad_dict else closest_tad_dict[g] for g in gene_df.index]

# ------------------------------------------------------------------------------------------------
# GTF-annotated features
# ------------------------------------------------------------------------------------------------
# Add the 5' TSS as column, and also all other annotated TSS as separate column.
gene_5tss = GeneFeature_Helpers.gene_window_bed(annotation, extend=200, gene_set=gene_set, tss_type='5', dict_only=True)
gene_alltss = GeneFeature_Helpers.gene_window_bed(annotation, extend=200, gene_set=gene_set, tss_type='all', dict_only=True)
gene_df['5prime TSS'] = [None if g not in gene_5tss else str(next(iter(gene_5tss[g]['tss']))) for g in gene_df.index]
gene_df['alternative TSS'] = [None if g not in gene_alltss else ','.join([str(t) for t in (gene_alltss[g]['tss'] - gene_5tss[g]['tss'])])
                              for g in gene_df.index]
gene_df['Strand'] = [None if g not in gene_5tss else gene_5tss[g]['strand'] for g in gene_df.index]

# Get the longest 3' and 5' UTR.
gene_utrs = GeneFeature_Helpers.gene_feature_bed(annotation, feature='UTR', gene_set=gene_set, dict_only=True, merge=False)
utr_map = {g: {'5': 0, '3': 0} for g in gene_set}
for gene in gene_set:
    if gene in gene_utrs:
        for utr in gene_utrs[gene]:
            if gene_alltss[gene]['strand'] == '+':
                utr_start = utr[1]
            else:  # On - strand have to check with the positional end of the UTR.
                utr_start = utr[2]
            if int(utr_start) in gene_alltss[gene]['tss']:
                utr_map[gene]['5'] = max([utr_map[gene]['5'], int(utr[2]) - int(utr[1])])
            else:
                utr_map[gene]['3'] = max([utr_map[gene]['3'], int(utr[2]) - int(utr[1])])
gene_df['longest 5prime UTR'] = [utr_map[g]['5'] for g in gene_df.index]
gene_df['longest 3prime UTR'] = [utr_map[g]['3'] for g in gene_df.index]

# ------------------------------------------------------------------------------------------------
# Number of CREs
# ------------------------------------------------------------------------------------------------
# Add the number of peCREs in the window.
gene_windows = GeneFeature_Helpers.gene_window_bed(annotation, extend=500000, gene_set=gene_set, tss_type='5')
encode_inter = gene_windows.intersect(encode_pecres_file, wa=True)
encode_hits = Counter([x.fields[3] for x in encode_inter])
gene_df['#ENCODE peCREs'] = [None if g not in encode_hits else encode_hits[g] for g in gene_df.index]

# ------------------------------------------------------------------------------------------------
# H3K27ac input features
# ------------------------------------------------------------------------------------------------
# Get the number of non-zero entries in the input matrices for the different feature setups.
cre_input_pattern = '*/IHEC_Activity_1MB_hg38/*.txt.gz'
binned_input_pattern = '*/*_w1MB_100bs_BinnedActivity.txt.gz'
stitchit_train_pattern = '*/Segmentation_*_Pearson_10.txt'
stitchit_test_pattern = '*/*_stitchit_feature_counts.txt.gz'
stitchit_losttest_pattern = '*/*_stitchit_feature_counts.txt.gz'


def input_count_zeros(args):
    """Read an input matrix and count the non-zero entries."""
    input_file, gene = args
    input_matrix = pd.read_table(input_file, sep='\t').set_index("Sample").iloc[:, :-1]
    train_matrix = input_matrix.loc[list(train_samples)]
    test_matrix = input_matrix.loc[list(test_samples)]
    train_nonzeros = (train_matrix != 0).sum().sum()
    test_nonzeros = (test_matrix != 0).sum().sum()
    return [gene, train_nonzeros, train_nonzeros / train_matrix.size, test_nonzeros, test_nonzeros / test_matrix.size]


def praise_stitchit(args):
    """STITCHIT being STITCHIT."""
    input_file, gene = args
    input_matrix = pd.read_table(input_file, sep='\t', index_col=0).iloc[:, :-1]
    nonzeros = (input_matrix != 0).sum().sum()
    return [gene, nonzeros, nonzeros / input_matrix.size]

def praise_stitchit_test(args):
    """STITCHIT being extra STITCHIT. For the test we have a separate file with an addition sample."""
    input_file, gene, separate_input_file = args
    input_matrix = pd.read_table(input_file, sep='\t', index_col=0).iloc[:, :-1]
    nonzeros = (input_matrix != 0).sum().sum()
    separate_matrix = pd.read_table(separate_input_file, sep='\t', index_col=0).iloc[:, :-1]
    separate_nonzeros = (separate_matrix != 0).sum().sum()
    return [gene, nonzeros + separate_nonzeros, (nonzeros + separate_nonzeros) / (input_matrix.size + separate_matrix.size)]


for feature_setup, input_folder in [['CRE', cre_input_pattern],
                                    ['Binned', binned_input_pattern]]:
    input_files = GeneFeature_Helpers.fn_patternmatch(input_folder)
    startw = clock()
    process_pool = Pool(processes=n_cores)
    file_nonzeros = process_pool.map(input_count_zeros, [[val, g] for val, g in input_files.items() if g in gene_set])
    process_pool.close()
    print(clock() - startw, feature_setup)
    nonzero_vals = pd.DataFrame(file_nonzeros, columns=['Ensembl ID', feature_setup + ': #Non-zero input entries in train',
                                                        feature_setup + ': Fraction non-zero input entries in train', 
                                                        feature_setup + ': #Non-zero input entries in test', 
                                                        feature_setup + ': Fraction non-zero input entries in test']).set_index("Ensembl ID")
    gene_df = gene_df.merge(nonzero_vals, left_index=True, right_index=True, how='left')

# Special treatment for STITCHIT.
for part_tag, input_folder in [['train', stitchit_train_pattern],
                                ['test', stitchit_test_pattern]]:
    input_files = GeneFeature_Helpers.fn_patternmatch(input_folder)
    startw = clock()
    process_pool = Pool(processes=n_cores)
    if part_tag == 'train':
        file_nonzeros = process_pool.map(praise_stitchit, [[val, g] for val, g in input_files.items() if g in gene_set])
    elif part_tag == 'test':
        separate_input_files = {g: val for val, g in GeneFeature_Helpers.fn_patternmatch(stitchit_losttest_pattern).items()}
        file_nonzeros = process_pool.map(praise_stitchit_test, [[val, g, separate_input_files[g]] for val, g in input_files.items() if g in gene_set])
    process_pool.close()
    print(clock() - startw, 'STITCHIT ' + part_tag)
    nonzero_vals = pd.DataFrame(file_nonzeros, columns=['Ensembl ID', 'STITCHIT: #Non-zero input entries in ' + part_tag,
                                                        'STITCHIT: Fraction non-zero input entries in ' + part_tag]).set_index("Ensembl ID")
    gene_df = gene_df.merge(nonzero_vals, left_index=True, right_index=True, how='left')


# ------------------------------------------------------------------------------------------------
# Put out
# ------------------------------------------------------------------------------------------------
gene_df.columns = [c.replace(' ', '_') for c in gene_df.columns]
gene_df.to_csv(out_path, sep='\t', header=True, index=True)

