import os
import gzip
import pandas as pd
import numpy as np
import copy
import sys
from timeit import default_timer as clock
#sys.path.append("*ISP/")  # Append the directory with the scripts if it's not the working directory.
import CRISPRi_Helpers

input_folder = "*/ValidateInteractions/"  # Directory with a folder written for each model, with the in-silico mutagenesis simulation.
out_dir = "*/"
joint_crispri_file = "*/JointValidatedInteractions_hg38_500kb_Gschwind.txt"
kept_genes_file = "*/kept_genes.txt"
kept_genes = set([x.strip().split('\t')[1] for x in open(kept_genes_file).readlines()[1:]])

k562_sample = 'IHECRE00001887'
models = ['ENCODE_MLP_ism_overlap', 'ENCODE_MLP_ism_all_reverse', 
          'ENCODE_RF_ism_overlap', 'ENCODE_RF_ism_all_reverse', 
          'STITCHIT_ism_overlap_backscaled_val_regions_exclusive', 'STITCHIT_ism_all_backscaled',
          'CNN_bestwarm_os_ism_overlap_expand_len10', 'binnedRF_os_ism_overlap_expand_len10']


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


inter_files = {}  # Add the effect of the screens themselves.

for tag in models:
    tag_out = out_dir + '/' + tag + "_CRISPRi_formatted.txt.gz"
    print(tag)
    start = clock()
    match_folder = [input_folder + '/' + x for x in os.listdir(input_folder) if x == tag][0] + '/' + k562_sample
    df_list = []  # Collecting a list is faster than concatenating a df per gene.
    for file in [match_folder+'/'+x for x in os.listdir(match_folder) if x.endswith('.bed') or x.endswith('.txt') or x.endswith('.gz') or x.endswith('.zip')]:
        this_gene = file.split('/')[-1].split('.')[0].split('_')[0]

        if this_gene not in kept_genes:
            continue

        header_rows = []
        if file.endswith('.gz'):
            file_opener = gzip.open(file, 'rt')
        else:
            file_opener = open(file)
        with file_opener as read_header:
            for row in read_header:
                if not row.startswith('#'):
                    break
                header_rows.append(row.strip())

        this_df = pd.read_table(file, header=0, sep='\t', skiprows=len(header_rows))
        this_df.columns = [match_colname(c) for c in this_df.columns]
        this_df['chr'] = [str(x).replace('chr', '') for x in this_df['chr'].values]
        this_df['Ensembl ID'] = this_gene
        this_df = this_df.astype({'start': 'int', 'end': 'int', "ISM": 'float64', "predicted": 'float64'})
        if this_df['ISM'].min() < 0 or this_df['predicted'].min() < 0:
            this_df['ISM'] = this_df['ISM'] + abs(min([this_df['ISM'].min(), this_df['predicted'].min()]))
            this_df['predicted'] = this_df['predicted'] + abs(min([this_df['ISM'].min(), this_df['predicted'].min()]))

        this_df['log fraction change'] = np.log2((this_df['ISM']+1) / (this_df['predicted']+1))
        this_df['abs log fraction change'] = this_df['log fraction change'].abs()
        if 'ism_all' in tag.lower() or 'ismall' in tag.lower():  # Divide the changes by the sum per gene.
            for c in ['log fraction change', 'abs log fraction change']:
                if this_df[c].sum() != 0:
                    this_df[c + " normG"] = this_df[c] / this_df[c].sum()
        df_list.append(this_df)
    
    if df_list:
        collected_df = pd.concat(df_list)
        collected_df.to_csv(tag_out, header=True, index=False, sep='\t', na_rep='0')
    else:
        print("No remaining genes", tag)
        continue
    print(clock() - start)

    inter_files[tag] = {"file": tag_out,
                        "score": ['abs log fraction change'],
                        'genome_version': "hg38",
                        'inter_f': 1e-09 if 'all' in tag else 1,
                        'inter_F': 1e-09 if 'all' in tag else 1}
    if 'ism_all' in tag.lower() or 'ismall' in tag.lower():
        inter_files[tag]['score'] = ["abs log fraction change normG"]

# List of groups of files to run on, either all, or the individual ones.
validation_dict = CRISPRi_Helpers.get_crispris_from_joint(inter_files, screens=['gschwind'], distance=False, meta_out=out_dir, meta_tag='',
                                                            joint_file=joint_crispri_file, all_required=True, tmp_dir='')

val_df = validation_dict['gschwind']['df']
val_df.to_csv(out_dir + "CRISPRiTable_Gschwind.tsv", sep='\t', header=True, index=False)

