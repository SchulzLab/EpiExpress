import pandas as pd
from itertools import chain
import sys
#sys.path.append("*ISP/")  # Append the directory with the scripts if it's not the working directory.
import eQTL_Helpers
from CRISPRi_Helpers import rename_scores

out_dir = "*/"
recall_file = "*/eQTLs_CollectedRecall.tsv"
recall_df = pd.read_table(recall_file, sep='\t', header=0)
tissue_order = ['Brain_Cortex', 'Nerve_Tibial', 'Thyroid', 'Whole_Blood']

models = ['ENCODE_RF_ism_all_reverse', 'ENCODE_MLP_ism_all_reverse', 'STITCHIT_ism_all_eqtl_backscaled', 'binnedRF_os_ism_all_expand_len10', 'CNN_bestwarm_os_ism_all_expand_len10']
colour_dict = {'models': {'CRE-RF': '#8accff',
                'CRE-MLP': '#0a63a5',
                'STITCHIT': '#b0b0b0',
                'Binned-RF': '#dbc430',
                'Binned-CNN': '#ab0303'}}
colour_list = [colour_dict['models'][m] for m in ['CRE-RF', 'CRE-MLP', 'STITCHIT', 'Binned-RF', 'Binned-CNN']]


for eqtl in ['CAVIAR', 'CaVEMaN', 'DAP-G']:
    to_plot_df = recall_df[recall_df['eqtl'] == eqtl]
    to_plot_df['model'] = [rename_scores(c).replace(' normG' ,'')+' normG'*('normG' in c) for c in to_plot_df['model'].values]
    # Get the better working one between normG and non-normG.
    to_keep = []
    for m in [rename_scores(m).replace(' normG', '') for m in models]:
        mean_nonNormG = to_plot_df[to_plot_df['model'] == m]['NES'].mean()
        mean_NormG = to_plot_df[to_plot_df['model'] == m+' normG']['NES'].mean()
        print(m, mean_nonNormG, mean_NormG)
        if mean_nonNormG > mean_NormG:
            to_keep.append(m)
        else:
            to_keep.append(m + ' normG')
    # Still remove the normG suffix.
    better_df = to_plot_df[to_plot_df['model'].isin(to_keep)]
    better_df['model'] = [m.replace(' normG', '') for m in better_df['model'].values]
    # And remove the underscores from the GTEx tissues.
    better_df['GTEx tissue'] = [t.replace('_', ' ') for t in better_df['GTEx tissue']]
    eQTL_Helpers.basic_violin(plot_df=better_df, y_col='NES', x_col='GTEx tissue', x_order=[t.replace('_', ' ') for t in tissue_order],
                                hue_col='model', hue_order=[m.replace(' normG', '') for m in to_keep], title=eqtl,
                                output_path=out_dir+eqtl+"_BestWithinModel", legend=False,
                                palette=colour_list, xsize=12, ysize=10, boxplot=True, jitter=True, jitter_colour=['black'], jitter_size=4, font_s=14)

    # A version for the supplements comparing with normG, having all 10 setups as boxplots.
    normG_colour_list = ['#8accff', '#8afffb', '#0a63a5', '#0aa5a0', '#b0b0b0', '#666666', '#dbc430', '#c1db30', '#ab0303', '#ab5103']
    # And remove the underscores from the GTEx tissues.
    to_plot_df['GTEx tissue'] = [t.replace('_', ' ') for t in to_plot_df['GTEx tissue']]
    norm_models = list(chain(*[[rename_scores(m).replace(' normG', ''), rename_scores(m)] for m in models]))
    eQTL_Helpers.basic_violin(plot_df=to_plot_df, y_col='NES', x_col='GTEx tissue', x_order=[t.replace('_', ' ') for t in tissue_order],
                                hue_col='model', hue_order=norm_models, title=eqtl,
                                output_path=out_dir+eqtl+"_BothNormModel", legend=True,
                                palette=normG_colour_list, xsize=20, ysize=10, boxplot=True, jitter=True, jitter_colour=['black'], jitter_size=4, font_s=14)

