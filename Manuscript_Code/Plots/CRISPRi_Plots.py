import copy
import pandas as pd
import sys
#sys.path.append("*ISP/")  # Append the directory with the scripts if it's not the working directory.
import CRISPRi_Helpers

"""Plots precision recall curves and subsequent barplots based on tables created with CRISPRi_Validation.py."""

tag = ''
plot_out = "*/" + tag + "_"
crispri_file = '*/CRISPRiTable_Gschwind.tsv'  # Written with CRISPRi_Validation.py.
model_order = ['CRE-RF', 'CRE-RF normG', 'CRE-MLP', 'CRE-MLP normG', 'STITCHIT', 'STITCHIT normG', 'Binned-RF', 'Binned-CNN']

colour_dict = {'models': {'CRE-RF': '#8accff',
                'CRE-MLP': '#0a63a5',
                'STITCHIT': '#b0b0b0',
                'Binned-RF': '#dbc430',
                'Binned-CNN': '#ab0303'}}
# Manually defining the order because of all the normG suffix handling.
colour_list = [colour_dict['models'][m] for m in ['CRE-RF', 'CRE-MLP', 'STITCHIT', 'Binned-RF', 'Binned-CNN']]


# Use the df from CRISPRi_Validation to plot PR curves for all setups. 
# The PR plots are admittedly ugly, we only need them for the AUPRC values.
crispri_df = pd.read_table(crispri_file, sep='\t', header=0)
auc_out, pr_dict = CRISPRi_Helpers.pre_rec_plotter(crispri_df, sig_col='Significant', legend_out=2.8,
                                                           steps=10000, mode='pr', plot_cols=crispri_df.columns[list(crispri_df.columns).index('Significant')+1:], 
                                                           output_path=plot_out)

auprc_df = pd.DataFrame(auc_out, columns=['model', 'AUPRC'])
auprc_df['model'] = [CRISPRi_Helpers.rename_scores(x) for x in auprc_df['model'].values]

# Pick better version for each model between normG and non normG.
to_keep = [] 
for m in ['CRE-RF', 'CRE-MLP', 'STITCHIT']:  # Manual as not all models have it.
    mean_nonNormG = auprc_df[auprc_df['model'] == m]['AUPRC'].values[0]
    mean_NormG = auprc_df[auprc_df['model'] == m+' normG']['AUPRC'].values[0]
    print(m, mean_nonNormG, mean_NormG)
    if mean_nonNormG > mean_NormG:
        to_keep.append(m)
    else:
        to_keep.append(m + ' normG')
to_keep += ['Binned-RF', 'Binned-CNN']  # For those we don't have normG.
# Make a separate df so that we can remove the normG from the names.
to_keep_df = auprc_df[auprc_df['model'].isin(to_keep)]
to_keep_df['model'] = [m.replace(' normG', '') for m in to_keep_df['model'].values]
to_keep_df['screen'] = 'Gschwind'  # To enable plotting with hue for the same x-position.
CRISPRi_Helpers.basic_bars(to_keep_df, x_col='screen', y_col='AUPRC', hue_col='model', hue_order=[m.replace(' normG', '') for m in to_keep], title="",
                        output_path=plot_out + "Gschwind_BestofAll",
                        x_size=7, y_size=5, rotation=None, palette=colour_list, legend_out=1.45)

# And one for the supplements where we show always normG and non-normG.
supp_df = copy.deepcopy(auprc_df)
supp_df['main model'] = [x.replace(' normG', '') for x in supp_df['model']]
supp_df['normalized'] = ['normalized per gene' if 'normG' in x else 'non-normalized' for x in supp_df['model']]
CRISPRi_Helpers.basic_bars(supp_df, x_col='main model', x_order=[x for x in model_order if 'normG' not in x], y_col='AUPRC', hue_col='normalized', 
                                      hue_order=['non-normalized', 'normalized per gene'], title="",
                        output_path=plot_out + "Gschwind_SuppsBothNorm",
                        x_size=8, y_size=5, rotation=None, palette=['#5B4E3F', '#BAA392'], legend_out=1.5)

# And a PR curve with only the better version of normG or non-normG.
sub_crispri_df = copy.deepcopy(crispri_df)
sub_crispri_df.columns = [CRISPRi_Helpers.rename_scores(c) for c in sub_crispri_df.columns]
# And again keeping the best of normG/non-normG and removing the suffix.
sub_crispri_df = sub_crispri_df[[c for c in sub_crispri_df.columns if c in to_keep+['chr', 'start', 'end', 'gene', 'effect', 'p-adj', 'Significant']]]
sub_crispri_df.columns = [c.replace(' normG', '') for c in sub_crispri_df.columns]
auc_out, pr_dict = CRISPRi_Helpers.pre_rec_plotter(sub_crispri_df, sig_col='Significant', colours=colour_list, legend_out=1.5, lines_solid=True,
                                                            steps=10000, mode='pr', plot_cols=['CRE-RF', 'CRE-MLP', 'STITCHIT', 'Binned-CNN', 'Binned-RF'], 
                                                            output_path=plot_out + "PRC_BestofAll_Gschwind")

