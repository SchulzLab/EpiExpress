import pandas as pd
import numpy as np
from itertools import chain
import scipy.stats
import plotnine as pn
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

        # A version for the supplements comparing with normG, having all 10 setups as boxplots.
    normG_colour_list = ['#8accff', '#8afffb', '#0a63a5', '#0aa5a0', '#b0b0b0', '#666666', '#dbc430', '#c1db30', '#ab0303', '#ab5103']
    # And remove the underscores from the GTEx tissues.
    to_plot_df['GTEx tissue'] = [t.replace('_', ' ') for t in to_plot_df['GTEx tissue']]
    norm_models = list(chain(*[[rename_scores(m).replace(' normG', ''), rename_scores(m)] for m in models]))
    to_plot_df['Significant'] = [str(x) for x in to_plot_df['NOM p-val'] <= 0.05]
    to_plot_df['forcex'] = to_plot_df['GTEx tissue'] + '#' + to_plot_df['model']
    forcex_order = sorted(set(to_plot_df['forcex'].values), key=lambda x: (tissue_order.index(x.split('#')[0]), norm_models.index(x.split('#')[1])))
    to_plot_df["forcex"] = pd.Categorical(to_plot_df['forcex'], categories=forcex_order, ordered=True)
    to_plot_df['model'] = pd.Categorical(to_plot_df['model'], categories=norm_models, ordered=True)

    # Use plotnine for this plot to allow different shapes for the jitter.
    custom_theme = pn.theme_matplotlib() + pn.theme(text=pn.element_text(size=14, colour='black'),
                         panel_background=pn.element_rect(fill='white', alpha=0.2))
    plotnine = (
                pn.ggplot(to_plot_df, pn.aes(x="forcex", y="NES", fill='model'))
                + pn.geom_boxplot(outlier_shape='')
                + pn.scale_fill_manual(values=normG_colour_list*len(tissue_order))
                + pn.geom_jitter(pn.aes(shape='Significant'), position=pn.position_jitterdodge(jitter_width=0.5, dodge_width=0.7))
                + pn.scale_shape_manual(values=['^', 'o'])
                + pn.guides(fill=pn.guide_legend(override_aes={'shape': 's', 'size': 12, 'color': 'white'}), linetype='')
                + pn.ggtitle(eqtl)
                + custom_theme
                )
    plotnine.save(out_dir+eqtl+"_BothNormModel_PlotNine.pdf", height=6, width=12)

    # And a heatmap of row versus column in how many samples the NES is greater, including a p-value.
    epirr_order = sorted(set(to_plot_df['EpiRR']))
    model_nes_vecs = {m: to_plot_df[to_plot_df['model'] == m].set_index("EpiRR").loc[epirr_order]['NES'].values for m in norm_models}
    epirr_comparison = np.zeros([len(norm_models), len(norm_models)])
    sigs = np.full([len(norm_models), len(norm_models)], fill_value='', dtype='object')
    for r_i, row_m in enumerate(norm_models):
        for c_i, col_m in enumerate(norm_models):
            if row_m == col_m:
                continue
            epirr_comparison[r_i][c_i] = (model_nes_vecs[row_m] > model_nes_vecs[col_m]).sum() / len(epirr_order)
            wilcox_test = scipy.stats.wilcoxon(x=model_nes_vecs[row_m], y=model_nes_vecs[col_m])
            if wilcox_test.pvalue <= 0.05:
                sigs[r_i][c_i] = '*'
        epirr_comparison[r_i][r_i] = np.nan
    epirr_comparison_df = pd.DataFrame(epirr_comparison, index=norm_models, columns=norm_models)
    sigs_df = pd.DataFrame(sigs, index=norm_models, columns=['sig-'+x for x in norm_models])
    joint_comparison = pd.concat([epirr_comparison_df, sigs_df], axis=1)
    cmap_cols = {0: {'cols': list(epirr_comparison_df.columns), 'cmap': 'bwr', 'centre': 0.5, 'cbar_label': 'NES of row > NES of column [%]'}}
    eQTL_Helpers.heatmap_cols(joint_comparison, cmap_cols, plot_out=out_dir+"PairwiseBestAcrossScores_"+eqtl, formats=['pdf'],
                        annot_cols={c: 'sig-'+c for c in epirr_comparison_df.columns}, annot_s=10, x_size=6, y_size=5, title=eqtl, ticksize=14, heat_ticksize=16, x_rotation=90)
