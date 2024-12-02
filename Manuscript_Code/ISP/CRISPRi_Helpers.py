import pandas as pd
import pybedtools
import gzip
import copy
import numpy as np
import json
from datetime import datetime
from matplotlib import pyplot as plt
import sklearn.metrics as metrics
import numpy as np
import gzip
import seaborn as sns

"""
Helper functions for CRISPRi validation. Take a list of enhancer-gene interactions with a score and intersects 
them with experimentally validated interactions from CRISPRi screens. Returns a df, which is then used 
to plot PR curves.
"""

tol_vibrant = ['#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB']


def rename_scores(score_name, keep_abs=False):
    """This beautiful piece of code translates the model names into something readable."""
    if 'abs' in score_name and keep_abs:
        abs_suffix = ' abs'
    else:
        abs_suffix = ''
    if 'CNN' in score_name:
        if 'all' in score_name:
            return 'Binned-CNN normG' + abs_suffix
        else:
            return 'Binned-CNN' + abs_suffix
    elif score_name.split(' ')[0] == 'ENCODE_MLP_ism_all_reverse':
        return 'CRE-MLP normG' + abs_suffix
    elif score_name.split(' ')[0] == 'ENCODE_MLP_ism_overlap':
        return 'CRE-MLP' + abs_suffix
    elif score_name.split(' ')[0] == 'ENCODE_MLP_all_reverse_NEW_ISM':
        return 'CorrISM CRE-MLP normG' + abs_suffix
    elif score_name.split(' ')[0] == 'ENCODE_RF_all_reverse_NEW_ISM':
        return 'CorrISM CRE-RF normG' + abs_suffix
    elif score_name.split(' ')[0] == 'ENCODE_RF_ism_all_reverse':
        return 'CRE-RF normG' + abs_suffix
    elif score_name.split(' ')[0] == 'ENCODE_RF_ism_overlap':
        return 'CRE-RF' + abs_suffix
    elif 'STITCHIT_ism_all' in score_name:
        return 'STITCHIT normG' + abs_suffix
    elif 'STITCHIT_ism_overlap' in score_name:
        return 'STITCHIT' + abs_suffix
    elif 'binnedRF' in score_name:
        if 'all' in score_name:
            return 'Binned-RF normG' + abs_suffix
        else:
            return 'Binned-RF' + abs_suffix
    elif ' effect' in score_name:  # Space important to not match the intersection dfs.
        return 'CRISPRi measured change'
    else:
        print("WARNING no matching name", score_name)
        return score_name
    

def get_crispris_from_joint(to_validate, joint_file, screens=['gschwind'],
                            all_required=True, distance=False, meta_out='', meta_tag='', tmp_dir=''):
    """
    Takes a dictionary with files of interactions, which should be bed-styled with a header holding the column names
    for the gene and a score column (e.g. ABC-Score). Will return a df for each
    CRISPRi screen and the score of the intersecting interactions from the file. If multiple enhancer intersect a
    validated interaction the sum is taken.
    :param to_validate: {label: {
                       'file': file_path,
                       'score': column name which score to take, can be an iterable with multiple columns,
                       'gene': column name to find the gene,
                       'avg_cols': a list of columns to add to the df, will be averaged instead of summed.}}
    :param screens: Names screens to validate on, can be used to reduce which CRISPRi data sets to include.
    :param all_required: If all score columns must have an intersection with the validated interactions, if False
    sets the score to 0 of the ones missing a score.
    :param meta_out: Directory where to write the JSON file to with the metadata of the function call.
    :param meta_tag: Prefix for the metadata file. The filename will also include the distance and all_required info.
    :return: A dictionary for each screen with a df with the validated interactions and the score from the given
    interaction file, also has entries for which columns were retrieved.
    """
    if tmp_dir:
        pybedtools.helpers.set_tempdir(tmp_dir)

    # Store where to find the files and which columns to grab. The order must be the same, it relies on the indices.
    validation_info = pd.read_table(joint_file, sep='\t', header=0)

    # We make a copy to be able to streamline some entries.
    score_files = copy.deepcopy(to_validate)
    validation_dfs = {v: {'df': None, 'score_cols': None, 'avg_cols': None} for v in screens}
    score_cols = []
    average_cols = []  # The columns from 'avg_cols' will be averaged instead of summed like the 'score' column.
    # Add default values to the score_files dictionary.
    for score in score_files:
        if 'gene' not in score_files[score]:
            score_files[score]['gene'] = 'Ensembl ID'
        if 'score' not in score_files[score]:
            score_files[score]['score'] = 'ABC-Score'
        if 'genome_version' not in score_files[score]:
            score_files[score]['genome_version'] = 'hg19'
        if 'aggregate_mode' not in score_files[score]:
            score_files[score]['aggregate_mode'] = 'sum'
        if 'inter_f' not in score_files[score]:
            score_files[score]['inter_f'] = 1e-09
        if 'inter_F' not in score_files[score]:
            score_files[score]['inter_F'] = 1e-09
        if 'inter_e' not in score_files[score]:
            score_files[score]['inter_e'] = False
    for score, vals in score_files.items():
        if not type(vals['score']) == list:
            score_files[score]['score'] = [score_files[score]['score']]
        score_cols += [score + ' ' + c for c in vals['score']]
        if 'avg_cols' in vals:
            average_cols += [score + ' ' + c for c in vals['avg_cols']]

    # After filling the defaults write a metadata file to trace the function call.
    json_dict = {'Date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                 "to_validate": to_validate,
                 'joint_file': joint_file,
                 'screens': screens,
                 'all_required': all_required,
                 'distance': distance,
                 'meta_out': meta_out,
                 'meta_tag': meta_tag}
    json.dump(json_dict, open(meta_out+'/'+meta_tag+"_reqAll"+str(all_required)+"_Dist"+str(distance)+"_Metadata.json", 'w'))

    for mode in screens:
        print(mode)
        screen_df = validation_info[validation_info['Screen'] == mode]
        if distance:
            screen_df = screen_df[screen_df['Distance'] <= distance]
        validation_scores = {
            '\t'.join([str(x[c]) for c in ['#chr', 'start', 'end', 'Ensembl ID']]): {'effect': float(x['effect']),
                                                                         'p-adj': float(x['p-adj']),
                                                                         'Significant': x['Significant']} for x in
            screen_df.to_dict(orient='records')}
        for entry in validation_scores:
            validation_scores[entry].update({c: [] for c in score_cols + average_cols})

        validation_bed = pybedtools.BedTool('\n'.join(validation_scores.keys()), from_string=True)

        def match_scores(abc_file, score_tag, match_cols=['ABC-Score'], gene_col='Ensembl ID', genome_version='hg19',
                         inter_f=1e-09, inter_F=1e-09, inter_e=False):
            """Takes an ABC-scoring file and validated interactions to fetch the matching scores."""
            if type(match_cols) != list:
                match_cols = list(values['score'])

            if abc_file.endswith('.gz'):
                with gzip.open(abc_file, 'rt') as abc_in:
                    abc_head = {x: i for i, x in enumerate(abc_in.readline().replace('#', '').strip().split('\t'))}
                    abc_rows = abc_in.readlines()
            else:
                abc_head = {x: i for i, x in
                            enumerate(open(abc_file).readline().replace('#', '').strip().split('\t'))}
                abc_rows = open(abc_file).readlines()[1:]
            abc_rows = [x.strip().split('\t') for x in abc_rows]
            abc_rows = [[str(x[0]), str(int(float(x[1]))), str(int(float(x[2])))] + x[3:] for x in abc_rows]

            # First find the enhancers that intersect a validated region at all.
            chr_prefix = "chr" if not abc_rows[0][0].startswith('chr') else ''

            peaks_bed = pybedtools.BedTool('\n'.join(set([chr_prefix + '\t'.join(x[:3]) for x in abc_rows])),
                                           from_string=True)
            hit_peaks = set(
                [str(x).strip() for x in peaks_bed.intersect(validation_bed, u=True, f=inter_f, F=inter_F, e=inter_e)])
            peak_preturb_map = {x: set() for x in hit_peaks}
            for inter in str(
                    peaks_bed.intersect(validation_bed, wo=True, f=inter_f, F=inter_F, e=inter_e)).strip().split('\n'):
                if '\t'.join(inter.split('\t')[:3]):  # For whatever reason there can be an empty intersection.
                    peak_preturb_map['\t'.join(inter.split('\t')[:3])].add('\t'.join(inter.split('\t')[3:6]))

            # Now go through the queried interactions again to find peaks that intersect a perturbed site.
            for inter in abc_rows:
                peak_str = chr_prefix + '\t'.join(inter[:3])
                if peak_str in hit_peaks:
                    inter_gene = inter[abc_head[gene_col]].split('.')[0]
                    for target in peak_preturb_map[peak_str]:
                        if target + '\t' + inter_gene in validation_scores:
                            for match_col in match_cols:
                                validation_scores[target + '\t' + inter_gene][score_tag + ' ' + match_col].append(
                                    float(inter[abc_head[match_col]]))

        score_cols_agg = {}
        for s_tag, values in score_files.items():
            print(s_tag)
            match_scores(values['file'], score_tag=s_tag, match_cols=values['score'], gene_col=values['gene'],
                         genome_version=values['genome_version'], inter_f=values['inter_f'],
                         inter_F=values['inter_F'], inter_e=values['inter_e'])
            for s_col in values['score']:
                score_cols_agg[s_tag + ' ' + s_col] = score_files[s_tag]['aggregate_mode']
            if 'avg_cols' in values:
                match_scores(values['file'], score_tag=s_tag, match_cols=values['avg_cols'], gene_col=values['gene'],
                             genome_version=values['genome_version'], inter_f=values['inter_f'],
                             inter_F=values['inter_F'], inter_e=values['inter_e'])

        compare_cols = ['chr', 'start', 'end', 'gene'] + list(list(validation_scores.values())[0].keys())
        col_types = ['str', 'int64', 'int64', 'str', 'float64', 'float64', bool] + ['float64'] * len(
            score_cols + average_cols)
        all_hits = []
        not_all = []
        for k, vals in validation_scores.items():
            if (all_required and np.all([1 if vals[score] else 0 for score in score_cols])) or \
                    (not all_required and np.any([1 if vals[score] else 0 for score in score_cols])):
                all_hits.append(k.split('\t') + [vals['effect'], vals['p-adj'], vals['Significant']] +
                                [(sum(vals[col]) if score_cols_agg[col] == 'sum' else np.mean(vals[col])) if vals[
                                    col] else None for col in score_cols] + 
                                [np.mean(vals[col]) if vals[col] else 0 for col in average_cols])
            else:
                not_all.append([k, vals])
        print("Interactions with only partial columns", len(not_all))

        compare_abcpp_df = pd.DataFrame(all_hits, columns=compare_cols)
        compare_abcpp_df = compare_abcpp_df.astype({c: col_types[i] for i, c in enumerate(compare_cols)})

        validation_dfs[mode]['df'] = compare_abcpp_df
        validation_dfs[mode]['score_cols'] = score_cols
        validation_dfs[mode]['avg_cols'] = average_cols

        pybedtools.helpers.cleanup(verbose=False, remove_all=False)

    return validation_dfs



def pre_rec_fetcher(binning, ele_gene_pairs, true_enhancers, score_col='ABC Score', lower_limit=0, upper_limit=1,
                    mark_thresh=None, mode='pr'):
    """Calculates precision and recall for a given range and number of steps. Is used by the function below to
    build Precision-Recall curves."""
    def fetch_prerec(thresholds):
        rec_prec_list = []
        for thresh in thresholds:
            true_positives = true_enhancers[pd.to_numeric(true_enhancers[score_col]) >= thresh]
            # How many of the positives did we find.
            recall = true_positives.shape[0] / true_enhancers.shape[0]
            if ele_gene_pairs.dtypes['Significant'] == bool:
                false_positives = ele_gene_pairs[(ele_gene_pairs['Significant'] == 0) & (pd.to_numeric(ele_gene_pairs[score_col]) >= thresh)]
                true_negatives = len(ele_gene_pairs[ele_gene_pairs['Significant'] == 0])
            else:
                false_positives = ele_gene_pairs[(ele_gene_pairs['Significant'] == 'False') & (pd.to_numeric(ele_gene_pairs[score_col]) >= thresh)]
                true_negatives = len(ele_gene_pairs[ele_gene_pairs['Significant'] == 'False'])
            if true_positives.shape[0] != 0:
                precision = true_positives.shape[0] / (true_positives.shape[0] + false_positives.shape[0])
            else:
                precision = 0
            if true_negatives > 0:
                fpr = len(false_positives) / true_negatives
            else:
                fpr = 0
            if mode == 'pr':
                rec_prec_list.append([recall, precision, thresh])
            else:  # ROC
                rec_prec_list.append([fpr, recall, thresh])
        return rec_prec_list

    binned_list = fetch_prerec([lower_limit + (upper_limit-lower_limit) * x / binning for x in range(0, binning+1)])
    if mark_thresh:
        singular_pos = fetch_prerec([mark_thresh])
    else:
        singular_pos = None
    return binned_list, singular_pos


def pre_rec_plotter(plot_df, sig_col, plot_cols, steps=10000, output_path='', colours=None, no_plot=False,
                    recall_start=0, zorder=None, legend_s=18, mode='pr', thresh=None, legend_out=False,
                    lines_solid=False):
    """Takes a df of validated interactions and iterates through the plot_cols to generate Precision-Recall curves and get
    their AUPRCs. The sig_col indicates if an interaction is true or not."""
    if not colours:
        colours = tol_vibrant*10
    lines = ['solid', 'dashed', 'solid', 'dashdot']*10
    pr_dict = {c: None for c in plot_cols}
    if plot_df.dtypes[sig_col] == bool:
        true_enhancers = plot_df[plot_df[sig_col] == 1]
    else:
        true_enhancers = plot_df[plot_df[sig_col].isin(['True', 'TRUE', 'true', 'T', '1'])]

    f, ax = plt.subplots(1, figsize=(8, 8))
    plt.subplots_adjust(left=0.08, right=0.98)
    ax.set_axisbelow(True)
    ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='both')
    auc_coll = []
    auc_output = []
    if not zorder:
        zorder = [i+10 for i in range(len(plot_cols))]
    for n, this_col in enumerate(plot_cols):
        lower_lim = min(plot_df[this_col])
        upper_lim = max(plot_df[this_col])
        rec_prec_list, thresh_pos = pre_rec_fetcher(steps, plot_df, true_enhancers, score_col=this_col,
                                                    lower_limit=lower_lim, upper_limit=upper_lim, mark_thresh=thresh, mode=mode)
        pr_dict[this_col] = rec_prec_list
        sorted_by_recall = sorted(rec_prec_list, key=lambda x: x[0])
        sorted_by_recall = [x for x in sorted_by_recall if round(x[0], 3) >= recall_start]
        auc = metrics.auc([x[0] for x in sorted_by_recall], [x[1] for x in sorted_by_recall])
        auc_coll.append(this_col + ' AUC: ' + str(round(auc, 3)))
        auc_output.append([this_col, auc])
        print(this_col, round(auc, 4))

        plt.plot([x[0] for x in rec_prec_list], [x[1] for x in rec_prec_list],
                 linestyle='solid' if lines_solid else lines[n], c=colours[n],
                 label=this_col, linewidth=4, zorder=zorder[n])

    if legend_out:
        ax.legend(prop={'size': legend_s, 'weight': 'bold'}, loc='upper right',
                  bbox_to_anchor=(2 if type(legend_out) == bool and not legend_out else legend_out, 1))
    else:
        ax.legend(prop={'size': legend_s, 'weight': 'bold'})
    plt.xlabel('Recall' if mode == 'pr' else 'FPR', size=22)
    plt.ylabel('Precision' if mode == 'pr' else 'TPR', size=22)
    ax.tick_params(axis='both', labelsize=18)
    plt.ylim(0, 1)
    plt.xlim(recall_start, 1)
    ax.set_title(str(len(true_enhancers))+' sig/ '+str(len(plot_df)) + ' interactions' + '\n' + '\n'.join(auc_coll), y=1.05)
    ax.set_facecolor('white')
    if not no_plot:
        f.savefig(output_path + '.pdf', bbox_inches='tight')
    plt.close()
    return auc_output, pr_dict


def basic_bars(plot_df, x_col, y_col, x_order=None, hue_col=None, hue_order=None, title=None, output_path='', y_label='',
               x_size=8, y_size=6, rotation=None, palette=None, legend=True, font_s=14, legend_out=False, ylim=None,
               formats=['pdf']):
    """Plots a basic barplot, allows to select hue levels.
    @param y_col: Can be list, then the df will be transformed long format and var_name set to hue_col. Use y_label to
    have an appropriate y-axis label."""
    if x_col not in plot_df.columns:  # Assumes the x_col is the index if the column doesn't exist.
        plot_df[x_col] = plot_df.index
    if type(y_col) == list:  # Assume that the columns should be stacked and hued.
        plot_df = pd.melt(plot_df, id_vars=x_col, value_vars=y_col, value_name=y_label, var_name=hue_col,
                          ignore_index=False)
        y_col = y_label

    f, ax = plt.subplots(figsize=(x_size, y_size))
    ax.set_axisbelow(True)
    ax.grid(True, axis='y', color='#f2f2f2', linewidth=1, which='major')
    sns.barplot(data=plot_df, x=x_col, y=y_col, order=x_order, hue=hue_col, hue_order=hue_order, ax=ax, alpha=1, edgecolor='k',
                linewidth=1, color='#2d63ad', palette='tab10' if hue_col and not palette else palette)
    if ylim:
        ax.set_ylim(ylim)
    ax.tick_params(axis='both', labelsize=font_s)
    ax.set_ylabel(y_col, fontsize=font_s+2)
    ax.set_xlabel(x_col, fontsize=font_s+2)
    if rotation:
        ax.tick_params(axis='x', rotation=rotation)
    if legend_out:
        ax.legend(prop={'size': 14, 'weight': 'bold'}, loc='upper right',
                  bbox_to_anchor=(2 if type(legend_out) == bool else legend_out, 1))
    else:
        ax.legend(prop={'size': 14, 'weight': 'bold'})
    if not legend and ax.get_legend():
        ax.get_legend().remove()
    plt.title(title, fontsize=font_s+4, fontweight='bold')
    if type(formats) != list:
        formats = list(formats)
    for form in formats:
        f.savefig((output_path + '_'.join([str(x_col), str(y_col)]) + '_Bars.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()



