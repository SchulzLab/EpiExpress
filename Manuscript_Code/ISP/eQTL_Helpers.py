import gzip
import numpy as np
import re
from pybedtools import BedTool
from itertools import chain
from matplotlib import pyplot as plt
import seaborn as sns
from collections import Counter
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec
from matplotlib import cm, colors, colorbar


def sanitize_path(path_string):
    """Function to remove unwanted characters from a file path before saving. Put here since we need it for a lot
    of functions and eases adjustment e.g. for OS."""
    return re.sub(r'[^A-Za-z0-9._\-\\\\////]+', '', path_string)


def gene_window_bed(gtf_file, extend=200, gene_set=set(), tss_type='5', dict_only=False, merge=False,
                    open_regions=False):
    """
    Based on a gtf file fetches all or the most 5' TSS for all genes, and returns a BedTool object with windows
    around the TSS, expanding by 'extend' in each direction, resulting in a total window size of 2*'extend'+1.
    Alternatively gives a dictionary with the TSS.
    The BedTools intervals will be 0-based, the TSS in the dictionary still 1-based like in the gtf-file.
    Care: removes the .-suffixes from all gene IDs.
    @param gtf_file: gtf-file in GENCODE's format, either .gz or .gtf
    @param extend: number of base pairs to extend the TSS in each direction
    @param gene_set: Limits the output to the given gene set, leave empty to get all.
    @param tss_type: "5" to get only the 5' TSS or "all" to get all unique TSS of all transcripts in the gtf-file
    @param dict_only: Returns a dictionary instead of a BedTool's object.
    @param merge: If True, merges all intersecting promoter of the same gene into one row in the BedTool's object.
    @param open_regions: Optional bed file or BedTools' object, only overlapping parts of promoters will be kept for the
                         BedTool's object.
    """
    if tss_type == '5':
        identifier = 'gene'
    elif tss_type == 'all':
        identifier = 'transcript'
    if gtf_file.endswith('.gz'):
        file_opener = gzip.open(gtf_file, 'rt')
    else:
        file_opener = open(gtf_file)

    if gene_set:
        gene_set = set([g.split('.')[0] for g in gene_set])

    tss_locs = {}
    with file_opener as gtf_in:
        for entry in gtf_in:
            if not entry.startswith('#') and entry.split('\t')[2] == identifier:
                line = entry.strip().split('\t')
                # Some gene IDs are non-unique if they have a _PAR_Y version.
                if not line[8].split('gene_id "')[-1].split('";')[0].endswith("_PAR_Y"):
                    this_gene = line[8].split('gene_id "')[-1].split('";')[0].split('.')[0]
                    gene_name = line[8].split('gene_name "')[-1].split('";')[0]
                    if not gene_set or this_gene in gene_set or gene_name in gene_set:
                        if this_gene not in tss_locs:
                            tss_locs[this_gene] = {'chr': None, 'tss': set(), '#transcripts': 0}

                        tss_locs[this_gene]['chr'] = line[0]
                        tss_locs[this_gene]['name'] = gene_name
                        if line[6] == '+':
                            if identifier == 'gene' and (not tss_locs[this_gene]['tss'] or list(tss_locs[this_gene]['tss'])[0] > int(line[3])):
                                tss_locs[this_gene]['tss'] = {int(line[3])}
                            elif identifier == 'transcript':
                                tss_locs[this_gene]['tss'].add(int(line[3]))
                                tss_locs[this_gene]['#transcripts'] += 1
                            tss_locs[this_gene]['strand'] = '+'
                        if line[6] == '-':
                            if identifier == 'gene' and (not tss_locs[this_gene]['tss'] or list(tss_locs[this_gene]['tss'])[0] < int(line[4])):
                                tss_locs[this_gene]['tss'] = {int(line[4])}
                            elif identifier == 'transcript':
                                tss_locs[this_gene]['tss'].add(int(line[4]))
                                tss_locs[this_gene]['#transcripts'] += 1
                            tss_locs[this_gene]['strand'] = '-'

    if dict_only:
        return tss_locs

    promoter_bed = BedTool('\n'.join(chain(*[[vals['chr'] + '\t' + str(max([0, tss - int(extend) - 1])) + '\t' +
                                             str(tss + int(extend)) + '\t' + g + '\t.\t' + vals['strand'] for tss in vals['tss']]
                                             for g, vals in tss_locs.items()])), from_string=True)

    if open_regions and str(open_regions).lower() != "false":
        if type(open_regions) == str:
            open_regions = BedTool('\n'.join(['\t'.join(x.strip().split('\t')[:3]) for x
                                              in open(open_regions).readlines() if not x.startswith('#')]), from_string=True)
        promoter_bed = promoter_bed.intersect(open_regions)

    if merge:  # Flip the chr and geneID column to merge promoter of the same gene, and afterwards flip again.
        promoter_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in promoter_bed]), from_string=True).sort().merge(c=[4, 5, 6], o='distinct')
        promoter_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in promoter_bed]), from_string=True)

    return promoter_bed



def basic_violin(plot_df, y_col, x_col, x_order=None, hue_col=None, hue_order=None, title=None, output_path='',
                 numerate=False, ylim=None, palette=None, xsize=12, ysize=8, boxplot=False, boxplot_meanonly=False,
                 rotation=None, numerate_break=True, jitter=False, colour='#2d63ad', font_s=14, saturation=0.75,
                 jitter_colour='black', jitter_size=5, vertical_grid=False, legend_title=True, legend=True, grid=True,
                 formats=['pdf']):
    """Plots a basic violin plot which allows for hue, whose order can be defined as well.
    Use y_col=None and x_col=None for seaborn to interpret the columns as separate plots on the x-asis.
    @param boxplot_meanonly: Remove all lines from the boxplot and show just the mean as horizontal line."""
    f, ax = plt.subplots(figsize=(xsize, ysize))
    ax.set_axisbelow(True)
    if grid:
        ax.grid(True, axis='both', color='#f2f2f2', linewidth=1, which='major')
    if not boxplot:
        vio = sns.violinplot(data=plot_df, y=y_col, x=x_col, order=x_order, hue=hue_col, ax=ax,
                             color=colour if not palette else None, saturation=saturation,
                             palette='tab10' if hue_col and not palette else palette, hue_order=hue_order)
    else:
        if boxplot_meanonly:
            vio = sns.boxplot(data=plot_df, y=y_col, x=x_col, order=x_order, hue=hue_col, ax=ax, hue_order=hue_order,
                              showfliers=False, showbox=False, showcaps=False, showmeans=True, meanline=True,
                              meanprops={'color': 'k', 'ls': '-', 'lw': 2}, medianprops={'visible': False},
                              whiskerprops={'visible': False}, zorder=1, saturation=saturation,
                              palette='tab10' if hue_col and not jitter_colour else jitter_colour)
        else:
            vio = sns.boxplot(data=plot_df, y=y_col, x=x_col, order=x_order, hue=hue_col, ax=ax, saturation=saturation,
                              color=colour if not palette else None, showfliers=False if jitter else True,
                              palette='tab10' if hue_col and not palette else palette, hue_order=hue_order)
    if jitter:
        sns.stripplot(data=plot_df, x=x_col, y=y_col, jitter=True, ax=ax, hue=hue_col, hue_order=hue_order, zorder=10,
                      order=x_order, palette=jitter_colour, dodge=True, legend=False, edgecolor='black', linewidth=1,
                      size=jitter_size)
    ax.tick_params(axis='both', labelsize=font_s+4)
    ax.set_ylabel(y_col, fontsize=font_s+8)
    ax.set_xlabel(x_col, fontsize=font_s+8)
    if ylim:
        ax.set_ylim(ylim)
    if numerate:
        if not x_col:
            ax.set_xticklabels(['(#' + str((~plot_df[y_col].isna()).sum()) + ')' for x in ax.get_xmajorticklabels()])
        else:
            count_df = plot_df[[x_col, y_col]][~plot_df[y_col].isna()]
            x_counts = Counter(count_df[x_col].values)
            ax.set_xticklabels([x._text+'\n'*numerate_break+'(#'+str(x_counts[x._text])+')' for x in ax.get_xmajorticklabels()])
    if hue_col:
        plt.setp(vio.get_legend().get_texts(), fontsize=font_s)
        plt.setp(vio.get_legend().get_title(), fontsize=font_s+2)
        if not legend_title:
            vio.get_legend().set_title('')
        sns.move_legend(vio, prop={'size': 14, 'weight': 'bold'}, loc='best')
    if rotation:
        plt.xticks(rotation=rotation, ha='center')
    if vertical_grid:  # Fun part is minor ticks are always x5.
        for x in range(len(set(plot_df[x_col]))):
            plt.axvline(x+0.5, color='#f2f2f2', linewidth=1, zorder=0)
    if not legend:
        ax.get_legend().remove()
    plt.title(title, fontsize=22, fontweight='bold', y=1.02)
    if type(formats) != list:
        formats = list(formats)
    for form in formats:
        f.savefig((output_path + str(x_col) + '_' + str(y_col) + '_' + str(hue_col) + '_Violin.'+form).replace(' ', ''),
                  bbox_inches='tight', format=form)
    plt.close()



def heatmap_cols(plot_df, cmap_cols, plot_out, row_label_col=None, column_labels=None, class_col=None,
                 x_size=20, y_size=40, title="", annot_cols=None, width_ratios=None, wspace=0.4, rasterized=True,
                 annot_s=10, ticksize=14, heat_ticksize=14, square=False, x_rotation=70, y_rotation=0,
                 ax_fontweight='normal', row_label_first=False, x_label_pos='top', formats=['pdf']):
    """
    Multiple heatmaps side-by-side but the same rows. Allows to show several metrics for the same rows with different
    colourmaps etc. E.g. for a list of top differential genes first a heatmap of baseline expression coloured by TPM,
    followed by a separate heatmap-block with the log2FC for the same genes.

    Args:
        cmap_cols: Dictionary with one entry for each block. The keys don't matter as long as they are unique. E.g.
            {0: {'cols': ['Mean_Control_FM', 'Mean_FM_Mock_Ctrl', 'Mean_Tcf15_FM', 'Mean_FM_Tcf15_OE'],
                     'centre': 0, (optional)
                     'cmap': 'mako',
                     'cbar_label': 'TPM',
                     'vmax': 200, (optional)
                     'vmin': 0, (optional)
                     }
        row_label_col: Column where to fetch the row-strings from. Set to None to use the index.
        column_labels: Alternative to using the column names as indicated in cmap_cols.
        class_col: Column that should be added as separate first heatmap-block, should be categorical.
        annot_cols: Dictionary of {"column": "column with annotation string"} to write the strings in the value into the cells of columns.
        width_ratios: Ratios of the widths of each heatmap-block.
        wspace: Additional horizontal space between blocks.
        rasterized: Whether to draw thin white lines around cells.
        square: Whether cells should be squares.
        row_label_first: Only write the row names for the first entry and skip for the others.
    """
    all_cols = list(chain(*[c['cols'] for c in cmap_cols.values()]))
    f, axes = plt.subplots(nrows=1, ncols=len(cmap_cols)+bool(class_col), figsize=(x_size, y_size),
                           gridspec_kw={'width_ratios': [0.1] * bool(class_col) + [0.9*len(c['cols'])/len(all_cols) for c in cmap_cols.values()] if not width_ratios else width_ratios})
    # Since we might have different colourmaps we define two matrices for each, one for the values
    # and one for the annotation.
    for n, c_attrs in enumerate(cmap_cols.values()):
        if class_col:
            n += 1
        if len(cmap_cols) == 1 and not class_col:  # If we only had one heatmap we can't index axes.
            this_ax = axes
        else:
            this_ax = axes[n]
        # this_cmap = cm.get_cmap(c_attrs['cmap'])

        value_mat = np.zeros([len(plot_df), len(c_attrs['cols'])])
        annot_mat = np.full([len(plot_df), len(c_attrs['cols'])], '', dtype=object)  # Numpy complains otherwise.
        for c, col in enumerate(c_attrs['cols']):
            value_mat[:, c] = plot_df[col].values
            if annot_cols:
                if col in annot_cols:
                    annot_mat[:, c] = plot_df[annot_cols[col]].values
        heat = sns.heatmap(value_mat, ax=this_ax,  rasterized=rasterized,
                           yticklabels=plot_df.index if not row_label_col else plot_df[row_label_col].values,
                           xticklabels=c_attrs['cols'] if not column_labels else column_labels, cbar=True,
                           cmap=c_attrs['cmap'], fmt='', annot=annot_mat,
                           annot_kws={'size': annot_s}, cbar_kws={'label': c_attrs['cbar_label'], 'shrink': 0.7},
                           center=None if 'centre' not in c_attrs or c_attrs['centre'] is False else c_attrs['centre'],
                           vmin=None if 'vmin' not in c_attrs else c_attrs['vmin'],
                           vmax=None if 'vmax' not in c_attrs else c_attrs['vmax'], square=square)
        if row_label_first and n > 0:
            heat.tick_params(left=False, labelleft=False)
        heat.set_xticklabels(heat.get_xmajorticklabels(), fontsize=ticksize, rotation=x_rotation, fontweight=ax_fontweight)
        heat.set_yticklabels(heat.get_ymajorticklabels(), fontsize=ticksize, rotation=y_rotation, fontweight=ax_fontweight)
        if 'row_labels' in c_attrs and not c_attrs['row_labels']:
                heat.axes.get_yaxis().set_visible(False)
        if x_label_pos == 'top':
            this_ax.tick_params(axis='x', labeltop=True, top=True, labelbottom=False, bottom=False)
        else:
            this_ax.tick_params(axis='x', labeltop=False, top=False, labelbottom=True, bottom=True)
        heat_cbar = heat.collections[0].colorbar
        heat_cbar.ax.tick_params(labelsize=heat_ticksize)
        heat_cbar.ax.yaxis.label.set_fontsize(heat_ticksize)

        if class_col and n == 1:
            # Add the class bar as separate one-column heatmap to the left.
            class_to_int = {c: i for i, c in enumerate(set(plot_df[class_col]))}
            if len(class_to_int) == 2:
                class_cmap = LinearSegmentedColormap.from_list("two_contrast", ['#E1BE6A', '#40B0A6'], N=2)
            else:
                class_cmap = cm.get_cmap("tab20", len(class_to_int))
            sns.heatmap([[class_to_int[c]] for c in plot_df[class_col].values], cmap=class_cmap, ax=axes[0],
                        xticklabels=False, yticklabels=False, square=square,
                        cbar_kws={'label': class_col, 'location': 'left', 'shrink': 1})  # Shrink is somehow ignored here.
            colorbar = axes[0].collections[0].colorbar
            r = colorbar.vmax - colorbar.vmin
            colorbar.set_ticks([colorbar.vmin + r / len(class_to_int) * (0.5 + i) for i in range(len(class_to_int))])
            colorbar.set_ticklabels(list(class_to_int.keys()))
            colorbar.ax.yaxis.set_label_position('left')
            colorbar.ax.yaxis.set_ticks_position('left')
            colorbar.ax.yaxis.label.set_fontsize(14)
            colorbar.ax.tick_params(labelsize=14)
    if title:
        plt.title(title, size=18, fontweight='bold')
    plt.subplots_adjust(wspace=0.4 if not wspace else wspace)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        plt.savefig(sanitize_path(plot_out + "_MultiColHeatmap."+form), bbox_inches='tight', format=form)
    plt.close()
