import gzip
from pybedtools import BedTool
from itertools import chain
from matplotlib import pyplot as plt
import seaborn as sns
from collections import Counter


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
