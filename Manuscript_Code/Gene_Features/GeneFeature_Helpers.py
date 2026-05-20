import pandas as pd
import os
import re
import fnmatch
from pybedtools import BedTool
import gzip
from itertools import chain
from collections import Counter


def fn_patternmatch(pattern):
    """
    Grabs all files in the file system that match the pattern and returns a dictionary with {file: wildcard}.
    Only works if the wildcard is in the file name and not in a directory name.
    Files that start with ._ (e.g. temporary or system files) will be skipped.
    E.g. for a directory BirdCollection/ that contains Bird_Kakapo.txt and Bird_Kea.txt:
    fn_patternmatch("BirdCollection/*.txt" = {"BirdCollection/Bird_Kakapo.txt": "Kakao", "BirdCollection/Bird_Kea.txt": "Kea"}
     """
    parent_folder = '/'.join(pattern.split('/')[:-1]) + '/'
    children_pattern = pattern.split('/')[-1]
    re_pattern = re.compile(children_pattern.replace('*', '(.*)'))
    matched_files = {parent_folder + x: re_pattern.search(x).group(1)
                     for x in os.listdir(parent_folder)
                     if fnmatch.fnmatch(x, children_pattern) and not x.startswith('._')}
    return matched_files


def bed_to_length_dict(bedobj):
    """
    Takes a bedtools object with an identifier (Ensembl ID) in the fourth column, and sums the length of all
    regions of that gene. Used by other functions.
    """
    length_dict = {}
    for entry in bedobj:
        if entry.fields[3] not in length_dict:
            length_dict[entry.fields[3]] = 0
        length_dict[entry.fields[3]] += entry.length
    return length_dict


def bed_to_feature_dict(bedobj):
    """
    Takes a bedtools object and returns a dictionary with the fourth column as key and the bedtools fields as value.
    Multiple occurrences of the identifier are appended.
    """
    feature_dict = {}
    for entry in bedobj:
        if entry.fields[3] not in feature_dict:
            feature_dict[entry.fields[3]] = []
        feature_dict[entry.fields[3]].append(entry.fields)
    return feature_dict



def gene_window_bed(gtf_file, extend=200, gene_set=set(), tss_type='5', dict_only=False, merge=False,
                    open_regions=False):
    """
    Based on a gtf file fetches all or the most 5' TSS for all genes, and returns a BedTool object with windows
    around the TSS, expanding by 'extend' in each direction, resulting in a total window size of 2*'extend'+1.
    Alternatively gives a dictionary with the TSS, also containing the number of transcripts, gene name and strand.
    The BedTools intervals will be 0-based, the TSS in the dictionary still 1-based like in the gtf-file.
    Care: removes the .-suffixes from all gene IDs.

    Args:
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        extend: Number of base pairs to extend the TSS in each direction. 200 means a window of size 401.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        tss_type: "5" to get only the 5' TSS or "all" to get all unique TSS of all transcripts in the gtf-file.
        dict_only: Returns a dictionary instead of a BedTool's object.
        merge: If True, merges all overlapping promoters of the same gene into one row in the BedTool's object.
        open_regions: Optional bed file or BedTools' object, only overlapping parts of promoters will be kept for the
            BedTool's object. Can therefore be used to find genes whose promoter overlap a set of peaks, for example to
            find genes that are accessible.
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

    if merge:  # Flip the chr and Ensembl ID column to merge promoter of the same gene, and afterwards flip again.
        # Due to using the Ensembl ID as the first coordinate we don't need to consider the strand for merging (relying
        # on the annotations putting a gene only on one strand exclusively).
        promoter_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in promoter_bed]), from_string=True).sort().merge(c=[4, 5, 6], o='distinct')
        promoter_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in promoter_bed]), from_string=True)

    return promoter_bed



def gene_body_bed(gtf_file, gene_set=set(), dict_only=False):
    """
    From a gtf-file fetches the gene bodies, meaning start-end for the entries labelled as 'gene'.

    Args:
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        dict_only: Returns a dict {Ensembl ID: [chr, start, end, Ensembl ID, '.', strand]. Else, returns a bedtool-
            object with the coordinates per gene.
    """
    if gene_set:
        gene_set = set([g.split('.')[0] for g in gene_set])

    hits = []
    if gtf_file.endswith('.gz'):
        file_opener = gzip.open(gtf_file, 'rt')
    else:
        file_opener = open(gtf_file)

    with file_opener as gtf_in:
        for line in gtf_in:
            if not line.startswith('#') and line.split('\t')[2] == 'gene':
                # Some gene IDs are non-unique if they have a _PAR_Y version.
                if not line.strip().split('\t')[8].split('gene_id "')[-1].split('";')[0].endswith("_PAR_Y"):
                    line = line.strip().split('\t')
                    gene = line[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]
                    gene_name = line[8].split('gene_name "')[-1].split('"; ')[0].split('.')[0]
                    if not gene_set or gene in gene_set or gene_name in gene_set:
                        hits.append([line[0], line[3], line[4], gene, line[5], line[6]])

    if dict_only:
        return {x[3]: x for x in hits}
    else:
        return BedTool('\n'.join(['\t'.join(x) for x in hits]), from_string=True)



def gene_biotypes(gtf_file, gene_set=set()):
    """
    Gives a dict of Ensembl IDs with their biotype {g: {gtf: as noted in the gtf, general: manually summarized}.
    For the manually summarized biotypes, all gtf-biotypes containing 'RNA' or equal 'ribozyme' are converted to
    'non-coding RNA'. All biotypes containing 'pseudogene' are assigned 'pseudogene'. 'TcR gene' and 'Ig gene' are
    assigned to 'protein-coding'. Removes Ensembl ID version suffixes.

    Args:
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
    """

    def biotyper(biotype):
        """Translates the specific gtf biotypes into more general ones."""
        if 'pseudogene' in biotype:
            return 'pseudogene'
        elif 'RNA' in biotype or biotype == 'ribozyme':
            return "non-coding RNA"
        elif 'TR_' in biotype:
            return 'protein-coding'  #"TcR gene"  # T-cell receptor gene
        elif 'IG_' in biotype:
            return 'protein-coding'  #'Ig gene'  # Immunoglobulin variable chain gene
        elif biotype == 'protein_coding':
            return 'protein-coding'
        else:
            return biotype

    if gtf_file.endswith('.gz'):
        gtf_opener = gzip.open(gtf_file, 'rt')
    else:
        gtf_opener = open(gtf_file)
    biotype_dict = {}
    with gtf_opener as gtf_in:
        for entry in gtf_in:
            if not entry.startswith('#') and entry.split('\t')[2] == 'gene':
                ensembl_id = entry.split('\t')[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]
                gene_name = entry[8].split('gene_name "')[-1].split('"; ')[0].split('.')[0]
                if not gene_set or ensembl_id in gene_set or gene_name in gene_set:
                    gtf_type = entry.split('\t')[8].split('gene_type "')[-1].split('"; ')[0]
                    biotype_dict[ensembl_id] = {'gtf': gtf_type,
                                                'general': biotyper(gtf_type)}
    return biotype_dict



def gene_feature_bed(gtf_file, feature, gene_set=set(), dict_only=False, merge=True, length_only=False,
                     keep_strand=False):
    """
    On gene-level fetch a specific feature of a gene, e.g. all exons, matching the 3rd column in the gtf-file.

    Args:
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        feature: one of gtf's annotated regions: 'CDS', 'Selenocysteine', 'UTR', 'exon', 'gene', 'start_codon',
            'stop_codon', 'transcript'.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        merge: Merge the features per gene.
        dict_only: Return a dict instead.
        length_only: Return a dict with {gene: feature_len}, only makes sense with merge=True.
        keep_strand: Whether the strand should be kept as information.
    """
    if not merge and length_only:
        print("WARNING: Getting length without merging")

    if gene_set:
        gene_set = set([g.split('.')[0] for g in gene_set])

    hits = []
    if gtf_file.endswith('.gz'):
        file_opener = gzip.open(gtf_file, 'rt')
    else:
        file_opener = open(gtf_file)

    with file_opener as gtf_in:
        for line in gtf_in:
            if not line.startswith('#') and line.split('\t')[2] == feature:
                # Some gene IDs are non-unique if they have a _PAR_Y version.
                if not line.strip().split('\t')[8].split('gene_id "')[-1].split('";')[0].endswith("_PAR_Y"):
                    line = line.strip().split('\t')
                    gene = line[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]
                    gene_name = line[8].split('gene_name "')[-1].split('"; ')[0].split('.')[0]
                    if not gene_set or gene in gene_set or gene_name in gene_set:
                        if keep_strand:
                            hits.append([line[0], line[3], line[4], gene, line[5], line[6]])
                        else:
                            hits.append([line[0], line[3], line[4], gene])

    feature_bed = BedTool('\n'.join(['\t'.join(x) for x in hits]), from_string=True)

    if merge:  # Flip the chr and Ensembl ID column to merge features of the same gene, and afterwards flip again.
        # We can ignore the strand for merging when using the Ensembl ID as first coordinate.
        if keep_strand:
            feature_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in feature_bed]), from_string=True).sort().merge(c=4, o='distinct')
            feature_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in feature_bed]), from_string=True)
        else:
            feature_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0]]) for x in feature_bed]), from_string=True).sort().merge(c=4, o='distinct')
            feature_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0]]) for x in feature_bed]), from_string=True)

    if length_only:
        return bed_to_length_dict(feature_bed)

    if dict_only:
        return bed_to_feature_dict(feature_bed)

    else:
        return feature_bed


def gene_feature_table(gtf_file, gene_set=(), extend=500000):
    """
    Based on a gtf-file fetch several features on gene-level into one DataFrame, relying heavily on other
    gtf-processing functions.

    Args:
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        extend: How much the gene windows should be extended to each side of the TSS to calculate gene density. 1000 means a window of size 2001.

    Returns:
        pandas DataFrame:
        A DataFrame with the features, with the indices being either the gene_set or the gtf-file genes. Gene density refers to the number of other genes within a window of ±extend.
    """
    print("Fetching transcripts, size and biotype")
    tss_annot = gene_window_bed(gtf_file=gtf_file, gene_set=gene_set, tss_type='all', dict_only=True)
    bodies = gene_body_bed(gtf_file=gtf_file, gene_set=gene_set, dict_only=True)
    biotypes = gene_biotypes(gtf_file=gtf_file, gene_set=gene_set)
    exons = gene_feature_bed(gtf_file, feature='exon', gene_set=set(), merge=True, length_only=True)

    # For the gene density we need the whole annotation.
    print("Getting gene density")
    gene_windows = gene_window_bed(gtf_file=gtf_file, tss_type='5', extend=extend, dict_only=False)
    tss_bed = gene_window_bed(gtf_file=gtf_file, tss_type='5', dict_only=False, extend=0)
    gene_inter_counter = Counter([x.fields[3] for x in gene_windows.intersect(tss_bed, wa=True)])
    gene_inter_counter = {g: val - 1 for g, val in gene_inter_counter.items()}  # To correct the hit with itself.

    if not gene_set:
        df_genes = list(tss_annot.keys())
    else:
        df_genes = list(gene_set)

    gene_df = pd.DataFrame()
    gene_df.index = df_genes
    gene_df.index.name = 'Ensembl ID'
    gene_df['Gene name'] = [None if g not in tss_annot else tss_annot[g]['name'] for g in df_genes]
    gene_df['#TSS'] = [None if g not in tss_annot else len(tss_annot[g]['tss']) for g in df_genes]
    gene_df['#Transcripts'] = [None if g not in tss_annot else tss_annot[g]['#transcripts'] for g in df_genes]
    gene_df['Gene length'] = [None if g not in bodies else int(bodies[g][2]) - int(bodies[g][1]) for g in df_genes]
    gene_df['Exons length'] = [None if g not in exons else exons[g] for g in df_genes]
    gene_df['Biotype gtf'] = [None if g not in biotypes else biotypes[g]['gtf'] for g in df_genes]
    gene_df['Biotype general'] = [None if g not in biotypes else biotypes[g]['general'] for g in df_genes]
    gene_df['Gene density'] = [None if g not in gene_inter_counter else gene_inter_counter[g] for g in df_genes]

    return gene_df