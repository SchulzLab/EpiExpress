import pyBigWig
from timeit import default_timer as clock
from multiprocess import Pool
import pandas as pd


def fetch_counts(args):
    """Gets a vector of counts for the bed_regions for one individual bw_file, so that we can parallelize over the
    bw files in the outer function."""
    bw_file, bed_regions = args
    collected_counts = []
    this_bw = pyBigWig.open(bw_file)

    for region in bed_regions:
        try:
            bw_count = this_bw.stats(region[0], int(region[1]), int(region[2]), type="mean")[0]
            if bw_count is None:
                bw_count = 0
                # errors.append([region, "count was NaN"])
            collected_counts.append(bw_count)
        except (RuntimeError, IndexError):
            collected_counts.append(0)
            # errors.append([region, "Invalid interval bounds"])

    if not len(collected_counts) == len(bed_regions):
        print("ERROR: Mismatch between fetched counts and number of bed regions", bw_file)
    return collected_counts


def bigwig_counts(bed_file, bigwigs, n_cores=1):
    """
    For a bed file gets the mean signal for each bigwig file. If pyBigWig can't retrieve a count it will be set to 0.
    @param bed_file: Either a BedTools object or a path to a bed-file.
    @param bigwigs: List of bigwigs for which pyBigWig will be used to extract the mean signal for.
    @return Will give a dataframe of the bed-file's regions with a column for each bigwig file and the file name
    as column. errors is a list of regions that are also part of the bed_regions but failed due to not
    returning a count or throwing an error of invalid interval bounds. Those have a count of 0.
    """
    start = clock()

    if type(bed_file) == str:
        bed_regions = [x.strip().split('\t') for x in open(bed_file).readlines() if not x.startswith('#')]
    else:
        bed_regions = [str(x).strip().split('\t') for x in str(bed_file).strip().split('\n')]
    if type(bigwigs) == set:
        bigwig_order = list(bigwigs)
    else:
        bigwig_order = bigwigs
    errors = []

    process_pool = Pool(processes=n_cores)
    bw_counts = process_pool.map(fetch_counts, map(lambda x: (x, bed_regions), bigwig_order))
    process_pool.close()
    region_counts = pd.DataFrame(bw_counts).T
    region_counts = pd.concat([pd.DataFrame([x[:3] for x in bed_regions]), region_counts], ignore_index=True, axis=1)
    region_counts.columns = ['#chr', 'start', 'end'] + bigwig_order

    print(clock() - start)

    return region_counts, errors


