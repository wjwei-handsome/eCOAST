from math import log
import pandas as pd
from rich.progress import Progress

from sklearn.cluster import AgglomerativeClustering
import sys

# in_file = sys.argv[1]
# out_file = sys.argv[2]

BEDGRAPH_DIR = "/public/home/wjwei/WORK/EXP_BIAS/999_all_cov/each_samples_covbed"
SAMPLE_LIST = sys.argv[1]
GENE_LIST = sys.argv[2]
OUT_FILE = sys.argv[3]
dist_threshold = float(sys.argv[4])


def read_fmt_df(in_file):
    df = pd.read_table(in_file, header=None)
    # add col name: chro; start; end; value; type
    df.columns = ["chro", "start", "end", "value", "type"]
    # drop type != exon
    df = df[df["type"] == "exon"]
    # convert value to int
    df["value"] = df["value"].astype(int)
    # compute log value
    df["log"] = df["value"].apply(lambda x: 0 if x == 0 else log(x))
    # convert log to float
    df["log"] = df["log"].astype(float)
    # add length
    df["length"] = df["end"] - df["start"]
    return df


def cluster(df):
    # perform cluster
    agg_dist = AgglomerativeClustering(
        n_clusters=None, linkage="average", distance_threshold=dist_threshold
    )
    data = df["log"].values.reshape(-1, 1)
    try:
        labels_dist = agg_dist.fit_predict(data)
    except ValueError:  # if just one sample, then return 0
        labels_dist = [0] * len(data)
    # df["labels_a2"] = labels_a2
    # df["labels_a3"] = labels_a3
    df["labels_dist"] = labels_dist
    return df


def convert_label(df):
    k = df["labels_dist"].nunique()
    mean_values = {}
    for i in range(k):
        mean_values[i] = df[df["labels_dist"] == i]["log"].mean()
    rank = sorted(mean_values, key=mean_values.get)
    # define rank
    df["rank"] = df["labels_dist"].apply(lambda x: f"rank_{rank.index(x)}")

    return df


def compute_high_ratio(df):
    lenth = df["length"].sum()
    # compute ratio in agg_dist by drop lowest
    drop_lowest = df[df["rank"] != "rank_0"]["length"].sum()
    drop_lowest_ratio = drop_lowest / lenth
    return high_ratio, highest_ratio, mid_high_ratio, drop_lowest_ratio


# main
def main(in_file):
    df = read_fmt_df(in_file)
    df = cluster(df)
    df = convert_label(df)
    high_ratio, highest_ratio, mid_high_ratio, drop_lowest_ratio = compute_high_ratio(
        df
    )
    return high_ratio, highest_ratio, mid_high_ratio, drop_lowest_ratio, df


# get sample_list
sample_list = []
with open(SAMPLE_LIST, "r") as f:
    for line in f:
        sample_list.append(line.strip())

# get gene_list
gene_list = []
with open(GENE_LIST, "r") as f:
    for line in f:
        gene_list.append(line.strip())

with Progress() as pb:
    t_sample = pb.add_task("[green]Sample", total=len(sample_list))
    t_gene = pb.add_task("[green]Gene", total=len(gene_list))

    with open(OUT_FILE, "w") as f:
        ## write header
        f.write(
            f"sample\tgene\thigh_ratio\thighest_ratio\tmid_high_ratio\tdrop_lowest_ratio_{dist_threshold}\n"
        )
        for i, sample in enumerate(sample_list):
            for j, gene in enumerate(gene_list):
                ## write sample and gene
                f.write(f"{sample}\t{gene}\t")
                in_file = f"{BEDGRAPH_DIR}/{sample}/{sample}.{gene}.bedgraph"
                high_ratio, highest_ratio, mid_high_ratio, drop_lowest_ratio, df = main(
                    in_file
                )
                ## write high_ratio, keep_highest_ratio, drop_lowest_ratio
                f.write(
                    f"{high_ratio}\t{highest_ratio}\t{mid_high_ratio}\t{drop_lowest_ratio}\n"
                )
                df.to_csv(f"{sample}.{gene}.low", sep="\t", index=False)

                pb.update(t_gene, completed=j + 1)
            pb.update(t_sample, completed=i + 1)
