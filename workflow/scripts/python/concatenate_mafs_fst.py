import polars as pl
from tqdm import tqdm
from functools import reduce

dfs = [pl.read_csv(x, separator="\t") for x in tqdm(snakemake.input)]
result = reduce(
    lambda left, right: left.join(right, on=['chrom', 'pos'], how='full', coalesce=True),
    dfs
)
result.write_csv(snakemake.output[0])
