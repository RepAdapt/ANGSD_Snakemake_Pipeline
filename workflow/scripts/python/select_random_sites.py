import random

lines = open(snakemake.input[0], "r").readlines()
chr_prop_length_dict = {}
total_length = 0
for l in lines:
    sl = l.split("\t")
    chr = sl[0]
    length = int(sl[1])
    if chr in snakemake.params['chroms']:
        total_length += int(length)
for l in lines:
    sl = l.split("\t")
    chr = sl[0]
    length = int(sl[1])
    if chr in snakemake.params['chroms']:
        chr_prop_length_dict[chr] = round(int(snakemake.params['total_sites'] * (length / total_length))

sites_dict = {}
random.seed(42)
for l in lines:
    sl = l.split("\t")
    chr = sl[0]
    length = int(sl[1])
    if chr in snakemake.params['chroms']:
        sites_dict[chr] = sorted(random.sample(range(1, length), chr_prop_length_dict[chr]))
with open(snakemake.output[0], "w") as fout:
    for k, v in sites_dict.items():
        for site in v:
            fout.write(f"{k}\t{site}\n")
