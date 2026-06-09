import random

# Read all verified sites from per-chromosome sites files
chr_sites_dict = {}
total_sites_available = 0

for chrom in snakemake.params['chroms']:
    sites_file = snakemake.input[0].replace('random.sites', f'MVP2026-armstrongoyster_{chrom}.sites')
    sites = []
    with open(sites_file, 'r') as f:
        for line in f:
            sl = line.strip().split()
            sites.append(int(sl[1]))
    chr_sites_dict[chrom] = sites
    total_sites_available += len(sites)

# Sample proportionally by chromosome based on available verified sites
random.seed(42)
selected_sites = {}
for chrom, sites in chr_sites_dict.items():
    prop = len(sites) / total_sites_available
    n_sites = round(snakemake.params['total_sites'] * prop)
    n_sites = min(n_sites, len(sites))  # can't sample more than available
    selected_sites[chrom] = sorted(random.sample(sites, n_sites))

with open(snakemake.output[0], 'w') as fout:
    for chrom, sites in selected_sites.items():
        for site in sites:
            fout.write(f"{chrom}\t{site}\n")
