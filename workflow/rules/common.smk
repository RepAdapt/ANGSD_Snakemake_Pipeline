def get_population_saf_files(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byPopulation.output.saf_idx, population=POPULATIONS, chrom=wildcards.chrom)
    first_hab = wildcards.pop_comb.split('_')[0]
    second_hab = wildcards.pop_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if f'{first_hab}' == os.path.basename(x).split("_")[1].split(".")[0]]
    saf2 = [x for x in all_saf_files if f'{second_hab}' == os.path.basename(x).split("_")[1].split(".")[0]]
    return saf1 + saf2
