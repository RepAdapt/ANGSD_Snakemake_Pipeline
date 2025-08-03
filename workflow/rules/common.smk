def get_population_saf_files(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byPopulation.output.saf_idx, population=POPULATIONS, chrom=wildcards.chrom)
    pop1 = wildcards.pop_comb.split('_')[0]
    pop2 = wildcards.pop_comb.split('_')[1]
    regex1 = r'(?<=_)' + re.escape(pop1) + '(?=\.saf\.idx)'
    regex2 = r'(?<=_)' + re.escape(pop2) + '(?=\.saf\.idx)'
    saf1 = [x for x in all_saf_files if re.search(regex1, os.path.basename(x), re.IGNORECASE)]
    saf2 = [x for x in all_saf_files if re.search(regex2, os.path.basename(x), re.IGNORECASE)]
    return saf1 + saf2

def get_population_saf_files_random100Mb(wildcards):
    all_saf_files = expand(rules.angsd_saf_random100Mb_byPopulation.output.saf_idx, population=POPULATIONS)
    pop1 = wildcards.pop_comb.split('_')[0]
    pop2 = wildcards.pop_comb.split('_')[1]
    regex1 = r'^' + re.escape(pop1) + r'(?=_random100Mb.saf\.idx)'
    regex2 = r'^' + re.escape(pop2) + r'(?=_random100Mb.saf\.idx)'
    saf1 = [x for x in all_saf_files if re.search(regex1, os.path.basename(x), re.IGNORECASE)]
    saf2 = [x for x in all_saf_files if re.search(regex2, os.path.basename(x), re.IGNORECASE)]
    return saf1 + saf2
