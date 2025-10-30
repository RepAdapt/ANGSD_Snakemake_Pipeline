rule angsd_snps_allSamples:
    """
    Identify SNPs across all samples using ANGSD
    """
    input:
        bams = f"{config['bam_lists']}/all_samples_bams.txt",
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        gls = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps.beagle.gz',
        mafs = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps.mafs.gz',
        snp_stats = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps.snpStat.gz',
        hwe = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps.hwe.gz',
        pos = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps.pos.gz',
        counts = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps.counts.gz'
    log: f"{LOG_DIR}/angsd_snps_allSamples/{{chrom}}_angsd_snps.log"
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 24000,
        runtime = 1440
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doGlf 2 \
            -baq 2 \
            -ref {input.ref} \
            -doCounts 1 \
            -dumpCounts 3 \
            -minQ 30 \
            -minMapQ 30 \
            -remove_bads 1 \
            -skipTriallelic 1 \
            -uniqueOnly 1 \
            -only_proper_pairs 1 \
            -dosnpstat 1 \
            -doHWE 1 \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule create_sites_file:
    """
    Create position file of ngsParalog
    """
    input:
        rules.angsd_snps_allSamples.output.pos
    output:
        f"{PROGRAM_RESOURCE_DIR}/angsd_sites/{PREFIX}_{{chrom}}.sites"
    shell:
        """
        zcat {input} | tail -n +2 | cut -f1,2 > {output}
        """

rule index_snps:
    """
    Index SNPs sites files for ANGSD
    """
    input:
        sites = rules.create_sites_file.output 
    output:
        idx = f"{PROGRAM_RESOURCE_DIR}/angsd_sites/{PREFIX}_{{chrom}}.sites.idx",
        bin = f"{PROGRAM_RESOURCE_DIR}/angsd_sites/{PREFIX}_{{chrom}}.sites.bin"
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    shell:
        """
        angsd sites index {input}
        """

rule angsd_alleleCounts_freq_byPopulation:
    input:
        bams = lambda w: f"{config['bam_lists']}/{w.population}_bams.txt",
        sites = rules.create_sites_file.output,
        sites_idx = rules.index_snps.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        mafs = f'{ANGSD_DIR}/snps/byPopulation/{{chrom}}/{PREFIX}_{{chrom}}_{{population}}_snps.mafs.gz',
        pos = f'{ANGSD_DIR}/snps/byPopulation/{{chrom}}/{PREFIX}_{{chrom}}_{{population}}_snps.pos.gz',
        counts = f'{ANGSD_DIR}/snps/byPopulation/{{chrom}}/{PREFIX}_{{chrom}}_{{population}}_snps.counts.gz',
    log: f'{LOG_DIR}/angsd_alleleCounts_freqs_byPopulation/{{chrom}}_{{population}}_snps.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/snps/byPopulation/{{chrom}}/{PREFIX}_{{chrom}}_{{population}}_snps'
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        runtime = lambda wildcards, attempt: attempt * 60
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -baq 2 \
            -ref {input.ref} \
            -doCounts 1 \
            -dumpCounts 3 \
            -doMaf 1 \
            -minQ 20 \
            -minMapQ 30 \
            -sites {input.sites} \
            -anc {input.ref} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule get_af_columns_only:
    input:
        rules.angsd_alleleCounts_freq_byPopulation.output.mafs
    output:
        temp(f'{PROGRAM_RESOURCE_DIR}/snps/{{population}}/{{chrom}}_{{population}}_mafs.txt')
    shell:
        """
        echo "chrom\tpos\t{wildcards.population}" > {output}
        zcat {input} | tail -n +2 | cut -f1,2,7 >> {output}
        """

rule concatenate_allele_frequencies:
    input:
        lambda w: expand(rules.get_af_columns_only.output, population=POPULATIONS, chrom=w.chrom)
    output:
        f"{ANGSD_DIR}/snps/{PREFIX}_all_populations_{{chrom}}_mafs.txt"
    conda: "../envs/python.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        runtime = 360
    script:
        "../scripts/python/concatenate_mafs_fst.py"

rule angsd_gea_allele_frequencies_done:
    input:
        expand(rules.concatenate_allele_frequencies.output, chrom=CHROMOSOMES)
    output:
        f"{ANGSD_DIR}/angsd_gea_allele_frequencies.done"
    shell:
        """
        touch {output}
        """
