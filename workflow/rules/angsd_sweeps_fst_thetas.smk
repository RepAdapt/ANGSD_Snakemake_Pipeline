################################
#### SAF AND SFS ESTIMATION ####
################################

rule angsd_saf_likelihood_byPopulation:
    """
    Generate Site Allele Frequency (SAF) likelihood file for each population using ANGSD. 
    """
    input:
        bams = lambda w: f"{config['bam_lists']}/{w.population}_bams.txt",
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        saf = f'{ANGSD_DIR}/saf/{{population}}/{{chrom}}/{{chrom}}_{{population}}.saf.gz',
        saf_idx = f'{ANGSD_DIR}/saf/{{population}}/{{chrom}}/{{chrom}}_{{population}}.saf.idx',
        saf_pos = f'{ANGSD_DIR}/saf/{{population}}/{{chrom}}/{{chrom}}_{{population}}.saf.pos.gz'
    log: f'{LOG_DIR}/angsd_saf_likelihood_byPopulation/{{population}}_{{chrom}}_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/saf/{{population}}/{{chrom}}/{{chrom}}_{{population}}'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        runtime = lambda wildcards, attempt: attempt * 180
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -baq 2 \
            -ref {input.ref} \
            -minQ 30 \
            -minMapQ 30 \
            -doSaf 1 \
            -anc {input.ref} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """
 
rule angsd_estimate_joint_population_sfs:
    """
    Estimated folded, two-dimensional SFS for each population pair using realSFS
    """
    input:
        safs = get_population_saf_files
    output:
        f'{ANGSD_DIR}/sfs/2dsfs/{{chrom}}/{{chrom}}_{{pop_comb}}.2dsfs'
    log: f'{LOG_DIR}/angsd_estimate_joint_population_sfs/{{chrom}}_{{pop_comb}}_2dsfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 12
    resources:
        mem_mb = 3000,
        runtime = lambda wildcards, attempt: attempt * 360 
    shell:
        """
        realSFS {input.safs} \
            -tole 1e-6 \
            -maxIter 30000 \
            -seed 42 \
            -fold 1 \
            -P {threads} > {output} 2> {log}
        """

rule angsd_estimate_sfs_byPopulation:
    """
    Estimate folded SFS separately for each population (i.e., 1D SFS) using realSFS. 
    """
    input:
        saf = rules.angsd_saf_likelihood_byPopulation.output.saf_idx
    output:
        f'{ANGSD_DIR}/sfs/1dsfs/{{chrom}}/{{chrom}}_{{population}}.sfs'
    log: f'{LOG_DIR}/angsd_estimate_sfs_byPopulation/{{chrom}}_{{population}}_1dsfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 6
    resources:
        mem_mb = 2000,
        runtime = lambda wildcards, attempt:  attempt * 180 
    shell:
        """
        realSFS {input.saf} \
            -P {threads} \
            -tole 1e-6 \
            -fold 1 \
            -maxIter 30000 \
            -seed 42 > {output} 2> {log}
        """

########################
### FST AND THETAS #####
########################

rule angsd_population_fst_index:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Fst estimation
    """
    input: 
        saf_idx = get_population_saf_files,
        joint_sfs = rules.angsd_estimate_joint_population_sfs.output
    output:
        fst = f'{ANGSD_DIR}/summary_stats/fst/{{chrom}}/{{chrom}}_{{pop_comb}}.fst.gz',
        idx = f'{ANGSD_DIR}/summary_stats/fst/{{chrom}}/{{chrom}}_{{pop_comb}}.fst.idx'
    log: f'{LOG_DIR}/angsd_habitat_fst_index/{{chrom}}_{{pop_comb}}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 4
    resources:
        mem_mb = 4000,
        runtime = 120
    params:
        fstout = f'{ANGSD_DIR}/summary_stats/fst/{{chrom}}/{{chrom}}_{{pop_comb}}'
    shell:
        """
        realSFS fst index {input.saf_idx} \
            -sfs {input.joint_sfs} \
            -fold 1 \
            -P {threads} \
            -fstout {params.fstout} 2> {log}
        """

rule angsd_fst_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        rules.angsd_population_fst_index.output.idx
    output:
        f'{ANGSD_DIR}/summary_stats/fst/{{chrom}}/{{chrom}}_{{pop_comb}}_readable.fst'
    log: f'{LOG_DIR}/angsd_fst_allSites_readable/{{chrom}}_{{pop_comb}}_readable_fst.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """

rule angsd_estimate_thetas_byPopulation:
    """
    Generate per-site thetas in each habitat from 1DSFS
    """
    input:
        saf_idx = rules.angsd_saf_likelihood_byPopulation.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_byPopulation.output
    output:
        idx = f'{ANGSD_DIR}/summary_stats/thetas/{{chrom}}/{{chrom}}_{{population}}.thetas.idx',
        thet = f'{ANGSD_DIR}/summary_stats/thetas/{{chrom}}/{{chrom}}_{{population}}.thetas.gz'
    log: f'{LOG_DIR}/angsd_estimate_thetas_byPopulation/{{chrom}}_{{population}}_thetas.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 4
    params:
        out = f'{ANGSD_DIR}/summary_stats/thetas/{{chrom}}/{{chrom}}_{{population}}'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        runtime = 60
    shell:
        """
        realSFS saf2theta {input.saf_idx} \
            -P {threads} \
            -fold 1 \
            -sfs {input.sfs} \
            -outname {params.out} 2> {log}
        """

rule angsd_thetas_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        rules.angsd_estimate_thetas_byPopulation.output.idx
    output:
        f'{ANGSD_DIR}/summary_stats/thetas/{{chrom}}/{{chrom}}_{{population}}_readable.thetas'
    log: f'{LOG_DIR}/angsd_thetas_allSites_readable/{{chrom}}_{{population}}_readable_thetas.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    shell:
        """
        thetaStat print {input} > {output} 2> {log}
        """

###########################
### WINDOWED ANALYSES #####
###########################

rule windowed_theta:
    input:
        rules.angsd_estimate_thetas_byPopulation.output.idx
    output:
        f"{ANGSD_DIR}/summary_stats/thetas/{{chrom}}/{{chrom}}_{{population}}_windowedThetas.gz.pestPG"
    log: f'{LOG_DIR}/windowed_theta/{{chrom}}_{{population}}_windowTheta.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f"{ANGSD_DIR}/summary_stats/thetas/{{chrom}}/{{chrom}}_{{population}}_windowedThetas.gz",
        win = 20000,
        step = 20000
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        runtime = 60
    shell:
        """
        thetaStat do_stat {input} -win {params.win} -step {params.step} -outnames {params.out} 2> {log}
        """

rule windowed_fst:
    input:
        rules.angsd_population_fst_index.output.idx
    output:
        f"{ANGSD_DIR}/summary_stats/fst/{{chrom}}/{{chrom}}_{{pop_comb}}_windowed.fst"
    log: f'{LOG_DIR}/windowed_fst/{{chrom}}_{{pop_comb}}_windowedFst.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        win = 20000,
        step = 20000
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        runtime = 60
    shell:
        """
        realSFS fst stats2 {input} -win {params.win} -step {params.step} > {output} 2> {log}
        """

##############
#### POST ####
##############

rule angsd_sweeps_fst_thetas_done:
    """
    Generate empty flag file signalling successful completion of SFS and summary stat for habitats
    """
    input:
        expand(rules.windowed_fst.output, pop_comb=POP_COMBS, chrom=CHROMOSOMES[0]), 
        expand(rules.windowed_theta.output, population=POPULATIONS, chrom=CHROMOSOMES[0]),
        expand(rules.angsd_fst_readable.output, pop_comb=POP_COMBS, chrom=CHROMOSOMES[0]), 
        expand(rules.angsd_thetas_readable.output, population=POPULATIONS, chrom=CHROMOSOMES[0]) 
    output:
        f'{ANGSD_DIR}/angsd_sweeps_fst_thetas.done'
    shell:
        """
        touch {output}
        """
