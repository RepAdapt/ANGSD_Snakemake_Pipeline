rule copy_ref:
    """
    Copy reference genome to different directory
    """
    input:
        REFERENCE_GENOME
    output:
        f"{REF_DIR}/{PREFIX}.fasta.gz"
    conda: '../envs/ref.yaml'
    log: f"{LOG_DIR}/copy_ref/copy_ref.log"
    shell:
        """
        cp {input} {output}
        """

rule samtools_index_ref:
    """
    Index reference with Samtools
    """
    input:
        rules.copy_ref.output
    output:
        f"{REF_DIR}/{PREFIX}.fasta.gz.fai"
    conda: '../envs/ref.yaml'
    log: f"{LOG_DIR}/samtools_index_ref/samtools_index_ref.log"
    shell:
        """
        sleep 10
        samtools faidx {input} 2> {log}
        """
