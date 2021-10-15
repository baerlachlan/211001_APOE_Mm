rule refs_downloadFa:
    output:
        temp(analysis.ref_path)
    params:
        ref_url = analysis.ref_url
    threads:
        1
    shell:
        """
        curl -L {params.ref_url} | gzip -d > {output}
        """

rule refs_downloadGtf:
    output:
        temp(analysis.gtf_path)
    params:
        gtf_url = analysis.gtf_url
    threads:
        1
    shell:
        """
        curl -L {params.gtf_url} > {output}
        """

rule refs_downloadDbsnp:
    output:
        temp(analysis.dbsnp_path)
    params:
        dbsnp_url = analysis.dbsnp_url
    threads:
        1
    shell:
        """
        curl -L {params.dbsnp_url} > {output}
        """

rule refs_refDict:
    input:
        rules.refs_downloadFa.output
    output:
        temp(analysis.ref_path.rstrip("fa") + "dict")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-00:10:00"
    shell:
        "gatk CreateSequenceDictionary -R {input}"

rule refs_refIndex:
    input:
        rules.refs_downloadFa.output
    output:
        temp(analysis.ref_path + ".fai")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-00:10:00"
    shell:
        "samtools faidx {input}"

rule refs_starIndex:
    input:
        ref_fa = rules.refs_downloadFa.output,
        gtf = rules.refs_downloadGtf.output
    output:
        temp(directory(os.path.join(analysis.refs_dir, "star")))
    params:
        overhang = analysis.read_length - 1
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 16,
        ntasks = 1,
        mem_mb = 32000,
        time = "00-00:30:00"
    shell:
        """
        zcat {input.gtf} > temp.gtf

        STAR \
            --runThreadN {resources.cpu} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.ref_fa} \
            --sjdbGTFfile temp.gtf \
            --sjdbOverhang {params.overhang}

        rm temp.gtf
        """

rule refs_dbsnpIndex:
    input:
        rules.refs_downloadDbsnp.output
    output:
        analysis.dbsnp_path + ".tbi"
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-00:10:00"
    shell:
        "tabix -p vcf {input}"