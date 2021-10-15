rule aseRC:
    input:
        bam = rules.wasp_merge.output.keep_sorted,
        bamIndex = rules.wasp_merge.output.keep_sortedIndex,
        vcf = rules.variants_select.output.vcf,
        refFa = rules.refs_downloadFa.output,
        intervals = rules.intervals.output
    output:
        tsv = os.path.join(analysis.aseRC_dir, "wasp/{SAMPLE}.tsv")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-04:00:00"
    shell:
        """
        gatk \
            ASEReadCounter \
            -I {input.bam} \
            -V {input.vcf} \
            -R {input.refFa} \
            -L {input.intervals} \
            -O {output.tsv} \
            --min-mapping-quality 10 \
            --min-base-quality 20
        """

rule aseRC_nowasp:
    input:
        bam = rules.bqsr_apply.output.bam,
        bamIndex = rules.bqsr_apply.output.bamIndex,
        vcf = rules.variants_select.output.vcf,
        refFa = rules.refs_downloadFa.output,
        intervals = rules.intervals.output
    output:
        tsv = os.path.join(analysis.aseRC_dir, "no_wasp/{SAMPLE}.tsv")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-04:00:00"
    shell:
        """
        gatk \
            ASEReadCounter \
            -I {input.bam} \
            -V {input.vcf} \
            -R {input.refFa} \
            -L {input.intervals} \
            -O {output.tsv} \
            --min-mapping-quality 10 \
            --min-base-quality 20
        """