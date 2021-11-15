rule variants_gvcf:
    input:
        bam = rules.bqsr_apply.output.bam,
        bamIndex = rules.bqsr_apply.output.bamIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
        refDbsnp = rules.refs_downloadDbsnp.output,
        intervals = rules.intervals.output
    output:
        gvcf = temp(os.path.join(analysis.variants_dir, "1_gvcf/{SAMPLE}.g.vcf.gz")),
        gvcfIndex = temp(os.path.join(analysis.variants_dir, "1_gvcf/{SAMPLE}.g.vcf.gz.tbi")),
        detailMetrics = os.path.join(analysis.variants_dir, "1_gvcf/log/{SAMPLE}.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join(analysis.variants_dir, "1_gvcf/log/{SAMPLE}.variant_calling_summary_metrics")
    params:
        metricsBname = os.path.join(analysis.variants_dir, "1_gvcf/log/{SAMPLE}")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "03-00:00:00"
    shell:
        """
        gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
            -R {input.refFa} \
            -I {input.bam} \
            -L {input.intervals} \
            -O {output.vcf} \
            -dont-use-soft-clipped-bases \
            --standard-min-confidence-threshold-for-calling 20 \
            --emit-ref-confidence GVCF

        gatk \
            CollectVariantCallingMetrics \
            --DBSNP {input.refDbsnp} \
            --INPUT {output.vcf} \
            --OUTPUT {params.metricsBname}
        """

rule variants_genomicsDB:
    input:
        sample_map = "misc/cohort.sample_map",
        intervals = rules.intervals.output
    output:
        touch("done")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 4,
        ntasks = 2,
        mem_mb = 32000,
        time = "01-00:00:00"
    shell:
        """
        gatk --java-options "-Xmx4g -Xms4g" \
            GenomicsDBImport \
                --genomicsdb-workspace-path 07_variants/2_genomicsDB \
                --intervals {input.intervals} \
                --sample-name-map {input.sample_map} \
                --tmp-dir . \
                --merge-input-intervals
        """

rule variants_genotype:
    input:
        gendb_dir = "07_variants/2_genomicsDB",
        refFa = rules.refs_downloadFa.output,
    output:
        vcf = "07_variants/3_jointGenotype/output.vcf.gz",
        vcfIndex = "07_variants/3_jointGenotype/output.vcf.gz.tbi"
    params:
        gendb = "gendb://07_variants/2_genomicsDB"
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 32000,
        time = "01-00:00:00"
    shell:
        """
        gatk --java-options "-Xmx16g" GenotypeGVCFs \
            -R {input.refFa} \
            -V {params.gendb} \
            -O {output.vcf}
        """

rule variants_extract:
    input:
        vcf = rules.variants_genotype.output.vcf,
        vcfIndex = rules.variants_genotype.output.vcfIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
        refDbsnp = rules.refs_downloadDbsnp.output
    output:
        vcf = temp(os.path.join(analysis.variants_dir, "4_extract/output.vcf.gz")),
        vcfIndex = temp(os.path.join(analysis.variants_dir, "4_extract/output.vcf.gz.tbi")),
        detailMetrics = os.path.join(analysis.variants_dir, "4_extract/log/output.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join(analysis.variants_dir, "4_extract/log/output.variant_calling_summary_metrics")
    params:
        metricsBname = os.path.join(analysis.variants_dir, "4_extract/log/output")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 16000,
        time = "00-01:00:00"
    shell:
        """
        gatk \
            SelectVariants \
            -R {input.refFa} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            -O {output.vcf}

        gatk \
            CollectVariantCallingMetrics \
            --DBSNP {input.refDbsnp} \
            --INPUT {output.vcf} \
            --OUTPUT {params.metricsBname}
        """

rule variants_filter:
    input:
        vcf = rules.variants_extract.output.vcf,
        vcfIndex = rules.variants_extract.output.vcfIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
        refDbsnp = rules.refs_downloadDbsnp.output
    output:
        vcf = temp(os.path.join(analysis.variants_dir, "5_filter/output.vcf.gz")),
        vcfIndex = temp(os.path.join(analysis.variants_dir, "5_filter/output.vcf.gz.tbi")),
        detailMetrics = os.path.join(analysis.variants_dir, "5_filter/log/output.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join(analysis.variants_dir, "5_filter/log/output.variant_calling_summary_metrics")
    params:
        metricsBname = os.path.join(analysis.variants_dir, "5_filter/log/output")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 8000,
        time = "00-03:00:00"
    shell:
        """
        gatk \
            VariantFiltration \
            --R {input.refFa} \
            --V {input.vcf} \
            --window 35 \
            --cluster 3 \
            --filter-name "FS" --filter "FS > 60.0" \
            --filter-name "QD" --filter "QD < 2.0" \
            --filter-name "MQ" --filter "MQ < 40.0" \
            --filter-name "SOR" --filter "SOR > 4.0" \
            --filter-name "MQRankSum" --filter "MQRankSum < -12.5" \
            --filter-name "ReadPosRankSum" --filter "ReadPosRankSum < -8.0" \
            -O {output.vcf}

        gatk \
            CollectVariantCallingMetrics \
            --DBSNP {input.refDbsnp} \
            --INPUT {output.vcf} \
            --OUTPUT {params.metricsBname}
        """

rule variants_select:
    input:
        vcf = rules.variants_filter.output.vcf,
        vcfIndex = rules.variants_filter.output.vcfIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
        refDbsnp = rules.refs_downloadDbsnp.output
    output:
        vcf = os.path.join(analysis.variants_dir, "6_select/output.vcf.gz"),
        vcfIndex = os.path.join(analysis.variants_dir, "6_select/output.vcf.gz.tbi"),
        detailMetrics = os.path.join(analysis.variants_dir, "6_select/log/output.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join(analysis.variants_dir, "6_select/log/output.variant_calling_summary_metrics")
    params:
        metricsBname = os.path.join(analysis.variants_dir, "6_select/log/output")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-03:00:00"
    shell:
        """
        gatk \
            SelectVariants \
            --exclude-filtered \
            -R {input.refFa} \
            -V {input.vcf} \
            -O {output.vcf}

        gatk \
            CollectVariantCallingMetrics \
            --DBSNP {input.refDbsnp} \
            --INPUT {output.vcf} \
            --OUTPUT {params.metricsBname}
        """