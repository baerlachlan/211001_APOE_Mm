rule variants_call:
    input:
        bam = rules.bqsr_apply.output.bam,
        bamIndex = rules.bqsr_apply.output.bamIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
        refDbsnp = rules.refs_downloadDbsnp.output,
        intervals = rules.intervals.output
    output:
        vcf = temp(os.path.join(analysis.variants_dir, "1_called/{SAMPLE}.vcf.gz")),
        vcfIndex = temp(os.path.join(analysis.variants_dir, "1_called/{SAMPLE}.vcf.gz.tbi")),
        detailMetrics = os.path.join(analysis.variants_dir, "1_called/log/{SAMPLE}.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join(analysis.variants_dir, "1_called/log/{SAMPLE}.variant_calling_summary_metrics")
    params:
        metricsBname = os.path.join(analysis.variants_dir, "1_called/log/{SAMPLE}")
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
            --standard-min-confidence-threshold-for-calling 20

        gatk \
            CollectVariantCallingMetrics \
            --DBSNP {input.refDbsnp} \
            --INPUT {output.vcf} \
            --OUTPUT {params.metricsBname}
        """

## We are only interested in Single Nucleotide Variants for ASE analysis
rule variants_extract:
    input:
        vcf = rules.variants_call.output.vcf,
        vcfIndex = rules.variants_call.output.vcfIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
        refDbsnp = rules.refs_downloadDbsnp.output
    output:
        vcf = temp(os.path.join(analysis.variants_dir, "2_extracted/{SAMPLE}.vcf.gz")),
        vcfIndex = temp(os.path.join(analysis.variants_dir, "2_extracted/{SAMPLE}.vcf.gz.tbi")),
        detailMetrics = os.path.join(analysis.variants_dir, "2_extracted/log/{SAMPLE}.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join(analysis.variants_dir, "2_extracted/log/{SAMPLE}.variant_calling_summary_metrics")
    params:
        metricsBname = os.path.join(analysis.variants_dir, "2_extracted/log/{SAMPLE}")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-12:00:00"
    shell:
        """
        gatk \
            SelectVariants \
            -R {input.refFa} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            --restrict-alleles-to BIALLELIC \
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
        vcf = temp(os.path.join(analysis.variants_dir, "3_filtered/{SAMPLE}.vcf.gz")),
        vcfIndex = temp(os.path.join(analysis.variants_dir, "3_filtered/{SAMPLE}.vcf.gz.tbi")),
        detailMetrics = os.path.join(analysis.variants_dir, "3_filtered/log/{SAMPLE}.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join(analysis.variants_dir, "3_filtered/log/{SAMPLE}.variant_calling_summary_metrics")
    params:
        metricsBname = os.path.join(analysis.variants_dir, "3_filtered/log/{SAMPLE}")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 8000,
        time = "00-00:30:00"
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
        vcf = os.path.join(analysis.variants_dir, "4_selected/{SAMPLE}.vcf.gz"),
        vcfIndex = os.path.join(analysis.variants_dir, "4_selected/{SAMPLE}.vcf.gz.tbi"),
        detailMetrics = os.path.join(analysis.variants_dir, "4_selected/log/{SAMPLE}.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join(analysis.variants_dir, "4_selected/log/{SAMPLE}.variant_calling_summary_metrics")
    params:
        metricsBname = os.path.join(analysis.variants_dir, "4_selected/log/{SAMPLE}")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-00:30:00"
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