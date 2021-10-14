rule bqsr_firstPass:
    input:
        bam = rules.splitNCigar.output.bam,
        bamIndex = rules.splitNCigar.output.bamIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
        refDbsnp = rules.refs_downloadDbsnp.output
    output:
        temp("06_bqsr/recal/{SAMPLE}.firstPass.table")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-04:00:00"
    shell:
        """
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R {input.refFa} \
            -I {input.bam} \
            -O {output} \
            --known-sites {input.refDbsnp}
        """

rule bqsr_apply:
    input:
        bam = rules.splitNCigar.output.bam,
        bamIndex = rules.splitNCigar.output.bamIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
        recal = rules.bqsr_firstPass.output
    output:
        bam = temp("06_bqsr/bam/{SAMPLE}.bam"),
        bamIndex = temp("06_bqsr/bam/{SAMPLE}.bai"),
        metrics = "06_bqsr/metrics/{SAMPLE}.tsv"
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-04:00:00"
    shell:
        """
        gatk \
            --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R {input.refFa} \
            -I {input.bam} \
            -O {output.bam} \
            --bqsr-recal-file {input.recal}

        samtools stats -d {output.bam} > {output.metrics}
        """

rule bqsr_secondPass:
    input:
        bam = rules.bqsr_apply.output.bam,
        bamIndex = rules.bqsr_apply.output.bamIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
        refDbsnp = rules.refs_downloadDbsnp.output
    output:
        temp("06_bqsr/recal/{SAMPLE}.secondPass.table")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-04:00:00"
    shell:
        """
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R {input.refFa} \
            -I {input.bam} \
            -O {output} \
            --known-sites {input.refDbsnp}
        """

rule bqsr_analyzeCovariates:
    input:
        firstPass = rules.bqsr_firstPass.output,
        secondPass = rules.bqsr_secondPass.output
    output:
        "06_bqsr/recal/{SAMPLE}.analyzeCovariates.csv"
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-01:00:00"
    shell:
        """
        gatk AnalyzeCovariates \
            -before {input.firstPass} \
             -after {input.secondPass} \
             -csv {output}
        """