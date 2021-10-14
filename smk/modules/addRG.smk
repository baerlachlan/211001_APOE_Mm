## gatk HaplotypeCaller requires read group (RG) tags
## An explanation of this is found at https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
rule addRG:
    input:
        bam = rules.align.output.bam,
        bamIndex = rules.align.output.bamIndex
    output:
        bam = temp("03_addRG/bam/{SAMPLE}.bam"),
        bamIndex = temp("03_addRG/bam/{SAMPLE}.bai"),
        samstats = "03_addRG/samstats/{SAMPLE}.tsv"
    params:
        read_group = lambda wildcard: settings.read_groups[wildcard.SAMPLE]
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 4000,
        time = "00-02:00:00"
    shell:
        """
        gatk \
            AddOrReplaceReadGroups \
            -I {input.bam} \
            -O {output.bam} \
            -SORT_ORDER coordinate \
            -RGID {params.read_group} \
            -RGPU null \
            -RGSM {wildcards.SAMPLE} \
            -RGPL Illumina \
            -RGLB null \
            -CREATE_INDEX True

        samtools stats -d {output.bam} > {output.samstats}
        """