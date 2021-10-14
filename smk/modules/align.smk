rule align:
    input:
        R1 = rules.trim.output.R1,
        R2 = rules.trim.output.R2,
        starIndex = rules.refs_starIndex.output
    output:
        bam = temp("02_align/bam/{SAMPLE}.bam"),
        bamIndex = temp("02_align/bam/{SAMPLE}.bam.bai"),
        STARgenome = temp(directory("02_align/bam/{SAMPLE}_STARgenome")),
        STARpass1 = temp(directory("02_align/bam/{SAMPLE}_STARpass1"))
    params:
        overhang = settings.read_length - 1,
        bname = "02_align/bam/{SAMPLE}",
        bamUnsorted = "02_align/bam/{SAMPLE}Aligned.out.bam"
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 16,
        ntasks = 1,
        mem_mb = 32000,
        time = "00-02:00:00"
    shell:
        """
        STAR \
            --genomeDir {input.starIndex}\
            --runThreadN {resources.cpu} \
            --readFilesIn {input.R1} {input.R2} \
            --readFilesCommand "gunzip -c" \
            --sjdbOverhang {params.overhang} \
            --outSAMtype BAM Unsorted \
            --twopassMode Basic \
            --outFileNamePrefix {params.bname}

        samtools sort {params.bamUnsorted} > {output.bam}
        samtools index {output.bam}
        rm {params.bamUnsorted}

        mkdir -p 02_align/log
        mv {params.bname}*out 02_align/log
        mv {params.bname}*tab 02_align/log
        """

rule align_fastqc:
    input:
        rules.align.output.bam
    output:
        "02_align/FastQC/{SAMPLE}_fastqc.zip",
        "02_align/FastQC/{SAMPLE}_fastqc.html"
    params:
        outDir = "02_align/FastQC/"
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-01:00:00"
    shell:
        "fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"