rule align:
    input:
        R1 = rules.trim.output.R1,
        R2 = rules.trim.output.R2,
        starIndex = rules.refs_starIndex.output
    output:
        bam = temp(os.path.join(analysis.align_dir, "bam/{SAMPLE}.bam")),
        bamIndex = temp(os.path.join(analysis.align_dir, "bam/{SAMPLE}.bam.bai")),
        STARgenome = temp(directory(os.path.join(analysis.align_dir, "bam/{SAMPLE}_STARgenome"))),
        STARpass1 = temp(directory(os.path.join(analysis.align_dir, "bam/{SAMPLE}_STARpass1")))
    params:
        overhang = analysis.read_length - 1,
        bname = os.path.join(analysis.align_dir, "bam/{SAMPLE}"),
        bamUnsorted = os.path.join(analysis.align_dir, "bam/{SAMPLE}Aligned.out.bam"),
        align_dir = analysis.align_dir
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

        mkdir -p {params.align_dir}/log
        mv {params.bname}*out {params.align_dir}/log
        mv {params.bname}*tab {params.align_dir}/log
        """

rule align_fastqc:
    input:
        rules.align.output.bam
    output:
        os.path.join(analysis.align_dir, "FastQC/{SAMPLE}_fastqc.zip"),
        os.path.join(analysis.align_dir, "FastQC/{SAMPLE}_fastqc.html")
    params:
        outDir = os.path.join(analysis.align_dir, "FastQC")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-01:00:00"
    shell:
        "fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"