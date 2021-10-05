rule trim:
    input:
        R1 = "00_rawData/fastq/{SAMPLE}" + config.read_pairs[0] + config.fastq_ext,
        R2 = "00_rawData/fastq/{SAMPLE}" + config.read_pairs[1] + config.fastq_ext
    output:
        R1 = temp("01_trim/fastq/{SAMPLE}" + config.read_pairs[0] + config.fastq_ext),
        R2 = temp("01_trim/fastq/{SAMPLE}" + config.read_pairs[1] + config.fastq_ext),
        html = "01_trim/log/{SAMPLE}.html"
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 2000,
        time = "00-02:00:00"
    shell:
        """
        fastp \
            -i {input.R1}  \
            -I {input.R2}  \
            -o {output.R1} \
            -O {output.R2} \
            --qualified_quality_phred 20 \
            --length_required 35 \
            --trim_poly_g \
            --thread 1 \
            --html {output.html} \
            --json /dev/null \
        """

rule trim_fastqc:
    input:
        "01_trim/fastq/{SAMPLE}" + config.fastq_ext
    output:
        "01_trim/FastQC/{SAMPLE}_fastqc.zip",
        "01_trim/FastQC/{SAMPLE}_fastqc.html"
    params:
        outDir = "01_trim/FastQC/"
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-01:00:00"
    shell:
        "fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"