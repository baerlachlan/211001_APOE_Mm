rule trim:
    input:
        R1 = os.path.join(analysis.raw_dir, "fastq/{SAMPLE}" + analysis.read_pairs[0] + analysis.fastq_ext),
        R2 = os.path.join(analysis.raw_dir, "fastq/{SAMPLE}" + analysis.read_pairs[1] + analysis.fastq_ext)
    output:
        R1 = temp(os.path.join(analysis.trim_dir, "fastq/{SAMPLE}" + analysis.read_pairs[0] + analysis.fastq_ext)),
        R2 = temp(os.path.join(analysis.trim_dir, "fastq/{SAMPLE}" + analysis.read_pairs[1] + analysis.fastq_ext)),
        html = os.path.join(analysis.trim_dir, "log/{SAMPLE}.html")
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
        os.path.join(analysis.trim_dir, "fastq/{SAMPLE}" + analysis.fastq_ext)
    output:
        os.path.join(analysis.trim_dir, "FastQC/{SAMPLE}_fastqc.zip"),
        os.path.join(analysis.trim_dir, "FastQC/{SAMPLE}_fastqc.html")
    params:
        outDir = os.path.join(analysis.trim_dir, "FastQC")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-01:00:00"
    shell:
        "fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"