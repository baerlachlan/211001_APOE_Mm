rule raw_fastqc:
    input:
        os.path.join(analysis.raw_dir, "fastq/{SAMPLE}" + analysis.fastq_ext)
    output:
        os.path.join(analysis.raw_dir, "FastQC/{SAMPLE}_fastqc.zip"),
        os.path.join(analysis.raw_dir, "00_rawData/FastQC/{SAMPLE}_fastqc.html")
    params:
        outDir = os.path.join(analysis.raw_dir, "FastQC/")
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-01:00:00"
    shell:
        "fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"