rule raw_fastqc:
    input:
        "00_rawData/fastq/{SAMPLE}" + settings.fastq_ext
    output:
        "00_rawData/FastQC/{SAMPLE}_fastqc.zip",
        "00_rawData/FastQC/{SAMPLE}_fastqc.html"
    params:
        outDir = "00_rawData/FastQC/"
    conda:
        "../envs/gatk.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-01:00:00"
    shell:
        "fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"