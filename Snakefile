from smk.modules.config import ASEAnalysis


config = ASEAnalysis(
    species="mus_musculus",
    assembly="GRCm39",
    ensembl_release=104,
    read_length=50,
    fastq_ext=".fastq.gz",
    read_pairs=[".bam_R1", ".bam_R2"]
)

localrules: refs_downloadFa, refs_downloadGtf, refs_downloadDbsnp

rule all:
    input:
        config.outputs

include: "smk/modules/refs.smk"
include: "smk/modules/intervals.smk"
include: "smk/modules/fastqc_raw.smk"
include: "smk/modules/trim.smk"
include: "smk/modules/align.smk"
include: "smk/modules/addRG.smk"
include: "smk/modules/markDuplicates.smk"
include: "smk/modules/splitNCigar.smk"
