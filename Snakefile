import os

from smk.modules.config import ASEAnalysis

## Do not name this object "config" as this is the default object name for snakemakes inbuilt configuration system
settings = ASEAnalysis(
    species="mus_musculus",
    assembly="GRCm39",
    ensembl_release=104,
    read_length=50,
    fastq_ext=".fastq.gz",
    read_pairs=[".bam_R1", ".bam_R2"],
    proj_root="/hpcfs/users/a1647910/211001_APOE_Mm_ASE/"
)

localrules: refs_downloadFa, refs_downloadGtf, refs_downloadDbsnp

rule all:
    input:
        settings.outputs

include: "smk/modules/refs.smk"
include: "smk/modules/intervals.smk"  ## Requires internet access so was run locally
include: "smk/modules/fastqc_raw.smk"
include: "smk/modules/trim.smk"
include: "smk/modules/align.smk"
include: "smk/modules/addRG.smk"
include: "smk/modules/markDuplicates.smk"
include: "smk/modules/splitNCigar.smk"
include: "smk/modules/bqsr.smk"
include: "smk/modules/variants.smk"
include: "smk/modules/wasp.smk"
