from smk.modules.config import ASEAnalysis


config = ASEAnalysis(
    species="mus_musculus",
    assembly="GRCm39",
    ensembl_release=104,
    read_length=50,
    fastq_ext=".fastq.gz"
)

rule all:
	input:
		config.outputs

include: "smk/modules/newrefs.smk"
