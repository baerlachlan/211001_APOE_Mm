from snakemake.io import expand
import os
import pandas as pd

class ASEAnalysis:

    species = None
    assembly = None
    ensembl_release = None
    read_length = None
    fastq_ext = None
    read_pairs = None
    read_groups = None

    ref_fa = None
    ref_path = None
    ref_url = None
    gtf = None
    gtf_path = None
    gtf_url = None

    def __init__(self, species, assembly, ensembl_release, read_length, fastq_ext, read_pairs, rg_tsv="config/read_groups.tsv"):

        self.species = species
        self.assembly = assembly
        self.ensembl_release = ensembl_release
        self.read_length = read_length
        self.fastq_ext = fastq_ext
        self.read_pairs = read_pairs
        self.rg_tsv = rg_tsv

        ## Samples
        self.samples = os.listdir("00_rawData/fastq/")
        self.samples = [sample.rstrip(self.fastq_ext) for sample in self.samples]
        for id in self.read_pairs:
            self.samples = [sample.replace(id, "") for sample in self.samples]
        self.samples = list(set(self.samples))  ## Remove duplicates

        ## Read groups
        rg_df = pd.read_csv(rg_tsv, sep="\t")
        self.read_groups = dict(zip(rg_df["sample"], rg_df["read_group"]))

        ## Reference files
        self.ref_fa = os.path.join(".".join([self.species.capitalize(), self.assembly, "dna.primary_assembly.fa.gz"]))
        self.ref_path = os.path.join("refs", self.ref_fa.rstrip(".gz"))  ## Remove .gz because file will be unzipped during download
        self.ref_url = os.path.join("http://ftp.ensembl.org/pub", "release-" + str(self.ensembl_release), "fasta", self.species, "dna", self.ref_fa)
        self.gtf = os.path.join(".".join([self.species.capitalize(), self.assembly, str(self.ensembl_release), "chr.gtf.gz"]))
        self.gtf_path = os.path.join("refs", self.gtf)
        self.gtf_url = os.path.join("http://ftp.ensembl.org/pub", "release-" + str(self.ensembl_release), "gtf", self.species, self.gtf)
        self.dbsnp = os.path.join(".".join([self.species, "vcf.gz"]))
        self.dbsnp_path = os.path.join("refs", self.dbsnp)
        self.dbsnp_url = os.path.join("http://ftp.ensembl.org/pub", "release-" + str(self.ensembl_release), "variation/vcf", self.species, self.dbsnp)

        ## FQC outputs
        self.fqcRaw = expand("00_rawData/FastQC/{SAMPLE}{READID}_fastqc.{EXT}", SAMPLE=self.samples, READID=self.read_pairs, EXT=["html", "zip"])
        self.fqcTrim = expand("01_trim/FastQC/{SAMPLE}{READID}_fastqc.{EXT}", SAMPLE=self.samples, READID=self.read_pairs, EXT=["html", "zip"])
        self.fqcAlign = expand("02_align/FastQC/{SAMPLE}_fastqc.{EXT}", SAMPLE=self.samples, EXT=["html", "zip"])

        ## tmp output
        self.tmp = expand("05_splitNCigar/bam/{SAMPLE}.bam", SAMPLE=self.samples)

        ## All outputs
        self.outputs = [self.fqcRaw, self.fqcTrim, self.fqcAlign, self.dbsnp_path, "refs/exons.intervals"]

config = ASEAnalysis(
    species="mus_musculus",
    assembly="GRCm39",
    ensembl_release=104,
    read_length=50,
    fastq_ext=".fastq.gz",
    read_pairs=[".bam_R1", ".bam_R2"]
)

print(config.ensembl_release)
