import os
import pandas as pd
from snakemake.io import expand

class ASEAnalysis:

    def __init__(self, species, assembly, ensembl_release, read_length, fastq_ext, read_pairs, proj_root, rg_tsv="config/read_groups.tsv"):

        ## User defined parameters
        self.species = species
        self.assembly = assembly
        self.ensembl_release = ensembl_release
        self.read_length = read_length
        self.fastq_ext = fastq_ext
        self.read_pairs = read_pairs
        self.proj_root = proj_root
        self.rg_tsv = rg_tsv

        ## Run methods required for setup
        self.dirs()
        self.samples()
        self.references()
        self.targets()

    def dirs(self):
        self.refs_dir = "refs"
        self.raw_dir = "00_rawData"
        self.trim_dir = "01_trim"
        self.align_dir = "02_align"
        self.addRG_dir = "03_addRG"
        self.markDuplicates_dir = "04_markDuplicates"
        self.splitNCigar_dir = "05_splitNCigar"
        self.bqsr_dir = "06_bqsr"
        self.variants_dir = "07_variants"
        self.wasp_dir = "08_wasp"
        self.aseRC_dir = "09_aseRC"
        self.geneiase_dir = "10_geneiase"

    def samples(self):
        self.samples = os.listdir("00_rawData/fastq/")
        self.samples = [sample.rstrip(self.fastq_ext) for sample in self.samples]
        for id in self.read_pairs:
            self.samples = [sample.replace(id, "") for sample in self.samples]
        self.samples = list(set(self.samples))  ## Remove duplicates

    def references(self):
        self.ref_fa = os.path.join(".".join([self.species.capitalize(), self.assembly, "dna.primary_assembly.fa.gz"]))
        self.ref_path = os.path.join(self.refs_dir, self.ref_fa.rstrip(".gz"))  ## Remove .gz because file will be unzipped during download
        self.ref_url = os.path.join("http://ftp.ensembl.org/pub", "release-" + str(self.ensembl_release), "fasta", self.species, "dna", self.ref_fa)
        self.gtf = os.path.join(".".join([self.species.capitalize(), self.assembly, str(self.ensembl_release), "chr.gtf.gz"]))
        self.gtf_path = os.path.join(self.refs_dir, self.gtf)
        self.gtf_url = os.path.join("http://ftp.ensembl.org/pub", "release-" + str(self.ensembl_release), "gtf", self.species, self.gtf)
        self.dbsnp = os.path.join(".".join([self.species, "vcf.gz"]))
        self.dbsnp_path = os.path.join(self.refs_dir, self.dbsnp)
        self.dbsnp_url = os.path.join("http://ftp.ensembl.org/pub", "release-" + str(self.ensembl_release), "variation/vcf", self.species, self.dbsnp)

    def read_groups(self):
        rg_df = pd.read_csv(self.rg_tsv, sep="\t")
        self.read_groups = dict(zip(rg_df["sample"], rg_df["read_group"]))

    def targets(self):
        ## FQC outputs
        self.fqc_raw = expand("00_rawData/FastQC/{SAMPLE}{READID}_fastqc.{EXT}", SAMPLE=self.samples, READID=self.read_pairs, EXT=["html", "zip"])
        self.fqc_trim = expand("01_trim/FastQC/{SAMPLE}{READID}_fastqc.{EXT}", SAMPLE=self.samples, READID=self.read_pairs, EXT=["html", "zip"])
        self.fqc_align = expand("02_align/FastQC/{SAMPLE}_fastqc.{EXT}", SAMPLE=self.samples, EXT=["html", "zip"])
        self.fqc = self.fqc_raw + self.fqc_trim + self.fqc_align
        ## BQSR tables
        self.bqsr_firstPass = expand("06_bqsr/recal/{SAMPLE}.firstPass.table", SAMPLE=self.samples)
        self.bqsr_secondPass = expand("06_bqsr/recal/{SAMPLE}.secondPass.table", SAMPLE=self.samples)
        self.bqsr_analyzeCovariates = expand("06_bqsr/recal/{SAMPLE}.analyzeCovariates.csv", SAMPLE=self.samples)
        self.bqsr = self.bqsr_firstPass + self.bqsr_secondPass + self.bqsr_analyzeCovariates
        ## tmp output
        self.tmp = expand("08_wasp/5_merge/{SAMPLE}.keep.merge.sort.bam", SAMPLE=self.samples)
        self.tmp = sorted(self.tmp)
        self.tmp = self.tmp[0:12]
        ## All outputs
        self.outputs = self.fqc + self.bqsr + self.tmp + ["refs/exons.intervals"]

# settings = ASEAnalysis(
#     species="mus_musculus",
#     assembly="GRCm39",
#     ensembl_release=104,
#     read_length=50,
#     fastq_ext=".fastq.gz",
#     read_pairs=[".bam_R1", ".bam_R2"]
# )

# settings.all_outputs()

# print(settings.fqc)
