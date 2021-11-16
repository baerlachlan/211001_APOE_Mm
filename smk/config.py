import os
import pandas as pd
from snakemake.io import expand

class ASEAnalysis:

    def __init__(
            self, species, assembly, ensembl_release, read_length, fastq_ext,
            read_pairs, proj_root, rg_tsv="config/read_groups.tsv"):

        ## User defined parameters
        self.species = species
        self.assembly = assembly
        self.ensembl_release = ensembl_release
        self.read_length = read_length
        self.fastq_ext = fastq_ext
        self.read_pairs = read_pairs
        self.proj_root = proj_root
        self.rg_tsv = rg_tsv

        ## Run methods
        self.dirs()
        self.samples()
        self.references()
        self.targets()

    def dirs(self):
        ## Build directory structure
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

    def samples(self):
        ## Define samples based on raw data files
        self.samples = os.listdir(os.path.join(self.raw_dir, "fastq"))
        self.samples = [
            sample.rstrip(self.fastq_ext) for sample in self.samples
        ]
        for id in self.read_pairs:
            self.samples = [sample.replace(id, "") for sample in self.samples]
        self.samples = list(set(self.samples))  ## Remove duplicates after removal of read pair tags

    def references(self):
        ## Define reference filenames, paths and urls for downloading
        self.ref_fa = os.path.join(
            ".".join([
                self.species.capitalize(),
                self.assembly,
                "dna.primary_assembly.fa.gz"
            ])
        )
        self.ref_path = os.path.join(self.refs_dir, self.ref_fa.rstrip(".gz"))  ## Remove .gz because file will be unzipped during download
        self.ref_url = os.path.join(
            "http://ftp.ensembl.org/pub",
            "release-" + str(self.ensembl_release),
            "fasta",
            self.species,
            "dna",
            self.ref_fa
        )
        self.gtf = os.path.join(
            ".".join([
                self.species.capitalize(),
                self.assembly,
                str(self.ensembl_release),
                "chr.gtf.gz"
            ])
        )
        self.gtf_path = os.path.join(self.refs_dir, self.gtf)
        self.gtf_url = os.path.join(
            "http://ftp.ensembl.org/pub",
            "release-" + str(self.ensembl_release),
            "gtf",
            self.species,
            self.gtf
        )
        self.dbsnp = os.path.join(".".join([self.species, "vcf.gz"]))
        self.dbsnp_path = os.path.join(self.refs_dir, self.dbsnp)
        self.dbsnp_url = os.path.join(
            "http://ftp.ensembl.org/pub",
            "release-" + str(self.ensembl_release),
            "variation/vcf",
            self.species,
            self.dbsnp
        )

    def read_groups(self):
        rg_df = pd.read_csv(self.rg_tsv, sep="\t")
        self.read_groups = dict(zip(rg_df["sample"], rg_df["read_group"]))

    def targets(self):
        ## Define all file endpoints

        ## Intervals
        self.intervals = [os.path.join(self.refs_dir, "exons.intervals")]  ## Single element list to allow concatenation with other lists

        ## FastQC reports
        self.fqc_raw = expand(
            os.path.join(self.raw_dir, "FastQC/{SAMPLE}{READID}_fastqc.{EXT}"),
            SAMPLE=self.samples,
            READID=self.read_pairs,
            EXT=["html", "zip"]
        )
        self.fqc_trim = expand(
            os.path.join(self.trim_dir, "FastQC/{SAMPLE}{READID}_fastqc.{EXT}"),
            SAMPLE=self.samples,
            READID=self.read_pairs,
            EXT=["html", "zip"]
        )
        self.fqc_align = expand(
            os.path.join(self.align_dir, "FastQC/{SAMPLE}_fastqc.{EXT}"),
            SAMPLE=self.samples,
            EXT=["html", "zip"]
        )
        self.fqc = self.fqc_raw + self.fqc_trim + self.fqc_align

        ## BQSR tables and summary
        self.bqsr_firstPass = expand(
            os.path.join(self.bqsr_dir, "recal/{SAMPLE}.firstPass.table"),
            SAMPLE=self.samples
        )
        self.bqsr_secondPass = expand(
            os.path.join(self.bqsr_dir, "recal/{SAMPLE}.secondPass.table"),
            SAMPLE=self.samples
        )
        self.bqsr_analyzeCovariates = expand(
            os.path.join(self.bqsr_dir, "recal/{SAMPLE}.analyzeCovariates.csv"),
            SAMPLE=self.samples
        )
        self.bqsr = self.bqsr_firstPass + self.bqsr_secondPass + self.bqsr_analyzeCovariates

        ## Final variants
        self.variants = [os.path.join(self.variants_dir, "6_select/output.vcf.gz")]

        ## All outputs
        self.outputs = self.intervals + self.fqc + self.bqsr + self.variants
