import os

class ASEAnalysis:

    species = None
    assembly = None
    ensembl_release = None
    read_length = None
    fastq_ext = None

    ref_fa = None
    ref_path = None
    ref_url = None
    gtf = None
    gtf_path = None
    gtf_url = None

    def __init__(self, species=None, assembly = None, ensembl_release=None, read_length=None, fastq_ext=None):

        self.species = species
        self.assembly = assembly
        self.ensembl_release = ensembl_release
        self.read_length = read_length
        self.fastq_ext = fastq_ext

        # ## Samples
        # self.samples = os.listdir("../../00_rawData/fastq")
        # self.samples.rstrip(self.fastq_ext)

        ## Reference data
        self.ref_fa = os.path.join(".".join([self.species.capitalize(), self.assembly, "dna.primary_assembly.fa.gz"]))
        self.ref_path = os.path.join("refs", self.ref_fa.rstrip(".gz"))  ## Remove .gz because file will be unzipped during download
        self.ref_url = os.path.join("http://ftp.ensembl.org/pub", "release-" + str(self.ensembl_release), "fasta", self.species, "dna", self.ref_fa)
        self.gtf = os.path.join(".".join([self.species.capitalize(), self.assembly, str(self.ensembl_release), "chr.gtf.gz"]))
        self.gtf_path = os.path.join("refs", self.gtf)
        self.gtf_url = os.path.join("http://ftp.ensembl.org/pub", "release-" + str(self.ensembl_release), "gtf", self.species, self.gtf)

        self.outputs = [self.ref_path]

# config = ASEAnalysis(
#     species="mus_musculus",
#     assembly="GRCm39",
#     ensembl_release=104,
#     read_length=50,
#     fastq_ext=".fastq.gz"
# )

# print(config.samples)