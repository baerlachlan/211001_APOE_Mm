import os

samples = os.listdir("../../a1647910/211001_APOE_Mm_ASE/00_rawData/fastq/")
samples = [sample.rstrip(".fastq.gz") for sample in samples]
# samples.rstrip(".fastq.gz")

print(samples)