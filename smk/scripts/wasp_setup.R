library(magrittr)
library(stringr)
library(dplyr)
library(readr)
library(SeqArray)

## If running interactively, edit the variables below
if (!exists("snakemake")) {
  proj_root <- "/hpcfs/users/a1647910/211001_APOE_Mm_ASE"
  variants_dir <- "07_variants"
  wasp_dir <- "08_wasp"
  ## Otherwise if running via snakemake, variables will be setup automatically
} else if (exists("snakemake")) {
  proj_root <- snakemake@params[["proj_root"]]
  variants_dir <- snakemake@params[["variants_dir"]]
  wasp_dir <- snakemake@params[["wasp_dir"]]
}

vcfPaths <- list.files(
  path = file.path(proj_root, variants_dir, "4_selected"),
  pattern = ".vcf.gz$",
  full.names = TRUE
)

lapply(vcfPaths, function(vcf){
  gdsPath <- paste0(
    str_remove(dirname(vcf), "4_selected"),
    "5_gds/",
    str_remove(basename(vcf), ".vcf.gz"),
    ".gds"
  )
  if (!dir.exists(dirname(gdsPath))) {
    dir.create(dirname(gdsPath), recursive = TRUE)
  }
  seqVCF2GDS(vcf, gdsPath)
})

gdsPaths <- list.files(
  path = file.path(proj_root, variants_dir, "5_gds"),
  pattern = ".gds$",
  full.names = TRUE
)

snvList <- lapply(seq_along(gdsPaths), function(x){
  gdsPath <- gdsPaths[x]
  sample <- basename(gdsPath) %>%
    str_remove(".gds")
  gds <- seqOpen(gdsPath, readonly = FALSE)
  snvs <- tibble(
    variant.id = seqGetData(gds, "variant.id"),
    chromosome = seqGetData(gds, "chromosome"),
    position = seqGetData(gds, "position"),
    allele = seqGetData(gds, "allele"),
    filter = seqGetData(gds, "annotation/filter")
  ) %>%
    mutate(
      sample = sample,
      refAllele = str_split(allele, ",", simplify = TRUE)[,1],
      altAllele = str_split(allele, ",", simplify = TRUE)[,2]
    )
  seqClose(gds)
  return(snvs)
})

lapply(snvList, function(x){
  chrList <- x %>%
    dplyr::select(chromosome, position, refAllele, altAllele, sample) %>%
    split(.[,'chromosome'])
  lapply(chrList, function(x){
    browser()
    chr <- unique(x$chromosome)
    sample <- unique(x$sample)
    path <- file.path(
      proj_root,
      wasp_dir,
      "1_snvs",
      sample,
      paste0(chr, ".snps.txt.gz")
    )
    if (!dir.exists(dirname(path))) {
      dir.create(dirname(path), recursive = TRUE)
    }
    tibble(
      position = x$position,
      refAllele = x$refAllele,
      altAllele = x$altAllele
    ) %>%
      write_delim(file = path, delim = " ", col_names = FALSE)
  })
})
