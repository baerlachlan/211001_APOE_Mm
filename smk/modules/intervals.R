library(magrittr)
library(AnnotationHub)
library(stringr)

species <- snakemake@params[["species"]] %>%
  str_replace("_", " ")
substr(species, 1, 1) <- str_to_upper(substr(species, 1, 1))

ensembl_release <- snakemake@params[["ensembl_release"]] %>%
  as.character()
ensembl_release <- 104

ah <- AnnotationHub() %>%
  subset(species == species) %>%
  subset(rdataclass == "EnsDb")
ahId <- ah$ah_id[str_detect(ah$title, ensembl_release)]
ensDb <- ah[[ahId]]
exons <- exonsBy(ensDb, by = "gene")

exonRanges <- exons %>%
  unlist() %>%
  GenomicRanges::reduce() %>%
  as.data.frame() %>%
  dplyr::filter(seqnames %in% drChrs) %>%
  dplyr::select(chromosome = seqnames, start, end)

intervalPath <- "../../refs/exons.intervals"
makeIntervals <- !file.exists(intervalPath)
if (makeIntervals) {
  if (!dir.exists(dirname(intervalPath))) {
    dir.create(dirname(intervalPath), recursive = TRUE)
  }
  tibble(
    interval = paste0(
      exonRanges$chromosome,
      ":",
      exonRanges$start,
      "-",
      exonRanges$end
    )
  ) %>%
    write.table(
      file = intervalPath,
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )
}