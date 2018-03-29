#how to read a .sam file into R.
#does this work? Untested.

file <- 'accepted_hits.sam'
colnames <- c(
    "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10", "v11", "v12",
    "v13", "v14", "v15", "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23",
    "v24"
)

sam <- read.table(
    file,
    header = FALSE,
    sep = "",
    quote = "'",
    dec = ".",
    col.names=colnames,
    fill=T
)

drops <- c(
    "v11","v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22",
    "v23","v24"
)

sam_2 <- sam[,!(names(sam) %in% drops)]

colnames(sam_2) <- c(
    "read_name", "sum_of_flags", "chr", "start", "mapping_quality", "CIGAR",
    "mate_alignment_chr", "mate_alignment_start", "inferred_frag_length",
    "sequence"
)
