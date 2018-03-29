##works with cufflinks transcripts.gtf
#was written awhile ago, format may have changed.
#only keeps transcript annotations; can be adjusted on line 23.
colnames <- c(
    "seqname", "source", "feature", "start", "end",
    "score", "strand", "frame", "attributes"
)

colclasses <- c(
    rep("character", 3),
    rep("numeric", 2),
    rep("character", 4)
)

gtf <- read.table(
    "transcripts.gtf",
    header=FALSE,
    sep="\t",
    col.names = colnames,
    colClasses = colclasses
)

gtf <- subset(gtf,feature=="transcript")

renameAttributes <- function(gtf, regex, name) {
    c <- which(colnames(gtf) == name)
    gtf[,c] <- gsub(regex, "\\1", gtf$attributes)
    return(gtf)
}

regex <- "gene_id ([A-Z]*.[0-9]*\\.[0-9]*).*"
gtf <- renameAttributes(gtf, regex, "gene_id")

regex <- ".*gene_type ([A-Za-z0-9]*_?[A-Za-z0-9]*_?[A-Za-z0-9]*)\\;.*"
gtf <- renameAttributes(gtf, regex, "gene_type")

regex <- ".*gene_status ([A-Za-z]*)\\;.*"
gtf <- renameAttributes(gtf, regex, "gene_status")


regex <- ".*gene_name ([A-Za-z0-9]*[-_\\.]*?[A-Za-z0-9]*[-_\\.]*?[A-Za-z0-9]*[-_\\.]*?[A-Za-z0-9]*)\\;.*"
gtf <- renameAttributes(gtf, regex, "gene_name")

regex <- ".*transcript_id ([A-Z]*.[0-9]*).*"
gtf <- renameAttributes(gtf, regex, "transcript_id")

regex <- ".*exon_number ([0-9]*);.*"
gtf <- renameAttributes(gtf, regex, "exon_number")

regex <- ".*FPKM ([0-9]*.[0-9]*);.*"
gtf <- renameAttributes(gtf, regex, "FPKM")

regex <- ".*frac ([0-9]*.[0-9]*);.*"
gtf <- renameAttributes(gtf, regex, "frac")

regex <- ".*conf_lo ([0-9]*.[0-9]*);.*"
gtf <- renameAttributes(gtf, regex, "conf_lo")

regex <- ".*conf_hi ([0-9]*.[0-9]*);.*"
gtf <- renameAttributes(gtf, regex, "conf_hi")

regex <- ".*cov ([0-9]*.[0-9]*);.*"
gtf <- renameAttributes(gtf, regex, "cov")

gtf$attributes <- NULL
