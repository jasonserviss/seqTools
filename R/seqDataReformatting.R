
#' move.genes.to.rownames
#'
#' Moves the first column of the counts data.frame to rownames and removes
#' the old column.
#'
#' @name move.genes.to.rownames
#' @rdname move.genes.to.rownames
#' @param counts data.frame; A data frame with scRNAseq counts data.
#' @return The counts data.frame is returned with the first column moved to
#' the rownames of the data.frame.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(LETTERS, a = runif(26, 1, 100))
#' move.genes.to.rownames(counts)
#'
NULL
#' @export

move.genes.to.rownames <- function(counts) {
  if(!length(dim(counts))) {
    stop("dim(counts) must have a positive length.")
  }
  if(!is.data.frame(counts)) {
    stop("Counts is not a data.frame")
  }
  if(!is.character(counts[[1]])) {
    stop("The first column of counts is not class character.")
  }
  rownames(counts) <- counts[[1]]
  counts[[1]] <- NULL
  return(counts)
}

#' convert.counts.to.matrix
#'
#' Coerces the counts data.frame into a matrix and runs some basic checks before doing so.
#'
#' @name convert.counts.to.matrix
#' @rdname convert.counts.to.matrix
#' @param counts data.frame; A data frame with scRNAseq counts data.
#' @return The counts data.frame coerced into a matrix.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(a = runif(26, 1, 100), b = runif(26, 1, 100))
#' convert.counts.to.matrix(counts)
#'
NULL
#' @export

convert.counts.to.matrix <- function(counts) {
  if(class(counts) != "data.frame") {
    stop("Counts is not a data.frame")
  }
  if(!all(sapply(counts, class) %in% c("numeric", "integer"))) {
    stop("Non-numeric columns detected. Should gene names be moved to rownames?")
  }
  #check that the data.frame has rownames? gene names should be in the rownames
  as.matrix(counts)
}

#' label.singlets.and.multiplets
#'
#' Adds the prefix "s." to colnames of columns that contain singlets and "m." to
#' columns that contain multiplets.
#'
#' @name label.singlets.and.multiplets
#' @rdname label.singlets.and.multiplets
#' @param data data.frame or character; A data frame with sample IDs as colnames
#'  or a vector of sample IDs.
#' @param ids character; a character vector of regex statments used to indicate
#'  colnames of columns that contain SINGLETS.
#' @return If data is of class data.frame colnames are modified in place with
#' the appropriate prefix attatched. If data is of class character, the modified
#' character vector is returned. Throws a warning if no singlets are detected.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(a = runif(26), b = runif(26), c = runif(26))
#' label.singlets.and.multiplets(counts, ids = c("a", "b"))
#'
NULL
#' @export

label.singlets.and.multiplets <- function(data, ids) {
  .inputChecks_labelSingletsAndMultiplets(data, ids)

  if(class(data) == "data.frame") {

    .dataframeChecks_labelSingletsAndMultiplets(data)
    bool <- sapply(ids, function(x) grepl(x, colnames(data)))
    .boolChecks_labelSingletsAndMultiplets(bool)
    colnames(data) <- ifelse(
      rowSums(bool) > 0,
      paste("s.", colnames(data), sep = ""),
      paste("m.", colnames(data), sep = "")
    )
    return(data)

  } else {

    bool <- sapply(ids, function(x) grepl(x, data))
    .boolChecks_labelSingletsAndMultiplets(bool)
    data <- ifelse(
      rowSums(bool) > 0,
      paste("s.", data, sep = ""),
      paste("m.", data, sep = "")
    )
    return(data)
  }
}

.inputChecks_labelSingletsAndMultiplets <- function(data, ids) {
  if(class(ids) != "character") {
    stop("The id argument must be a character vector.")
  }
  if(!class(data) %in% c("data.frame", "character")) {
    stop("The data argument must be a character vector or data.frame.")
  }
  if(any(is.na(ids))) {
    stop("The ids argument contains NA values.")
  }
}

.dataframeChecks_labelSingletsAndMultiplets <- function(data) {
  if(is.null(colnames(data))) {
    stop("is.null(colnames(data)) returned TRUE.")
  }
  if(!all(sapply(data, class) %in% c("numeric", "integer"))) {
    stop("Non-numeric columns detected. Should gene names be moved to rownames?")
  }
}

.boolChecks_labelSingletsAndMultiplets <- function(bool) {
  if(sum(bool) == 0) {
    warning("No singlets were detected. Is your id argument valid?")
  }
}

#' remove.HTSEQ.suffix
#'
#' HTSeq adds the suffix ".htseq" to column names when it reports counts. This
#' function removes that suffix from the column names.
#'
#'
#' @name remove.HTSEQ.suffix
#' @rdname remove.HTSEQ.suffix
#' @param data data.frame or character; A data frame with sample IDs as colnames
#'  or a vector of sample IDs.
#' @return The counts data.frame with the ".htseq" suffix removed from the
#'  column names.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(a.htseq = runif(26), b.htseq = runif(26))
#' remove.HTSEQ.suffix(counts)
#'
NULL
#' @export

remove.HTSEQ.suffix <- function(data) {
  .inputChecks_removeHTSEQsuffix(data)
  if(class(data) == "data.frame") {
    .dataframeChecks_removeHTSEQsuffix(data)
    colnames(data) <- gsub("(.*)\\.htseq$", "\\1", colnames(data))
    return(data)
  }
  if(class(data) == "character") {
    data <- gsub("(.*)\\.htseq$", "\\1", data)
    return(data)
  }
}

.inputChecks_removeHTSEQsuffix <- function(data) {
  if(!class(data) %in% c("data.frame", "character")) {
    stop("The data argument must be a character vector or data.frame.")
  }
}

.dataframeChecks_removeHTSEQsuffix <- function(data) {
  if(is.null(colnames(data))) {
    stop("is.null(colnames(data)) returned TRUE.")
  }
  if(all(colnames(data) == paste0("V", 1:ncol(data)))) {
    warning("Your colnames are V1..Vi. Are these sample names?" )
  }
  if(!any(grepl("\\.htseq", colnames(data)))) {
    warning("Could not find the .htseq suffix in colnames(data)")
  }
}

#' detect.ERCC.reads
#'
#' Detects the 92 ERCC reads in the rownames of the counts data.frame. ERCC
#' reads must be named with the standard naming convention that matches the
#' regex "^ERCC\\-\\[0-9\\]*$".
#'
#' @name detect.ERCC.reads
#' @rdname detect.ERCC.reads
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row contains an ERCC read. A warning is issued if all 92
#' ERCC reads are not detected.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(runif(100), row.names = c(1:8, paste0("ERCC-", 1:92)))
#' detect.ERCC.reads(counts)
#'
NULL
#' @export

detect.ERCC.reads <- function(counts) {
  if(class(rownames(counts)) != "character") {
    stop("rownames(counts) is not of class character.")
  }
  if(all(rownames(counts) == as.character(1:nrow(counts)))) {
    m <- "rownames(counts) = 1:nrow(counts). Are gene names in rownames counts?"
    stop(m)
  }
  ercc <- grepl("^ERCC\\-[0-9]*$", rownames(counts))

  if(sum(ercc) != 92) {
    warning("Couldn't detect all 92 ERCC reads.")
  }

  return(ercc)
}

#' detect.non.genes
#'
#' HTSeq outputs several "non-genes" in the counts output. These include:
#' "__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned",
#' "__alignment_not_unique". Since these are measurements of the quantification
#' quality and are often undesirable in downstream analysis, this function
#' detects them in the submitted counts data.frame to allow easy removal.
#'
#' @name detect.non.genes
#' @rdname detect.non.genes
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row contains a non-gene.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(runif(10), row.names = c(1:9, "__no_feature"))
#' detect.non.genes(counts)
#'
NULL
#' @export

detect.non.genes <- function(counts) {
  if(class(rownames(counts)) != "character") {
    stop("rownames(counts) is not of class character.")
  }
  if(all(rownames(counts) == as.character(1:nrow(counts)))) {
    m <- "rownames(counts) = 1:nrow(counts). Are gene names in rownames counts?"
    stop(m)
  }

  nonGenes <- c(
    "__no_feature", "__ambiguous", "__too_low_aQual",
    "__not_aligned", "__alignment_not_unique"
  )
  rownames(counts) %in% nonGenes
}

#' Reformat wells
#'
#' Reformats wells named e.g. E2 to E02.
#'
#' @name reformat.well
#' @rdname reformat.well
#' @param names character; A character vector of the old names. Names that are
#' already in the new nomenclature will not be effected.
#' @return A character vector of the updated names.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom stringr str_replace str_extract

reformat.well <- function(names){
  base <- str_replace(names, "(.*)\\.[^\\.]*$", "\\1")
  pos <- paste0(".", str_extract(names, "[^.]+$"))
  posidx <- which(nchar(pos) == 3)
  pos[posidx] <- str_replace(pos[posidx], "(..)(.)", "\\10\\2")
  paste0(base, pos)
}
