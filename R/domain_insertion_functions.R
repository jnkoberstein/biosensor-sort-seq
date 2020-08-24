#' Generate a metadata file to be used for reading in DIP count files
#' 
#' This function is used to generate a metadata file for annotating the 
#' insertion sites in a transposon-mediate domain-insertion library. Each 
#' insertion site is described by the 
#' @param seq.length Integer corresponding to the to total number of rows in 
#' the DIP count files. This is approximately two times the length of DNA 
#' sequence that was transposed to account for forward and reverse insertions.
#' @param AA.start Integer indication the insertion site that corresponds to an 
#' insertion of cpGFP that is forward and in-frame in the first amino acid of 
#' ligand-binding domain.
#' @param AA.length Integer indicating the total number of amino acids in the
#' ligand-binding domain.
#' @param fname Filename for saving the metadata file produced as output.
#' @return Dataframe containing metadata describing mapping between insertion 
#' sites and amino acid sequence. The file will also be written out to filename
#' indicated by \code{fname}
#' @details metadata columns:
#' * insertion_site_name: Character string indicating the site of insertion
#' and direction.
#' * site: Integer indicating the number of DNA bases before the site of 
#' cpFP insertion.
#' * forward_insertion: Logical indicating whether the insertion is in the
#' forward direction (TRUE) or the reverse direction (FALSE).
#' * in_frame: Logical indicating whether the insertion site is in-frame (TRUE) 
#' or out-of-frame (FALSE). Only insertion after second nucleotide in codon are
#' result in an in-frame translation of inserted domain. 
#' * productive: Logical indicating whether the insertion is in the forward
#' direction and in-frame (TRUE) or either reversed or out-of-frame (FALSE).
#' * AA: Integer indicating which amino acid in the LBD the insertion of cpFP
#' occurs after. All unproductive insertions are set to AA=0.
#' @md
#' @examples
#' \dontrun{
#' # generating metadata for MBP domain-insertion profiling experiment
#' metaFile(seq.length = 1120, AA.start = 8, AA.length = 370, 
#' fname = "data/meta.txt")
#' }
#' @export
metaFile <- function(seq.length, AA.start, AA.length, fname) {
  insertion_site_name <- c(paste0(1:seq.length, "_fwd"),
                            paste0(1:seq.length, "_rev"))
  site <- rep(1:seq.length, 2)
  forward_insertion <- c(rep(TRUE, seq.length), rep(FALSE, seq.length))
  in_frame <- rep(FALSE, seq.length*2)
  in_frame[site%%3==AA.start%%3 & site >= AA.start] <- TRUE
  productive <- forward_insertion & in_frame
  AA <- rep(0, seq.length*2)
  metadata <- data.frame(insertion_site_name, site, forward_insertion, in_frame,
                     productive, AA)

  metadata[seq(AA.start, (AA.length*3)+AA.start-3, by = 3), "AA"] <- 1:AA.length

  write.table(metadata, fname, quote = F, col.names = T, row.names = F, sep = "\t")
  return(metadata)
}

#' Read domain-insertion profiling count files into a dataframe
#'
#' This function provides a means to read the count files produced by the python
#' dipseq package (https://github.com/SavageLab/dipseq) into a single dataframe
#' for sort-seq analysis.
#' @param fnames List of filenames for read count tables to be read
#' @param meta.file Filename for metadata file describing mapping between insertion 
#' sites and amino acid sequence.
#' @return Matrix with the first n columns corresponding to 
#' the n files/samples provided by \code{fnames}, rows corresponding to 
#' insertion sites and entries containing read counts. The last 6 columns 
#' provide insertion site annotations from the metadata file.
#' @examples
#' \dontrun{
#' read.DIP.Counts(fnames = "example.csv", meta.file = "meta.txt")
#' }
#' @export
read.DIP.Counts <- function(fnames, meta.file) {
  metadata <- read.table(meta.file, header = T, sep = "\t")
  numrow <- nrow(metadata)
  cdf <- as.data.frame(matrix(data = 0, nrow = numrow, ncol = length(fnames)))
  tmp <- unlist(lapply(fnames, function(x)
  {strsplit(basename(x), split = "_", fixed = T)[[1]][1]}))
  colnames(cdf) <- tmp
  rownames(cdf) <- metadata$insertion_site_name

  for (i in 1:length(fnames)) {
    dat <- read.table(fnames[i], sep = ",", header = T)
    dat <- dat[dat$insertion_site_name %in% metadata$insertion_site_name,]
    dat$sum <- rowSums(dat[,2:3])
    cdf[match(dat$insertion_site_name, row.names(cdf)),i] <- dat$sum
  }

  cdf <- cbind(cdf, metadata)

  return(cdf)
}
