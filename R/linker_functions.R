#' Convert IUPAC degenerate base codes to the full set of encoded bases
#'
#' This function is used internally to generate the rows of count matrix for
#' linker library sort-seq experiments. Each row will correspond to a unique
#' sequence encoded by the degenerate bases provided in \code{deg_codons} 
#' argument.
#' @param deg_codons List of characters indicating the bases used at each 
#' position including degenerate bases that encode multiple possible nucleotides. 
#' For example S corresponds to G or C.
#' @param parent Character string of the parent sequence used for cloning 
#' library. This argument allows for a single sequence not included in 
#' degenerate library to be included in processing steps in the event there is 
#' carryover of a single sequence during cloning process.
#' @return Dataframe containing single column with characters for all sequences 
#' encoded by degenerate bases.
#' @examples
#' \dontrun{
#' deg_codons <- c(rep(c("S", "N", "A"), 2), rep(c("V", "S", "T"), 2))
#' codonMap(deg_codons = deg_codons, parent = NA)
#' }
#' @keywords internal
codonMap <- function(deg_codons, parent = NA) {
  tmp <- Biostrings::IUPAC_CODE_MAP[deg_codons]

  tmp2 <- lapply(tmp, function(v) {
    as.vector(stringr::str_split_fixed(v, pattern = "", n = nchar(v)))
  })

  dat <- expand.grid(tmp2)
  seq <- apply(dat, 1, paste, collapse='')

  if(!is.na(parent)) {
    seq[length(seq)+1] <- parent
  }

  seq <- seq[order(seq)]
  out <- data.frame(seq,
                    stringsAsFactors = F)

  return(out)
}

#' Trim the constant sequences from the sequencing reads leaving only the 
#' variable regions of interest.
#'
#' This function is used internally to process sequencing reads by removing the
#' constant regions in the DNA that were not mutated in order to produce only
#' the variable regions that are of interest. Cutadapt 
#' (https://cutadapt.readthedocs.io/) must be installed and the full path 
#' supplied.
#' @param fl Filename indicating the location of sequencing reads. For 
#' paired-end reads only the first (R1) file needs to be included.
#' @param leftseq Character indicating the constant sequence on the 5' side of 
#' variable sequence. N can be used to indicate any base is allowed at certain
#' positions.
#' @param rightseq Character indicating the constant sequence on the 3' side 
#' of variable sequence. N can be used to indicate any base is allowed at 
#' certain positions.
#' @param linklen Integer indicating the number of basepairs that comprise the
#' mutated linker sequence (or mutated sequence in general).
#' @param cutadapt Character indicating the full path to installation of the 
#' software cutadapt (https://cutadapt.readthedocs.io/)
#' @return No return. Writes of trimmed fastq files that are deleted after use.
#' @examples
#' \dontrun{
#' fl <- "inst/extdata/naive_R1.fastq.gz"
#' leftseq <- "^NNNNNACGTGAGACAGAATTTTGAGCTCCTG"
#' rightseq <- "^NNNNNTAGACACGAGTGGCAGCATTTCGCGCCGAGAATA"
#' cutadapt <- "/home/jnk/software/anaconda3/bin/cutadapt"
#' trimReads(fl, leftseq, rightseq, linklen = 6, cutadapt)
#' }
#' @keywords internal
#' @export
trimReads <- function(fl, leftseq, rightseq, linklen, cutadapt) {
  #fl1 <- gsub("#", "1", file)
  fl1 <- fl
  fl2 <- gsub("R1", "R2", fl)

  outfile <- gsub("R1.*", "trimmed.fastq", fl)

  cmd <- paste(cutadapt = cutadapt,
               "-o", outfile,
               "-j 0",
               "--interleaved",
               "--no-indels",
               "-e 0.2",
               "-g", leftseq,
               "-G", rightseq,
               "-l", linklen,
               fl1,
               fl2)

  system(cmd)
}

#' Count the number of times each unique sequence is present in each fastq file.
#'
#' This function is used internally to count the number of times each unique
#' sequence is present in the trimmed read files.
#' @param fl Filename indicating the location of sequencing reads. For 
#' paired-end reads only the first (R1) file needs to be included.
#' @param deg_codons List of characters indicating the bases used at each 
#' position of sequence. Only sequences encoded by the provided bases including
#' degenerate bases will be counted.
#' @param parent Character string of the parent sequence used for cloning 
#' library. This argument allows for a single sequence not included in 
#' degenerate library to be included in processing steps in the event there is 
#' carryover of a single sequence during cloning process.
#' @return No return. Writes files containing read counts for each sequence in
#' each file that can later be concatenated using read.Linker.Counts.
#' @examples
#' \dontrun{
#' fl <- "inst/extdata/PyronicSF_example.fastq.gz"
#' deg_codons <- c(rep(c("S", "N", "A"), 2), rep(c("V", "S", "T"), 2))
#' countReads(fl = fl, deg_codons = deg_codons)
#' }
#' @keywords internal
countReads <- function(fl, deg_codons, parent = NA) {
  fl <- gsub("R1.*", "trimmed.fastq", fl)
  tmp <- gsub("trimmed.*", "counts.txt", fl)
  out <- gsub("/fastq/", "/counts/", tmp)

  counts <- codonMap(deg_codons, parent)
  counts$counts <- rep(0, nrow(counts))
  stream <- ShortRead::FastqStreamer(fl)

  repeat {
    fq <- ShortRead::yield(stream)

    if (length(fq) == 0)
      break

    ind <- seq(1, length(fq)-1, by = 2)

    left <- fq@sread[ind,]
    right <- Biostrings::reverseComplement(fq@sread[ind + 1,])

    linker <- paste0(left, right)

    ldat <- plyr::count(linker)
    colnames(ldat) <- c("seq", "counts")

    ldat <- ldat[ldat$seq %in% counts$seq,]
    ldat <- ldat[order(ldat$seq),]

    ind <- counts$seq %in% ldat$seq
    counts[ind, 2] <- counts[ind, 2] + ldat$counts

    }
  close(stream)
  invisible(file.remove(fl))
  write.table(counts, out, quote = F, sep = "\t", row.names = F)
}

#' Trim and count the variable sequences in fastq read files
#'
#' This function is used to process raw sequencing reads to trim the constant
#' regions. The bases encoded by indicated degenerate bases should remain
#' after trimming which are then counted and saved to a new counts file. This 
#' function should only need to be run one time to produce saved fils can be
#' reused by read.Linker.Counts().
#' @param fnames Filename indicating the location of sequencing reads. For 
#' paired-end reads only the first (R1) file needs to be included.
#' @param deg_codons List of characters indicating the bases used at each 
#' position of sequence. Only sequences encoded by the provided bases including
#' degenerate bases will be counted.
#' @param leftseq Character indicating the constant sequence on the 5' side of 
#' variable sequence. N can be used to indicate any base is allowed at certain
#' positions.
#' @param rightseq Character indicating the constant sequence on the 3' side 
#' of variable sequence. N can be used to indicate any base is allowed at 
#' certain positions.
#' @param linklen Integer indicating the number of basepairs that comprise the
#' mutated linker sequence (or mutated sequence in general).
#' @param cutadapt Character indicating the full path to installation of the 
#' software cutadapt (https://cutadapt.readthedocs.io/)
#' @param parent Character string of the parent sequence used for cloning 
#' library. This argument allows for a single sequence not included in 
#' degenerate library to be included in processing steps in the event there is 
#' carryover of a single sequence during cloning process.
#' @return No return. Writes files containing trimmed reads and read counts for 
#' each sequence in each file that can later be concatenated using 
#' the function read.Linker.Counts().
#' @examples
#' \dontrun{
#' fnames <- "inst/extdata/PyronicSF_example_R1.fastq.gz"
#' leftseq <- "^NNNNNACGTGAGACAGAATTTTGAGCTCCTG"
#' rightseq <- "^NNNNNTAGACACGAGTGGCAGCATTTCGCGCCGAGAATA"
#' deg_codons <- c(rep(c("S", "N", "A"), 2), rep(c("V", "S", "T"), 2))
#' cutadapt <- "/home/jnk/software/anaconda3/bin/cutadapt"
#' countLinkers(fnames = fnames, deg_codons = deg_codons, 
#' leftseq = leftseq, rightseq = rightseq, linklen = 6, cutadapt = cutadapt)
#' }
#' @export
countLinkers <- function(fnames, deg_codons, leftseq, rightseq, linklen, 
                          cutadapt, parent = NA) {
  
  invisible(lapply(fnames, function(fl) {
    trimReads(fl, leftseq, rightseq, linklen, cutadapt = cutadapt)
    countReads(fl, deg_codons, parent)
  }))
}

#' Read linker counts from multiple files into a single dataframe
#' 
#' This function is used to concatenate linker read count files for each sample
#' into a single dataframe for further analysis.
#' @param count.dir Character indicating the directory containing linker read 
#' count files.
#' @return Dataframe containing columns corresponding to samples, rows 
#' corresponding to linker sequences and entries corresponding to read counts.
#' @export
read.Linker.Counts <- function(count.dir) {
  fullnames <- list.files(count.dir, full.names = T)
  fnames <- tools::file_path_sans_ext(basename(fullnames))

  fnames <- vapply(strsplit(fnames, "_"), `[`, 1, FUN.VALUE=character(1))

  seqs <- Reduce(intersect, lapply(fullnames,
                                   function(f) {read.table(f, sep = "\t",
                                                           header = F)[-1,1]}))

  cdf <- as.data.frame(matrix(data = 0, nrow=length(seqs), ncol=length(fnames)),
                    row.names = seqs,
                    stringsAsFactors = F)
  colnames(cdf) <- fnames

  for (fname in fnames){
    fl <- fullnames[grep(fname, fnames)]
    ind <- grep(fname, colnames(cdf))
    tmp <- read.table(fl, sep = "\t", header = F, stringsAsFactors = F)
    tmp <- tmp[tmp[,1] %in% seqs,]
    tmp <- tmp[match(tmp$V1, rownames(cdf)),]
    cdf[,ind] <- as.numeric(tmp$V2)
  }
  return(cdf)
}
