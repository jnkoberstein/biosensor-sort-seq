#' Log-likelihood function for binned sort-seq
#' 
#'
#' @param par parameters to be optimized
#' @param ri read counts for variant i
#' @param lj lower boundaries for j gates
#' @param uj upper boundaries for j gates
#' @param dj proportionality factor describing the total number of cell sorted 
#' to the total number of read counts in each bin. 
#' @param eps error rate for sorting
#' @return integer representing function evaluated at tested parameters
#' @keywords internal  
ll_func <- function(par, ri, lj, uj, dj, eps) {
  mu.t <- par[1]
  sig.t <- par[2]
  if (sig.t < 0) return(-Inf)
  u.prob <- eps + (1 - sum(eps)) * pnorm(uj, mu.t, sig.t)
  l.prob <- eps + (1 - sum(eps)) * pnorm(lj, mu.t, sig.t)
  tmp <- dj * ri * log10(u.prob - l.prob)
  tmp[is.nan(tmp)] <- 0
  out <- sum(unlist(tmp))
  return(out)
}

#' Maximum-likelihood estimates for binned sort-seq using 
#' log-normal distribution
#' 
#'
#' @param ri read counts for variant i
#' @param mu mean of fluorescence distribution
#' @param sig standard deviation of fluorescence distribution
#' @param dj proportionality factor for each bin describing the total number 
#' of cells sorted in bin j to the total number of read counts in bin j. 
#' @param eps error rates for sorting
#' @return Dataframe containing maximum likelihood estimates of mean and 
#' standard deviation for each variant i.
#' @keywords internal
lnorm_mle <- function(ri, mu, sig, gates, dj, eps) {
  lj <- gates[-length(gates)]
  uj <- gates[-1]
  dat <- as.data.frame(t(apply(ri, 1,
                               function(row) {optim(par = c(mu, sig),
                                                    fn = ll_func,
                                                    ri = row,
                                                    lj = lj, uj = uj, dj = dj,
                                                    eps = eps,
                                                    control = list(
                                                      fnscale = -1))$par })))
  colnames(dat) <- c("mean", "sd")
  return(dat)
}

#' Maximum likelihood estimation of mean and standard deviation using sort-seq
#' read count data
#' 
#' This function is used to obtain maximum likelihood estimates of mean and 
#' standard deviation for each variant in a sort-seq experiment. The read count
#' data for each variantis used with the parameters describing sorting gates and
#' total cells sorted to mininize the log-likelihood function.
#' @param gates List of integers indicating the gate boundaries
#' @param read.counts Dataframe with rows corresponding to variants and columns
#' to gates with elements indicating read counts
#' @param sorted.cells List of integers indicating the total number of cells 
#' sorted into each gate
#' @param eps Integer describing the error rate of sorting and sequencing.
#' @return Dataframe with columns indicating the estimated cells for each 
#' variant sorted into each bin, the total number of cells sorted, and the MLE 
#' mean and s.d. for each variant.
#' @export
sortseqMLE <- function(gates, read.counts, sorted.cells, eps = 0) {
  dj <- sorted.cells / colSums(read.counts)

  cell.counts <- reads2cells(read.counts, sorted.cells)

  f.params <- lnorm_mle(ri = read.counts, 2, 0.5, gates, dj, eps)

  out <- data.frame(cell.counts, totalCells = rowSums(cell.counts), f.params)

  return(out)
}

#' Convert read counts for each variant to estimates of cells sorted containing
#' each variant in each bin using the proportionality constant dj
#' 
#' @param read.counts Dataframe with rows corresponding to variants and columns
#' to gates with elements indicating read counts
#' @param sorted.cells List of integers indicating the total number of cells 
#' sorted into each gate
reads2cells <- function(read.counts, sorted.cells) {
  dj <- sorted.cells / colSums(read.counts)
  bj <- sweep(read.counts, MARGIN=2, dj, `*`)
  return(bj)
}

#' Convert log-normal parameters μ and σ to the log-normal mean
#' 
#' @param x list of log-normal parameters
#' @return integer idicating log-normal mean
#' @keywords internal
lnorm_mean <- function(x) {
  10^(x$mean + (x$sd^2 / 2))
}

#' Combine parameter estimates for both F0 and Fl datasets and calculate the 
#' dynamic range for each variant
#' 
#' This function is used to combine the maximum likelihood estimates for the 
#' two datasets F0 and Fl in order to calculate the dynamic range for each 
#' variant.
#' @param F0 Dataframe containing maximum likelihood estimates for the 
#' ligand-free sample.
#' @param Fl Dataframe containing maximum likelihood estimates for the 
#' ligand-bound sample.
#' @return Dataframe combining the two dataframes with an additional column
#' with the calculated dynamic range.
#' @export
summarizeAssay <- function(F0, Fl) {
  j <- ncol(F0) - 3

  delta <- lnorm_mean(Fl) - lnorm_mean(F0)
  DR <- delta / lnorm_mean(F0)

  out <- data.frame(row.names(F0), F0, Fl, DR)
  colnames(out) <- c("variant",
                     paste0("F0.bin", 1:j), "F0.n", "F0.mean", "F0.sd",
                     paste0("Fl.bin", 1:j),  "Fl.n", "Fl.mean", "Fl.sd",
                     "DynamicRange")
  row.names(out) <- NULL
  return(out)
}

