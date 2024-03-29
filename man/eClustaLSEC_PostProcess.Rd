% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eClustaLSEC_PostProcess.R
\name{eClustaLSEC_PostProcess}
\alias{eClustaLSEC_PostProcess}
\title{Post-process Hamiltonian Monte Carlo}
\usage{
eClustaLSEC_PostProcess(
  HMC_obj,
  burnin = 1,
  thin = 1,
  summary_function = mean,
  A = NULL,
  showPB = TRUE
)
}
\arguments{
\item{HMC_obj}{Object of class 'eClustaLSEC_HMC'.}

\item{burnin}{integer. Number of leading posterior draws to remove.}

\item{thin}{integer. Only each \code{Thin}-th sample will be processed.}

\item{summary_function}{function, e.g. mean or median. Function used to summarize
the posterior samples.}

\item{A}{igraph.}

\item{showPB}{logical. Whether or not to show the progress bar.}
}
\value{
Named List.
\itemize{
\item HMC_obj object of class 'eClustaLSEC_HMC'
\item keep_seq integer vector giving the index of the original HMC draws kept in HMC_obj
\item HMC_summary object of class 'eClustaLSEC_HMC_summary'
\item logPosterior vector giving the log posterior values at each HMC draw.
}
}
\description{
'eClustaLSEC_PostProcess' post processes the HMC output for the aLSEC model to
make summaries of the posterior draws possible.
}
