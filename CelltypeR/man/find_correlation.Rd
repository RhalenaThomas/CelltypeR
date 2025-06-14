% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CelltypeR.R
\name{find_correlation}
\alias{find_correlation}
\title{CAM predict cell types in Flow Cytometry data by correlation to an existing reference matrix}
\usage{
find_correlation(test, reference, min_corr = 0.1, min_diff = 0.05)
}
\description{
The function requires an reference matix already created with the same markers
(antibodies) as used in the Flow Cytometry expression data frame. Correlation of
expression across markers for each cell type in the reference matrix and each cell in the
sample data is calculated. A cutoff for the R value (pearson correlation coefficient) is set,
where an R value must be above this threshold for a cell type to be assigned.  Cells with only
correlations below the threshold are set as "unknown". The threshold is defined by "min_corr".
The best and second best R value is recorded for each cell. If the difference between R max 1
and R max 2 is less than or equal to "min_diff" that cell is assigned as a joint cell type. For
example "neuron-npc". Test is a data frame of the expression for each marker and cell. Reference
is the reference matrix, cell type by marker.
}
