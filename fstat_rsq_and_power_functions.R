library(TwoSampleMR)

# Calculating instrument F-statistic for Two-sample MR using beta and se
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5446088/
# F = b^2/se^2
get_fval_bowden_method = function(b, se) {
    F = b^2/se^2
    return(F)
}

## Calculating F-statistic and r-square value from QTL p-value and N, using the
## Two-sample MR R package functions.
## This should give almost identical values to the Bowden method
get_fval_2SMR=function(p,n) {
    Fval <- suppressWarnings(qf(p, 1, n - 1, low = FALSE))
    return(Fval)
}

## Use function get_r_from_pn(p=d$p.value, n=d$N) to get r
get_rsq_2SMR=function(p,n) {
    r=get_r_from_pn(p, n)
    rsq=r^2
    return(rsq)
}

## Calculate power
## The power is estimated using the method by Burgess et al.
## https://www.ncbi.nlm.nih.gov/pubmed/24608958
## and can also be estimated using the online power calculator:
## https://sb452.shinyapps.io/power/

# ratio is 1:X (case:ctrl)
# For SCZ: n = 105318; ratio = 1.59  (40675cases,64643 ctrls)
# For BIP: n = 51710; ratio = 1.54 (20352 cases, 31358 ctrls)
# For MDD: n = 480359; ratio = 2.55 (40675 cases, 64643 ctrls)
# n=sample size,
# rsq is proportion of variance in exposure explained by SNP

calc_power <- function(b1, n, ratio, rsq) {
    power=pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))
    return(power)
}

