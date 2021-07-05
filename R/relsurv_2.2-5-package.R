

#' Relative Survival Data
#' 
#' Survival of patients with colon and rectal cancer diagnosed in 1994-2000.
#' 
#' 
#' @name colrec
#' @docType data
#' @format A data frame with 5971 observations on the following 7 variables:
#' \describe{ \item{sex}{sex (1=male, 2=female).} \item{age}{age (in days).}
#' \item{diag}{date of diagnosis (in date format).} \item{time}{survival time
#' (in days).} \item{stat}{censoring indicator (0=censoring, 1=death).}
#' \item{stage}{cancer stage. Values 1-3, code \code{99} stands for unknown.}
#' \item{site}{cancer site. } }
#' @references Provided by Slovene Cancer Registry. The \code{age}, \code{time}
#' and \code{diag} variables are randomly perturbed to make the identification
#' of patients impossible.
#' @keywords datasets
NULL





#' Survival Data
#' 
#' Survival data.
#' 
#' 
#' @name rdata
#' @docType data
#' @format A data frame with 1040 observations on the following 6 variables:
#' \describe{ \item{time}{survival time (in days).} \item{cens}{censoring
#' indicator (0=censoring, 1=death).} \item{age}{age (in years).}
#' \item{sex}{sex (1=male, 2=female).} \item{year}{date of diagnosis (in date
#' format).} \item{agegr}{age group.} }
#' @references Pohar M., Stare J. (2006) "Relative survival analysis in R."
#' Computer Methods and Programs in Biomedicine, \bold{81}: 272-278.
#' @keywords datasets
NULL





#' Census Data Set for the Slovene Population
#' 
#' Census data set for the Slovene population.
#' 
#' 
#' @name slopop
#' @docType data
#' @keywords datasets
#' @examples
#' 
#' data(slopop)
#' 
NULL



