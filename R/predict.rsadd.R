#' Subject-specific prediction from rsadd
#'
#' Function
#' @param object An rsadd object
#' @param newdata A data.frame with one row (add covariate values in columns)
#' @param ... Not used for now
#' @return A data.frame with times, excess and population hazard.
#'
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}
#' @export
predict.rsadd <- function(object, newdata, ...){
  # object: a rsadd object
  # newdata: data.frame with one row, in columns add covariate values.

  # subject-specific prediction.

  # Zaenkrat delamo z object$times. Treba bo dodati se add.times. Ta add.times bo verjetno v rsadd

  lin.pred <- 0
  for(cov in names(object$coefficients)){
    lin.pred <- lin.pred + sum(newdata[1,cov]*object$coefficients[cov])
  }

  # haz_function(object$formula, rdata, object$ratetable, rmap=list(age=age*365.241), add.times=0, include.all.times = FALSE)

  kar <- deparse(object$formula[[2]])
  autkam <- gsub(' ', '', strsplit(substr(kar, start = 6, stop=nchar(kar)-1), ',')[[1]])
  newdata_2 <- data.frame(matrix(1, nrow=1, ncol=length(autkam)))
  if(length(autkam)==3) newdata_2[,1] <- 0
  colnames(newdata_2) <- autkam
  newdata_2 <- cbind(newdata, newdata_2)

  if(is.null(object$rmap)){
    rform <- suppressWarnings(rformulate(object$formula,
                                         newdata_2,object$ratetable))
  } else{
    rform <- suppressWarnings(rformulate(object$formula,
                                         newdata_2,object$ratetable, rmap = object$rmap))
  }

  fk <- (attributes(rform$ratetable)$factor != 1)
  nfk <- length(fk)
  pop.surv <- sapply(1:length(object$times), function(i) exp_prep(rform$data[1, 4:(nfk + 3), drop = FALSE], object$times[i], object$ratetable))
  Haz.p <- -log(pop.surv)

  data.frame(time=object$times, Haz.e=object$Lambda0*exp(lin.pred), Haz.p=Haz.p)
}
