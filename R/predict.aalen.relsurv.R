predict.aalen.relsurv <- function(object, newdata, ...){
  # object: an aalen.model object
  # newdata: data.frame with one row, in columns add covariate values.

  # subject-specific prediction.

  # Zaenkrat delamo z object$times. Treba bo dodati se add.times. Ta add.times bo verjetno v rsadd

  lin.pred <- object$coefficients[,2]

  if(ncol(object$coefficients) >= 3){
    for(cov in colnames(object$coefficients)[3:length(colnames(object$coefficients))]){
      lin.pred <- lin.pred + newdata[1,cov]*object$coefficients[,cov]
    }
  }

  # haz_function(object$formula, rdata, object$ratetable, rmap=list(age=age*365.241), add.times=0, include.all.times = FALSE)

  kar <- deparse(object$formula[[2]])
  autkam <- gsub(' ', '', strsplit(substr(kar, start = 6, stop=nchar(kar)-1), ',')[[1]])
  newdata_2 <- data.frame(matrix(1, nrow=1, ncol=length(autkam)))
  if(length(autkam)==3) newdata_2[,1] <- 0
  colnames(newdata_2) <- autkam
  newdata_2 <- cbind(newdata, newdata_2)

  rform <- suppressWarnings(rformulate(object$formula,
                                       newdata_2, object$ratetable, stats::na.omit(), rmap = object$rmap))


  fk <- (attributes(rform$ratetable)$factor != 1)
  nfk <- length(fk)
  pop.surv <- sapply(1:nrow(object$coefficients), function(i) exp_prep(rform$data[1, 4:(nfk + 3), drop = FALSE], object$coefficients[i,1], object$ratetable))
  Haz.p <- -log(pop.surv)

  data.frame(time=object$coefficients[,1], Haz.e=lin.pred, Haz.p=Haz.p)
}
