"maxlik" <-
function (rform,interval,method, subset,init, control) 
{

	data <- data.frame(start = rform$start, Y = rform$Y, stat = rform$status,rform$R)
	if(rform$m!=0) data <- cbind(data,rform$X)

	
	
	if(any(data$start>data$Y))stop("Negative follow-up time")
	data$stat[data$Y >= max(interval) * 365.24] <- 0
	data$Y <- pmin(data$Y, max(interval) * 365.24)
    	
	

	max.time <- max(data$Y)/365.24
	if(max.time<max(interval))    interval<-interval[1:(sum(max.time>interval)+1)]

	
	 
	fk <- (attributes(rform$ratetable)$factor != 1)
	nfk <- length(fk)

	data <- cbind(data,offset=rform$offset)
	data <- survsplit(data, cut = interval[-1] * 365.24, end = "Y", event = "stat", 
		start = "start", episode = "epi",interval=interval)
	offset <- data$offset
	data$offset <- NULL
	
	
	d.int <- diff(interval)
	
	
	data[, 4:(nfk + 3)] <- data[, 4:(nfk + 3)] + data$start %*% t(fk)
	data$lambda<-rep(0,nrow(data))
	
	
	nsk<-nrow(data[data$stat==1,])
	xx <- srvxp.fit(data[data$stat==1, 4:(nfk + 3)]+(data[data$stat==1,]$Y-data[data$stat==1,]$start)%*% t(fk), 
			rep(1,nsk), rform$ratetable)
	
	data$lambda[data$stat==1]<- -log(xx)*365.24
	
	
	xx <- srvxp.fit(data[, 4:(nfk + 3)],data$Y-data$start, rform$ratetable)
	data$epi<-NULL
	
	data$ds <- -log(xx)
	
	data$Y <- data$Y/365.24
	data$start <- data$start/365.24
	
	data <- data[,-(4:(3+nfk))]

	fit <- lik.fit(data,rform$m,length(interval[-1]),init,control,offset)
	
	fit$int <- interval
	class(fit) <- "rsadd"
	fit
}

