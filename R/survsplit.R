"survsplit" <-
function (data, cut, end, event, start, id = NULL, zero = 0, 
    episode = NULL,interval=NULL) 
{
	ntimes <- length(cut)
	n <- nrow(data)

	p <- ncol(data)
	if(length(interval)>0){
		ntimes <- ntimes -1
		sttime <- c(rep(0, n), rep(cut[-length(cut)], each = n))
		endtime <- rep(cut, each = n)
	}
	else{
		endtime <- rep(c(cut, Inf), each = n)
		sttime <- c(rep(0, n), rep(cut, each = n))	
	}
	newdata <- lapply(data, rep, ntimes + 1)
	
	eventtime <- newdata[[end]]
	if (start %in% names(data)) 
	starttime <- newdata[[start]]
	else starttime <- rep(zero, length = (ntimes + 1) * n)
	
	starttime <- pmax(sttime, starttime)
	epi <- rep(0:ntimes, each = n)

	if(length(interval)>0) status <- ifelse(eventtime <= endtime & eventtime >= starttime, newdata[[event]], 0)
	else status <- ifelse(eventtime <= endtime & eventtime > starttime, newdata[[event]], 0)
	endtime <- pmin(endtime, eventtime)
	if(length(interval)>0)	drop <- (starttime > endtime)|(starttime==endtime&status==0)
	else drop <- starttime >= endtime
	newdata <- do.call("data.frame", newdata)

	newdata <- newdata[!drop, ]
	newdata[, start] <- starttime[!drop]
	newdata[, end] <- endtime[!drop]
	newdata[, event] <- status[!drop]
	if (!is.null(id)) 
		newdata[, id] <- rep(rownames(data), ntimes + 1)[!drop]
  fu<-NULL
        if(length(interval)>2){
		for(it in 1:length(interval[-1])){
			drop1 <- sum(!drop[1:(it*n-n)])
			drop2 <- sum(!drop[(it*n-n+1):(it*n)])
			drop3 <- sum(!drop[(it*n+1):(length(interval[-1])*n)])
			if(it==1)  fu<-cbind(fu, c( rep(1,drop2) , rep(0,drop3) ) )
			else if(it==length(interval[-1])) fu<-cbind(fu, c( rep(0,drop1) , rep(1,drop2) ) )
			else fu<-cbind(fu,c( rep(0,drop1),rep(1,drop2),rep(0,drop3)))		
		}
		fu<-as.data.frame(fu)
		names(fu)<-c(paste("fu [",interval[-length(interval)],",",interval[-1],")",sep=""))
		newdata <- cbind(newdata,fu)
	}
	else if(length(interval)==2){
		fu <- rep(1,sum(!drop))
	        newdata <- cbind(newdata,fu)
	       	names(newdata)[ncol(newdata)] <- paste("fu [",interval[1],",",interval[2],")",sep="")
         }
        if (!is.null(episode)) 
        newdata[, episode] <- epi[!drop]
	newdata
}
