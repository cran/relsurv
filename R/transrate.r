transrate<-function(men,women,yearlim,int.length=1){
	
	if(any(dim(men)!=dim(women)))stop("The men and women matrices must be of the same size. \n In case of missing values at the end carry the last value forward")
	if((yearlim[2]-yearlim[1])/int.length+1!=dim(men)[2])stop("'yearlim' cannot be divided into intervals of equal length")
	if(!is.matrix(men)|!is.matrix(women))stop("input tables must be of class matrix")

	dimi<-dim(men)
	temp<-array(c(men,women),dim=c(dimi,2))

	temp<- -log(temp)/365.24
	temp<-aperm(temp,c(1,3,2))

	cp<-as.date(apply(matrix(yearlim[1]+int.length*(0:(dimi[2]-1)),ncol=1),1,function(x){paste("1jan",x,sep="")}))

	attributes(temp)<-list(
		dim=c(dimi[1],2,dimi[2]),
		dimnames=list(as.character(0:(dimi[1]-1)),c("male","female"),as.character(yearlim[1]+int.length*(0:(dimi[2]-1)))),
		dimid=c("age","sex","year"),
		factor=c(0,1,0),
		cutpoints=list((0:(dimi[1]-1))*(365.24),NULL,cp),
		class="ratetable"
	)
	temp
}