
###########
#ROC CURVES
###########

roc.data <- function(truth,scores) {
	out <- do.call(rbind,lapply(c(-Inf,sort(scores),Inf), function(t) {
		calls <- scores >= t
		tp <- sum(calls & truth)
		tn <- sum(!calls & !truth)
		fp <- sum(calls & !truth)
		fn <- sum(!calls & truth)
		ppv.prec <- tp/(tp+fp)
		tpr.sens <- tp/(tp+fn)
		fpr.fall <- fp/(tn+fp)
		cbind(thresh=t,ppv.prec=ppv.prec,tpr.sens=tpr.sens,fpr.fall=fpr.fall)
	}))
	# rbind(
	# 	c(thresh=-Inf,ppv.prec=0,tpr.sens=1,fpr.fall=1),
	# 	out,
	# 	c(thresh=Inf,ppv.prec=1,tpr.sens=0,fpr.fall=0)
	# )
	out
}
draw.roc <- function(data,col="steelblue3",main="",lwd=2,add=FALSE) {
	# data <- roc.data(truth,scores)
	# data <- na.omit(data)
	if (add) {
		lines(
			100*data[,"fpr.fall"],100*data[,"tpr.sens"],
			col=col,lwd=lwd
		)
	} else {
		plot(
			100*data[,"fpr.fall"],100*data[,"tpr.sens"],
			type="l",
			xlab="False positive rate (%)\n(= 100%-specificity)", ylab="Sensitivity or True positive rate (%)",
			main=main,
			xlim=c(0,100),ylim=c(0,100),col=col,lwd=lwd
		)
	}
}
draw.prc <- function(data,col="firebrick3",main="",lwd=2,add=FALSE) {
	# data <- roc.data(truth,scores)
	# data <- na.omit(data)

	for (i in 2:nrow(data)) {
		if (!is.na(data[i,"ppv.prec"]) && data[i,"ppv.prec"] < data[i-1,"ppv.prec"]) {
			data[i,"ppv.prec"] <- data[i-1,"ppv.prec"]
		}
	}

	if (add) {
		lines(
			100*data[,"tpr.sens"],100*data[,"ppv.prec"],
			col=col,lwd=lwd
		)
	} else {
		plot(
			100*data[,"tpr.sens"],100*data[,"ppv.prec"],
			type="l",
			xlab="Recall (%)", ylab="Precision (%)",
			main=main,
			xlim=c(0,100),ylim=c(0,100),col=col,lwd=lwd
		)
	}
}
auroc <- function(data) {
	# data <- roc.data(truth,scores)
	# data <- na.omit(data)
	sum(sapply(1:(nrow(data)-1),function(i) {
		delta.x <- data[i,"fpr.fall"]-data[i+1,"fpr.fall"]
		y <- (data[i,"tpr.sens"]+data[i+1,"tpr.sens"])/2
		delta.x * y
	}))
}
auprc <- function(data) {
	data <- na.omit(data)
	sum(sapply(1:(nrow(data)-1),function(i) {
		delta.x <- data[i,"tpr.sens"]-data[i+1,"tpr.sens"]
		y <- (data[i,"ppv.prec"]+data[i+1,"ppv.prec"])/2
		delta.x * y
	}))
}
