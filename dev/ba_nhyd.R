






ba_nhyd <- function(probs, verbose = T){


	N <- length(probs)
	
	





}

if(!exists("tprob")){
	tprob = c(0.2,0.4,0.7,0.9)
}

df <- data.frame(
	prob=double(),
	cbrob=double(),
	nhyd=integer(),
	depth=integer()
)


for(i in 1:length(tprob)){

	message(sprintf("Prob %d is %0.2f, nrow(df) is %d",i,tprob[i],nrow(df)))
	
	ht <- data.frame(prob = tprob[i], cbprob = tprob[i], nhyd = 1, depth = i)
	ot <- data.frame(prob = 1-tprob[i], cbprob = 1-tprob[i], nhyd = 0, depth = i)
	
	newdf <- data.frame(
		prob=double(),
		cbrob=double(),
		nhyd=integer(),
		depth=integer()
	)
	
	if(nrow(df)==0){
	
		df <- rbind(df,ht)
		df <- rbind(df,ot)
	}
	else{
		for(j in 1:nrow(df)){df
		
			#another hydroxylation
			htn <- data.frame(
				prob = tprob[i], 
				cbprob = tprob[i] * df$cbprob[j],
				nhyd = df$nhyd[j] + 1,
				depth = i
				)
			newdf <- rbind(newdf,htn)	
	
			#no new hydroxylation
			otn <- data.frame(
				prob = (1-tprob[i]), 
				cbprob = (1-tprob[i]) * df$cbprob[j],
				nhyd = df$nhyd[j] + 0,
				depth = i
				)
			newdf <- rbind(newdf,otn)	
		}
		
		df = newdf
		rm(newdf)
	}
	
	
	for(j in 1:nrow(df)){
		#message(sprintf("beep %d",j))
		message(sprintf("%0.3f\t%0.3f\t%d\t%d",df$prob[j],df$cbprob[j],df$nhyd[j],df$depth[j]))
		#message(sprintf("%0.3f\t%0.3f\t%d\t%d",df$prob[j],df$dbprob[j],df$nhyd[j],df$depth[j]))
	}
	
	message(" ")


	result <- data.frame(
		prob = double(),
		nhyd = integer()
	)
	#finally, let's get the probabilities for each hydroxylation level.
	for( i in 0:length(tprob)){
		r <- data.frame( prob = sum(df$cbprob[df$nhyd == i]) , nhyd=i) 	
		result <- rbind(result,r)
	}
}
