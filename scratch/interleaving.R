df <- data.frame(matrix(ncol=3,nrow=1000))
names(df) <- c("twos","threes","fours")
set <- c(1,1,2,2)
for(i in 1:1000){
  trial_ord <- c()
  for(j in 1:7){
    trial_ord <- c(trial_ord,sample(set))
  }
  runs <- rle(trial_ord)
  df$ones[i] <- length(which(runs$lengths == 1))
  df$twos[i] <- length(which(runs$lengths == 2))
  df$threes[i] <- length(which(runs$lengths == 3))
  df$fours[i] <- length(which(runs$lengths == 4))
}



df <- data.frame(matrix(ncol=3,nrow=1000))
names(df) <- c("twos","threes","fours")
setodd <- c(1,1,2)
seteven <- c(2,2,1)
for(i in 1:1000){
  trial_ord <- c()
  for(j in 1:9){
    if(j %% 2 == 0){
      toadd <- sample(seteven)
    }else{
      toadd <- sample(setodd)
    }
    trial_ord <- c(trial_ord,toadd)
  }
  runs <- rle(trial_ord)
  df$ones[i] <- length(which(runs$lengths == 1))
  df$twos[i] <- length(which(runs$lengths == 2))
  df$threes[i] <- length(which(runs$lengths == 3))
  df$fours[i] <- length(which(runs$lengths == 4))
}