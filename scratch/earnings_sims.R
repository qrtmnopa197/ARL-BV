#This code gives distributions of expected earnings given the reward distributions and the assumption of
#the (highly unrealistic) worst-case scenario: that the subject chooses the better or worse shape every trial. 

#This loop assumes that subjects will either select the better or worse shape every time. Using central
#limit theorem and the mean and standard deviations of the reward distributions (means and SDs can be
#found in rewards.docx), I calculated distributions of total earnings given 20 draws from each shape. Those distributions, for the initial six
#pairs (EVs 1, .68, .325), are represented in rnorm below (they are all normal distributions, according
#to central limit theorem). This loop simulates twenty thousand subjects selecting the better/worse
#shape every time, creating a distribution (b) of their earnings.
b <- c()
sample(c(-5:5),prob=rep(1/11,11),size=56,replace=T)
sample(c(-2:5),prob=rep(1/8,8),size=56,replace=T)
for(i in 1:20000){
  sum_unipair <- sum(sample(c(-5:5),prob=rep(1/11,11),size=56,replace=T))
  sum_learnpair <- sum(sample(c(-2:5),prob=rep(1/8,8),size=18,replace=T))
  b[i]<- sum_unipair + sum_learnpair
}
hist(b)
mean(b)
sd(b)


#Now, I'm simulating a distribution of the number of earnings-correcting trials it'd take to get X amount of
#change in earnings, assuming that the reward dist. for each shape on those trials is just a random draw
#from gains or losses.
trials_necessary <- c() #trials needed to get $180 in movement
dist <- c(1,2,3,5)
for(i in 1:10000){
  total <- 0 #total movement. Keep doing more trials until you get $180 in movement
  trials_needed <- 0 #trials needed to get there
  while(total < 180){
    added <- sample(dist,1) #randomly draw from dist
    total <- total + added #add that to the total change
    trials_needed <- trials_needed + 1 #one draw and add is one trial
  }
  trials_necessary[i] <- trials_needed #all the trials needed for all hypothetical subjects
}
mean(trials_necessary)
sd(trials_necessary)