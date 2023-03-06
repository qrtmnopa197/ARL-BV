probs <- c(0.04,0.04,0.04,0.04,0.04,0.05,0.075,0.075,0.1,0.2,0.3)
values <- c(-5:5)

gt <- c()
eq <- 0
for(i in 1:10000){
  fo <- sample(values,prob=probs,size=1)
  ro <- sample(values,prob=rep(1/11,11),size=1)
  if(fo > ro){
    gt <- c(gt,1)
  } else{
    gt <- c(gt,0)
  } 
  if(fo == ro){
    eq <- eq + 1
  }
}

#probs <- c(0.05,0.06,0.06,0.06,0.07,0.07,0.08,0.09,0.11,0.15,0.2)
probs <- c(0.04,0.05,0.05,0.0,0.06,0.07,0.08,0.09,0.12,0.16,0.23)
values <- c(-5:5)
probs %*% values
gt <- c()
for(i in 1:25000){
  fo <- sample(values,prob=probs,size=1)
  ro <- sample(values,prob=rep(1/11,11),size=1)
  if(fo > ro){
    gt <- c(gt,1)
  }else if (ro > fo){
    gt <- c(gt,0)
  }
}
length(which(gt==0))/length(gt)

probs <- c(0.12,0.12,0.13,0.16,0.2,0.27)
values <- c(-1:4)
probs %*% values
gt <- c()
for(i in 1:25000){
  fo <- sample(values,prob=probs,size=1)
  ro <- sample(values,prob=rep(1/6,6),size=1)
  if(fo > ro){
    gt <- c(gt,1)
  }else if (ro > fo){
    gt <- c(gt,0)
  }
}
length(which(gt==0))/length(gt)


