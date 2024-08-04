b <- c()
for(i in 1:20000){
  sum_bvs <- sum(sample(c(-2,2),prob=c(0.5,0.5),size=72,replace=T))
  sum_shapes <- sum(sample(c(-.25,0,.25),prob=c(.33,.33,.33),size=72,replace=T))
  b[i]<- sum_bvs + sum_shapes
}
hist(b)
mean(b)
sd(b)