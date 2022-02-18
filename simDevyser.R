#######################################
#
# Setting 4:
# MAFs from Devyser assay
# Diploid (i.e. no X/Y)
#
#######################################

# initialize matrices to store simulation results
res.related.devyser <- matrix(-1, nrow = nsim, ncol = 24)
res.unrelated.devyser <- matrix(-1, nrow = nsim, ncol = 24)
res.pc.devyser <- matrix(-1, nrow = nsim, ncol = 24)
res.related.devyser.homo <- matrix(-1, nrow = nsim, ncol = 24)
res.unrelated.devyser.homo <- matrix(-1, nrow = nsim, ncol = 24)
res.pc.devyser.homo <- matrix(-1, nrow = nsim, ncol = 24)
res.related.devyser.allinf <- matrix(-1, nrow = nsim, ncol = 24)
res.unrelated.devyser.allinf <- matrix(-1, nrow = nsim, ncol = 24)
res.pc.devyser.allinf <- matrix(-1, nrow = nsim, ncol = 24)

# calculate MAFs from empirical data
hardy <- read.csv("~/Publicaties/Chimerisme_nmarkers/Data/hardy.csv")
p <- hardy[,1]*2+hardy[,2]
q <- hardy[,3]*2+hardy[,2]
pop.freq <- 0.5 - abs(p/(p+q)-0.5)


for(sim in 1:nsim){
  cat("sim ", sim, " of ", nsim,"\n")
  for(marker in 1:24){
    # generate parental and offspring genotypes
    p1 <- sample(c("A", "B"), size = 2, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
    p2 <- sample(c("A", "B"), size = 2, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
    o1 <- c(sample(p1, size = 1), sample(p2, size = 1))
    o2 <- c(sample(p1, size = 1), sample(p2, size = 1))
    
    # true/false informative
    if(sum(sort(o1) != sort(o2)) == 2 || (sum(sort(o1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 ))){
      res.related.devyser[sim, marker] <- 1
    } else {
      res.related.devyser[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(p2)) == 2 || (sum(sort(p1) == c("A", "B")) == 2 && (sum(p2 == c("A", "A")) == 2 || sum(p2 == c("B", "B")) == 2 ))){
      res.unrelated.devyser[sim, marker] <- 1
    } else {
      res.unrelated.devyser[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(o2)) == 2 || (sum(sort(p1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 ))){
      res.pc.devyser[sim, marker] <- 1
    } else {
      res.pc.devyser[sim, marker] <- 0
    }
    
    # true/false informative homozygous
    if(sum(sort(o1) != sort(o2)) == 2){
      res.related.devyser.homo[sim, marker] <- 1
    } else {
      res.related.devyser.homo[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(p2)) == 2){
      res.unrelated.devyser.homo[sim, marker] <- 1
    } else {
      res.unrelated.devyser.homo[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(o2)) == 2){
      res.pc.devyser.homo[sim, marker] <- 1
    } else {
      res.pc.devyser.homo[sim, marker] <- 0
    }
    
    # true/false informative + potentially informative
    if(sum(sort(o1) != sort(o2)) == 2 || (sum(sort(o1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 )) || (sum(sort(o2) == c("A", "B")) == 2 && (sum(o1 == c("A", "A")) == 2 || sum(o1 == c("B", "B")) == 2 ))){
      res.related.devyser.allinf[sim, marker] <- 1
    } else {
      res.related.devyser.allinf[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(p2)) == 2|| (sum(sort(p1) == c("A", "B")) == 2 && (sum(p2 == c("A", "A")) == 2 || sum(p2 == c("B", "B")) == 2 )) || (sum(sort(p2) == c("A", "B")) == 2 && (sum(p1 == c("A", "A")) == 2 || sum(p1 == c("B", "B")) == 2 ))){
      res.unrelated.devyser.allinf[sim, marker] <- 1
    } else {
      res.unrelated.devyser.allinf[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(o2)) == 2|| (sum(sort(p1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 )) || (sum(sort(o2) == c("A", "B")) == 2 && (sum(p1 == c("A", "A")) == 2 || sum(p1 == c("B", "B")) == 2 ))){
      res.pc.devyser.allinf[sim, marker] <- 1
    } else {
      res.pc.devyser.allinf[sim, marker] <- 0
    }
  }
}

# generate average informative values across 
#  simulations, as function of number of markers
means.vec.related.devyser <- vector("numeric", 24-2)
means.vec.unrelated.devyser <- vector("numeric", 24-2)
means.vec.pc.devyser <- vector("numeric", 24-2)
means.vec.related.devyser.homo <- vector("numeric", 24-2)
means.vec.unrelated.devyser.homo <- vector("numeric", 24-2)
means.vec.pc.devyser.homo <- vector("numeric", 24-2)
means.vec.related.devyser.allinf <- vector("numeric", 24-2)
means.vec.unrelated.devyser.allinf <- vector("numeric", 24-2)
means.vec.pc.devyser.allinf <- vector("numeric", 24-2)
for(calcmean in 3:24){
  means.vec.related.devyser[calcmean-2] <- mean(rowSums(res.related.devyser[,1:calcmean]) >= 3)
  means.vec.unrelated.devyser[calcmean-2] <- mean(rowSums(res.unrelated.devyser[,1:calcmean]) >= 3)
  means.vec.pc.devyser[calcmean-2] <- mean(rowSums(res.pc.devyser[,1:calcmean]) >= 3)
  means.vec.related.devyser.homo[calcmean-2] <- mean(rowSums(res.related.devyser.homo[,1:calcmean]) >= 3)
  means.vec.unrelated.devyser.homo[calcmean-2] <- mean(rowSums(res.unrelated.devyser.homo[,1:calcmean]) >= 3)
  means.vec.pc.devyser.homo[calcmean-2] <- mean(rowSums(res.pc.devyser.homo[,1:calcmean]) >= 3)
  means.vec.related.devyser.allinf[calcmean-2] <- mean(rowSums(res.related.devyser.allinf[,1:calcmean]) >= 3)
  means.vec.unrelated.devyser.allinf[calcmean-2] <- mean(rowSums(res.unrelated.devyser.allinf[,1:calcmean]) >= 3)
  means.vec.pc.devyser.allinf[calcmean-2] <- mean(rowSums(res.pc.devyser.allinf[,1:calcmean]) >= 3)
}
