#######################################
#
# Setting 2:
# MAFs all 0.5 ("ideal scenario")
# Diploid (i.e. no X/Y)
#
#######################################

# initialize matrices to store simulation results
res.related.ideal <- matrix(-1, nrow = nsim, ncol = max.markers)
res.unrelated.ideal <- matrix(-1, nrow = nsim, ncol = max.markers)
res.pc.ideal <- matrix(-1, nrow = nsim, ncol = max.markers)
res.related.ideal.homo <- matrix(-1, nrow = nsim, ncol = max.markers)
res.unrelated.ideal.homo <- matrix(-1, nrow = nsim, ncol = max.markers)
res.pc.ideal.homo <- matrix(-1, nrow = nsim, ncol = max.markers)
res.related.ideal.allinf <- matrix(-1, nrow = nsim, ncol = max.markers)
res.unrelated.ideal.allinf <- matrix(-1, nrow = nsim, ncol = max.markers)
res.pc.ideal.allinf <- matrix(-1, nrow = nsim, ncol = max.markers)

pop.freq <- rep(0.5, max.markers)
for(sim in 1:nsim){
  cat("sim ", sim, " of ", nsim,"\n")
  for(marker in 1:max.markers){
    # generate parental and offspring genotypes
    p1 <- sample(c("A", "B"), size = 2, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
    p2 <- sample(c("A", "B"), size = 2, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
    o1 <- c(sample(p1, size = 1), sample(p2, size = 1))
    o2 <- c(sample(p1, size = 1), sample(p2, size = 1))
    
    # true/false informative
    if(sum(sort(o1) != sort(o2)) == 2 || (sum(sort(o1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 ))){
      res.related.ideal[sim, marker] <- 1
    } else {
      res.related.ideal[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(p2)) == 2 || (sum(sort(p1) == c("A", "B")) == 2 && (sum(p2 == c("A", "A")) == 2 || sum(p2 == c("B", "B")) == 2 ))){
      res.unrelated.ideal[sim, marker] <- 1
    } else {
      res.unrelated.ideal[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(o2)) == 2 || (sum(sort(p1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 ))){
      res.pc.ideal[sim, marker] <- 1
    } else {
      res.pc.ideal[sim, marker] <- 0
    }
    
    # true/false informative homozygous
    if(sum(sort(o1) != sort(o2)) == 2){
      res.related.ideal.homo[sim, marker] <- 1
    } else {
      res.related.ideal.homo[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(p2)) == 2){
      res.unrelated.ideal.homo[sim, marker] <- 1
    } else {
      res.unrelated.ideal.homo[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(o2)) == 2){
      res.pc.ideal.homo[sim, marker] <- 1
    } else {
      res.pc.ideal.homo[sim, marker] <- 0
    }
    
    # true/false informative + potentially informative
    if(sum(sort(o1) != sort(o2)) == 2 || (sum(sort(o1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 )) || (sum(sort(o2) == c("A", "B")) == 2 && (sum(o1 == c("A", "A")) == 2 || sum(o1 == c("B", "B")) == 2 ))){
      res.related.ideal.allinf[sim, marker] <- 1
    } else {
      res.related.ideal.allinf[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(p2)) == 2|| (sum(sort(p1) == c("A", "B")) == 2 && (sum(p2 == c("A", "A")) == 2 || sum(p2 == c("B", "B")) == 2 )) || (sum(sort(p2) == c("A", "B")) == 2 && (sum(p1 == c("A", "A")) == 2 || sum(p1 == c("B", "B")) == 2 ))){
      res.unrelated.ideal.allinf[sim, marker] <- 1
    } else {
      res.unrelated.ideal.allinf[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(o2)) == 2|| (sum(sort(p1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 )) || (sum(sort(o2) == c("A", "B")) == 2 && (sum(p1 == c("A", "A")) == 2 || sum(p1 == c("B", "B")) == 2 ))){
      res.pc.ideal.allinf[sim, marker] <- 1
    } else {
      res.pc.ideal.allinf[sim, marker] <- 0
    }
  }
}

# generate average informative values across 
#  simulations, as function of number of markers
means.vec.related.ideal <- vector("numeric", max.markers-2)
means.vec.unrelated.ideal <- vector("numeric", max.markers-2)
means.vec.pc.ideal <- vector("numeric", max.markers-2)
means.vec.related.ideal.homo <- vector("numeric", max.markers-2)
means.vec.unrelated.ideal.homo <- vector("numeric", max.markers-2)
means.vec.pc.ideal.homo <- vector("numeric", max.markers-2)
means.vec.related.ideal.allinf <- vector("numeric", max.markers-2)
means.vec.unrelated.ideal.allinf <- vector("numeric", max.markers-2)
means.vec.pc.ideal.allinf <- vector("numeric", max.markers-2)
for(calcmean in 3:max.markers){
  means.vec.related.ideal[calcmean-2] <- mean(rowSums(res.related.ideal[,1:calcmean]) >= 3)
  means.vec.unrelated.ideal[calcmean-2] <- mean(rowSums(res.unrelated.ideal[,1:calcmean]) >= 3)
  means.vec.pc.ideal[calcmean-2] <- mean(rowSums(res.pc.ideal[,1:calcmean]) >= 3)
  means.vec.related.ideal.homo[calcmean-2] <- mean(rowSums(res.related.ideal.homo[,1:calcmean]) >= 3)
  means.vec.unrelated.ideal.homo[calcmean-2] <- mean(rowSums(res.unrelated.ideal.homo[,1:calcmean]) >= 3)
  means.vec.pc.ideal.homo[calcmean-2] <- mean(rowSums(res.pc.ideal.homo[,1:calcmean]) >= 3)
  means.vec.related.ideal.allinf[calcmean-2] <- mean(rowSums(res.related.ideal.allinf[,1:calcmean]) >= 3)
  means.vec.unrelated.ideal.allinf[calcmean-2] <- mean(rowSums(res.unrelated.ideal.allinf[,1:calcmean]) >= 3)
  means.vec.pc.ideal.allinf[calcmean-2] <- mean(rowSums(res.pc.ideal.allinf[,1:calcmean]) >= 3)
}

