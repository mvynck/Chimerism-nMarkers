#######################################
#
# Setting 3:
# MAFs sampled from U[0.3, 0.5]
# Diploid (i.e. no X/Y)
#
#######################################

# initialize matrices to store simulation results
res.related.real1 <- matrix(-1, nrow = nsim, ncol = max.markers)
res.unrelated.real1 <- matrix(-1, nrow = nsim, ncol = max.markers)
res.pc.real1 <- matrix(-1, nrow = nsim, ncol = max.markers)
res.related.real1.homo <- matrix(-1, nrow = nsim, ncol = max.markers)
res.unrelated.real1.homo <- matrix(-1, nrow = nsim, ncol = max.markers)
res.pc.real1.homo <- matrix(-1, nrow = nsim, ncol = max.markers)
res.related.real1.allinf <- matrix(-1, nrow = nsim, ncol = max.markers)
res.unrelated.real1.allinf <- matrix(-1, nrow = nsim, ncol = max.markers)
res.pc.real1.allinf <- matrix(-1, nrow = nsim, ncol = max.markers)

set.seed(1)
pop.freq <- runif(max.markers, 0.3, 0.5)
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
      res.related.real1[sim, marker] <- 1
    } else {
      res.related.real1[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(p2)) == 2 || (sum(sort(p1) == c("A", "B")) == 2 && (sum(p2 == c("A", "A")) == 2 || sum(p2 == c("B", "B")) == 2 ))){
      res.unrelated.real1[sim, marker] <- 1
    } else {
      res.unrelated.real1[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(o2)) == 2 || (sum(sort(p1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 ))){
      res.pc.real1[sim, marker] <- 1
    } else {
      res.pc.real1[sim, marker] <- 0
    }
    
    # true/false informative homozygous
    if(sum(sort(o1) != sort(o2)) == 2){
      res.related.real1.homo[sim, marker] <- 1
    } else {
      res.related.real1.homo[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(p2)) == 2){
      res.unrelated.real1.homo[sim, marker] <- 1
    } else {
      res.unrelated.real1.homo[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(o2)) == 2){
      res.pc.real1.homo[sim, marker] <- 1
    } else {
      res.pc.real1.homo[sim, marker] <- 0
    }
    
    # true/false informative + potentially informative
    if(sum(sort(o1) != sort(o2)) == 2 || (sum(sort(o1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 )) || (sum(sort(o2) == c("A", "B")) == 2 && (sum(o1 == c("A", "A")) == 2 || sum(o1 == c("B", "B")) == 2 ))){
      res.related.real1.allinf[sim, marker] <- 1
    } else {
      res.related.real1.allinf[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(p2)) == 2|| (sum(sort(p1) == c("A", "B")) == 2 && (sum(p2 == c("A", "A")) == 2 || sum(p2 == c("B", "B")) == 2 )) || (sum(sort(p2) == c("A", "B")) == 2 && (sum(p1 == c("A", "A")) == 2 || sum(p1 == c("B", "B")) == 2 ))){
      res.unrelated.real1.allinf[sim, marker] <- 1
    } else {
      res.unrelated.real1.allinf[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(o2)) == 2|| (sum(sort(p1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 )) || (sum(sort(o2) == c("A", "B")) == 2 && (sum(p1 == c("A", "A")) == 2 || sum(p1 == c("B", "B")) == 2 ))){
      res.pc.real1.allinf[sim, marker] <- 1
    } else {
      res.pc.real1.allinf[sim, marker] <- 0
    }
  }
}

# generate average informative values across 
#  simulations, as function of number of markers
means.vec.related.real1 <- vector("numeric", max.markers-2)
means.vec.unrelated.real1 <- vector("numeric", max.markers-2)
means.vec.pc.real1 <- vector("numeric", max.markers-2)
means.vec.related.real1.homo <- vector("numeric", max.markers-2)
means.vec.unrelated.real1.homo <- vector("numeric", max.markers-2)
means.vec.pc.real1.homo <- vector("numeric", max.markers-2)
means.vec.related.real1.allinf <- vector("numeric", max.markers-2)
means.vec.unrelated.real1.allinf <- vector("numeric", max.markers-2)
means.vec.pc.real1.allinf <- vector("numeric", max.markers-2)
for(calcmean in 3:max.markers){
  means.vec.related.real1[calcmean-2] <- mean(rowSums(res.related.real1[,1:calcmean]) >= 3)
  means.vec.unrelated.real1[calcmean-2] <- mean(rowSums(res.unrelated.real1[,1:calcmean]) >= 3)
  means.vec.pc.real1[calcmean-2] <- mean(rowSums(res.pc.real1[,1:calcmean]) >= 3)
  means.vec.related.real1.homo[calcmean-2] <- mean(rowSums(res.related.real1.homo[,1:calcmean]) >= 3)
  means.vec.unrelated.real1.homo[calcmean-2] <- mean(rowSums(res.unrelated.real1.homo[,1:calcmean]) >= 3)
  means.vec.pc.real1.homo[calcmean-2] <- mean(rowSums(res.pc.real1.homo[,1:calcmean]) >= 3)
  means.vec.related.real1.allinf[calcmean-2] <- mean(rowSums(res.related.real1.allinf[,1:calcmean]) >= 3)
  means.vec.unrelated.real1.allinf[calcmean-2] <- mean(rowSums(res.unrelated.real1.allinf[,1:calcmean]) >= 3)
  means.vec.pc.real1.allinf[calcmean-2] <- mean(rowSums(res.pc.real1.allinf[,1:calcmean]) >= 3)
}
