#######################################
#
# Setting 1:
# MAFs sampled from U[0.2, 0.4]
# Autosomes (i.e. no X/Y)
#
#######################################

# initialize matrices to store simulation results
res.related <- matrix(-1, nrow = nsim, ncol = max.markers)
res.unrelated <- matrix(-1, nrow = nsim, ncol = max.markers)
res.pc <- matrix(-1, nrow = nsim, ncol = max.markers)
res.related.homo <- matrix(-1, nrow = nsim, ncol = max.markers)
res.unrelated.homo <- matrix(-1, nrow = nsim, ncol = max.markers)
res.pc.homo <- matrix(-1, nrow = nsim, ncol = max.markers)
res.related.allinf <- matrix(-1, nrow = nsim, ncol = max.markers)
res.unrelated.allinf <- matrix(-1, nrow = nsim, ncol = max.markers)
res.pc.allinf <- matrix(-1, nrow = nsim, ncol = max.markers)

# sample MAFs
set.seed(1)
pop.freq <- runif(max.markers, 0.2, 0.4)

# loop simulations
for(sim in 1:nsim){
  cat("sim ", sim, " of ", nsim,"\n")
  
  # loop markers in assay
  for(marker in 1:max.markers){
    # generate parental and offspring genotypes
    p1 <- sample(c("A", "B"), size = 2, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
    p2 <- sample(c("A", "B"), size = 2, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
    o1 <- c(sample(p1, size = 1), sample(p2, size = 1))
    o2 <- c(sample(p1, size = 1), sample(p2, size = 1))
    
    # true/false informative
    if(sum(sort(o1) != sort(o2)) == 2 || (sum(sort(o1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 ))){
      res.related[sim, marker] <- 1
    } else {
      res.related[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(p2)) == 2 || (sum(sort(p1) == c("A", "B")) == 2 && (sum(p2 == c("A", "A")) == 2 || sum(p2 == c("B", "B")) == 2 ))){
      res.unrelated[sim, marker] <- 1
    } else {
      res.unrelated[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(o2)) == 2 || (sum(sort(p1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 ))){
      res.pc[sim, marker] <- 1
    } else {
      res.pc[sim, marker] <- 0
    }
    
    # true/false informative homozygous
    if(sum(sort(o1) != sort(o2)) == 2){
      res.related.homo[sim, marker] <- 1
    } else {
      res.related.homo[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(p2)) == 2){
      res.unrelated.homo[sim, marker] <- 1
    } else {
      res.unrelated.homo[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(o2)) == 2){
      res.pc.homo[sim, marker] <- 1
    } else {
      res.pc.homo[sim, marker] <- 0
    }
    
    # true/false informative + potentially informative
    if(sum(sort(o1) != sort(o2)) == 2 || (sum(sort(o1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 )) || (sum(sort(o2) == c("A", "B")) == 2 && (sum(o1 == c("A", "A")) == 2 || sum(o1 == c("B", "B")) == 2 ))){
      res.related.allinf[sim, marker] <- 1
    } else {
      res.related.allinf[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(p2)) == 2|| (sum(sort(p1) == c("A", "B")) == 2 && (sum(p2 == c("A", "A")) == 2 || sum(p2 == c("B", "B")) == 2 )) || (sum(sort(p2) == c("A", "B")) == 2 && (sum(p1 == c("A", "A")) == 2 || sum(p1 == c("B", "B")) == 2 ))){
      res.unrelated.allinf[sim, marker] <- 1
    } else {
      res.unrelated.allinf[sim, marker] <- 0
    }
    if(sum(sort(p1) != sort(o2)) == 2|| (sum(sort(p1) == c("A", "B")) == 2 && (sum(o2 == c("A", "A")) == 2 || sum(o2 == c("B", "B")) == 2 )) || (sum(sort(o2) == c("A", "B")) == 2 && (sum(p1 == c("A", "A")) == 2 || sum(p1 == c("B", "B")) == 2 ))){
      res.pc.allinf[sim, marker] <- 1
    } else {
      res.pc.allinf[sim, marker] <- 0
    }
  }
}

# generate average informative values across 
#  simulations, as function of number of markers
means.vec.related <- vector("numeric", max.markers-2)
means.vec.unrelated <- vector("numeric", max.markers-2)
means.vec.pc <- vector("numeric", max.markers-2)
means.vec.related.homo <- vector("numeric", max.markers-2)
means.vec.unrelated.homo <- vector("numeric", max.markers-2)
means.vec.pc.homo <- vector("numeric", max.markers-2)
means.vec.related.allinf <- vector("numeric", max.markers-2)
means.vec.unrelated.allinf <- vector("numeric", max.markers-2)
means.vec.pc.allinf <- vector("numeric", max.markers-2)
for(calcmean in 3:max.markers){
  means.vec.related[calcmean-2] <- mean(rowSums(res.related[,1:calcmean]) >= 3)
  means.vec.unrelated[calcmean-2] <- mean(rowSums(res.unrelated[,1:calcmean]) >= 3)
  means.vec.pc[calcmean-2] <- mean(rowSums(res.pc[,1:calcmean]) >= 3)
  means.vec.related.homo[calcmean-2] <- mean(rowSums(res.related.homo[,1:calcmean]) >= 3)
  means.vec.unrelated.homo[calcmean-2] <- mean(rowSums(res.unrelated.homo[,1:calcmean]) >= 3)
  means.vec.pc.homo[calcmean-2] <- mean(rowSums(res.pc.homo[,1:calcmean]) >= 3)
  means.vec.related.allinf[calcmean-2] <- mean(rowSums(res.related.allinf[,1:calcmean]) >= 3)
  means.vec.unrelated.allinf[calcmean-2] <- mean(rowSums(res.unrelated.allinf[,1:calcmean]) >= 3)
  means.vec.pc.allinf[calcmean-2] <- mean(rowSums(res.pc.allinf[,1:calcmean]) >= 3)
}
