#######################################
#
# Setting 6:
# MAFs all 0.5, double transplant
# Diploid (i.e. no X/Y)
#
#######################################

# initialize matrices to store simulation results
res.sibling.sibling <- matrix(-1, nrow = nsim, ncol = max.markers)
res.mud.sibling <- matrix(-1, nrow = nsim, ncol = max.markers)
res.mud.mud <- matrix(-1, nrow = nsim, ncol = max.markers)

pop.freq <- rep(0.5, max.markers)
for(sim in 1:nsim){
  cat("sim ", sim, " of ", nsim,"\n")
  for(marker in 1:max.markers){
    # generate parental and offspring genotypes, and p3 (third unrelated genotype)
    p1 <- sample(c("A", "B"), size = 2, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
    p2 <- sample(c("A", "B"), size = 2, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
    p3 <- sample(c("A", "B"), size = 2, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
    o1 <- c(sample(p1, size = 1), sample(p2, size = 1))
    o2 <- c(sample(p1, size = 1), sample(p2, size = 1))
    o3 <- c(sample(p1, size = 1), sample(p2, size = 1))
   
    # true/false informative: donors are identical and homozygous and the recipient
    #  carries at least one copy of the other allele
    if((identical(o1, c("A", "A")) && identical(o2, o3) && identical(o2, c("B", "B"))) ||
       (identical(o1, c("B", "B")) && identical(o2, o3) && identical(o2, c("A", "A"))) ||
       (identical(sort(o1), c("A", "B")) && identical(o2, o3) && identical(o2, c("A", "A"))) ||
       (identical(sort(o1), c("A", "B")) && identical(o2, o3) && identical(o2, c("B", "B")))){
      res.sibling.sibling[sim, marker] <- 1
    } else {
      res.sibling.sibling[sim, marker] <- 0
    }
    
    if((identical(o1, c("A", "A")) && identical(o2, p3) && identical(o2, c("B", "B"))) ||
       (identical(o1, c("B", "B")) && identical(o2, p3) && identical(o2, c("A", "A"))) ||
       (identical(sort(o1), c("A", "B")) && identical(o2, p3) && identical(o2, c("A", "A"))) ||
       (identical(sort(o1), c("A", "B")) && identical(o2, p3) && identical(o2, c("B", "B")))){
      res.mud.sibling[sim, marker] <- 1
    } else {
      res.mud.sibling[sim, marker] <- 0
    }
    if((identical(p1, c("A", "A")) && identical(p2, p3) && identical(p2, c("B", "B"))) ||
       (identical(p1, c("B", "B")) && identical(p2, p3) && identical(p2, c("A", "A"))) ||
       (identical(sort(p1), c("A", "B")) && identical(p2, p3) && identical(p2, c("A", "A"))) ||
       (identical(sort(p1), c("A", "B")) && identical(p2, p3) && identical(p2, c("B", "B")))){
      res.mud.mud[sim, marker] <- 1
    } else {
      res.mud.mud[sim, marker] <- 0
    }
  }
}

# generate average informative values across 
#  simulations, as function of number of markers
means.vec.sib.sib <- vector("numeric", max.markers-2)
means.vec.mud.sib <- vector("numeric", max.markers-2)
means.vec.mud.mud <- vector("numeric", max.markers-2)
for(calcmean in 3:max.markers){
  means.vec.sib.sib[calcmean-2] <- mean(rowSums(res.sibling.sibling[,1:calcmean]) >= 3)
  means.vec.mud.sib[calcmean-2] <- mean(rowSums(res.mud.sibling[,1:calcmean]) >= 3)
  means.vec.mud.mud[calcmean-2] <- mean(rowSums(res.mud.mud[,1:calcmean]) >= 3)
}

