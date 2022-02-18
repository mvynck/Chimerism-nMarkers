#######################################
#
# Setting 5:
# MAF 0.5
# X chromosome (donor-recipient)
#  sibling: brother-brother, sister-sister, brother-sister
#  or unrelated: female-female, female-male, male-male
# Y chromosome, male-male
#   i.e. male-male (unrelated)
#   
# combined fully informative and potentially informative
#   markers considered
#
#######################################

# compare to diploid marker
# related informativity (MAF 0.5): 41%
mean(res.related.ideal.allinf)
# unrelated informativity (MAF 0.5): 63%
mean(res.unrelated.ideal.allinf)

res.related.X <- matrix(-1, nrow = nsim, ncol = 3)
res.unrelated.X <- matrix(-1, nrow = nsim, ncol = 3)
res.unrelated.Y <- matrix(-1, nrow = nsim, ncol = 1)
set.seed(1)
pop.freq = c(0.5, 0.5)

for(sim in 1:nsim){
  cat("sim ", sim, " of ", nsim,"\n")
  for(marker in 1:2){
    # generate parental and offspring genotypes
    if(marker == 1){
      # X chromosome
      # mother, XX
      p1 <- sample(c("A", "B"), size = 2, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
      # female 2, XX
      p1b <- sample(c("A", "B"), size = 2, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
      # father, XY
      p2 <- sample(c("A", "B"), size = 1, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
      # male 2, XY
      p2b <- sample(c("A", "B"), size = 1, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
      # daughter 1, XX
      o1 <- c(sample(p1, size = 1), sample(p2, size = 1))
      # daughter 2, XX
      o1b <- c(sample(p1, size = 1), sample(p2, size = 1))
      # son 1, XY
      o2 <- c(sample(p1, size = 1))
      # son 2, XY
      o2b <- c(sample(p1, size = 1))
      
      # true/false informative brother-sister, XX-XY
      if(o1 == c("A", "A") && o2 == c("B") ||
         sort(o1) == c("A", "B") && o2 == c("A") ||
         sort(o1) == c("A", "B") && o2 == c("B") ||
         o1 == c("B", "B") && o2 == c("A")){
        res.related.X[sim, 1] <- 1
      }
      else {
        res.related.X[sim, 1] <- 0
      }
      # true/false informative brother-brother, XY-XY
      if(o2 != o2b){
        res.related.X[sim, 2] <- 1
      } else {
        res.related.X[sim, 2] <- 0
      }
      # true/false informative sister-sister, XX-XX
      if(sum(sort(o1) != sort(o1b)) == 2 ||
         (sum(sort(o1) == c("A", "B")) == 2 && (sum(o1b == c("A", "A")) == 2 ||
                                                sum(o1b == c("B", "B")) == 2 )) ||
         (sum(sort(o1b) == c("A", "B")) == 2 && (sum(o1 == c("A", "A")) == 2 ||
                                                 sum(o1 == c("B", "B")) == 2 ))){
        
        res.related.X[sim, 3] <- 1
      } else {
        res.related.X[sim, 3] <- 0
      }
      
      # true/false informative female-male, XX-XY
      if(p1 == c("A", "A") && p2 == c("B") ||
         sort(p1) == c("A", "B") && p2 == c("A") ||
         sort(p1) == c("A", "B") && p2 == c("B") ||
         p1 == c("B", "B") && p2 == c("A")){
        res.unrelated.X[sim, 1] <- 1
      } else {
        res.unrelated.X[sim, 1] <- 0
      }
      # true/false informative male-male, XY-XY
      if(p2 != p2b){
        res.unrelated.X[sim, 2] <- 1
      } else {
        res.unrelated.X[sim, 2] <- 0
      }
      # true/false informative female-female, XX-XX
      if(sum(sort(p1) != sort(p1b)) == 2 ||
         (sum(sort(p1) == c("A", "B")) == 2 && (sum(p1b == c("A", "A")) == 2 ||
                                                sum(p1b == c("B", "B")) == 2 )) ||
         (sum(sort(p1b) == c("A", "B")) == 2 && (sum(p1 == c("A", "A")) == 2 ||
                                                 sum(p1 == c("B", "B")) == 2 ))){
        res.unrelated.X[sim, 3] <- 1
      } else {
        res.unrelated.X[sim, 3] <- 0
      }
    }
    
    if(marker == 2){
      # Y chromosome
      # male 1
      p1 <- sample(c("A", "B"), size = 1, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
      # male 2
      p2 <- sample(c("A", "B"), size = 1, replace = TRUE, prob = c(pop.freq[marker], 1-pop.freq[marker]))
      
      # true/false informative
      if(p1 != p2){
        res.unrelated.Y[sim, 1] <- 1
      } else {
        res.unrelated.Y[sim, 1] <- 0
      }
    }   
    
  }
}

colMeans(res.related.X)
colMeans(res.unrelated.X)
colMeans(res.unrelated.Y)
# X chr:
#   brother-sister: 81%
#   brother-brother: 25%
#   sister-sister: 25%
#
#   male-female: 87%
#   male-male: 51%
#   female-female: 63%
#
# Y chr:
#   male-male: 50%
