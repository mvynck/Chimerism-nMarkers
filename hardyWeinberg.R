library(qvalue)
hardy <- read.csv("Data/hardy.csv")
pval <- rep(1, 24)
for(marker in 1:nrow(hardy)){
  p <- hardy[marker ,1]*2/(218*2)+hardy[marker ,2]/(218*2)
  q <- hardy[marker ,3]*2/(218*2)+hardy[marker ,2]/(218*2)
  chisq.val <- (hardy[marker, 1]-p^2*218)^2/(p^2*218) + (hardy[marker, 2]-2*p*q*218)^2/(2*p*q*218) + (hardy[marker, 3]-q^2*218)^2/(q^2*218)
  pval[marker] <- pchisq(chisq.val, 1, lower.tail = FALSE)
}
round(qvalue(pval)$qvalues, 2)