library(ggplot2);

cr <- function(n2,w1) {
  bc0 <- w1 * (1-1/w1)^n2;
  bc1 <- n2*(1-1/w1)^(n2-1);
  bc2 <- w1 - bc0 - bc1;
  return(bc2/(bc2+bc1));
}


min.rate <- function(n2,w1,rate) {
  bc0 <- w1 * (1-1/w1)^n2;
  bc1 <- n2*(1-1/w1)^(n2-1);
  bc2 <- w1 - bc0 - bc1;
  return(abs(bc2/(bc2+bc1)-rate));
}

outputs.sci <- NULL;
w1s <- seq(1,1500,1);

for(w1 in w1s) {
  outputs.sci <- c(outputs.sci, w1*log10(optim(10, min.rate, w1=w1, rate=0.05)$par)*0.95);
}


min.rate.cid <- function(lambda,rate) {
  n2<-rpois(50000,lambda);
  bc0 <- w1 * (1-1/w1)^n2;
  bc1 <- n2*(1-1/w1)^(n2-1);
  bc2 <- w1 - bc0 - bc1;
  tmp <- bc2/(bc1+bc2)*n2;
  cr <- sum(na.omit(tmp))/sum(n2[!is.na(tmp)]);
  return(abs(cr-rate));
}

outputs.cid <- NULL;

for(lambda in seq(0.01, 2, 0.01)) {
  print(lambda);
  lambda <-
  outputs.cid <- c(outputs.cid, log10(optim(10, min.rate.cid, rate=0.05)$par)*0.95);

}

cells <- NULL;

for(lambda in seq(0.01, 2, 0.01)) {
  cells <- c(cells, sum(rpois(50000,lambda)))
}
