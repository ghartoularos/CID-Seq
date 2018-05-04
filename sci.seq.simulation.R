w1 <- 384;
w2 <- 384;
n1 <- 2000;
n2 <- 300;

r0 <- NULL;
r1 <- NULL;
r2 <- NULL;
rm <- NULL;
cr <- NULL;

for(i in 1:1) {
 print(i)
  bc1 <- round(runif(n1*w1,1,w1))
  bc2 <- round(runif(n2*w2,1,w2))
  bc1.small <- bc1[1:(w2*n2)]

  all.bc <- paste(bc1.small, bc2,sep=".")
  a.table <- table(as.factor(all.bc))

  for(bc2.i in bc2) {
    c0 <- w1-length(unique(bc1.small[which(bc2.i==bc2)]));
    c1 <- length(which(table(bc1.small[which(bc2.i==bc2)])==1));
    c2 <- length(which(table(bc1.small[which(bc2.i==bc2)])==2));

    r0 <- c(r0, c0/w1);
    r1 <- c(r1, c1/w1);
    rm <- c(rm, (w1-c0-c1)/w1);
    r2 <- c(r2, c2/w1);
  }
}

cr <- rm*w1/((rm+r1)*w1);

ec0 = w1 * (1.0 - 1.0 / w1)^n2;
ec1 = n2 * (1.0 - 1.0 / w1)^(n2 - 1)
ecm =  w1 - ec0 - ec1
ecr = ecm / (ec1 + ecm)

## what i'm calculating is the probability any given w1 barcode is not presented once
## which could be represented 0 times or more than once
jr1 = (1-1/w1)^(n2-1);
jr0 = (1-1/w1)^(n2);
jrno01 = 1-jr1-jr0;

jcr <- jrno01/jr1;
