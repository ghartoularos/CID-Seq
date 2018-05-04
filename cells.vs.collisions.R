library(ggplot2);

## simulate cid-seq
w1 <- 16;

cells <- NULL;
crs <- NULL;
crs.mux <- NULL;
cr2s <- NULL;

for(lambda in c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2)) {
  n2<-rpois(50000,lambda);
  cells <- c(cells, sum(n2));
  bc0 <- w1 * (1-1/w1)^n2;
  bc1 <- n2*(1-1/w1)^(n2-1);
  bc2 <- w1 - bc0 - bc1;
  tmp <- bc2/(bc1+bc2)*n2;
  cr <- sum(na.omit(tmp))/sum(n2[!is.na(tmp)]);
  crs <- c(crs, cr);
}


## simluate sci-seq
w1 <- 384;

cells.sci <- NULL;
crs.sci <- NULL;
crs2.sci <- NULL;

for(n2 in c(1,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200)) {
  n2<-rep(n2, w1);
  cells.sci <- c(cells.sci, sum(n2));
  bc0 <- w1 * (1-1/w1)^n2;
  bc1 <- n2*(1-1/w1)^(n2-1);
  bc2 <- w1 - bc0 - bc1;
  tmp <- bc2/(bc1+bc2)*n2;
  cr <- sum(na.omit(tmp))/sum(n2[!is.na(tmp)]);
  crs.sci <- c(crs.sci, cr);
  cr2 <- sum(bc2)/sum(bc1+bc2);
  crs2.sci <- c(crs2.sci, cr2);

}

## simluate drop-seq
cells.drop <- NULL;
crs.drop <- NULL;

cells2.drop <- NULL;
crs2.drop <- NULL;

#
for(lambda in c(0.01,0.02,0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,1.1,1.2,1.3,2,3,4,5)) {
   n2<-rpois(50000,lambda);

# for(lambda in c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,0.2,0.3,0.4)) {
#     n2<-rpois(50000,lambda);

  cells.drop <- c(cells.drop, length(which(n2>0)));
  cr <- length(which(n2>1))/length(which(n2>0));
  crs.drop <- c(crs.drop, cr);

  cells2.drop <- c(cells2.drop, sum(n2))
  cr2 <- sum(n2[which(n2>1)])/sum(n2);
  crs2.drop <- c(crs2.drop, cr2);

}

## plot everything
plot.df <- rbind(data.frame(cells=cells.drop, collision.rate=crs.drop, method="drop-seq"),
data.frame(cells=cells.sci, collision.rate=crs.sci, method="sci-seq"),
data.frame(cells=cells, collision.rate=crs, method="cid-seq"))

ggsave("collision.vs.cells.pdf",ggplot(aes(cells, collision.rate, color=method),data=plot.df)+geom_line()+geom_point()+scale_y_continuous(limits=c(0,0.1),breaks=seq(0,0.1,0.01))+theme_bw())
#+scale_y_log10(breaks=seq(0,0.2,0.02));
