library(ggplot2);

## goal is to sequence 1000 cells per person

inds <- seq(1,100,1);

## simulate cid-seq
## assume 500 (tn5) + 934 rt + 1150 master mix + 300 = 2884
## cost per drop = 657/50000+(300)/50000 = 0.03
w1 <- 16;

cells <- NULL;
crs <- NULL;
cost.prep <- NULL;
cost.seq <- NULL;
singlets <- NULL;
collisions <- NULL;

##lambda.cid <- 1.6;

for(ind in inds) {
  lambda.cid <- ind*1000/50000;
  n2<-rpois(50000,lambda.cid);
  cells <- c(cells, sum(n2));
  bc0 <- w1 * (1-1/w1)^n2;
  bc1 <- n2*(1-1/w1)^(n2-1);
  bc2 <- w1 - bc0 - bc1;
  tmp <- bc2/(bc1+bc2)*n2;
  cr <- sum(na.omit(tmp))/sum(n2[!is.na(tmp)])*(1-1/ind);
  crs <- c(crs, cr);

  singlets <- c(singlets, sum(n2)*(1-cr));
  collisions <- c(collisions, sum(n2)*cr);

  ##cost.prep <- c(cost.prep,0.15*singlet(657*sum(n2)/10000+500+300)/ind/(1-cr));
  ##cost.seq <- c(cost.seq, sum(n2)*50000/350e6*2000/(1-cr)/ind);
}

cost.prep <- 1000/(singlets/inds)*50000*((657+300)/50000)/inds;
cost.seq <- 1000/(singlets/inds)*1000*50000/350e6*2000;

## simluate sci-seq
## let's assume some costs
## $8334 illumina tn5
## $3834 own tn5
## cost estimate is 3834/384 = $10/well
w1 <- 384;

cells.sci <- NULL;
crs.sci <- NULL;
cost.prep.sci <- NULL;
cost.seq.sci <- NULL;
singlets.sci <- NULL;
collisions.sci <- NULL;

##n2 <- 40;

for(ind in inds) {
  n2<-rep(1000*ind/w1, w1);
  cells.sci <- c(cells.sci, sum(n2));
  bc0 <- w1 * (1-1/w1)^n2;
  bc1 <- n2*(1-1/w1)^(n2-1);
  bc2 <- w1 - bc0 - bc1;
  tmp <- bc2/(bc1+bc2)*n2;
  cr <- sum(na.omit(tmp))/sum(n2[!is.na(tmp)]);
  crs.sci <- c(crs.sci, cr);

  singlets.sci <- c(singlets.sci, sum(n2)*(1-cr));
  collisions.sci <- c(collisions.sci, sum(n2)*cr);

  ##cost.prep.sci <- c(cost.prep.sci, 3834/ind/(1-cr));
  ##cost.seq.sci <- c(cost.seq.sci, (sum(n2)/(1-cr))*50000/350e6*2000/ind);
}

cost.prep.sci <- 1000/(singlets.sci/inds)*384*(3834/384)/inds;
cost.seq.sci <- 1000/(singlets.sci/inds)*1000*50000/350e6*2000;

## simluate drop-seq
## let's assume 10x costs
## cost per prep is $657/10k cells
## cost per prep is $657/50k drops = 0.02
##

cells.drop <- NULL;
crs.drop <- NULL;
cost.prep.drop <- NULL;
cost.seq.drop <- NULL;
singlets.drop <- NULL;
collisions.drop <- NULL;

for(ind in inds) {
  lambda <- ind*1000/50000;
  n2<-rpois(50000,lambda);

# for(lambda in c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,0.2,0.3,0.4)) {
#     n2<-rpois(50000,lambda);

  cells.drop <- c(cells.drop, length(which(n2>0)));
  cr <- length(which(n2>1))/length(which(n2>0));
  crs.drop <- c(crs.drop, cr);

  singlets.drop <- c(singlets.drop, length(which(n2>0))*(1-cr));
  collisions.drop <- c(collisions.drop, length(which(n2>0))*cr);

  ##cost.prep.drop <- c(cost.prep.drop, 657/10000*length(which(n2>0))/ind/(1-cr)*(ind*1000/length(which(n2>0))));
  ##cost.seq.drop <- c(cost.seq.drop, (length(which(n2>0))/(1-cr))*50000/350e6*2000/ind*(ind*1000/length(which(n2>0))));

}

cost.prep.drop <- 1000/(singlets.drop/inds)*50000*(1500*2/50000)/inds;
cost.seq.drop <- 1000/(singlets.drop/inds)*1000*50000/350e6*2000;


## plot everything
plot.df <- rbind(data.frame(cost=cost.prep.drop+cost.seq.drop, inds=inds, method="drop-seq",type="total"),
data.frame(cost=cost.prep.sci+cost.seq.sci, inds=inds, method="sci-seq",type="total"),
data.frame(cost=cost.prep+cost.seq, inds=inds, method="cid-seq",type="total"),
data.frame(cost=cost.prep.drop, inds=inds, method="drop-seq", type="prep"),
data.frame(cost=cost.seq.drop, inds=inds, method="drop-seq", type="seq"),
data.frame(cost=cost.prep.sci, inds=inds, method="sci-seq", type="prep"),
data.frame(cost=cost.seq.sci, inds=inds, method="sci-seq", type="seq"),
data.frame(cost=cost.prep, inds=inds, method="cid-seq", type="prep"),
data.frame(cost=cost.seq, inds=inds, method="cid-seq", type="seq"))

plot.total.df <- rbind(data.frame(cost=cost.prep.drop+cost.seq.drop, inds=inds, method="drop-seq"),
data.frame(cost=cost.prep.sci+cost.seq.sci, inds=inds, method="sci-seq"),
data.frame(cost=cost.prep+cost.seq, inds=inds, method="cid-seq"))

ggsave("cost.vs.ind.pdf",ggplot(aes(inds, cost, linetype=type,color=method),data=plot.df)+geom_line()+theme_bw()+scale_y_log10(breaks=c(0,1,2,5,10,20,50,100,200,500,1000,2000,5000))+scale_x_continuous(breaks=seq(0,100,10)));

##geom_bar(stat="identity",position="dodge")+theme_bw());

#+scale_y_log10(breaks=seq(0,0.2,0.02));
