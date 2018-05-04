library(ggplot2);
## power simulation
beta <- 0.8; ## mean beta for effect size

##error <- 0.1; ## mean sd
## number of tests
M <- 5000;
num.sims <- 100;

##N <- seq(1:500);
N <- seq(20, 250, 5);

##maf <- 02;

power.df = c();

## fix maf
##for(beta in c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1))

##r.squared <- 0.1;

for(maf in c(0.1, 0.2, 0.3, 0.4)) {
        cat("maf: ",maf,"\n");
    for(N.i in N) {
            cat("\t","N: ",N.i,"\n");
        geno <- c(rep(0,N.i*maf^2), rep(1, 2*N.i*maf*(1-maf)), rep(2, N.i*(1-maf)^2));

        power = 0;
        rsquare = 0;

        for (i in 1:num.sims) {
          ##y <- sqrt(r.squared/var(geno))*geno + rnorm(length(geno), 0, sqrt(1-var(sqrt(r.squared/var(geno))*geno)));
          y <- beta*geno + rnorm(length(geno), 0, sqrt(1-var(beta*geno)));

            ##browser();
            if(anova(lm(y~geno))$"Pr(>F)"[[1]]<0.05/100000) {
                power = power+ 1;
            }
            ##print(summary(lm(y~geno))$adj.r.squared);
            rsquare <- rsquare+summary(lm(y~geno))$adj.r.squared;
        }

        power.df <- rbind(power.df, c(power=power/num.sims, rsquare=rsquare/num.sims, maf=maf, N=N.i));

    }
}

power.df <- data.frame(power.df)

ggsave("power.pdf", ggplot(aes(N,power, color=factor(maf)), data=data.frame(power.df))+geom_line()+geom_point()+theme_bw());

## beta^2/r^2 = (sd(y)/sd(x))^2

## beta/r = sd(y)/sd(x)

## r = beta/(sd(y)/sd(x))

## beta = r*(sd(y)/sd(x))

## beta = 0.1*(sd(y)/sd(x))

## beta = 0.1*sd(y)
