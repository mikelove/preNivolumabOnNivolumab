library(DESeq2)
library(pbapply)
library(ggplot2)
library(tidyr)

mult <- seq(from=0, to=.4, length=50)
out <- pblapply(mult, function(m) {
  n <- 500
  z <- numeric(3)
  dds <- makeExampleDESeqDataSet(n=1000,m=n,interceptMean=10)
  dds$x <- as.numeric(dds$condition) - 1
  dds$batch <- rnorm(n) + m * dds$x
  z[1] <- cor(dds$x, dds$batch)
  new_counts <- round( t( t(counts(dds)[1:100,]) * exp(dds$batch * rep(rnorm(100),each=n)) ) )
  mode(new_counts) <- "integer"
  counts(dds)[1:100,] <- new_counts
  dds <- DESeq(dds, quiet=TRUE, fitType="mean")
  z[2] <- sum(results(dds)$padj < .1, na.rm=TRUE)
  design(dds) <- ~batch + condition
  dds <- DESeq(dds, quiet=TRUE, fitType="mean")
  z[3] <- sum(results(dds)$padj < .1, na.rm=TRUE)
  return(z)
})

tab <- do.call(rbind, out)
tab <- as.data.frame(tab)
colnames(tab) <- c("cor","condition","batch+condition")
tab$cor <- abs(tab$cor)
dat <- pivot_longer(tab, cols=c("condition","batch+condition"), names_to="design", values_to="FP")
ggplot(dat, aes(cor, FP, color=design)) +
  geom_point() +
  stat_smooth(method="gam") +
  ggtitle("Batch correlated with condition induces FP in mis-specified model (n=500 samples)") +
  scale_color_discrete(breaks=c("condition","batch+condition")) +
  xlab("correlation of batch variable with condition") +
  ylab("number of false positives (out of 100)")

write.csv(tab, file="batch_sim.csv")
