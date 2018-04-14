# one - sided -----
threshold <- c(-Inf, 1.96)
y <- seq(from = 2, to = 5, by = 0.1)
est <- sapply(y, univNormMLE, sd = 1, threshold = threshold)
oneside <- data.frame(y = y, est = est, threshold = threshold[2],
                      selection = "one sided")

# two sided ----
threshold <- c(-1.96, 1.96)
est <- sapply(y, univNormMLE, sd = 1, threshold = threshold)
twoside <- data.frame(y = y, est = est, threshold = threshold[2],
                      selection = "two sided")

library(ggplot2)
forplot <- rbind(oneside, twoside)
ggplot(forplot) +
  geom_line(aes(x = y, y = est), col = "red") +
  facet_wrap(~ selection, labeller = "label_both") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_vline(xintercept = threshold[2], col = "grey", linetype = 1) +
  ylim(-5, 5) +
  geom_hline(yintercept = 0) +
  ylab("conditional MLE") + xlab("observed")
ggsave(file = "figures/univest.pdf", units = "cm",
       width = 12, height = 6)


