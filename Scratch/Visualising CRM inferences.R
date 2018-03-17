
# Construct crm_samp, and then...

library(tidyverse)
library(ggridges)

prob_tox_samp <- as.data.frame(crm_samp, 'prob_tox')
prob_tox_samp_tall <- prob_tox_samp %>%
  gather(Label, ProbTox) %>%
  mutate(
    DoseLevel = rep(1:ncol(prob_tox_samp), each = nrow(prob_tox_samp)),
    Draw = rep(1:nrow(prob_tox_samp), times = ncol(prob_tox_samp))
  )
prob_tox_samp_tall %>% head()

(p1 <- ggplot(prob_tox_samp_tall, aes(x = DoseLevel, y = ProbTox, group = DoseLevel)) +
    geom_boxplot() + ylim(0, 1) +
    labs(title = 'boxplot of Pr(DLT) under CRM'))
(p2 <- ggplot(prob_tox_samp_tall, aes(x = DoseLevel, y = ProbTox, group = DoseLevel)) +
    geom_violin(fill = 'orange') + ylim(0, 1) +
    labs(title = 'violin plot of Pr(DLT) under CRM'))
(p3 <- ggplot(prob_tox_samp_tall %>% mutate(DoseLevel = factor(DoseLevel)),
              aes(x = ProbTox, y = DoseLevel, fill = DoseLevel)) +
    geom_density_ridges() + theme(legend.position = 'none') +
    labs(title = 'joyplot of Pr(DLT) under CRM'))
(p4 <- ggplot(prob_tox_samp_tall %>% arrange(Draw) %>% head(5 * 10^3),
              aes(x = DoseLevel, y = ProbTox, group = Draw)) +
  geom_line(alpha = 0.03, col = 'orange') + ylim(0, 1) +
    labs(title = 'overplot of Pr(DLT) under CRM'))

gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2,
                        top = "Four beautiful visualisations of the posterior dose-toxicity curve in a CRM model analysis")


