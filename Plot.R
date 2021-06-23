library(here)
setwd(here())

library(tidyverse)
library(viridis)
library(ggrepel)

tests = read_tsv('./test_results/Test_results_2.tsv')
tests = tests %>% mutate(Size = Variants * Samples)

tests1 = read_tsv('./test_results/Test_results_fixed_lr.tsv')
tests1 = tests1 %>% mutate(Size = Variants * Samples)

tests_long = tests %>% 
  filter(MSE_SQRT <= 1) %>%
  select(Size, Exec_Time, F1_score, MSE_SQRT, Precision, Recall) %>%
  pivot_longer(-Size, names_to='Test', values_to='Score')

plot_1 = tests_long %>%
  ggplot(aes(Size, Score)) +
  geom_boxplot(aes(group = cut_width(Size, 25)), outlier.alpha = 0) +
  geom_smooth(color = "blue") +
  geom_point(alpha = 0.2) +
  facet_grid(Test ~ ., scales = "free_y") +
  theme_bw()

ggsave('./img/Size_effect_dlr.png', height = 8, width = 6)

plot_2 = tests %>% ggplot(aes(Variants, Samples, color=F1_score, label=F1_score)) + 
  geom_point(size = 4, alpha = 0.5) +
  geom_text_repel(size = 2) +
  scale_color_viridis(name = "F1 score") +
  scale_y_continuous(breaks = 1:10, minor_breaks = F) +
  xlab('Number of variants') +
  ylab('Number of samples') +
  theme_bw()

ggsave('./img/F1_Size.png', height = 6, width = 6)

plot_3 = tests %>% filter(MSE_SQRT <= 1) %>% ggplot(aes(F1_score, MSE_SQRT)) +
  geom_density2d_filled() +
  geom_point(shape = 1, color = "white") +
  theme_minimal() +
  xlab('F1 Score') +
  ylab('Square root of MSE')

ggsave('./img/F1_MSE.png', height = 6, width = 6)

tests_long_2 = tests %>% 
  filter(MSE_SQRT <= 1) %>%
  select(Variants, Samples, Subpopulations, Precision, Recall, F1_score, MSE_SQRT) %>%
  pivot_longer(-c(Variants, Samples, Subpopulations),names_to='Score', values_to='Value') %>%
  pivot_longer(-c(Score, Value), names_to='Dimention', values_to='Size')
tests_long_2$LR = 'MIP'

tests_long_1 = tests1 %>% 
  filter(MSE_SQRT <= 1) %>%
  select(Variants, Samples, Subpopulations, Precision, Recall, F1_score, MSE_SQRT) %>%
  pivot_longer(-c(Variants, Samples, Subpopulations),names_to='Score', values_to='Value') %>%
  pivot_longer(-c(Score, Value), names_to='Dimention', values_to='Size')
tests_long_1$LR = 'Gradient descend'

tests_long = rbind(tests_long_1, tests_long_2)

plot_4 = tests_long %>% ggplot(aes(Size, Value, color = LR)) +
  geom_point(alpha=0.3, size = 0.5) +
  geom_smooth(alpha = 0.1) +
  facet_grid(Score ~ Dimention, scales = 'free')

ggsave('./img/GD_vs_MIP.png', height = 6, width = 6)

plot_5 = tests %>% ggplot(aes(Variants, Subpopulations)) + 
  geom_point(alpha = 0.3) +
  geom_abline(slope = 2)
