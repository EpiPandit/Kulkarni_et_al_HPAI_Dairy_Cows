library(tidyverse)

beta_dm <-read.csv("./summary_daily_100sims_optiters_dm.csv")
r_dm <- read.csv("./summary_farm_100sims_optiters_dm.csv")

per_farm_seir_dm <- beta_dm %>% group_by(farm, day) %>%
  summarise(mean_S = mean(S), sd_S = sd(S),
            mean_E = mean(E), sd_E = sd(E),
            mean_I = mean(I), sd_I = sd(I),
            mean_R = mean(R_count), sd_R = sd(R_count)) %>%
  mutate_if(is.numeric, floor)



ggplot(per_farm_seir_dm, aes(x = day)) + 
  geom_line(aes(y = mean_E), linewidth = 1.2, color = "orange") +
  geom_line(aes(y = mean_I), linewidth = 1.2, color = "palegreen") +
  geom_line(aes(y = mean_R), linewidth = 1.2, color = "plum") +
  geom_ribbon(aes(ymin = mean_E -  sd_E, y = mean_E,
                ymax = mean_E +  sd_E), alpha = 0.4, fill = "orange") +
  geom_ribbon(aes(ymin = mean_I -  sd_I, y = mean_I,
                ymax = mean_I +  sd_I), alpha = 0.4, fill = "palegreen") +
  geom_ribbon(aes(ymin = mean_R -  sd_R, y = mean_R,
                ymax = mean_R +  sd_R), alpha = 0.4, fill = "plum") +
  facet_wrap(~farm, scales = "free") +
  labs(subtitle = "*** Not showing Susceptible counts (S)",
       y = "E/ I/ R cow counts",
       x = "Day in Simulated Outbreak") +
  theme_classic() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 18, face = "bold"))


beta_plf <-read.csv("./summary_re_100sims_optiters_plf.csv")
r_plf <- read.csv("./summary_farm_100sims_optiters_plf.csv")

per_farm_seir_plf <- beta_plf %>% group_by(farm, day) %>%
  summarise(mean_S = mean(S), sd_S = sd(S),
            mean_E = mean(E), sd_E = sd(E),
            mean_I = mean(I), sd_I = sd(I),
            mean_R = mean(R_count), sd_R = sd(R_count)) %>%
  mutate_if(is.numeric, floor)



ggplot(per_farm_seir_plf, aes(x = day)) +
  geom_line(aes(y = mean_E), linewidth = 1.2, color = "orange") +
  geom_line(aes(y = mean_I), linewidth = 1.2, color = "palegreen") +
  geom_line(aes(y = mean_R), linewidth = 1.2, color = "plum") +
  geom_ribbon(aes(ymin = mean_E -  sd_E, y = mean_E,
                  ymax = mean_E +  sd_E), alpha = 0.4, fill = "orange") +
  geom_ribbon(aes(ymin = mean_I -  sd_I, y = mean_I,
                  ymax = mean_I +  sd_I), alpha = 0.4, fill = "palegreen") +
  geom_ribbon(aes(ymin = mean_R -  sd_R, y = mean_R,
                  ymax = mean_R +  sd_R), alpha = 0.4, fill = "plum") +
  facet_wrap(~farm, scales = "free") +
  labs(subtitle = "*** Not showing Susceptible counts (S)",
       y = "E/ I/ R cow counts",
       x = "Day in Simulated Outbreak") +
  theme_classic() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 18, face = "bold"))
ggsave("./SEIR_sim_curve_plf.png",
       dpi = 300, height = 12, width = 12)


r_all$low_beta_eff = ifelse(r_all$mean_beta_eff - 1.96*r_all$sd_beta_eff >= 0,r_all$mean_beta_eff - 1.96*r_all$sd_beta_eff, 0) 
r_all$hi_beta_eff = ifelse(r_all$mean_beta_eff + 1.96*r_all$sd_beta_eff <= 0.05,r_all$mean_beta_eff + 1.96*r_all$sd_beta_eff, 0.05) 

ggplot(r_all, aes(x = farm, y = mean_beta_eff, color = farm, group = method)) +
  geom_point(aes(shape = method), size = 5, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = low_beta_eff, ymax = hi_beta_eff), lineend = "square", linewidth = 1.75,
                width = 0.25, alpha = 0.6,
                position = position_dodge(width = 0.6),) +
  geom_text(aes(label = format(round(mean_beta_eff, 4), scientific = F)), vjust = -3, position = position_dodge(width = 0.9),
            check_overlap = T, show.legend = F, fontface = "bold", size = 5) +
  scale_color_brewer(type = "div") +
  scale_y_continuous(limits = c(0, 0.05), labels = seq(0, 0.05, 0.01), breaks = seq(0, 0.05, 0.01)) +
  theme_classic() +
  labs(x = "Farm", y = expression(paste("Mean ", beta[eff])), color = "Farm", shape = "Loss Function") +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = "bottom")
ggsave("./Figures/eff_beta_dm_plf_100iters.png", 
       dpi = 300, height = 6, width = 10)

ggplot(r_all, aes(x = farm, y = mean_peak_I, color = farm, group = method)) +
  geom_point(aes(shape = method), size = 5, position = position_dodge(width = 0.6)) +
  geom_text(aes(label = format(round(mean_peak_I, 0), scientific = F)), vjust = -1, position = position_dodge(width = 0.6),
            check_overlap = T, show.legend = F, fontface = "bold", size = 5) +
  scale_color_brewer(type = "div") +
  # scale_y_continuous(limits = c(0, 0.08), labels = seq(0, 0.08, 0.005), breaks = seq(0, 0.08, 0.005)) +
  theme_classic() +
  labs(x = "Farm", y = expression(paste("Mean ", beta[eff])), color = "Farm", shape = "Loss Function") +
  theme(axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "bottom")


#####################################################################

beta_dm <-read.csv("./summary_daily_100sims_optiters_dm.csv")
beta_plf <-read.csv("./summary_re_100sims_optiters_plf.csv")
library(tidyverse)
################################15days#####################################
max_period_df <- tibble(
  farm = c("Farm1",
           "Farm2",
           "Farm3",
           "Farm4"),
#   max_period = c(33,
#                  53,
#                  23,
#                  18)
# )
  max_period = rep(15, 4))


gamma_df_dm <- beta_dm %>%
  group_by(farm, sim_id) %>%
  arrange(day) %>%
  mutate(
    daily_recov = R_count - lag(R_count),
    daily_recov = ifelse(is.na(daily_recov), 0, pmax(daily_recov, 0))
  ) %>%
  summarise(
    gamma = sum(daily_recov, na.rm = TRUE) / sum(I, na.rm = TRUE),
    .groups = "drop"
  )

rt_daily_dm <- beta_dm %>%
  left_join(max_period_df, by = "farm") %>%
  filter(day <= max_period) %>%
  group_by(farm, sim_id, day) %>%
  summarise(
    S = first(S),
    E = first(E),
    I = first(I),
    R_count = first(R_count),
    new_E = sum(new_E, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(gamma_df_dm, by = c("farm", "sim_id")) %>%
  group_by(farm, sim_id) %>%
  mutate(
    N = S + E + I + R_count,
    Rt = ifelse(I > 0 & gamma > 0, (new_E / I) / gamma * (N / S), NA_real_)
  ) %>%
  ungroup()

gamma_df_plf <- beta_plf %>%
  group_by(farm, sim_id) %>%
  arrange(day) %>%
  mutate(
    daily_recov = R_count - lag(R_count),
    daily_recov = ifelse(is.na(daily_recov), 0, pmax(daily_recov, 0))
  ) %>%
  summarise(
    gamma = sum(daily_recov, na.rm = TRUE) / sum(I, na.rm = TRUE),
    .groups = "drop"
  )

rt_daily_plf <- beta_plf %>%
  left_join(max_period_df, by = "farm") %>%
  filter(day <= max_period) %>%
  group_by(farm, sim_id, day) %>%
  summarise(
    S = first(S),
    E = first(E),
    I = first(I),
    R_count = first(R_count),
    new_E = sum(new_E, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(gamma_df_plf, by = c("farm", "sim_id")) %>%
  group_by(farm, sim_id) %>%
  mutate(
    N = S + E + I + R_count,
    Rt = ifelse(I > 0 & gamma > 0, (new_E / I) / gamma * (N / S), NA_real_)
  ) %>%
  ungroup()
rt_daily_dm$method = rep("dm", nrow(rt_daily_dm))
rt_daily_plf$method = rep("plf", nrow(rt_daily_plf))

rt_daily <- do.call("rbind", list(rt_daily_dm, rt_daily_plf))

rt_summary <- rt_daily %>%
  group_by(farm, day, method) %>%
  summarise(
    n = sum(!is.na(Rt)),
    Rt_mean = mean(Rt, na.rm = TRUE),
    Rt_sd = sd(Rt, na.rm = TRUE),
    Rt_se = Rt_sd / sqrt(n),
    Rt_low = ifelse(n > 1, Rt_mean - qt(0.975, df = n - 1) * Rt_se, NA_real_),
    Rt_high = ifelse(n > 1, Rt_mean + qt(0.975, df = n - 1) * Rt_se, NA_real_),
    .groups = "drop"
  )

ggplot(rt_summary, aes(x = day, y = Rt_mean, color = farm, group = farm)) +
  geom_line(stat = "identity") +
  geom_ribbon(aes(ymin = Rt_low, ymax = Rt_high, x = day, fill = farm, group = farm), alpha = 0.3, stat = "identity", linewidth = 0) +
  facet_grid(cols = vars(farm), rows = vars(method), scales = "free") +
  labs(x = "Days of observed Outbreak", y = expression(paste("Daily Reproduction Number (", R[daily],")")),
       fill = "Farms", color = "Farms") +
  scale_fill_brewer(type = "div") +
  scale_color_brewer(type = "div") +
  theme_classic() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major.y = element_line(size = 0.4, linetype = 2, color = "grey50"))
ggsave("./Figures/R0_dm_plf_100iters.png",
       dpi = 300, height = 6, width = 10)

R0_15_day <- rt_summary %>% group_by(farm, method) %>% summarise(mean_R0 = mean(Rt_mean), sd_R0 = mean(Rt_sd)) %>%
  pivot_wider(names_from = method, values_from = mean_R0:sd_R0) %>% 
  select(farm, mean_R0_dm, sd_R0_dm, mean_R0_plf, sd_R0_plf)


################################Outbreak_period#####################################


max_period_df <- tibble(
  farm = c("Farm1",
           "Farm2",
           "Farm3",
           "Farm4"),
    max_period = c(33,
                   53,
                   23,
                   18)
  )
  # max_period = rep(10, 4))


gamma_df_dm <- beta_dm %>%
  group_by(farm, sim_id) %>%
  arrange(day) %>%
  mutate(
    daily_recov = R_count - lag(R_count),
    daily_recov = ifelse(is.na(daily_recov), 0, pmax(daily_recov, 0))
  ) %>%
  summarise(
    gamma = sum(daily_recov, na.rm = TRUE) / sum(I, na.rm = TRUE),
    .groups = "drop"
  )

rt_daily_dm <- beta_dm %>%
  left_join(max_period_df, by = "farm") %>%
  filter(day <= max_period) %>%
  group_by(farm, sim_id, day) %>%
  summarise(
    S = first(S),
    E = first(E),
    I = first(I),
    R_count = first(R_count),
    new_E = sum(new_E, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(gamma_df_dm, by = c("farm", "sim_id")) %>%
  group_by(farm, sim_id) %>%
  mutate(
    N = S + E + I + R_count,
    Rt = ifelse(I > 0 & gamma > 0, (new_E / I) / gamma * (N / S), NA_real_)
  ) %>%
  ungroup()

gamma_df_plf <- beta_plf %>%
  group_by(farm, sim_id) %>%
  arrange(day) %>%
  mutate(
    daily_recov = R_count - lag(R_count),
    daily_recov = ifelse(is.na(daily_recov), 0, pmax(daily_recov, 0))
  ) %>%
  summarise(
    gamma = sum(daily_recov, na.rm = TRUE) / sum(I, na.rm = TRUE),
    .groups = "drop"
  )

rt_daily_plf <- beta_plf %>%
  left_join(max_period_df, by = "farm") %>%
  filter(day <= max_period) %>%
  group_by(farm, sim_id, day) %>%
  summarise(
    S = first(S),
    E = first(E),
    I = first(I),
    R_count = first(R_count),
    new_E = sum(new_E, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(gamma_df_plf, by = c("farm", "sim_id")) %>%
  group_by(farm, sim_id) %>%
  mutate(
    N = S + E + I + R_count,
    Rt = ifelse(I > 0 & gamma > 0, (new_E / I) / gamma * (N / S), NA_real_)
  ) %>%
  ungroup()
rt_daily_dm$method = rep("dm", nrow(rt_daily_dm))
rt_daily_plf$method = rep("plf", nrow(rt_daily_plf))

rt_daily <- do.call("rbind", list(rt_daily_dm, rt_daily_plf))

rt_summary <- rt_daily %>%
  group_by(farm, day, method) %>%
  summarise(
    n = sum(!is.na(Rt)),
    Rt_mean = mean(Rt, na.rm = TRUE),
    Rt_sd = sd(Rt, na.rm = TRUE),
    Rt_se = Rt_sd / sqrt(n),
    Rt_low = ifelse(n > 1, Rt_mean - qt(0.975, df = n - 1) * Rt_se, NA_real_),
    Rt_high = ifelse(n > 1, Rt_mean + qt(0.975, df = n - 1) * Rt_se, NA_real_),
    .groups = "drop"
  )

R0_outbreak_day <- rt_summary %>% group_by(farm, method) %>% summarise(mean_R0 = mean(Rt_mean), sd_R0 = mean(Rt_sd)) %>%
  pivot_wider(names_from = method, values_from = mean_R0:sd_R0) %>% 
  select(farm, mean_R0_dm, sd_R0_dm, mean_R0_plf, sd_R0_plf)


##########################beta_eff###################################

a <- beta_dm %>% group_by(farm) %>% summarise(beta = mean(T_nominal), beta_sd = sd(T_nominal))
b <- beta_plf%>% group_by(farm) %>% summarise(beta = mean(T_nominal), beta_sd = sd(T_nominal))

a <- do.call("rbind", list(a,b))
a$method <- rep(c("dm", "plf"), each = 4)

ggplot(a, aes(y = farm, color = method)) +
  # geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = beta - 1.96*beta_sd, xmax = beta + 1.96*beta_sd), 
                height = 0.2, show.legend = F, position = position_dodge(width = 0.6)) +
  geom_label(aes(x = beta, label = round(beta, 3)), 
             vjust = 0.5, hjust = 0.5, size = 4.5, show.legend = F,
             position = position_dodge(width = 0.6),fontface = "bold") +
  labs(x = expression(beta[eff]), y = NULL, color = "Method") +
  scale_fill_brewer(type = "div") +
  theme_classic() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = "top",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major.x = element_line(size = 0.4, linetype = 2, color = "grey50"))
ggsave("./Figures/beta_dm_plf.png",
       dpi = 300, height = 6, width = 6)
mean(a$beta)
sd(a$beta_sd)
