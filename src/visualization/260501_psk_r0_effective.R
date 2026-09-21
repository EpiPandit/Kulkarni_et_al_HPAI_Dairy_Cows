library(tidyverse)

beta_dm <-read.csv("./summary_daily_100sims_optiters_dm.csv")
r_dm <- read.csv("./summary_farm_100sims_optiters_dm.csv")

per_farm_seir_dm <- beta_dm %>% group_by(farm, day) %>%
  summarise(mean_S = mean(S), sd_S = sd(S),
            mean_E = mean(E), sd_E = sd(E),
            mean_I = mean(I), sd_I = sd(I),
            mean_R = mean(R_count), sd_R = sd(R_count)) %>%
  mutate_if(is.numeric, floor)

beta_plf <-read.csv("./summary_re_100sims_optiters_plf.csv")
r_plf <- read.csv("./summary_farm_100sims_optiters_plf.csv")

per_farm_seir_plf <- beta_plf %>% group_by(farm, day) %>%
  summarise(mean_S = mean(S), sd_S = sd(S),
            mean_E = mean(E), sd_E = sd(E),
            mean_I = mean(I), sd_I = sd(I),
            mean_R = mean(R_count), sd_R = sd(R_count)) %>%
  mutate_if(is.numeric, floor)

outbreak_length <- c("Farm1" = 33,
                     "Farm2" = 53,
                     "Farm3" = 23,
                     "Farm4" = 18)

# est_R <- function(df, early_days) {
#   inc <- df %>% 
#     group_by(day) %>%
#     summarise(incidence = sum(new_E, na.rm = TRUE), .groups = "drop") %>%
#     arrange(day) %>%
#     filter(day <= min(day) + early_days - 1)
#   
#   if (nrow(inc) < 3 || sum(inc$incidence, na.rm = TRUE) == 0) return(NA_real_)
#   
#   fit <- lm(log(pmax(incidence, 1e-8)) ~ day, data = inc)
#   unname(coef(fit)[2])
# }
# r_to_R0_seir <- function(r, sigma, gamma) {
#   beta <- ((r + sigma) * (r + gamma)) / sigma
#   beta / gamma
# }
# sigma <- 1/latency_period 
# 
# r0_estimates <- beta_dm %>%
#   group_by(farm, sim_id) %>%
#   group_split() %>%
#   map_dfr(function(df, outbreak) {
#     r <- estimate_r_from_early_phase(df, early_days = outbreak)
#     
#     gamma_vals <- df$recovery_rate
#     gamma_vals <- gamma_vals[is.finite(gamma_vals)]
#     gamma <- if (length(gamma_vals) == 0) NA_real_ else median(gamma_vals)
#     
#     tibble(
#       farm = unique(df$farm),
#       sim_id = unique(df$sim_id),
#       r = r,
#       gamma = gamma,
#       R0 = if (is.na(r) || is.na(gamma) || is.na(sigma)) NA_real_ else r_to_R0_seir(r, sigma, gamma)
#     )
#   })
# 
# r0_estimates

latent_period <- c("Farm1" = 3,"Farm2" = 4,"Farm3"= 2,"Farm4" = 2)

estimate_r0_final_size <- function(df) {
  # Get initial and final susceptible counts
  S0 <- df$S[1]  # Initial susceptibles
  S_final <- df$S[nrow(df)]  # Final susceptibles
  
  # Check if outbreak occurred
  total_infected <- S0 - S_final
  if (total_infected <= 2) {
    return(tibble(S0 = S0, S_final = S_final, attack_rate = 0, R0 = 0))
  }
  
  # Calculate attack rate
  attack_rate <- total_infected / S0
  
  # Solve for R0 using final size equation
  # R0 = -ln(S_final/S0) / (1 - S_final/S0)
  s_ratio <- S_final / S0
  
  if (s_ratio >= 1 || s_ratio <= 0) {
    return(tibble(S0 = S0, S_final = S_final, attack_rate = attack_rate, R0 = NA_real_))
  }
  
  R0 <- -log(s_ratio) / (1 - s_ratio)
  
  tibble(S0 = S0, S_final = S_final, attack_rate = attack_rate, R0 = R0)
}

r0_final_size <- map_dfr(farm_ids, function(fm) {
  df_farm <- beta_dm %>% filter(farm == fm)
  
  df_farm %>%
    group_by(sim_id) %>%
    group_split() %>%
    map_dfr(function(df_sim) {
      estimate_r0_final_size(df_sim) %>%
        mutate(farm = fm, sim_id = unique(df_sim$sim_id))
    })
})

r0_summary <- r0_final_size %>%
  group_by(farm) %>%
  summarise(mean_R0 = mean(R0, na.rm = T),
            median_R0 = median(R0, na.rm = T),
            sd_r0 = sd(R0, na.rm = T),
            q25_R0 = quantile(R0,0.25),
            q75_R0 = quantile(R0, 0.75))

r0_final_size_plf <- map_dfr(farm_ids, function(fm) {
  df_farm <- beta_plf %>% filter(farm == fm)
  
  df_farm %>%
    group_by(sim_id) %>%
    group_split() %>%
    map_dfr(function(df_sim) {
      estimate_r0_final_size(df_sim) %>%
        mutate(farm = fm, sim_id = unique(df_sim$sim_id))
    })
})

r0_summary_plf <- r0_final_size %>%
  group_by(farm) %>%
  summarise(mean_R0 = mean(R0, na.rm = T),
            median_R0 = median(R0, na.rm = T),
            sd_r0 = sd(R0, na.rm = T),
            q25_R0 = quantile(R0,0.25),
            q75_R0 = quantile(R0, 0.75))

r0_summary_merged <- do.call("rbind", list(r0_summary, r0_summary_plf))
r0_summary_merged$method <- rep(c("dm", "plf"), each = 4)
write.csv(r0_summary_merged, "./summary_R0_all.csv", row.names = F)

r0_summary_merged <- read.csv("./summary_R0_all.csv")
r0_sum <- r0_summary_merged %>% filter(days == "1-15 days")
r0_sum <- r0_sum %>% pivot_longer(c(mean_R0_dm,mean_R0_plf), names_to = "metric", values_to = "mean")
r0_sum <- r0_sum %>% pivot_longer(c(sd_R0_dm,sd_R0_plf), names_to = "metric_sd", values_to = "sd")
r0_sum$metric <- gsub("mean_R0_", "", r0_sum$metric)
r0_sum$metric_sd <- gsub("sd_R0_", "", r0_sum$metric_sd)
r0_sum <- r0_sum %>% filter(metric == metric_sd) %>% select(-c(days, metric_sd))
ggplot(r0_sum, aes(y = farm, x = mean, color = metric)) +
  geom_errorbarh(aes(xmin = mean - 1.96*sd, xmax = mean + 1.96*sd), 
                 position = position_dodge(width = 0.9),
                 width = 0.2) +
  geom_label(aes(label = round(mean, digits = 2)), vjust = 0.5, hjust = 0.5, size = 4.5, 
             show.legend = F, position = position_dodge(width = 0.9),fontface = "bold", fill = "white") +
  scale_fill_brewer(type = "div")+
  theme_classic() +
  labs(y = "Farm", x = expression(paste("Effective ", R[0])), color = "Loss Func.") +
  theme_classic() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = "top",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major.x = element_line(size = 0.4, linetype = 2, color = "grey50"))
ggsave("R0_summary_all.png", width = 6, height = 6, dpi = 300)
