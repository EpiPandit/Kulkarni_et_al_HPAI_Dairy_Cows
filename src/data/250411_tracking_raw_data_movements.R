library(tidyverse)

farm1a <- readxl::read_excel("pone.0203190.s007.xlsx", sheet = 1)
farm1b <- readxl::read_excel("pone.0203190.s007.xlsx", sheet = 2)
farm2 <- readxl::read_excel("pone.0203190.s007.xlsx", sheet = 3)
farm3 <- readxl::read_excel("pone.0203190.s007.xlsx", sheet = 4)
farm4 <- readxl::read_excel("pone.0203190.s007.xlsx", sheet = 5)

par(mfrow = c(2,3))
farm1a$modelstate <- factor(farm1a$modelstate, levels = paste0("i", 1:14), ordered = T)
plot(farm1a$modelstate, main = "farm1a")
farm1b$modelstate <- factor(farm1b$modelstate, levels = paste0("i", 1:14), ordered = T)
plot(farm1b$modelstate, main = "farm1b")
farm2$modelstate <- factor(farm2$modelstate, levels = paste0("i", 1:14), ordered = T)
plot(farm2$modelstate, main = "farm2")
farm3$modelstate <- factor(farm3$modelstate, levels = paste0("i", 1:14), ordered = T)
plot(farm3$modelstate, main = "farm3")
farm4$modelstate <- factor(farm4$modelstate, levels = paste0("i", 1:14), ordered = T)
plot(farm4$modelstate, main = "farm4")


summ_1a <- farm1a %>% group_by(year, month, modelstate) %>% summarise(counts = n())
table(summ_1a)
summ_1a$time <- paste(substr(summ_1a$month,1, 3), summ_1a$year, sep = "-")
plot1a <- ggplot(summ_1a, aes(x = time, y = counts, fill = modelstate)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6))+
  scale_fill_viridis_d(option = "cividis") +
  theme_minimal() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) + 
  labs(x = "month-year", y = "cows in pen", fill = "Pen code", subtitle = "A. farm 1a") + 
  coord_flip()


summ_3 <- farm3 %>% group_by(year, month, modelstate) %>% summarise(counts = n())
table(summ_3)
summ_3$time <- paste(substr(summ_3$month,1, 3), summ_3$year, sep = "-")
plot3 <- ggplot(summ_3, aes(x = time, y = counts, fill = modelstate)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6))+
  scale_fill_viridis_d(option = "cividis") +
  theme_minimal() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) + 
  labs(x = "month-year", y = "cows in pen", fill = "Pen code", subtitle = "D. farm 3") +
  coord_flip()


summ_2 <- farm2 %>% group_by(year, month, modelstate) %>% summarise(counts = n())
table(summ_2)
summ_2$time <- paste(substr(summ_2$month,1, 3), summ_2$year, sep = "-")
plot2 <- ggplot(summ_2, aes(x = time, y = counts, fill = modelstate)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6))+
  scale_fill_viridis_d(option = "cividis") +
  theme_minimal() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) + 
  labs(x = "month-year", y = "cows in pen", fill = "Pen code", subtitle = "C. farm 2")+
  coord_flip()


summ_1b <- farm1b %>% group_by(year, month, modelstate) %>% summarise(counts = n())
table(summ_1b)
summ_1b$time <- paste(substr(summ_1b$month,1, 3), summ_1b$year, sep = "-")
plot1b <- ggplot(summ_1b, aes(x = time, y = counts, fill = modelstate)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6))+
  scale_fill_viridis_d(option = "cividis") +
  theme_minimal() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) + 
  labs(x = "month-year", y = "cows in pen", fill = "Pen code", subtitle = "B. farm 1b")+
  coord_flip()

summ_4 <- farm4 %>% group_by(year, month, modelstate) %>% summarise(counts = n())
table(summ_4)
summ_4$time <- paste(substr(summ_4$month,1, 3), summ_4$year, sep = "-")
plot4 <- ggplot(summ_4, aes(x = time, y = counts, fill = modelstate)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6))+
  scale_fill_viridis_d(option = "cividis") +
  theme_minimal() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) + 
  labs(x = "month-year", y = "cows in pen", fill = "Pen code", subtitle = "E. farm 4")+
  coord_flip()



ggsave(filename = "farm1a_movements.png", plot1a, dpi = 300, height = 8, width = 6)
ggsave(filename = "farm1b_movements.png", plot1b, dpi = 300, height = 8, width = 6)
ggsave(filename = "farm2_movements.png", plot2, dpi = 300, height = 10, width = 6)
ggsave(filename = "farm3_movements.png", plot3, dpi = 300, height = 8, width = 6)
ggsave(filename = "farm4_movements.png", plot4, dpi = 300, height = 8, width = 6)

write.csv(summ_1a, "farm1a_movements.csv", row.names = F)
write.csv(summ_1b, "farm1b_movements.csv", row.names = F)
write.csv(summ_2, "farm2_movements.csv", row.names = F)
write.csv(summ_3, "farm3_movements.csv", row.names = F)
write.csv(summ_4, "farm4_movements.csv", row.names = F)

