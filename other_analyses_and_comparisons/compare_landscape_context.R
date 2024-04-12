# these are the imp surface cover and canopy cover in all census blocks
# within 500m radius of the site. The data include both 2014 and 2020 but
# we will use 2020.

library(tidyverse)
library(gridExtra)

df <- read.csv("./data/land_cover_by_site.csv")

df1 <- df %>%
  group_by(site) %>%
  mutate(mean_imp = mean(Imp2020_pe),
         mean_canopy = mean(CC2020_per)) %>%
  slice(1) %>%
  select(site, latitude, longitude, category, 
         mean_herb_scaled, mean_woody_scaled, mean_herb, mean_woody,
         mean_imp, mean_canopy) %>%
  ungroup()

control <- df1 %>% filter(category == "comparison") %>% select(mean_imp, mean_canopy)
enhanced <- df1 %>% filter(category == "meadow") %>% select(mean_imp, mean_canopy)

t.test(control$mean_imp, enhanced$mean_imp)
t.test(control$mean_canopy, enhanced$mean_canopy)

df2 <- df %>%
  mutate(area_imp = Imp2020_pe * 0.01 * area,
         area_canopy = CC2020_per * 0.01 * area) %>% # multiply by 0.01 to get into a decimal
  group_by(site) %>%
  mutate(total_area = sum(area),
         total_imp = sum(area_imp),
         imp_standardized = total_imp/total_area,
         total_CC = sum(area_canopy),
         canopy_standardized = total_CC/total_area) %>%
  slice(1) %>%
  select(site, latitude, longitude, category, 
         mean_herb_scaled, mean_woody_scaled, mean_herb, mean_woody,
         imp_standardized, canopy_standardized) %>%
  ungroup()
  
control <- df2 %>% filter(category == "comparison") %>% select(imp_standardized, canopy_standardized)
enhanced <- df2 %>% filter(category == "meadow") %>% select(imp_standardized, canopy_standardized)

ttest1 <- t.test(control$imp_standardized, enhanced$imp_standardized)
ttest2 <- t.test(control$canopy_standardized, enhanced$canopy_standardized)

mean_control1 <- signif(ttest1$estimate[1], digits = 4) * 100
mean_enhanced1 <- signif(ttest1$estimate[2], digits = 4) * 100
t1 <- signif(ttest1$statistic, digits = 4)
pvalue1 <- signif(ttest1$p.value, digits = 4)
sd_control1 <- sd(control$imp_standardized)
sd_enhanced1 <- sd(enhanced$imp_standardized) 
 
(p <- ggplot(data = df2, aes(x=as.factor(category), y=imp_standardized)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
    annotate("text", x = 1, y = 0.25, size=6, label = paste0("mean = ", mean_control1, "%")) +
    annotate("text", x = 2, y = 0.25, size=6, label = paste0("mean = ", mean_enhanced1, "%")) +
    annotate("text", x = 1.5, y = 0.15, size=6, label = paste0("t = ", t1)) +
    annotate("text", x = 1.5, y = 0.07, size=6, label = paste0("p-value = ", pvalue1)) +
    scale_x_discrete(name="", breaks = c("comparison", "meadow"),
                     labels=c("control", "enhanced"
                     )) +
    ylab("(%) impervious surface cover \nin 500m radius") +
    scale_y_continuous(limits = c(0,1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = scales::percent 
    ) +
    theme_bw() +
    theme(legend.text=element_text(size=16),
          legend.position = c(0.75, 0.125),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
)  

mean_control2 <- signif(ttest2$estimate[1], digits = 4) * 100
mean_enhanced2 <- signif(ttest2$estimate[2], digits = 4) * 100
t2 <- signif(ttest2$statistic, digits = 4)
pvalue2 <- signif(ttest2$p.value, digits = 4)
sd_control2 <- sd(control$canopy_standardized)
sd_enhanced2 <- sd(enhanced$canopy_standardized) 

(q <- ggplot(data = df2, aes(x=as.factor(category), y=canopy_standardized)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
    annotate("text", x = 1, y = 0.75, size=6, label = paste0("mean = ", mean_control2, "%")) +
    annotate("text", x = 2, y = 0.75, size=6, label = paste0("mean = ", mean_enhanced2, "%")) +
    annotate("text", x = 1.5, y = 0.97, size=6, label = paste0("t = ", t2)) +
    annotate("text", x = 1.5, y = 0.89, size=6, label = paste0("p-value = ", pvalue2)) +
    scale_x_discrete(name="", breaks = c("comparison", "meadow"),
                     labels=c("control", "enhanced"
                     )) +
    ylab("(%) tree canopy cover \nin 500m radius") +
    scale_y_continuous(limits = c(0,1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = scales::percent 
    ) +
    theme_bw() +
    theme(legend.text=element_text(size=16),
          legend.position = c(0.75, 0.125),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
)  

grid.arrange(p, q, ncol=2)

write.csv(df2, "./data/land_cover_by_site_reduced.csv")
