# Manhatten plot
setwd("/Users/halowin/Desktop/UoN/thesis")
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
leaf_As75_Global <- read.csv("filtered1/Global/leaf_ionome_As75_GWAS.csv")
# Align the five chromosomes
leaf_As75_Global_mutate <- 
  leaf_As75_Global %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(leaf_As75_Global,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

# Alpha level adjusted with Bonferroni
sig <- 0.05 / nrow(leaf_As75_Global)

# Obtain the midpoint of every chromosome for drawing scales on x axis
axisdf <- leaf_As75_Global_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)

# Plot
p1 <- ggplot(leaf_As75_Global_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome", label = axisdf$chromosomes, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "As leaf") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") 
        
seed_As75_Global <- read.csv("filtered1/Global/seed_ionome_As75_GWAS.csv")
seed_As75_Global_mutate <- 
  seed_As75_Global %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_As75_Global,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

sig <- 0.05 / nrow(seed_As75_Global)
axisdf <- seed_As75_Global_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)

p2 <- ggplot(seed_As75_Global_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome", label = axisdf$chromosomes, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "As seed") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") 

leaf_P31_Global <- read.csv("filtered1/Global/leaf_ionome_P31_GWAS.csv")
leaf_P31_Global_mutate <- 
  leaf_P31_Global %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(leaf_P31_Global,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

axisdf2 <- leaf_P31_Global_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig <- 0.05 / nrow(leaf_P31_Global)

p3 <- ggplot(leaf_P31_Global_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome", label = axisdf2$chromosomes, breaks = axisdf2$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "P leaf") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") 

seed_P31_Global <- read.csv("filtered1/Global/seed_ionome_P31_GWAS.csv")
seed_P31_Global_mutate <- 
  seed_P31_Global %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_P31_Global,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

axisdf2 <- seed_P31_Global_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig <- 0.05 / nrow(seed_P31_Global)

p4 <- ggplot(seed_P31_Global_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome", label = axisdf2$chromosomes, breaks = axisdf2$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "P seed") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") 

leaf_Zn66_Global <- read.csv("filtered1/Global/leaf_ionome_Zn66_GWAS.csv")
leaf_Zn66_Global_mutate <- 
  leaf_Zn66_Global %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(leaf_Zn66_Global,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

axisdf3 <- leaf_Zn66_Global_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig <- 0.05 / nrow(leaf_Zn66_Global)
p5 <- ggplot(leaf_Zn66_Global_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome",label = axisdf3$chromosomes, breaks = axisdf3$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "Zn leaf") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") 

seed_Zn66_Global <- read.csv("filtered1/Global/seed_ionome_Zn66_GWAS.csv")
seed_Zn66_Global_mutate <- 
  seed_Zn66_Global %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_Zn66_Global,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

axisdf3 <- seed_Zn66_Global_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig <- 0.05 / nrow(seed_Zn66_Global)
p6 <- ggplot(seed_Zn66_Global_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome",label = axisdf3$chromosomes, breaks = axisdf3$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "Zn seed") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") 

leaf_Mn55_Global <- read.csv("filtered2/Global/leaf_ionome_Mn55_GWAS.csv")
leaf_Mn55_Global_mutate <- 
  leaf_Mn55_Global %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(leaf_Mn55_Global,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

axisdf5 <- leaf_Mn55_Global_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig5 <- 0.05 / nrow(leaf_Mn55_Global)
p7 <- ggplot(leaf_Mn55_Global_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(label = axisdf5$chromosomes, breaks = axisdf5$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "Mn leaf") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(sig5), color = "red", linetype = "dashed") 

seed_Mn55_Global <- read.csv("filtered2/Global/seed_ionome_Mn55_GWAS.csv")
seed_Mn55_Global_mutate <- 
  seed_Mn55_Global %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_Mn55_Global,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

axisdf5 <- seed_Mn55_Global_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig5 <- 0.05 / nrow(seed_Mn55_Global)
p8 <- ggplot(seed_Mn55_Global_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(label = axisdf5$chromosomes, breaks = axisdf5$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "Mn seed") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(sig5), color = "red", linetype = "dashed") 

leaf_K39_Global <- read.csv("filtered2/Global/leaf_ionome_K39_GWAS.csv")
leaf_K39_Global_mutate <- 
  leaf_K39_Global %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(leaf_K39_Global,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

axisdf5 <- leaf_K39_Global_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig5 <- 0.05 / nrow(leaf_K39_Global)
p9 <- ggplot(leaf_K39_Global_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(label = axisdf5$chromosomes, breaks = axisdf5$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "K leaf") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(sig5), color = "red", linetype = "dashed") 

seed_K39_Global <- read.csv("filtered2/Global/seed_ionome_K39_GWAS.csv")
seed_K39_Global_mutate <- 
  seed_K39_Global %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_K39_Global,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

axisdf5 <- seed_K39_Global_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig5 <- 0.05 / nrow(seed_K39_Global)
p10 <- ggplot(seed_K39_Global_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(label = axisdf5$chromosomes, breaks = axisdf5$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "K seed") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(sig5), color = "red", linetype = "dashed") 

leaf_Ni60_Global <- read.csv("filtered2/Global/leaf_ionome_Ni60_GWAS.csv")
leaf_Ni60_Global_mutate <- 
  leaf_Ni60_Global %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(leaf_Ni60_Global,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

axisdf5 <- leaf_Ni60_Global_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig5 <- 0.05 / nrow(leaf_Ni60_Global)
p11 <- ggplot(leaf_Ni60_Global_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(label = axisdf5$chromosomes, breaks = axisdf5$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "Ni leaf") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(sig5), color = "red", linetype = "dashed") 

seed_Ni60_Global <- read.csv("filtered2/Global/seed_ionome_Ni60_GWAS.csv")
seed_Ni60_Global_mutate <- 
  seed_Ni60_Global %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_Ni60_Global,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

axisdf5 <- seed_Ni60_Global_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig5 <- 0.05 / nrow(seed_Ni60_Global)
p12 <- ggplot(seed_Ni60_Global_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(label = axisdf5$chromosomes, breaks = axisdf5$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "Ni seed") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(sig5), color = "red", linetype = "dashed") 

# 
global_p1 <- grid.arrange(p3,p4,p5,p6,p1,p2,ncol=2)
global_p2 <- grid.arrange(p9,p10,p7,p8,p11,p12,ncol=2)
