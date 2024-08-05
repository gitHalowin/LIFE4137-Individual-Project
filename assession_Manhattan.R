library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
leaf_Zn66_7411 <- read.csv("filtered1/Anchored/leaf_Zn66/7411_leaf_ionome_Zn66_GWAS.csv")
leaf_Zn66_7411_mutate <- 
  leaf_Zn66_7411 %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(leaf_Zn66_7411,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

axisdf4 <- leaf_Zn66_7411_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig <- 0.05 / nrow(leaf_Zn66_7411)
p7 <- ggplot(leaf_Zn66_7411_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome", label = axisdf4$chromosomes, breaks = axisdf4$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "Zn leaf ionome anchored on accession 7411") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed")

leaf_As75_1158 <- read.csv("filtered1/Anchored/leaf_As75/1158_leaf_ionome_As75_GWAS.csv")
leaf_As75_1158_mutate <- 
  leaf_As75_1158 %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(leaf_As75_1158,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)
axisdf7 <- leaf_As75_1158_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig <- 0.05 / nrow(leaf_As75_1158)

p8 <- ggplot(leaf_As75_1158_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome", label = axisdf7$chromosomes, breaks = axisdf7$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "As leaf ionome anchored on accession 1158")+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") 

seed_Ni60_9558 <- read.csv("filtered2/Anchored/seed_Ni60/9558_seed_ionome_Ni60_GWAS.csv")
seed_Ni60_9558_mutate <- 
  seed_Ni60_9558 %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_Ni60_9558,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)
axisdf10 <- seed_Ni60_9558_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig10 <- 0.05 / nrow(seed_Ni60_9558)
p9 <- ggplot(seed_Ni60_9558_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome",label = axisdf10$chromosomes, breaks = axisdf10$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "Ni seed ionome anchored on accession 9558") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = -log10(sig10), color = "red", linetype = "dashed")
# combine  p7, p8, p9, arrange them vertically
fig3 <- grid.arrange(p7,p8,p9,ncol=1)

leaf_As75_7206 <- read.csv("filtered1/Anchored/leaf_As75/7206_leaf_ionome_As75_GWAS.csv")
leaf_As75_7206_mutate <- 
  leaf_As75_7206 %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(leaf_As75_7206,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)
axisdf7 <- leaf_As75_7206_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig <- 0.05 / nrow(leaf_As75_7206)

p10 <- ggplot(leaf_As75_7206_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome", label = axisdf7$chromosomes, breaks = axisdf7$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "As leaf ionome anchored on accession 7206")+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") 

seed_As75_1137 <- read.csv("filtered1/Anchored/seed_As75/1137_seed_ionome_As75_GWAS.csv")
seed_As75_1137_mutate <- 
  seed_As75_1137 %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_As75_1137,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)

axisdf <- seed_As75_1137_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig <- 0.05 / nrow(seed_As75_1137)

p1 <- ggplot(seed_As75_1137_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome",label = axisdf$chromosomes, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "As seed ionome anchored on accession 1137") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") 
 
fig4 <- grid.arrange(p10,p1,p11,ncol=1)
seed_As75_9391 <- read.csv("filtered1/Anchored/seed_As75/9391_seed_ionome_As75_GWAS.csv")
seed_As75_9391_mutate <- 
  seed_As75_9391 %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_As75_9391 ,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot) 

axisdf <- seed_As75_9391_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig <- 0.05 / nrow(seed_As75_9391)

p2 <- ggplot(seed_As75_9391_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome",label = axisdf$chromosomes, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "As seed ionome anchored on accession 9391") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") 
# combine p1 and p2, arrange them vertically
p2A <- grid.arrange(p1,p2,ncol=1)

leaf_K39_7010 <- read.csv("filtered2/Anchored/leaf_K39/7010_leaf_ionome_K39_GWAS.csv")
leaf_K39_7010_mutate <- 
  leaf_K39_7010 %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(leaf_K39_7010,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)
axisdf10 <- leaf_K39_7010_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig10 <- 0.05 / nrow(leaf_K39_7010)
p3 <- ggplot(leaf_K39_7010_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(label = axisdf10$chromosomes, breaks = axisdf10$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "K leaf ionome anchored on accession 7010") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = -log10(sig10), color = "red", linetype = "dashed")

seed_Mn55_6074 <- read.csv("filtered2/Anchored/seed_Mn55/6074_seed_ionome_Mn55_GWAS.csv")
seed_Mn55_6074_mutate <- 
  seed_Mn55_6074 %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_Mn55_6074,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)
axisdf10 <- seed_Mn55_6074_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig10 <- 0.05 / nrow(seed_Mn55_6074)
p4 <- ggplot(seed_Mn55_6074_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(label = axisdf10$chromosomes, breaks = axisdf10$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "Mn seed ionome anchored on accession 6074") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = -log10(sig10), color = "red", linetype = "dashed")

# combine p3 and p4, arrange them vertically
p2B <- grid.arrange(p3,p4,ncol=1)

seed_Mn55_7002 <- read.csv("filtered2/Anchored/seed_Mn55/7002_seed_ionome_Mn55_GWAS.csv")
seed_Mn55_7002_mutate <- 
  seed_Mn55_7002 %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_Mn55_7002,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)
axisdf10 <- seed_Mn55_7002_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig10 <- 0.05 / nrow(seed_Mn55_7002)
p4 <- ggplot(seed_Mn55_7002_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(label = axisdf10$chromosomes, breaks = axisdf10$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "Mn seed ionome anchored on accession 7002") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = -log10(sig10), color = "red", linetype = "dashed")

seed_Mn55_7344 <- read.csv("filtered2/Anchored/seed_Mn55/7344_seed_ionome_Mn55_GWAS.csv")
seed_Mn55_7344_mutate <- 
  seed_Mn55_7344 %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_Mn55_7344,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)
axisdf10 <- seed_Mn55_7344_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig10 <- 0.05 / nrow(seed_Mn55_7344)
p5 <- ggplot(seed_Mn55_7344_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(label = axisdf10$chromosomes, breaks = axisdf10$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "Mn seed ionome anchored on accession 7344") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = -log10(sig10), color = "red", linetype = "dashed")

leaf_K39_7133 <- read.csv("filtered2/Anchored/leaf_K39/7133_leaf_ionome_K39_GWAS.csv")
leaf_K39_7133_mutate <- 
  leaf_K39_7133 %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(leaf_K39_7133,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)
axisdf10 <- leaf_K39_7133_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig10 <- 0.05 / nrow(leaf_K39_7133)
p6 <- ggplot(leaf_K39_7133_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(label = axisdf10$chromosomes, breaks = axisdf10$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "K leaf ionome anchored on accession 7133") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = -log10(sig10), color = "red", linetype = "dashed")
# combine p4, p5 , p6, arrange them vertically
p2C <- grid.arrange(p4,p5,p6,ncol=1)

seed_K39_7333 <- read.csv("filtered2/Anchored/seed_K39/7333_seed_ionome_K39_GWAS.csv")
seed_K39_7333_mutate <- 
  seed_K39_7333 %>%
  # Compute chromosome size
  group_by(chromosomes) %>%
  summarise(chr_len=max(positions)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add to the initial dataset
  left_join(seed_K39_7333,.,by=c("chromosomes"="chromosomes")) %>%
  arrange(chromosomes,positions) %>%
  mutate(BPcum = positions + tot)
axisdf10 <- seed_K39_7333_mutate %>% group_by(chromosomes) %>% summarize(center = ( max(BPcum) + min(BPcum) )/2)
sig10 <- 0.05 / nrow(seed_K39_7333)
p11 <- ggplot(seed_K39_7333_mutate,aes(x=BPcum,y=-log10(pvals))) +
  geom_point(aes(color=as.factor(chromosomes)), alpha=0.5, size=2) +
  scale_color_manual(values = rep(c(blues9[7],blues9[9]),5)) +
  scale_x_continuous(name = "Chromosome",label = axisdf10$chromosomes, breaks = axisdf10$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(title = "K seed ionome anchored on accession 7333") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = -log10(sig10), color = "red", linetype = "dashed")
