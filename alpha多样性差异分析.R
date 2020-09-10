library(openxlsx)
library(ggplot2)
library(vegan)
library(devtools)

suppressWarnings(suppressMessages(library(amplicon)))
if (!requireNamespace("amplicon", quietly = TRUE))install_github("microbiota/amplicon")


setwd("C:\\Users\\xhx\\Desktop\\新增图片")
alpha_index <- read.xlsx("alpha.xlsx",rowNames = T)



#alpha_index$Group <- group$Group
alpha_index <- data.frame(alpha_index)
#alpha_index$Group <- as.factor(alpha_index$Group)

# altitude
group_altitude <- read.xlsx("altitude.xlsx", rowNames = T)
group_altitude$Group <- as.factor(group_altitude$Group)
group_altitude$sample <- rownames(group_altitude)
p1 <- alpha_boxplot(alpha_index, group_altitude, "Shannon", groupID="Group")
ggsave(paste0("alpha_altitude_shannon_boxplot.png"), p1, width=89, height=56, units="mm")

p2 <- alpha_boxplot(alpha_index, group_altitude, "Simpson", groupID="Group")
ggsave(paste0("alpha_altitude_simpson_boxplot.png"), p2, width=89, height=56, units="mm")

p3 <- alpha_boxplot(alpha_index, group_altitude, "Chao1", groupID="Group")
ggsave(paste0("alpha_altitude_chao1_boxplot.png"), p3, width=89, height=56, units="mm")


# population area
group_area <- read.xlsx("area.xlsx", rowNames = T)
group_area$Group <- as.factor(group_area$Group)
group_area$sample <- rownames(group_area)
p1 <- alpha_boxplot(alpha_index, group_area, "Shannon", groupID="Group")
ggsave(paste0("alpha_area_shannon_boxplot.png"), p1, width=89, height=56, units="mm")

p2 <- alpha_boxplot(alpha_index, group_area, "Simpson", groupID="Group")
ggsave(paste0("alpha_area_simpson_boxplot.png"), p2, width=89, height=56, units="mm")

p3 <- alpha_boxplot(alpha_index, group_area, "Chao1", groupID="Group")
ggsave(paste0("alpha_area_chao1_boxplot.png"), p3, width=89, height=56, units="mm")





unweighted <- read.xlsx("差异性矩阵unweighted.unifrac.distance.matrix.xlsx",rowNames=T)
unweighted.matrix <- as.matrix(unweighted)

p4 <- beta_pcoa(unweighted.matrix, group, groupID="Group", ellipse=T, label=T, PCo=13)
ggsave(paste0("unweighted.unifrac.PCoA.png"), p4, width=89, height=56, units="mm")