library(openxlsx)
library(tidyverse)

setwd("C:\\Users\\xhx\\Desktop\\新增图片")

## 导入数据
data <- read.xlsx("genus.xlsx", rowNames=T)
head(data)
y <- function(x) {x/sum(x)}
data <- apply(data, 2, y)
data <- data.frame(data*100)
data <- data.frame(data[-1,])
data$sum <- rowSums(data)

group_altitude <- read.xlsx("altitude.xlsx")
group_area <- read.xlsx("area.xlsx")


STAMP <- function(group_info, picture_name, whole=FALSE){
  
  if (whole==F){
  top20 <- data[order(data$sum, decreasing = T),][1:20,-10]
  top20.t <- t(top20)
  data1 <- data.frame(top20.t,group_info$Group)
  colnames(data1) <- c(colnames(top20.t),"Group")
  data1$Group <- as.factor(data1$Group)
  }
 
  else{
    data.t <- t(data[,-10])
    data1 <- data.frame(data.t,group_info$Group)
    colnames(data1) <- c(colnames(data.t),"Group")
    data1$Group <- as.factor(data1$Group)
    head(data1)
  }
 
  
  
  ## t-test
  
  diff <- data1 %>% 
    select_if(is.numeric) %>%
    map_df(~ broom::tidy(t.test(. ~ Group,data = data1)), .id = 'var')
  
  diff$p.value <- p.adjust(diff$p.value,"none")
  diff$p.value
  diff <- diff %>% filter(p.value < 0.05)
  
  ## 绘图数据构建
  ## 左侧条形图
  abun.bar <- data1[,c(diff$var,"Group")] %>% 
    gather(variable,value,-Group) %>% 
    group_by(variable,Group) %>% 
    summarise(Mean = mean(value))
  
  ## 右侧散点图
  diff.mean <- diff[,c("var","estimate","conf.low","conf.high","p.value")]
  diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data1$Group)[1],
                              levels(data1$Group)[2]))
  diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]
  
  ## 左侧条形图
  library(ggplot2)
  cbbPalette <- c("#E69F00", "#56B4E9")
  abun.bar$variable <- factor(abun.bar$variable,levels = rev(diff.mean$var))
  p1 <- ggplot(abun.bar,aes(variable,Mean,fill = Group)) +
    scale_x_discrete(limits = levels(diff.mean$var)) +
    coord_flip() +
    xlab("") +
    ylab("Mean proportion (%)") +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          legend.title=element_blank(),
          legend.text=element_text(size=12,face = "bold",colour = "black",
                                   margin = margin(r = 20)),
          legend.position = c(-0.35,-0.15),
          legend.direction = "horizontal",
          legend.key.width = unit(0.5,"cm"),
          legend.key.height = unit(0.5,"cm"))
  
  
  for (i in 1:(nrow(diff.mean) - 1)) 
    p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, 
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
  
  p1 <- p1 + 
    geom_bar(stat = "identity",position = "dodge",width = 0.7,colour = "black") +
    scale_fill_manual(values=cbbPalette)
  
  
  ## 右侧散点图
  diff.mean$var <- factor(diff.mean$var,levels = levels(abun.bar$variable))
  diff.mean$p.value <- signif(diff.mean$p.value,3)
  diff.mean$p.value <- as.character(diff.mean$p.value)
  p2 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          axis.text.y = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +
    scale_x_discrete(limits = levels(diff.mean$var)) +
    coord_flip() +
    xlab("") +
    ylab("Difference in mean proportions (%)") +
    labs(title="95% confidence intervals") 
  
  for (i in 1:(nrow(diff.mean) - 1)) 
    p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, 
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
  
  p2 <- p2 +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                  position = position_dodge(0.8), width = 0.5, size = 0.5) +
    geom_point(shape = 21,size = 3) +
    scale_fill_manual(values=cbbPalette) +
    geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')
  
  
  p3 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
    geom_text(aes(y = 0,x = var),label = diff.mean$p.value,
              hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
    geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value",
              srt = 90,fontface = "bold",size = 5) +
    coord_flip() +
    ylim(c(0,1)) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
  
  ## 图像拼接
  library(patchwork)
  p <- p1 + p2 + p3 + plot_layout(widths = c(4,6,2))
  
  ggsave(paste0(picture_name,".pdf"),p,width = 10,height = 4)
  ggsave(paste0(picture_name,".png"),p,width = 10,height = 4)
}


STAMP(group_altitude, "stamp_altitude")
STAMP(group_area, "stamp_area", whole=T)




