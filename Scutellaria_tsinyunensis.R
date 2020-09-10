library(vegan)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(plyr)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(GUniFrac)
library(ape)
library(BiocManager)
# BiocManager::install("phyloseq")
library(phyloseq)
library(ggalt)
library(ggvegan)
library(packfor)
library(dendextend)
library(MASS)
# 2.0各样本OTU分布及物种信息表格
# *******************************
work_path <- "C:\\Users\\xhx\\Desktop\\生信分析\\2.0各样本OTU分布及物种信息表格"
setwd(work_path)

# 原始数据
raw_data <- read.xlsx("1.各样本OTU分布及物种信息表格.xlsx",rowNames=T)

# 对数据进行均一化，转化为相对丰度,并保存相对丰度表格
# 均一化函数
y <- function(x){x/sum(x)}

relative_abundance <- apply(raw_data[,1:10],2,y)
relative_abundance.file <- cbind(relative_abundance,raw_data[,-(1:10)])
rownames(relative_abundance.file) <- row.names(raw_data)
write.xlsx(relative_abundance.file,
           file="各物种在样地内的相对丰度表格.xlsx",rowNames=T)
write.xlsx(relative_abundance,file="去除物种信息的相对丰度表格.xlsx",rowNames=T)


# 各分类级别物种柱状图图
data <- read.xlsx("物种统计.xlsx", rowNames=T)
data$sample <- rownames(data)
data.melt <- melt(data, id.vars=c("sample"))
colnames(data.melt)[2:3] <- c("taxonomy", "count")
data.melt$taxonomy <- as.factor(data.melt$taxonomy)
ggplot(data.melt, aes(x=taxonomy, y=count, group=sample, col=sample))+
  geom_line(size=1)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  ggsave("物种统计.pdf", width =8.0, height = 6.0)+
  ggsave("物种统计.png", width =8.0, height = 6.0)




# 2.1稀释性曲线图Rarefaction-curve
# *******************
# 将工作目录转至稀释性曲线图文件夹
setwd("C:\\Users\\xhx\\Desktop\\生信分析\\2.1稀释性曲线图Rarefaction-curve")

# 获取稀释性曲线图绘图数据
data <- raw_data[,1:9]
raremax <- colSums(data)
rarefaction.data <- data.frame()

# 这里的add变量是由于稀释时间隔100造成raremax时其它组别的表格数据空缺，
# 因此需要获取这些空缺数据来填补表格
for (i in colnames(data)){
  add <- subset(raremax,raremax[i]>=raremax)
  print(add)
  x <- rarefy(data[i],c(1,seq(100,raremax[i],100),add))
  x <- as.data.frame(x)
  rarefaction.data <- rbind.fill(rarefaction.data,x)
}
rarefaction.data.t <- as.data.frame(t(rarefaction.data))
colnames(rarefaction.data.t) <- colnames(data)
rownames(rarefaction.data.t) <- str_replace(rownames(rarefaction.data.t),"N","")

# 将表格重新排序使之更加美观
rarefaction.data.sort <- 
  rarefaction.data.t[order(as.integer(rownames(rarefaction.data.t))),]

# 保存数据
write.xlsx(rarefaction.data.sort,"稀释性曲线图数据.xlsx",rowNames=T)

#绘制稀释性曲线图
plot_data <- as.data.frame(rarefaction.data.t)
plot_data$numsample <- as.integer(rownames(plot_data))
plot_data.melt <- melt(plot_data,id.vars=c("numsample"))
names(plot_data.melt)[2] <- "sample"

ggplot(plot_data.melt,aes(x=numsample,y=value,color=sample))+
  geom_line(size=1)+
  theme_classic()+
  labs(x="Number of Sequences",y="Number of OTUs",title="Multy samples Rarefaction Curves")+
  theme(plot.title = element_text(hjust = 0.5))+
  annotate(geom="text", x=4000, y=600,label="B77")
  ggsave("稀释性曲线图.pdf", width =8.0, height = 6.0)+
  ggsave("稀释性曲线图.png", width =8.0, height = 6.0)



# 2.2丰度等级图
# **************
setwd("C:\\Users\\xhx\\Desktop\\生信分析\\2.2丰度等级图")

ra <- relative_abundance.file[,1:9]

# 排序
s <- function(x){sort(x,decreasing=T)}   # 降序
ra.sort <- apply(ra,2,s)
head(ra.sort)

# 保存排序后的文件
write.xlsx(ra.sort,file="丰度等级图数据.xlsx",rowNames=T)

# 绘制丰度等级图

ra.sort.melt <- melt(ra.sort)
colnames(ra.sort.melt)[2] <- "sample"
plot_data2 <- subset(ra.sort.melt,ra.sort.melt$value!=0)
max(plot_data2$value)
min(plot_data2$value)

ggplot(plot_data2,aes(x=Var1,y=log10(value),color=sample))+
  geom_line()+
  theme_classic()+
  labs(x="OTU Rank",y="Relative Abundance(log10)")+
  ggtitle("Rank abundance distribution curve")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks = -1:-5,
                     limits = c(-5, 0))+
  ggsave("丰度等级图.pdf", width =8.0, height = 6.0)



# 2.3多样本相似性树状图
# *******************
# 分别基于Bray-Curtis和Jaccard来构造相异矩阵
setwd("C:\\Users\\xhx\\Desktop\\生信分析\\2.3多样本相似性树状图")

# 采用之前保存的相对丰度表格
# 仅包含数据的相对丰度表格
ra.t <- t(relative_abundance.file[,1:9])
head(ra)

# 基于Bray-Curtis构建数据的相异矩阵
data.ra.bc <- vegdist(ra.t)
data.ra.bc.matrix <- as.matrix(data.ra.bc)

# 基于Jaccard构建数据(二元)相异矩阵
# 设置binary=T来使多度数据转化为二元数据
data.ra.jac <- vegdist(ra.t,"jaccard",binary=T)
data.ra.jac.matrix <- as.matrix(data.ra.jac)

# 保存这两组数据
write.xlsx(data.ra.bc.matrix,"多样本差异性矩阵Bray_cutis.matrix.xlsx")
write.xlsx(data.ra.bc.matrix,"多样本差异性矩阵Jaccard.matrix.xlsx")

# 分别绘制热图
pheatmap(data.ra.bc.matrix,main="Heatmap(Bray-Curtis)")
pheatmap(data.ra.jac.matrix,main="Heatmap(Jaccard)")

# 绘制聚类树
data.ra.bc.hc <- hclust(data.ra.bc,"complete")
dend <-as.dendrogram(data.ra.bc.hc)
dend %>% plot(horiz = TRUE)
abline(v = heights_per_k.dendrogram(dend)["3"],lwd = 2,lty = 2,col = "green")
abline(v = heights_per_k.dendrogram(dend)["4"],lwd = 2,lty = 2,col = "blue")
abline(v = heights_per_k.dendrogram(dend)["5"],lwd = 2,lty = 2,col = "red")
abline(v = heights_per_k.dendrogram(dend)["6"],lwd = 2,lty = 2,col = "grey")
abline(v = heights_per_k.dendrogram(dend)["8"],lwd = 2,lty = 2,col = "orange")
data.ra.jac.hc <- hclust(data.ra.jac,"average")
plot(data.ra.jac.hc,hang=-1)



# 2.4主成分分析PCA
# ******************
setwd("C:\\Users\\xhx\\Desktop\\生信分析\\2.4主成分分析PCA")

# 物种数据Hellinger转化以避免双零问题
data.hel <- decostand(raw_data[,1:9],"hellinger")

pca <- rda(data.hel)
axis <- pca$CA$v
percentage <- round(100*pca$CA$eig/sum(pca$CA$eig),1)
axis.dataframe <- as.data.frame(axis)

# 绘图
ggplot(axis.dataframe,aes(x=PC1,y=PC2))+
  geom_point(color="blue",size=2) +
  geom_text_repel(aes(label=rownames(axis)), size=3)+
  geom_hline(yintercept=0,colour="gray65")+
  geom_vline(xintercept=0,colour="gray65")+
  xlab(paste("PC1","(",percentage[1],"%)"))+
  ylab(paste("PC2","(",percentage[2],"%)"))+
  ggtitle("PCA plot")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  ggsave("PCA.pdf", width = 8, height = 6)


# 3.0各级别物种信息(包含各级别物种序列统计和相对丰度计算)

setwd("C:\\Users\\xhx\\Desktop\\生信分析\\3.0.各级别物种信息统计表格")

# 各级别物种序列统计
level_info <- data.frame()
classification.type <- c("phylum","class","order","family","genus","species")
for (i in 1:6){
  info <- read.xlsx(paste(classification.type[i],".xlsx",sep=""),rowNames=T)
  info$OTUsize <- rowSums(info[,1:9])
  info$number <- i
  info$classification <- classification.type[i]
  info$tax <- rownames(info)
  info.order <- info[order(info$OTUsize),]
  level_info <- rbind(level_info,info.order)

}

level_info <- level_info[c(11:13,1:10)]

write.xlsx(level_info,"1.各样本物种序列统计.xlsx",rowNames=F)







# 相对丰度计算
level_info.percentage <- data.frame()
y <- function(x){x/sum(x)}
for (i in 1:6){
  level <- subset(level_info,level_info$number==i)
  level.percentage <- apply(level[,4:13],2,y)
  level_info.percentage <- rbind(level_info.percentage,level.percentage)
}
level_info.percentage <- cbind(level_info[,1:3],level_info.percentage)
write.xlsx(level_info.percentage,"2.各样本物种相对丰度统计.xlsx",rowNames=F)




# 3.1热图Heatmap
# ***************
# 分别绘制门纲目科属种的热图

setwd("C:\\Users\\xhx\\Desktop\\生信分析\\3.1热图Heatmap\\绘图数据")

# 计数unclassfied的OTU的函数
unclassfied <- function(x){sum(is.na(x))}

# 用一个嵌套循环分别生成各级别绘图所需的数据，然后分别进行绘图(取前50个数据)
for (i in c("phylum","class","order","family","genus","species")){
    type <- names(table(raw_data[i]))
    count.df <- data.frame()
    for (j in type){
      count.sub <- subset(raw_data,raw_data[i]==j)
      count <- apply(count.sub[,1:9],2,sum)
      count.df <- rbind(count.df,count)
    }
    rownames(count.df) <- type
    colnames(count.df) <- colnames(raw_data[,1:9])
    count.unclassified.sub <- subset(raw_data,is.na(raw_data[i]))
    count.unclassified <- apply(count.unclassified.sub[,1:9],2,sum)
    count.df <- rbind(count.unclassified, count.df)
    rownames(count.df)[1] <- "unclassified"
    write.xlsx(count.df,file=paste(i,".xlsx",sep=""),rowNames=T)
}

  # 彩虹配色和红黑配色热图的绘制函数
painting <- function(x){
  pheatmap(x,scale="row",
           color=colorRampPalette(c("blue","cyan","green","yellow",
                                    "orange","red"))(50),
           cluster_cols=F,main="Relative Abundance")
  pheatmap(x,scale="row",
           color=colorRampPalette(c("black","red"))(50),
           cluster_cols=F,main="Relative Abundance")
}

# 一共十四张图，依次为门纲目科属种的彩虹配色和红黑配色
y <- function(x){x/sum(x)}
for (i in c("phylum","class","order","family","genus","species")){
  count <- read.xlsx(paste(i,".xlsx",sep=""),rowNames=T)
  count <- as.data.frame(apply(count,2,y))
  if (nrow(count)<=50){
    painting(count)
  }
  else {
    count.sort <-
      subset(count,rowSums(count)>=(sort(rowSums(count),decreasing=T)[50]))
    painting(count.sort)
  }}


# 去掉unclassified后的属水平热图
genus <- read.xlsx("genus.xlsx",rowNames=T)
genus <- subset(genus, rownames(genus)!="unclassified")
genus.sort <-
  subset(genus,rowSums(genus)>=(sort(rowSums(genus),decreasing=T)[20]))
y <- function(x){x/sum(x)}
genus <- as.data.frame(apply(genus.sort,2,y))
pheatmap(genus,
         color=colorRampPalette(c("blue","cyan","green","yellow",
                                  "orange","red"))(50),
         cluster_cols=F,main="Relative Abundance",filename="rm_un.png")
pheatmap(genus,display_numbers = T,
         color=colorRampPalette(c("black","red"))(50),
         cluster_cols=F,main="Relative Abundance")


# 3.2群落组成分析

setwd("C:\\Users\\xhx\\Desktop\\生信分析\\3.2群落组成分析")

barplot_names <- c("phylum","class","order","family","genus","species")
y <- function(x){x/sum(x)}
# 循环绘制六个水平的堆叠柱状图
for (i in 1:6){
  # level <- subset(level_info.percentage,level_info.percentage$number==i)
  level <- read.xlsx(paste0(barplot_names[i],".xlsx"),rowNames=T)
  level <- as.data.frame(apply(level,2,y))
  level$sum <- rowSums(level)
  level.order <- level[order(level$sum),]
  level_0.95 <- subset(level.order,level.order$sum>=0.005)[,-10]
  level_0.05 <- subset(level.order,level.order$sum<0.005)[,-10]
  other <- colSums(level_0.05)
  barplot_data <- rbind(other,level_0.95)
  rownames(barplot_data)[1] <- "other(<0.5%)"
  barplot_data$tax <- 
    factor(rownames(barplot_data), levels = rev(rownames(barplot_data)))
  barplot_data.melt <- melt(barplot_data,id.vars="tax")

  # 绘图
  # 堆叠图的颜色填充类型用fill参数，position='stack'代表堆叠图
  # 为使颜色区分明显，不采用默认绘图颜色，改用RColorBrewer，使用Set1调色板，
  # 由于RColorBrewer调色板最多只包含12种不同颜色，因此使用
  # scale_fill_manual(values = getPalette(colourCount))来自定义颜色数量，使
  # Set1插入现有的调色板来构建具有任意数量颜色的调色板，详情见
  # https://www.cnblogs.com/shaocf/p/9600340.html
  colourCount = length(unique(barplot_data$tax))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  ggplot(barplot_data.melt,aes(x=variable,y=value*100,fill=tax))+
    geom_col(position='stack', width=0.6) +
    labs(x='', y='Relative Abundance(%)') +
    theme(axis.text=element_text(size=12), axis.title = element_text(size = 13)) +
    theme(legend.text=element_text(size=12))+
    theme(panel.grid=element_blank(), panel.background=element_rect(color = 'black', fill = 'transparent')) +
    theme(legend.title=element_blank())+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(values = getPalette(colourCount))+
    ggtitle(barplot_names[i])+
    ggsave(paste(i,".",barplot_names[i],".barplot.pdf",sep=""), width = 10, height = 6)+
    ggsave(paste(i,".",barplot_names[i],".barplot.png",sep=""), width = 10, height = 6)
}


# 绘制门水平饼图
phylum <- read.xlsx("phylum.xlsx",rowNames=T)
y <- function(x){x/sum(x)}
phylum.per <- as.data.frame(apply(phylum,2,y))
phylum.per$sum <- rowSums(phylum.per)
phylum.per.order <- phylum.per[order(phylum.per$sum, decreasing=T),]
top7 <- phylum.per.order[1:7,-10]
rest <- phylum.per.order[-(1:7),-10]
other <- colSums(rest)
plot_data <- rbind(other,top7)
rownames(plot_data)[1] <- "other"
for (i in 1:ncol(plot_data)){
  data <- plot_data[,i]
  names(data) <- rownames(plot_data)
  pie(data,labels="",col=c("skyblue","lightgreen","red","blue","lightyellow",
             "orange","purple","yellow"))
}




# 去除unclassified的属水平群落结构组成柱状图
y <- function(x){x/sum(x)}
genus <- read.xlsx("genus.xlsx",rowNames=T)
genus$sum <- rowSums(genus)
genus.rm_un <- subset(genus, rownames(genus)!="unclassified")
genus.order.rm_un <- genus.rm_un[order(genus.rm_un$sum),]
genus.order.rm_un <- as.data.frame(apply(genus.order.rm_un,2,y))
write.xlsx(genus.order.rm_un,"rm_un_genus.xlsx",rowNames=T)
genus_0.95 <- subset(genus.order.rm_un,genus.order.rm_un$sum>=0.005)[,-10]
genus_0.05 <- subset(genus.order.rm_un,genus.order.rm_un$sum<0.005)[,-10]
other <- colSums(genus_0.05)
barplot_data <- rbind(other,genus_0.95)
rownames(barplot_data)[1] <- "other(<0.5%)"
barplot_data$tax <- 
  factor(rownames(barplot_data), levels = rev(rownames(barplot_data)))
barplot_data.melt <- melt(barplot_data,id.vars="tax")

colourCount = length(unique(barplot_data$tax))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

ggplot(barplot_data.melt,aes(x=variable,y=value*100,fill=tax))+
  geom_col(position='stack', width=0.6) +
  labs(x='', y='Relative Abundance(%)') +
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 13)) +
  theme(legend.text=element_text(size=8))+
  theme(panel.grid=element_blank(), panel.background=element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = getPalette(colourCount))+
  ggtitle("genus")+
  ggsave("genus.remove_unclassified.pdf", width = 10, height = 6)+
  ggsave("genus.remove_unclassified.png", width = 10, height = 6)








# 4.Alpha多样性分析
setwd("C:\\Users\\xhx\\Desktop\\生信分析\\4.alpha多样性分析")

data.t <- t(data)
write.xlsx(data.t,"sample.shared",rowNames=T)
# coverage
# 代表样本文库的覆盖率
data.coverage <- 1 - rowSums(data.t == 1) / rowSums(data.t)
# Shannon指数
data.shannon <- diversity(data.t,"shannon")

# Simpson指数
data.simpson <- diversity(data.t,"simpson")

# sob、chao、ace
est <- estimateR(data.t)
data.sob <- est[1,]
data.chao1 <- est[2,]
data.ace <- est[4,]

# α多样性指数原始数据表
alpha_ori  <- 
  t(rbind(data.coverage,data.sob,data.chao1,data.ace,data.shannon,data.simpson))
colnames(alpha_ori) <- c("coverage","sob","chao1","ace","shannon","simpson")
write.xlsx(alpha_ori,"α多样性指数原始数据.xlsx",rowNames=T)
# 获取稀释数据以计算各个梯度下的指数
# 该数据之前绘制稀释曲线时已经生成

# 各个抽样梯度下的数据,并保存
rare <- function(data,method){
  curve.df <- data.frame()
  raremax <- colSums(data)
  for (i in colnames(data)){
    curve.sample <- data.frame()
    for (j in c(1,seq(100,raremax[i],100),raremax[i])){
      rare.df <- rrarefy(t(data[i]),j)
      if (method %in% c("shannon","simpson")){result <- diversity(rare.df,method)}
      else if (method %in% c("sob","chao1","ace")){
        est <- estimateR(rare.df)
        if (method=="sob"){result <- est[1,]}
        else if(method=="chao"){result <- est[2,]}
        else if(method=="ace"){result <- est[4,]}
        
      }
      
      curve.sample <- rbind(curve.sample,result)
    }
    curve.df <- rbind.fill(curve.df,as.data.frame(t(curve.sample)))
  }
  curve.df <- as.data.frame(t(curve.df))
  colnames(curve.df) <- colnames(data)
  rownames(curve.df) <- c(1,seq(100,max(raremax),100),max(raremax))
  
  write.xlsx(curve.df,paste0(method,"指数数据.xlsx"),rowNames=T)
  
  return(curve.df)
}

# 将数据框转化为ggplot绘图格式并绘图
curve_paint <- function(curve.df,method){
  # 将数据框转化为ggplot绘图格式
  curve.df$numsample <- rownames(curve.df)
  curve.df.melt <- melt(curve.df)
  colnames(curve.df.melt)[2] <- "group"
  
  if (method=="simpson"){
    ggplot(curve.df.melt,aes(x=as.integer(numsample),y=value,col=group))+
      geom_line()+
      theme_classic()+
      #scale_y_continuous(breaks=0.989:0.996:0.001,limits=c(0.98,1))
      ylim(0.98,1)+
      labs(x="Number of Sequences",y="index",title=paste(method,"curve"))+
      theme(plot.title = element_text(hjust = 0.5))+
      ggsave(paste0(method,"指数图.pdf"), width =8.0, height = 6.0)
  }
  else {
    ggplot(curve.df.melt,aes(x=as.integer(numsample),y=value,col=group))+
      geom_line()+
      theme_classic()+
      labs(x="Number of Sequences",y="index",title=paste(method,"curve"))+
      theme(plot.title = element_text(hjust = 0.5))+
      ggsave(paste0(method,"指数图.pdf"), width =8.0, height = 6.0)
  }
}

curve_output <- function(data,method){
  curve_paint(rare(data,method),method)
}
# shannon曲线
curve_output(data,"shannon")
curve_output(data,"simpson")
curve_output(data,"chao")
curve_output(data,"sob")
curve_output(data,"ace")



# 5.Beta多样性分析-基于物种分类分析

setwd("C:\\Users\\xhx\\Desktop\\生信分析\\5.Beta多样性分析-基于物种分类分析")

# 物种各水平相对丰度表
level.ra <- level_info.percentage
colnames(level.ra)
for (i in c("phylum","class","order","family","genus","species")){
  data.sub <- subset(level.ra,level.ra$classification==i)
  data.t <- t(data.sub[,4:12])
  data.bc <- vegdist(data.t)
  data.bc.matrix <- as.matrix(data.bc)
  write.xlsx(data.bc.matrix,paste(i,"Bray-Curtis.matrix.xlsx"),rowNames=T)
  pheatmap(data.bc.matrix,main=paste("Heatmap(",i,")",sep=""),
           filename=paste0(i,".heatmap.pdf"))
}





# 5.Beta多样性分析-基于系统发育分析

setwd("C:\\Users\\xhx\\Desktop\\生信分析\\5.Beta多样性分析-基于系统发育分析")

wei_unif_dis <-
  read.xlsx("差异性矩阵(weighted).matrix.xlsx",rowNames=T)
pheatmap(wei_unif_dis)
weighted <- read.xlsx("差异性矩阵(weighted.unifrac).matrix.xlsx",rowNames=T)
unweighted <- read.xlsx("差异性矩阵unweighted.unifrac.distance.matrix.xlsx",rowNames=T)
unweighted.matrix <- as.matrix(unweighted)



# PCoA

# p_list = c("FactoMineR", "dplyr", "factoextra", "ggpubr", "pca3d")
# 
# for(p in p_list){
#   
#   if (!requireNamespace(p, quietly = TRUE))
#     
#     install.packages(p)
#   
# }
# 
# 
# 
# # 安装github来源R包
# 
# suppressWarnings(suppressMessages(library(devtools)))
# 
# if (!requireNamespace("ggbiplot", quietly = TRUE))
#   
#  install_github("vqv/ggbiplot")
library(factoextra)
pcoa <- cmdscale(unweighted.matrix,k=nrow(unweighted.matrix)-1,eig=T)
percentage <- round(100*pcoa$eig/sum(pcoa$eig),1)
axis.pcoa <- as.data.frame(pcoa$points)


ggplot(data=axis.pcoa,aes(x=V1,y=V2))+
  geom_point(color="blue",size=4) +
  geom_text_repel(aes(label=rownames(axis.pcoa)), size=4)+
  geom_hline(yintercept=0,colour="gray65")+
  geom_vline(xintercept=0,colour="gray65")+
  xlab(paste("PC1","(",percentage[1],"%)"))+
  ylab(paste("PC2","(",percentage[2],"%)"))+
  ggtitle("PCoA plot")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  ggsave("unweighted.unifrac.PCoA.pdf", width = 8, height = 6)+
  ggsave("unweighted.unifrac.PCoA.png", width = 8, height = 6)
# 聚类树(UPGMA)
tree <- hclust(as.dist(weighted),"average")
dend <-as.dendrogram(tree)
dend %>% plot(horiz = TRUE)

abline(v = heights_per_k.dendrogram(dend)["3"],lwd = 2,lty = 2,col = "green")
abline(v = heights_per_k.dendrogram(dend)["4"],lwd = 2,lty = 2,col = "blue")
abline(v = heights_per_k.dendrogram(dend)["5"],lwd = 2,lty = 2,col = "red")
abline(v = heights_per_k.dendrogram(dend)["6"],lwd = 2,lty = 2,col = "grey")
abline(v = heights_per_k.dendrogram(dend)["8"],lwd = 2,lty = 2,col = "orange")
pdf("3.weighted.unifrac.tre.pdf")
dev.off()


unweighted <- read.xlsx("差异性矩阵unweighted.unifrac.distance.matrix.xlsx",rowNames=T)
unweighted.dist <- as.dist(unweighted)
tree <- hclust(unweighted.dist,"average")
dend <-as.dendrogram(tree)
dend %>% plot(horiz = TRUE)

abline(v = heights_per_k.dendrogram(dend)["3"],lwd = 2,lty = 2,col = "green")
abline(v = heights_per_k.dendrogram(dend)["4"],lwd = 2,lty = 2,col = "blue")
abline(v = heights_per_k.dendrogram(dend)["5"],lwd = 2,lty = 2,col = "red")
abline(v = heights_per_k.dendrogram(dend)["6"],lwd = 2,lty = 2,col = "grey")
abline(v = heights_per_k.dendrogram(dend)["8"],lwd = 2,lty = 2,col = "orange")

# NMDS
nmds <- metaMDS(as.dist(weighted),k=nrow(weighted)-1)
axis.nmds <- as.data.frame(nmds$points)
ggplot(data=axis.nmds,aes(x=MDS1,y=MDS2))+
  labs(x="NMDS1",y="NMDS2")+
  geom_point(color="blue",size=2) +
  geom_text_repel(aes(label=rownames(axis.nmds)), size=3)+
  geom_hline(yintercept=0,colour="gray65")+
  geom_vline(xintercept=0,colour="gray65")+
  ggtitle("NMDS plot")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  ggsave("4.weighted.unifrac.NMDS.pdf", width = 8, height = 6)



# 6.RDA
setwd("C:\\Users\\xhx\\Desktop\\生信分析\\6.冗余分析")


genus <- read.xlsx("genus.xlsx",rowNames=T)
genus$sum <- rowSums(genus)
genus.rm_un <- subset(genus, rownames(genus)!="unclassified")
genus.top10 <- genus.rm_un[order(genus.rm_un$sum,decreasing=T),][1:10,1:9]
top10_names <- rownames(genus.top10)
env <- read.xlsx("各样地理化信息表格.xlsx",rowNames=T)
env <- log1p(env)
colnames(env) <- c("altitude","slope","aspect","crown density","area",
                   "pH","soil moisture","organic matter","available N","available P")
env.phy <- env[,1:5]
env.chem <- env[,6:10]

genus.h <- decostand(t(genus[,-10]),"hellinger")


dca_otu <- decorana(genus.h, ira = 0)
phy.rda <- rda(genus.h~.,env.phy)
phy.rda$CCA$eig/sum(phy.rda$CCA$eig)



# 方差膨胀检验
vif.cca(phy.rda)
ordiR2step(rda(genus.h~1,env),scope=formula(phy.rda),direction="forward")
R2 <- RsquareAdj(phy.rda)$r.squared
R2.adj <- RsquareAdj(phy.rda)$adj.r.squared

anova.cca(phy.rda,step = 1000)
anova.cca(phy.rda,by="axis",step=1000)

rda.summary <- summary(phy.rda)
# 由于物种太多，这里只显示相对丰度前十的物种
top10_axis <- subset(rda.summary$species,rownames(rda.summary$species) %in% top10_names)
spe.axis <- as.data.frame(top10_axis[,1:2])
sites.axis <- as.data.frame(rda.summary$sites[,1:2])
env.axis <- as.data.frame(rda.summary$biplot[,1:2])

ggplot() +
  geom_text_repel(data = sites.axis,
                  aes(RDA1,RDA2,label=row.names(sites.axis)),size=4)+
  geom_segment(data = spe.axis,aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "red")+
  geom_text_repel(data = spe.axis,aes(RDA1,RDA2,label=row.names(spe.axis)),col="red")+
  geom_segment(data = env.axis,aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "blue")+
  geom_text_repel(data = env.axis,aes(RDA1,RDA2,label=row.names(env.axis)),col="blue")+
  labs(x="RDA1 42.72%",y="RDA2 26.10%")+
  geom_hline(yintercept=0,linetype=3,size=1) +
  geom_vline(xintercept=0,linetype=3,size=1)+
  theme_bw()+theme(panel.grid=element_blank())







chem.rda <- rda(genus.h ~ pH + `soil moisture`+
                  `available N`,env)
chem.rda$CCA$eig/sum(chem.rda$CCA$eig)


# 方差膨胀检验
vif.cca(chem.rda)
R2 <- RsquareAdj(chem.rda)$r.squared
R2.adj <- RsquareAdj(chem.rda)$adj.r.squared
ordiR2step(rda(genus.h~1,env.chem,scale= FALSE), scope = formula(chem.rda),
                             direction = 'forward', permutations = 999)

anova.cca(chem.rda,step = 1000)
anova.cca(chem.rda,by="axis",step=1000)
rda.summary <- summary(chem.rda)
  # 由于物种太多，这里只显示相对丰度前十的物种
top10_axis <- subset(rda.summary$species,rownames(rda.summary$species) %in% top10_names)
spe.axis <- as.data.frame(top10_axis[,1:2])
sites.axis <- as.data.frame(rda.summary$sites[,1:2])
env.axis <- as.data.frame(rda.summary$biplot[,1:2])

ggplot() +
  geom_text_repel(data = sites.axis,
                  aes(RDA1,RDA2,label=row.names(sites.axis)),size=4)+
  geom_segment(data = spe.axis,aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "red")+
  geom_text_repel(data = spe.axis,aes(RDA1,RDA2,label=row.names(spe.axis)),col="red")+
  geom_segment(data = env.axis,aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "blue")+
  geom_text_repel(data = env.axis,aes(RDA1,RDA2,label=row.names(env.axis)),col="blue")+
  labs(x="RDA1 42.72%",y="RDA2 26.10%")+
  geom_hline(yintercept=0,linetype=3,size=1) +
  geom_vline(xintercept=0,linetype=3,size=1)+
  theme_bw()+theme(panel.grid=element_blank())




setwd("C:\\Users\\xhx\\Desktop\\生信分析\\6.冗余分析")

# spearman相关
library(corrplot)
library(psych)
library(pheatmap)
env <- read.xlsx("各样地理化信息表格.xlsx",rowNames=T)[,-8]
env <- log1p(env)
colnames(env) <- c("Altitude","Slope","Aspect","Crown density","Population area",
                   "pH","Soil moisture","Available N","Available P")
env.phy <- env[,1:5]
env.chem <- env[,6:9]
genus <- read.xlsx("genus.xlsx",rowNames=T)
genus.rb <- as.data.frame(apply(genus,2,y))
genus.rb$sum <- rowSums(genus.rb)
genus.rm_un <- subset(genus.rb, rownames(genus)!="unclassified")
genus.top10 <- genus.rm_un[order(genus.rm_un$sum,decreasing=T),][1:10,-10]

phylum_env_corr <- corr.test(t(genus.top10),env.phy, 
                             method = 'spearman',adjust="none")
phylum_env_corr$r
phylum_env_corr$p

pheatmap(phylum_env_corr$r, 
         display_numbers = matrix(ifelse(phylum_env_corr$p<0.05,"*",""),
                                  nrow(phylum_env_corr$p)))




phylum_env_corr <- corr.test(t(genus.top10),env.chem, 
                             method = 'spearman',adjust="none")
phylum_env_corr$r
phylum_env_corr$p

pheatmap(phylum_env_corr$r, 
         display_numbers = matrix(ifelse(phylum_env_corr$p<0.05,"*",""),
                                  nrow(phylum_env_corr$p)))






library(igraph)
library(psych)
library(openxlsx)

setwd("C:\\Users\\xhx\\Desktop\\缙云黄芩\\生信分析\\3.2群落组成分析")
data <- read.xlsx("rm_un_genus.xlsx",rowNames=T)[-1,]
data <- data[which(rowSums(data) >= 0.005), ]
data.t <- t(data)
occor <- corr.test(data.t,use="pairwise",method="spearman",adjust="fdr",alpha=.05)
occor.r <- occor$r # 取相关性矩阵R值
occor.p <- occor$p # 取相关性矩阵p值
p <- p.adjust(occor.p, method = 'BH') 
occor.r[occor.p>0.01|abs(occor.r)<0.6] <- 0 

diag(occor.r) <- 0 
write.csv(occor.r,file='genus_ccor.csv')
setwd("C:\\Users\\xhx\\Desktop")

edge <- read.csv('边.csv')
edge[which(edge$Weight > 0),'pn'] <- 'p'
edge[which(edge$Weight < 0),'pn'] <- 'n'
head(edge)
write.csv(edge,'edge.csv')

node <- read.csv("节点.csv")
head(node)
info <- read.xlsx("1.各样本OTU分布及物种信息表格.xlsx", rowNames=T)
phy <- data.frame()
for (i in node$Id) {
  phylum <- subset(info, info$genus==i)[1,]
  phy <- rbind(phy,phylum)
}
node$phylum <- phy$phylum
write.csv(node, "node_phy.csv")






# SEM
library(sem)
library(vegan)
library(openxlsx)
library(semPlot)

setwd("C:\\Users\\xhx\\Desktop\\生信分析\\6.冗余分析")

# 环境数据
env_data <- read.xlsx("各样地理化信息表格.xlsx", rowNames=T)

# alpha、beta多样性数据
alpha <- read.xlsx("α多样性指数原始数据.xlsx", rowNames=T)
spe <- read.xlsx("去除物种信息的相对丰度表格.xlsx", rowNames=T)[,-10]
Alpha_diversity <- alpha[,5]
unweighted.matrix <- as.matrix(
  read.xlsx("差异性矩阵unweighted.unifrac.distance.matrix.xlsx", rowNames=T))
pcoa <- cmdscale(unweighted.matrix,k=nrow(unweighted.matrix)-1,eig=T)
Beta_diversity <- pcoa$points[,1]

env_data2 <- cbind(env_data, Alpha_diversity, Beta_diversity)
env_data2 <- decostand(env_data2, "standardize")

# 检验显著相关
summary(lm(available.P ~ altitude, data=env_data2))
summary(lm(organic.matters ~ altitude, data=env_data2))  # p-value: 0.01902
summary(lm(alkaline.N ~ altitude, data=env_data2))
summary(lm(Alpha_diversity ~ altitude, data=env_data2))
summary(lm(soil.moisture ~ altitude, data=env_data2))
summary(lm(population.area ~ altitude, data=env_data2))

summary(lm(available.P ~ pH, data=env_data2))
summary(lm(organic.matters ~ pH, data=env_data2))
summary(lm(alkaline.N ~ pH, data=env_data2))
summary(lm(Alpha_diversity ~ pH, data=env_data2))  # p-value: 0.05056
summary(lm(soil.moisture ~ pH, data=env_data2))
summary(lm(canopy.density ~ pH, data=env_data2))  # p-value: 0.03229

summary(lm(organic.matters ~ canopy.density, data=env_data2))
summary(lm(alkaline.N ~ canopy.density, data=env_data2))
summary(lm(available.P ~ canopy.density, data=env_data2))

summary(lm(available.P ~ population.area, data=env_data2))
summary(lm(alkaline.N ~ population.area, data=env_data2))
summary(lm(organic.matters ~ population.area, data=env_data2))

summary(lm(Alpha_diversity ~ organic.matters, data=env_data2))
summary(lm(Alpha_diversity ~ available.P, data=env_data2))
summary(lm(Alpha_diversity ~ alkaline.N, data=env_data2))
summary(lm(Alpha_diversity ~ population.area, data=env_data2))
summary(lm(Alpha_diversity ~ soil.moisture, data=env_data2))
summary(lm(Alpha_diversity ~ canopy.density, data=env_data2))
summary(lm(Alpha_diversity ~ slope, data=env_data2))

summary(lm(Beta_diversity ~ available.P, data=env_data2))  # p-value: 0.07919
summary(lm(Beta_diversity ~ alkaline.N, data=env_data2))
summary(lm(Beta_diversity ~ population.area, data=env_data2))
summary(lm(Beta_diversity ~ pH, data=env_data2))
summary(lm(Beta_diversity ~ altitude, data=env_data2))
summary(lm(Beta_diversity ~ organic.matters, data=env_data2))  #  p-value: 0.006335
summary(lm(Beta_diversity ~ canopy.density, data=env_data2))
summary(lm(Beta_diversity ~ soil.moisture, data=env_data2))

# sem数据
sem_data <- env_data2[,-c(2,3,5,7,9)]

colnames(sem_data) <- c("Altitude","Canopy_density","pH", "Organic_matters",
                          "Available_P", "Alpha_diversity","Beta_diversity")

cor_num <- cor(sem_data)

sem_model <- specifyModel(text="
                           Altitude -> Organic_matters, ao,  NA
                           pH -> Alpha_diversity, pA,  NA
                           Canopy_density -> pH, cp,  NA
                           Organic_matters -> Beta_diversity, oA,  NA
                           Available_P -> Beta_diversity, aB,  NA
                           
                           Alpha_diversity <-> Alpha_diversity, AA,  NA
                           Altitude <-> Altitude, aa,  NA
                           pH <-> pH, pp,  NA
                           Available_P <-> Available_P, l,  NA
                           Canopy_density <-> Canopy_density, cc,  NA
                           Organic_matters <-> Organic_matters, oo,  NA
                           Beta_diversity <-> Beta_diversity, BB,  NA
                           ")

out_sem <- sem::sem(sem_model, cor_num, nrow(sem_data))

# 回归系数
coef <- out_sem$coeff

# 系数名
coeff_name <- out_sem$semmod[,1]
summary(out_sem)

library(DiagrammeR)
library(semPlot)

pathDiagram(out_sem, edge.labels="values",
            edge.color="red",edge.weight = "proportional",
            edge.font = 20, standardize = T)

semPaths(out_sem, what = "stand",
          style = "lisrel",
          layout = "circle2",
          curveAdjacent = TRUE,
          edge.label.cex = 1.5,
          # exoVar = FALSE,
          #exoCov = FALSE,
          esize  = 5,
          nDigits = 3,
          residuals = F,
          intercepts = TRUE, sizeMan = 10,
         )



