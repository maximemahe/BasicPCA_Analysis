#Set directory
setwd("~\Data")

library("factoextra")
library(ggplot2)
library(ggpubr)

#Import your data into R
#We use tab delimited .txt files
data <- read.table("data.txt", sep = "\t", header = T, row.names = 1)

# Standardize the data sets
df <- data[1:191, 1:13]
df <- df.scaled <- scale(df)

#Parse numeric data matrix ([samples, gene expression])
#Data per region and per condition

#Regional data for WAS only
was.jej <- data[1:24, 1:15]
was.ile <- data[49:72, 1:15]
was.pc <- data[96:119, 1:15]
was.dc <- data[144:167, 1:15]

#Scaled data per region for WAS only
df.was.jej <- df[1:24, 1:13]
df.was.ile <- df[49:72, 1:13]
df.was.pc <- df[96:119, 1:13]
df.was.dc <- df[144:167, 1:13]

#Compute PCA with prcomp function
p1was <- pca.was.jej <- fviz_pca_ind(prcomp(df.was.jej), title = "PCA - Jejunum",
             habillage = was.jej$Condition,
             palette = c("Blue3", "Red3"), # color by groups
             geom.ind = "point", # show points only (nbut not "text")
             ggtheme = theme_classic(),
             legend = "right",
             addEllipses = TRUE)

p2was <- pca.was.ile <- fviz_pca_ind(prcomp(df.was.ile), title = "PCA - Ileum",
              habillage = was.ile$Condition,
              palette = c("Blue3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p3was <- pca.was.pc <- fviz_pca_ind(prcomp(df.was.pc), title = "PCA - Proximal colon",
              habillage = was.pc$Condition,
              palette = c("Blue3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p4was <- pca.was.dc <- fviz_pca_ind(prcomp(df.was.dc), title = "PCA - Distal colon",
              habillage = was.dc$Condition,
              palette = c("Blue3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

#Visual inspection of the data
#plot(pca.was.jej) #Plot individual PCA
ggexport(plotlist = list(p1was,p2was,p3was,p4was), nrow=2, ncol=2, filename = "PCA_WAS.pdf", width = 800, height = 800)
ggexport(plotlist = list(p1was,p2was,p3was,p4was), nrow=2, ncol=2, filename = "PCA_WAS.png", width = 800, height = 800, res = 300)

##############################################################################################

#Parse numeric data matrix ([samples, gene expression])
#Data per region per condition per function

#Scaled data per region for WAS only
df.was.jej.ax <- df[1:24, c(1,2,3,4)]
df.was.ile.ax <- df[49:72, c(1,2,3,4)]
df.was.pc.ax <- df[96:119, c(1,2,3,4)]
df.was.dc.ax <- df[144:167, c(1,2,3,4)]

#Compute PCA with prcomp function
p1was.ax <- pca.was.jej.ax <- fviz_pca_ind(prcomp(df.was.jej.ax), title = "PCA - Jejunum",
             habillage = was.jej$Condition,
             palette = c("Blue3", "Red3"), # color by groups
             geom.ind = "point", # show points only (nbut not "text")
             ggtheme = theme_classic(),
             legend = "right",
             addEllipses = TRUE)

p2was.ax <- pca.was.ile.ax <- fviz_pca_ind(prcomp(df.was.ile.ax), title = "PCA - Ileum",
              habillage = was.ile$Condition,
              palette = c("Blue3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p3was.ax <- pca.was.pc.ax <- fviz_pca_ind(prcomp(df.was.pc.ax), title = "PCA - Proximal colon",
              habillage = was.pc$Condition,
              palette = c("Blue3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p4was.ax <- pca.was.dc.ax <- fviz_pca_ind(prcomp(df.was.dc.ax), title = "PCA - Distal colon",
              habillage = was.dc$Condition,
              palette = c("Blue3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

#Visual inspection of the data
#plot(pca.was.jej.ax) #Plot individual PCA
ggexport(plotlist = list(p1was.ax,p2was.ax,p3was.ax,p4was.ax), nrow=2, ncol=2, filename = "PCA_WAS_AntiOx.pdf", width = 800, height = 800)
ggexport(plotlist = list(p1was.ax,p2was.ax,p3was.ax,p4was.ax), nrow=2, ncol=2, filename = "PCA_WAS_AntiOx.png", width = 800, height = 800, res = 300)

##############################################################################################

#Parse numeric data matrix ([samples, gene expression])
#Data per region per condition per function

#Scaled data per region for WAS and inflamation genes only
df.was.jej.inf <- df[1:24, c(5,6,7,8,9)]
df.was.ile.inf  <- df[49:72, c(5,6,7,8,9)]
df.was.pc.inf <- df[96:119, c(5,6,7,8,9)]
df.was.dc.inf <- df[144:167, c(5,6,7,8,9)]

#Compute PCA with prcomp function
p1was.inf <- pca.was.jej.inf <- fviz_pca_ind(prcomp(df.was.jej.inf), title = "PCA - Jejunum",
             habillage = was.jej$Condition,
             palette = c("Blue3", "Red3"), # color by groups
             geom.ind = "point", # show points only (nbut not "text")
             ggtheme = theme_classic(),
             legend = "right",
             addEllipses = TRUE)

p2was.inf <- pca.was.ile.inf <- fviz_pca_ind(prcomp(df.was.ile.inf), title = "PCA - Ileum",
              habillage = was.ile$Condition,
              palette = c("Blue3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p3was.inf <- pca.was.pc.inf <- fviz_pca_ind(prcomp(df.was.pc.inf), title = "PCA - Proximal colon",
              habillage = was.pc$Condition,
              palette = c("Blue3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p4was.inf <- pca.was.dc.inf <- fviz_pca_ind(prcomp(df.was.dc.inf), title = "PCA - Distal colon",
              habillage = was.dc$Condition,
              palette = c("Blue3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

#Visual inspection of the data
#plot(pca.was.jej.ax) #Plot individual PCA
ggexport(plotlist = list(p1was.inf,p2was.inf,p3was.inf,p4was.inf), nrow=2, ncol=2, filename = "PCA_WAS_Inflammation.pdf", width = 800, height = 800)
ggexport(plotlist = list(p1was.inf,p2was.inf,p3was.inf,p4was.inf), nrow=2, ncol=2, filename = "PCA_WAS_Inflammation.png", width = 800, height = 800, res = 300)

##############################################################################################

#Parse numeric data matrix ([samples, gene expression])
#Data per region per condition per function

#Scaled data per region for WAS and tight junction genes only
df.was.jej.tj <- df[1:24, c(10,11,12,13)]
df.was.ile.tj  <- df[49:72, c(10,11,12,13)]
df.was.pc.tj <- df[96:119, c(10,11,12,13)]
df.was.dc.tj <- df[144:167, c(10,11,12,13)]

#Compute PCA with prcomp function
p1was.tj <- pca.was.jej.tj <- fviz_pca_ind(prcomp(df.was.jej.tj), title = "PCA - Jejunum",
             habillage = was.jej$Condition,
             palette = c("Blue3", "Red3"), # color by groups
             geom.ind = "point", # show points only (nbut not "text")
             ggtheme = theme_classic(),
             legend = "right",
             addEllipses = TRUE)

p2was.tj <- pca.was.ile.tj <- fviz_pca_ind(prcomp(df.was.ile.tj), title = "PCA - Ileum",
              habillage = was.ile$Condition,
              palette = c("Blue3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p3was.tj <- pca.was.pc.tj <- fviz_pca_ind(prcomp(df.was.pc.tj), title = "PCA - Proximal colon",
              habillage = was.pc$Condition,
              palette = c("Blue3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p4was.tj <- pca.was.dc.tj <- fviz_pca_ind(prcomp(df.was.dc.tj), title = "PCA - Distal colon",
              habillage = was.dc$Condition,
              palette = c("Blue3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

#Visual inspection of the data
#plot(pca.was.jej.tj) #Plot individual PCA
ggexport(plotlist = list(p1was.tj,p2was.tj,p3was.tj,p4was.tj), nrow=2, ncol=2, filename = "PCA_WAS_TightJunction.pdf", width = 800, height = 800)
ggexport(plotlist = list(p1was.tj,p2was.tj,p3was.tj,p4was.tj), nrow=2, ncol=2, filename = "PCA_WAS_TightJunction.png", width = 800, height = 800, res = 300)

##########################################################################################
##########################################################################################

#Parse numeric data matrix ([samples, gene expression])
#Data per region and per condition

#Regional data for WAS + treatment
ib.jej <- data[1:48, 1:15]
ib.ile <- data[49:95, 1:15]
ib.pc <- data[96:143, 1:15]
ib.dc <- data[144:191, 1:15]

#Scaled data per region for WAS only
df.ib.jej <- df[1:48, 1:13]
df.ib.ile <- df[49:95, 1:13]
df.ib.pc <- df[96:143, 1:13]
df.ib.dc <- df[144:191, 1:13]

#Compute PCA with prcomp function
p1ib <- pca.ib.jej <- fviz_pca_ind(prcomp(df.ib.jej), title = "PCA - Jejunum",
             habillage = ib.jej$Condition,
             palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
             geom.ind = "point", # show points only (nbut not "text")
             ggtheme = theme_classic(),
             legend = "right",
             addEllipses = TRUE)

p2ib <- pca.ib.ile <- fviz_pca_ind(prcomp(df.ib.ile), title = "PCA - Ileum",
              habillage = ib.ile$Condition,
              palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p3ib <- pca.ib.pc <- fviz_pca_ind(prcomp(df.ib.pc), title = "PCA - Proximal colon",
              habillage = ib.pc$Condition,
              palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p4ib <- pca.ib.dc <- fviz_pca_ind(prcomp(df.ib.dc), title = "PCA - Distal colon",
              habillage = ib.dc$Condition,
              palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

#Visual inspection of the data
#plot(pca.ib.jej) #Plot individual PCA
ggexport(plotlist = list(p1ib,p2ib,p3ib,p4ib), nrow=2, ncol=2, filename = "PCA_IB.pdf", width = 800, height = 800)
ggexport(plotlist = list(p1ib,p2ib,p3ib,p4ib), nrow=2, ncol=2, filename = "PCA_IB.png", width = 800, height = 800, res = 300)

##############################################################################################

#Parse numeric data matrix ([samples, gene expression])
#Data per region per condition per function

#Scaled data per region for WAS and antioxydant genes only
df.ib.jej.ax <- df[1:48, c(1,2,3,4)]
df.ib.ile.ax <- df[49:95, c(1,2,3,4)]
df.ib.pc.ax <- df[96:143, c(1,2,3,4)]
df.ib.dc.ax <- df[144:191, c(1,2,3,4)]

p1ib.ax <- pca.ib.jej.ax <- fviz_pca_ind(prcomp(df.ib.jej.ax), title = "PCA - Jejunum",
             habillage = ib.jej$Condition,
             palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
             geom.ind = "point", # show points only (nbut not "text")
             ggtheme = theme_classic(),
             legend = "right",
             addEllipses = TRUE)

p2ib.ax <- pca.ib.ile.ax <- fviz_pca_ind(prcomp(df.ib.ile.ax), title = "PCA - Ileum",
              habillage = ib.ile$Condition,
              palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p3ib.ax <- pca.ib.pc.ax <- fviz_pca_ind(prcomp(df.ib.pc.ax), title = "PCA - Proximal colon",
              habillage = ib.pc$Condition,
              palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p4ib.ax <- pca.ib.dc.ax <- fviz_pca_ind(prcomp(df.ib.dc.ax), title = "PCA - Distal colon",
              habillage = ib.dc$Condition,
              palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

#Visual inspection of the data
#plot(pca.ib.jej.ax) #Plot individual PCA
ggexport(plotlist = list(p1ib.ax,p2ib.ax,p3ib.ax,p4ib.ax), nrow=2, ncol=2, filename = "PCA_IB_Antioxydant.pdf", width = 800, height = 800)
ggexport(plotlist = list(p1ib.ax,p2ib.ax,p3ib.ax,p4ib.ax), nrow=2, ncol=2, filename = "PCA_IB_Antioxydant.png", width = 800, height = 800, res = 300)

##############################################################################################

#Parse numeric data matrix ([samples, gene expression])
#Data per region per condition per function

#Scaled data per region for WAS and inflamation genes only
df.ib.jej.inf <- df[1:48, c(5,6,7,8,9)]
df.ib.ile.inf  <- df[49:95, c(5,6,7,8,9)]
df.ib.pc.inf <- df[96:143, c(5,6,7,8,9)]
df.ib.dc.inf <- df[144:191, c(5,6,7,8,9)]

p1ib.inf <- pca.ib.jej.inf <- fviz_pca_ind(prcomp(df.ib.jej.inf), title = "PCA - Jejunum",
             habillage = ib.jej$Condition,
             palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
             geom.ind = "point", # show points only (nbut not "text")
             ggtheme = theme_classic(),
             legend = "right",
             addEllipses = TRUE)

p2ib.inf <- pca.ib.ile.inf <- fviz_pca_ind(prcomp(df.ib.ile.inf), title = "PCA - Ileum",
              habillage = ib.ile$Condition,
              palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p3ib.inf <- pca.ib.pc.inf <- fviz_pca_ind(prcomp(df.ib.pc.inf), title = "PCA - Proximal colon",
              habillage = ib.pc$Condition,
              palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p4ib.inf <- pca.ib.dc.inf <- fviz_pca_ind(prcomp(df.ib.dc.inf), title = "PCA - Distal colon",
              habillage = ib.dc$Condition,
              palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

#Visual inspection of the data
#plot(pca.ib.jej.inf) #Plot individual PCA
ggexport(plotlist = list(p1ib.inf,p2ib.inf,p3ib.inf,p4ib.inf), nrow=2, ncol=2, filename = "PCA_IB_Inflammation.pdf", width = 800, height = 800)
ggexport(plotlist = list(p1ib.inf,p2ib.inf,p3ib.inf,p4ib.inf), nrow=2, ncol=2, filename = "PCA_IB_Inflammation.png", width = 800, height = 800, res = 300)

##############################################################################################

#Parse numeric data matrix ([samples, gene expression])
#Data per region per condition per function

#Scaled data per region for WAS and tight junction genes only
df.ib.jej.tj <- df[1:48, c(10,11,12,13)]
df.ib.ile.tj  <- df[49:95, c(10,11,12,13)]
df.ib.pc.tj <- df[96:143, c(10,11,12,13)]
df.ib.dc.tj <- df[144:191, c(10,11,12,13)]

p1ib.tj <- pca.ib.jej.tj <- fviz_pca_ind(prcomp(df.ib.jej.tj), title = "PCA - Jejunum",
             habillage = ib.jej$Condition,
             palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
             geom.ind = "point", # show points only (nbut not "text")
             ggtheme = theme_classic(),
             legend = "right",
             addEllipses = TRUE)

p2ib.tj <- pca.ib.ile.tj <- fviz_pca_ind(prcomp(df.ib.ile.tj), title = "PCA - Ileum",
              habillage = ib.ile$Condition,
              palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p3ib.tj <- pca.ib.pc.tj <- fviz_pca_ind(prcomp(df.ib.pc.tj), title = "PCA - Proximal colon",
              habillage = ib.pc$Condition,
              palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

p4ib.tj <- pca.ib.dc.tj <- fviz_pca_ind(prcomp(df.ib.dc.tj), title = "PCA - Distal colon",
              habillage = ib.dc$Condition,
              palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
              geom.ind = "point", # show points only (nbut not "text")
              ggtheme = theme_classic(),
              legend = "right",
              addEllipses = TRUE)

#Visual inspection of the data
#plot(pca.ib.jej.tj) #Plot individual PCA
ggexport(plotlist = list(p1ib.tj,p2ib.tj,p3ib.tj,p4ib.tj), nrow=2, ncol=2, filename = "PCA_IB_TightJunction.pdf", width = 800, height = 800)
ggexport(plotlist = list(p1ib.tj,p2ib.tj,p3ib.tj,p4ib.tj), nrow=2, ncol=2, filename = "PCA_IB_TightJunction.png", width = 800, height = 800, res = 300)

##########################################################################################
##########################################################################################

#FAMD - Factor Analysis of Mixed Data

library("FactoMineR")
library("factoextra")

data <- read.table("data.txt", sep = "\t", header = T, row.names = 1)
res.famd <- FAMD(data, graph = FALSE) # Standardize the data sets

famd1 <- fviz_mfa_ind(res.famd,
             habillage = "Condition", # color by groups
             palette = c("Blue3", "Green3", "Purple3", "Red3"), # color by groups
             geom = c("point"),
             addEllipses = TRUE, ellipse.type = "confidence",
             repel = TRUE) # Avoid text overlapping

famd2 <- fviz_ellipses(res.famd, c("Organ", "Condition"),
             palette = c("Blue3", "coral", "Green3", "Purple3", "gold1", "orange2", "orangered2", "Red3"),
             geom = c("point"),
             repel = TRUE)

ggexport(plotlist = list(famd1, famd2), nrow=2, filename = "FAMD.pdf", width = 800, height = 800)
ggexport(plotlist = list(famd1, famd2), nrow=2, filename = "FAMD.png", width = 800, height = 800, res = 300)
