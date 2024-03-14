
library("devEMF")
library("readxl")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("ComplexHeatmap")
library("viridis")

#set working directory

df1<- read.csv("T720.csv") #import csv file of data
df2<- read.csv("Z5871.csv")
df3<- read.csv("8178.csv")
df4<- read.csv("Z5761.csv")
df5<- read.csv("T704.csv")
df6<- read.csv("T766.csv")
df7<- read.csv("Z1331.csv")

WISH<-df1[,1] #store names of mutated genes in vector

df1 <- df1[,-1] %>% #change data to matrix
  as.matrix()
df2 <- df2[,-1] %>% 
  as.matrix()
df3 <- df3[,-1] %>% 
  as.matrix()
df4 <- df4[,-1] %>% 
  as.matrix()
df5 <- df5[,-1] %>% 
  as.matrix()
df6 <- df6[,-1] %>% 
  as.matrix()
df7 <- df7[,-1] %>% 
  as.matrix()

rownames(df1) <- WISH #add gene names to matrix
rownames(df2) <- WISH
rownames(df3) <- WISH
rownames(df4) <- WISH
rownames(df5) <- WISH
rownames(df6) <- WISH
rownames(df7) <- WISH

#define colors of heatmap
#setting colors
my_palette <- colorRampPalette(c("navy", "white", "red","darkred"))(n = 302)

#setting breaks to 0 to 0.2 then 5 to 50 then 50 to 210
breaks1 = c(seq(0, 0.2, length.out = 101),seq(5, 50, length.out = 101),seq(50, 210, length.out = 101))
breaks1[length(breaks1)] <- max(max(df1, na.rm=T),max(breaks1))
breaks1[1] <- min(min(df1, na.rm=T),min(breaks1))

breaks2 = c(seq(0, 0.2, length.out = 101),seq(5, 50, length.out = 101),seq(50, 210, length.out = 101))
breaks2[length(breaks2)] <- max(max(df2, na.rm=T),max(breaks2))
breaks2[1] <- min(min(df2, na.rm=T),min(breaks2))

breaks3 = c(seq(0, 0.2, length.out = 101),seq(5, 50, length.out = 101),seq(50, 210, length.out = 101))
breaks3[length(breaks3)] <- max(max(df3, na.rm=T),max(breaks3))
breaks3[1] <- min(min(df3, na.rm=T),min(breaks3))

breaks4 = c(seq(0, 0.2, length.out = 101),seq(5, 50, length.out = 101),seq(50, 210, length.out = 101))
breaks4[length(breaks4)] <- max(max(df4, na.rm=T),max(breaks4))
breaks4[1] <- min(min(df4, na.rm=T),min(breaks4))

breaks5 = c(seq(0, 0.2, length.out = 101),seq(5, 50, length.out = 101),seq(50, 210, length.out = 101))
breaks5[length(breaks5)] <- max(max(df5, na.rm=T),max(breaks5))
breaks5[1] <- min(min(df5, na.rm=T),min(breaks5))

breaks6 = c(seq(0, 0.2, length.out = 101),seq(5, 50, length.out = 101),seq(50, 210, length.out = 101))
breaks6[length(breaks6)] <- max(max(df6, na.rm=T),max(breaks6))
breaks6[1] <- min(min(df6, na.rm=T),min(breaks6))

breaks7 = c(seq(0, 0.2, length.out = 101),seq(5, 50, length.out = 101),seq(50, 210, length.out = 101))
breaks7[length(breaks7)] <- max(max(df7, na.rm=T),max(breaks7))
breaks7[1] <- min(min(df7, na.rm=T),min(breaks7))

#heatmap
heatmap1<-ComplexHeatmap::pheatmap(df1, cellheight = 20, cellwidth = 30, scale="none", cluster_cols=FALSE,
                    main="T720/Wt", border_color =  "white", color = my_palette, breaks = breaks1,
                    cluster_rows = FALSE,display_numbers = TRUE, number_format = "%.2f", number_color = "darkgrey", fontsize_number = 9,
                    column_names_side = c("top"),row_names_side = c("left"))

heatmap2<-ComplexHeatmap::pheatmap(df2, cellheight = 20, cellwidth = 30, scale="none", cluster_cols=FALSE,
                                   main="Z5871/Wt", border_color =  "white", color = my_palette, breaks = breaks1,
                                   cluster_rows = FALSE,legend = F,display_numbers = TRUE, number_format = "%.2f", number_color = "darkgrey", fontsize_number = 9,
                                   column_names_side = c("top"),row_names_side = c("left"))

heatmap3<-ComplexHeatmap::pheatmap(df3, cellheight = 20, cellwidth = 30, scale="none", cluster_cols=FALSE,
                                   main="8178/Wt", border_color =  "white", color = my_palette, breaks = breaks1,
                                   cluster_rows = FALSE,legend = F,display_numbers = TRUE, number_format = "%.2f", number_color = "darkgrey", fontsize_number = 9,
                                   column_names_side = c("top"),row_names_side = c("left"))

heatmap4<-ComplexHeatmap::pheatmap(df4, cellheight = 20, cellwidth = 30, scale="none", cluster_cols=FALSE,
                                   main="Z5761/Wt", border_color =  "white", color = my_palette, breaks = breaks1,
                                   cluster_rows = FALSE,legend = F,display_numbers = TRUE, number_format = "%.2f", number_color = "darkgrey", fontsize_number = 9,
                                   column_names_side = c("top"),row_names_side = c("left"))

heatmap5<-ComplexHeatmap::pheatmap(df5, cellheight = 20, cellwidth = 30, scale="none", cluster_cols=FALSE,
                                   main="T704/Wt", border_color =  "white", color = my_palette, breaks = breaks1,
                                   cluster_rows = FALSE,legend = F,display_numbers = TRUE, number_format = "%.2f", number_color = "darkgrey", fontsize_number = 9,
                                   column_names_side = c("top"),row_names_side = c("left"))

heatmap6<-ComplexHeatmap::pheatmap(df6, cellheight = 20, cellwidth = 30, scale="none", cluster_cols=FALSE,
                                   main="T766/Wt", border_color =  "white", color = my_palette, breaks = breaks1,
                                   cluster_rows = FALSE,legend = F,display_numbers = TRUE, number_format = "%.2f", number_color = "darkgrey", fontsize_number = 9,
                                   column_names_side = c("top"),row_names_side = c("left"))

heatmap7<-ComplexHeatmap::pheatmap(df7, cellheight = 20, cellwidth = 30, scale="none", cluster_cols=FALSE,
                                   main="Z1331/Wt", border_color =  "white", color = my_palette, breaks = breaks1,
                                   cluster_rows = FALSE,legend = F,display_numbers = TRUE, number_format = "%.2f", number_color = "darkgrey", fontsize_number = 9,
                                   column_names_side = c("top"),row_names_side = c("left"))

#export to EMF #adjust width and height according to plot size
emf(file = "MP_heatmap.emf", width=17, height=23, emfPlus = F)
heatmap1 + heatmap2 + heatmap3 + heatmap4 + heatmap5 + heatmap6 + heatmap7
dev.off()

#export to pdf
pdf("pdf.pdf", width=18, height=23)
heatmap1 + heatmap2 + heatmap3 + heatmap4 + heatmap5 + heatmap6 + heatmap7
dev.off() 
graphics.off() #only if dev.off does not work

