####################################AUC with growthcurver
library("readxl")
library("tidyverse")
library("ggplot2")
library("dplyr")
library("plater")
library("DescTools")
library("scales")
library("lubridate")
library("growthcurver")
library("RColorBrewer")
library("pheatmap")
library("devEMF")

#############################normalize auc 
#set working directory

df1<- read_excel("231020_PM1_500.xlsx") #import file making sure the first sheet includes a table with time, then all the wells
df2<- read_excel("231025_PM1_500.xlsx")
df3<- read_excel("231028_PM1_500.xlsx")
df4<- read_excel("231108_PM1_500.xlsx")
df5<- read_excel("231104_PM1_500.xlsx")
df6<- read_excel("231110_PM1_500.xlsx")
df7<- read_excel("231114_PM1_500.xlsx")
df8<- read_excel("231122_PM1_500.xlsx")
df9<- read_excel("231123_PM1_500.xlsx")

#transform time and table into correct format
df1 <- df1 %>%
  mutate(time = (hour(Time)*60 + minute(Time) + second(Time)/60)) %>%
  dplyr::select(-Time) %>% 
  relocate(time, .before=NULL) %>% 
  mutate_at(vars(-time),funs(.-first(.)))

df2 <- df2 %>%
  mutate(time = (hour(Time)*60 + minute(Time) + second(Time)/60)) %>%
  dplyr::select(-Time) %>% 
  relocate(time, .before=NULL) %>% 
  mutate_at(vars(-time),funs(.-first(.)))

df3 <- df3 %>%
  mutate(time = (hour(Time)*60 + minute(Time) + second(Time)/60)) %>%
  dplyr::select(-Time) %>% 
  relocate(time, .before=NULL) %>% 
  mutate_at(vars(-time),funs(.-first(.)))

df4 <- df4 %>%
  mutate(time = (hour(Time)*60 + minute(Time) + second(Time)/60)) %>%
  dplyr::select(-Time) %>% 
  relocate(time, .before=NULL) %>% 
  mutate_at(vars(-time),funs(.-first(.)))

df5 <- df5 %>%
  mutate(time = (hour(Time)*60 + minute(Time) + second(Time)/60)) %>%
  dplyr::select(-Time) %>% 
  relocate(time, .before=NULL) %>% 
  mutate_at(vars(-time),funs(.-first(.)))

df6 <- df6 %>%
  mutate(time = (hour(Time)*60 + minute(Time) + second(Time)/60)) %>%
  dplyr::select(-Time) %>% 
  relocate(time, .before=NULL) %>% 
  mutate_at(vars(-time),funs(.-first(.)))

df7 <- df7 %>%
  mutate(time = (hour(Time)*60 + minute(Time) + second(Time)/60)) %>%
  dplyr::select(-Time) %>% 
  relocate(time, .before=NULL) %>% 
  mutate_at(vars(-time),funs(.-first(.)))

df8 <- df8 %>%
  mutate(time = (hour(Time)*60 + minute(Time) + second(Time)/60)) %>%
  dplyr::select(-Time) %>% 
  relocate(time, .before=NULL) %>% 
  mutate_at(vars(-time),funs(.-first(.)))

df9 <- df9 %>%
  mutate(time = (hour(Time)*60 + minute(Time) + second(Time)/60)) %>%
  dplyr::select(-Time) %>% 
  relocate(time, .before=NULL) %>% 
  mutate_at(vars(-time),funs(.-first(.)))

#fit the growth curves to a function, which allows calculations 
fit1 <- SummarizeGrowthByPlate(df1)
fit2 <- SummarizeGrowthByPlate(df2)
fit3 <- SummarizeGrowthByPlate(df3)
fit4 <- SummarizeGrowthByPlate(df4)
fit5 <- SummarizeGrowthByPlate(df5)
fit6 <- SummarizeGrowthByPlate(df6)
fit7 <- SummarizeGrowthByPlate(df7)
fit8 <- SummarizeGrowthByPlate(df8)
fit9 <- SummarizeGrowthByPlate(df9)

#add medium and bacteria to the samples
fit1 <- add_plate(fit1, "Layout1.csv",
                 well_ids_column = "sample")
fit2 <- add_plate(fit2, "Layout2.csv",
                 well_ids_column = "sample")
fit3 <- add_plate(fit3, "Layout3.csv",
                  well_ids_column = "sample")
fit4 <- add_plate(fit4, "Layout4.csv",
                  well_ids_column = "sample")
fit5 <- add_plate(fit5, "Layout5.csv",
                  well_ids_column = "sample")
fit6 <- add_plate(fit6, "Layout6.csv",
                  well_ids_column = "sample")
fit7 <- add_plate(fit7, "Layout7.csv",
                  well_ids_column = "sample")
fit8 <- add_plate(fit8, "Layout8.csv",
                  well_ids_column = "sample")
fit9 <- add_plate(fit9, "Layout9.csv",
                  well_ids_column = "sample")

fit <- fit1 %>% 
  bind_rows(fit2) %>% 
  bind_rows(fit3) %>%
  bind_rows(fit4) %>%
  bind_rows(fit5) %>%
  bind_rows(fit6) %>%
  bind_rows(fit7) %>%
  bind_rows(fit8) %>%
  bind_rows(fit9) %>%  
  arrange(Medium, Bacteria) 

notes <- fit %>% filter(note!="")

#chose the colors of the heatmap #two options are added here,  but this can be personalized
#col <- colorRampPalette(brewer.pal(9, "Reds"))(256)
col <- colorRampPalette(brewer.pal(9, "RdBu"))(256)

summary_fit <- fit %>%
  mutate(t_gen_min = t_gen * 60) %>% 
  filter(Bacteria!="Blank") %>% 
  group_by(Medium, Bacteria) %>%
  summarise_if(is.numeric, funs(mean, sd)) %>% 
  arrange(Medium, Bacteria)

#the heatmap will be sorted based on one strain
SL1344 <- summary_fit %>%
  select(Bacteria, Medium, auc_e_mean) %>%
  filter(Bacteria == "S. Tm SB300") %>%
  arrange(auc_e_mean)

target <- c(unique(SL1344$Medium))
target <- factor(target)

#important to add the names of the bacteria, should be the same as in layout file
bacteria <- c("S. Tm SB300", "S. Tm ATCC 14028", "E. coli 8178", "E. coli Z5761", "E. coli T704", "E. coli T720", "E. coli Nissle T766", "E. coli Z5871", "E. coli T711")
bacteria <- factor(bacteria)

summary <- summary_fit %>%
  select(Bacteria, Medium, auc_e_mean) %>%
  spread(Bacteria, auc_e_mean) %>% 
  relocate(any_of(c(bacteria))) %>%
  relocate(Medium, .before = "S. Tm SB300")

summary$Medium <- factor(summary$Medium, levels=c(target))
summary <- summary[order(summary$Medium),]


summary <- summary %>% 
  column_to_rownames(var="Medium") %>%
  as.matrix() %>% 
  t()

auc_500 <- pheatmap(summary, cellheight = 20, cellwidth = 30, scale="none", cluster_cols=FALSE,
                main="auc_e_mean", border_color =  "white", color = col,
                cluster_rows = FALSE, display_numbers = TRUE, number_color = "darkgrey", fontsize_number = 9)

#export as EMF #adjust width and height according to plot size
emf(file = "auc500.emf", width=42, height=5, emfPlus = F)
auc_500
dev.off()

pdf("auc500.pdf", width=45, height=10)
auc_500
dev.off()
