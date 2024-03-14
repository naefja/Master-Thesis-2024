####################################plot growth curves on a single pannel 

library("devEMF")
library("readxl")
library("tidyverse")
library("ggplot2")
library("dplyr")
library("plater")
library("DescTools")
library("lubridate")

#set working directory

df<- read_excel("231020_PM1.xlsx") #import file making sure the first sheet includes a table with time, then all the wells

df <- df %>% #time in minutes and merge the data into a single column with ODs and their respective id
  pivot_longer(-c(Time), names_to = "Well", values_to = "OD") %>%
  mutate(Time = (hour(Time)*60 + minute(Time) + second(Time)/60)) 

#general plot to check if all went well
ggplot(df, aes(x=Time, y=OD, color = Well)) + geom_point() + geom_smooth() +
  scale_y_continuous()

#instead of letters and numbers, change the id to their respective carbon source
growth <- add_plate(df, "Layout.csv",
                    well_ids_column = "Well")

growth <- growth %>% 
  arrange(Well, Time) %>% 
  na.omit(OD)

#plot each growth curve individually
ggplot(growth, aes(x=Time, y=OD, color = Bacteria)) + geom_line() +
  scale_y_continuous() +
  facet_wrap(~Medium) 

#put the files nicely in a pdf or emf
final <- ggplot(growth, aes(x=Time, y=OD, color = Bacteria)) + geom_line() +
  scale_y_continuous() +
  facet_wrap(~Medium) 

#export to EMF #adjust width and height according to plot size
emf(file = "GrowthCurve.emf", width=11.69, height=8.27, emfPlus = F)
final
dev.off()

pdf("GrowthCurves.pdf", width=11.69, height=8.27)
final
dev.off()

##################################################using ipolygrowth package

library("readxl")
library("tidyverse")
library("ggplot2")
library("dplyr")
library("plater")
library("DescTools")
library("lubridate")
library("growthcurver")
library("RColorBrewer")
library("pheatmap")
library("gdata")
library("ipolygrowth")
library("kableExtra")

df<- read_excel("231020_PM1.xlsx") 

df <- df %>%
  pivot_longer(-c(Time), names_to = "Well", values_to = "OD") %>%
  mutate(Time = (hour(Time)*60 + minute(Time) + second(Time)/60)) 

df <- add_plate(df, "Layout.csv",
                 well_ids_column = "Well")

out.multi.f <- ipg_multisample(data = df, id = "Medium", time.name = "Time", y.name = "OD")

#fancy table
summary <- out.multi.f$estimates %>%
  kable() %>% 
  kable_styling("striped", full_width = F) %>% 
  scroll_box(width = "800px", height = "300px")  # table formatting for rmarkdown

#excel table
summary <- out.multi.f$estimates
write.table(summary, file = "summary.csv", dec = ".", sep = ",", row.names = F)
dev.off()
dev.set(dev.next())

#plots with fitted line
ggplot()+
  geom_point(data = df, aes(x = Time, y = OD, color = Bacteria))+ 
  geom_line(data = out.multi.f$fitted, aes(x = time, y = fit))+ 
  facet_wrap(~ Medium)+
  theme_bw()


################################using growthcurver package

df<- read_excel("231020_PM1.xlsx") 

df <- df %>%
  mutate(time = (hour(Time)*60 + minute(Time) + second(Time)/60)) %>%
  dplyr::select(-Time) %>% 
  relocate(time, .before=NULL) %>% 
  mutate_at(vars(-time),funs(.-first(.)))

fit <- SummarizeGrowthByPlate(df)

fit <- add_plate(fit, "Layout.csv",
                    well_ids_column = "sample")

logistic_growth <- function(t, n0, r, k) {
  OD <- k / (1 +((k - n0)/n0) * exp((-r) * t))
}

summary <- fit %>%
  mutate(OD_t_mid = logistic_growth(t_mid, n0, r, k)) %>%
  mutate(slope_derivative = ((k*n0*r*(k-n0)*exp(r*t_mid))/((k+n0*(exp(r*t_mid)-1))^2))) %>% 
  mutate(c = OD_t_mid - (slope_derivative * t_mid)) %>% 
  mutate(lagtime = -c/slope_derivative) %>%
  mutate(t_gen_min = t_gen * 60) %>% 
  group_by(Medium, Bacteria) %>%
  summarise_if(is.numeric, funs(mean, sd)) %>% 
  arrange(Medium, Bacteria, auc_e_mean)


write.table(summary, file = "summary.csv", dec = ".", sep = ",", row.names = F)

dev.off()

dev.set(dev.next())


###########################################################