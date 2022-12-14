library(tidyverse)
library(gridExtra)

temp <- read.csv("data/termocoupleData.csv")

temp$Date <- temp$ï..Date
#temp$Time <- as.POSIXct(temp$Time, format = "%H:%M:%S")
#temp$Time<- format(as.POSIXct(temp$Time,format='%Y/%m/%d %H:%M:%S'),format='%H:%M:%S')
temp$MugA_Time <- as.numeric(temp$MugA_Time)
temp$MugB_Time <- as.numeric(temp$MugB_Time)
temp$MugC_Time <- as.numeric(temp$MugC_Time)
str(temp)

Bare_5_450 <- ggplot(data = temp[210:308,], aes(x = MugA_Time, y = MugA_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugA_Time, y = MugA_5cm), color = "orange") +
  geom_point(aes(x = MugA_Time, y = MugA_8cm), color = "green") + 
  theme_bw()  + ylab("Temperature (C)") + xlab("Time (minutes)") + ggtitle("5 Bare 450")

Bare_12_450 <- ggplot(data = temp[240:350,], aes(x = MugB_Time, y = MugB_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugB_Time, y = MugB_5cm), color = "orange") +
  geom_point(aes(x = MugB_Time, y = MugB_8cm), color = "green") + 
  theme_bw() +  ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("12 Bare 450")

Bare_11_450 <- ggplot(data = temp[282:394,], aes(x = MugC_Time, y = MugC_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugC_Time, y = MugC_5cm), color = "orange")+
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("11 Bare 450")

Bare_7_450 <- ggplot(data = temp[342:560,], aes(x = MugA_Time, y = MugA_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugA_Time, y = MugA_5cm), color = "orange") +
  geom_point(aes(x = MugA_Time, y = MugA_8cm), color = "green") + 
  theme_bw() +ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("7 Bare 450")


Bare_17_450 <- ggplot(data = temp[382:460,], aes(x = MugB_Time, y = MugB_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugB_Time, y = MugB_5cm), color = "orange") +
  geom_point(aes(x = MugB_Time, y = MugB_8cm), color = "green") + 
  theme_bw() +  ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("17 Bare 450")


Bare_16_450 <- ggplot(data = temp[428:526,], aes(x = MugC_Time, y = MugC_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugC_Time, y = MugC_5cm), color = "orange")+
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("16 Bare 450")


grid.arrange(Bare_5_450, Bare_12_450, Bare_11_450, Bare_7_450, Bare_17_450, Bare_16_450, nrow =2)


Bare_9_650 <- ggplot(data = temp[614:689,], aes(x = MugA_Time, y = MugA_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugA_Time, y = MugA_5cm), color = "orange") +
  geom_point(aes(x = MugA_Time, y = MugA_8cm), color = "green") + 
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("9 Bare 600")

Bare_14_650 <- ggplot(data = temp[640:718,], aes(x = MugB_Time, y = MugB_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugB_Time, y = MugB_8cm), color = "orange") +
  geom_point(aes(x = MugB_Time, y = MugB_5cm), color = "green") + 
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("14 Bare 600")

Bare_21_650 <- ggplot(data = temp[678:756,], aes(x = MugC_Time, y = MugC_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugC_Time, y = MugC_5cm), color = "orange")+
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("21 Bare 600")

Bare_2_650 <- ggplot(data = temp[728:800,], aes(x = MugA_Time, y = MugA_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugA_Time, y = MugA_5cm), color = "orange") +
  geom_point(aes(x = MugA_Time, y = MugA_8cm), color = "green") + 
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("2 Bare 600")


Bare_10_650 <- ggplot(data = temp[766:844,], aes(x = MugB_Time, y = MugB_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugB_Time, y = MugB_5cm), color = "orange") +
  geom_point(aes(x = MugB_Time, y = MugB_8cm), color = "green") + 
  theme_bw()+ ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("10 Bare 600")


Bare_19_650 <- ggplot(data = temp[796:877,], aes(x = MugC_Time, y = MugC_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugC_Time, y = MugC_5cm), color = "orange")+
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("19 Bare 600")

grid.arrange(Bare_9_650, Bare_14_650, Bare_21_650, Bare_2_650, Bare_10_650, Bare_19_650, nrow =2)

BSC_5_450 <- ggplot(data = temp[886:968,], aes(x = MugA_Time, y = MugA_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugA_Time, y = MugA_5cm), color = "orange") +
  geom_point(aes(x = MugA_Time, y = MugA_8cm), color = "green") + 
  theme_bw()  + ylab("Temperature (C)") + xlab("Time (minutes)")+ ggtitle("5 BSC 450")

BSC_7_450 <- ggplot(data = temp[922:1004,], aes(x = MugB_Time, y = MugB_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugB_Time, y = MugB_5cm), color = "orange") +
  geom_point(aes(x = MugB_Time, y = MugB_8cm), color = "green") + 
  theme_bw() +  ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("7 BSC 450")

BSC_11_450 <- ggplot(data = temp[962:1038,], aes(x = MugC_Time, y = MugC_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugC_Time, y = MugC_5cm), color = "orange")+
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("11 BSC 450")

BSC_12_450 <- ggplot(data = temp[992:1070,], aes(x = MugA_Time, y = MugA_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugA_Time, y = MugA_5cm), color = "orange") +
  geom_point(aes(x = MugA_Time, y = MugA_8cm), color = "green") + 
  theme_bw() +ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("12 BSC 450")


BSC_16_450 <- ggplot(data = temp[1026:1112,], aes(x = MugB_Time, y = MugB_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugB_Time, y = MugB_5cm), color = "orange") +
  geom_point(aes(x = MugB_Time, y = MugB_8cm), color = "green") + 
  theme_bw() +  ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("16 BSC 450")


BSC_17_450 <- ggplot(data = temp[1064:1144,], aes(x = MugC_Time, y = MugC_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugC_Time, y = MugC_5cm), color = "orange")+
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("17 BSC 450")


grid.arrange(BSC_5_450, BSC_12_450, BSC_11_450, BSC_7_450, BSC_17_450, BSC_16_450, nrow =2)


BSC_2_650 <- ggplot(data = temp[1111:1191,], aes(x = MugA_Time, y = MugA_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugA_Time, y = MugA_5cm), color = "orange") +
  geom_point(aes(x = MugA_Time, y = MugA_8cm), color = "green") + 
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("2 BSC 600")

BSC_9_650 <- ggplot(data = temp[1160:1242,], aes(x = MugB_Time, y = MugB_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugB_Time, y = MugB_5cm), color = "orange") +
  geom_point(aes(x = MugB_Time, y = MugB_8cm), color = "green") + 
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("9 BSC 600")

BSC_10_650 <- ggplot(data = temp[1195:1275,], aes(x = MugC_Time, y = MugC_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugC_Time, y = MugC_5cm), color = "orange")+
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("10 BSC 600")

BSC_14_650 <- ggplot(data = temp[1225:1303,], aes(x = MugA_Time, y = MugA_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugA_Time, y = MugA_5cm), color = "orange") +
  geom_point(aes(x = MugA_Time, y = MugA_8cm), color = "green") + 
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("14 BSC 600")


BSC_19_650 <- ggplot(data = temp[1256:1334,], aes(x = MugB_Time, y = MugB_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugB_Time, y = MugB_5cm), color = "orange") +
  geom_point(aes(x = MugB_Time, y = MugB_8cm), color = "green") + 
  theme_bw()+ ylab("Temperature (C)")+ xlab("Time (minutes)")+ ggtitle("19 BSC 600")


BSC_21_650 <- ggplot(data = temp[1298:1392,], aes(x = MugC_Time, y = MugC_2cm)) + geom_point(color = "red")  +
  geom_point(aes(x = MugC_Time, y = MugC_5cm), color = "orange")+
  theme_bw() + ylab("Temperature (C)")+ xlab("Time (minutes)")+ylim(0,500)+ ggtitle("21 BSC 600")

grid.arrange(BSC_2_650, BSC_9_650, BSC_10_650, BSC_14_650, BSC_19_650, BSC_21_650, nrow =2)

library(patchwork)

pdf("figures/thermocouple data.pdf", width = 8, height = 12)
(Bare_5_450 + Bare_12_450 + Bare_11_450)/ (Bare_7_450 + Bare_17_450 + Bare_16_450) /(Bare_9_650 + Bare_14_650 + Bare_21_650)/ (Bare_2_650 + Bare_10_650 + Bare_19_650) / (BSC_5_450 + BSC_12_450 + BSC_11_450) / (BSC_7_450 + BSC_17_450 + BSC_16_450) /(BSC_9_650 + BSC_14_650 + BSC_21_650)/ (BSC_2_650 + BSC_10_650 + BSC_19_650)
dev.off()