library(tidyverse)

temp700 <- read.csv("data/temp700model.csv")
# temp 700
# length at peak 480 seconds 

depth <- (seq(0,20, by=0.1))
colnames(temp700) <- c("time", depth)


temp700.long <-  temp700 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp700.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none")+ geom_hline(yintercept=107, linetype="dashed")

# top 0.5 cm exposed to lethal heat

ggplot(temp700.long, aes(x = as.numeric(depth), y = value)) + geom_line(aes(color = depth)) + theme_bw() + theme(legend.position = "none") + xlim(0,10)



temp450 <- read.csv("data/temp450model.csv")
colnames(temp450) <- c("time", depth)
# temp 450
# length at peak 480 seconds 

temp450.long <-  temp450 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp450.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none") + geom_hline(yintercept=107, linetype="dashed")
# top 0.4 cm exposed to lethal heat 

ggplot(temp450.long, aes(x = as.numeric(depth), y = value)) + geom_line(aes(color = depth)) + theme_bw() + theme(legend.position = "none") + xlim(0,10) + ylim(20,500)+ geom_hline(yintercept=107, linetype="dashed")


temp450.long <- read.csv("data/temp450model_long.csv")
colnames(temp450.long) <- c("time", depth)
# temp 450
# length at peak 720 seconds 

temp450_long.long <-  temp450.long %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp450_long.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none")+ geom_hline(yintercept=107, linetype="dashed")
# top 0.5 cm experience lethal heat, more than 450C for 480 seconds 

ggplot(temp450_long.long, aes(x = as.numeric(depth), y = value)) + geom_line(aes(color = depth)) + theme_bw() + theme(legend.position = "none") + xlim(0,10) + geom_hline(yintercept=107, linetype="dashed")


temp700.long <- read.csv("data/temp700model_long.csv")
colnames(temp700.long) <- c("time", depth)
# temp 700
# length at peak 720 seconds

temp700_long.long <-  temp700 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp700_long.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none")+ geom_hline(yintercept=107, linetype="dashed")
# 0.5 cm experience lethal heat 

ggplot(temp700_long.long, aes(x = as.numeric(depth), y = value)) + geom_line(aes(color = depth)) + theme_bw() + theme(legend.position = "none") + xlim(0,10) + geom_hline(yintercept=107, linetype="dashed")

# find temperature at each time where top cm is burned 
temp1000_time480 <- read.csv("data/temp1000_time480.csv")
colnames(temp1000_time480) <- c("time", depth)
# temp 1000
# length at peak 480 seconds 

temp1000_time480.long <-  temp1000_time480 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp1000_time480.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none") + geom_hline(yintercept=107, linetype="dashed")
# top 0.7 cm exposed to lethal heat 

temp1500_time480 <- read.csv("data/temp1500_time480.csv")
colnames(temp1500_time480) <- c("time", depth)
# temp 1500
# length at peak 480 seconds 

temp1500_time480.long <-  temp1500_time480 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp1500_time480.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none") + geom_hline(yintercept=107, linetype="dashed")
# top 0.8 cm exposed to lethal heat 

temp2000_time480 <- read.csv("data/temperature2000_time480.csv")
colnames(temp2000_time480) <- c("time", depth)
# temp 2000
# length at peak 480 seconds 

temp2000_time480.long <-  temp2000_time480 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp2000_time480.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none") + geom_hline(yintercept=107, linetype="dashed")
# top 1 cm exposed to lethal heat 

temp1500_time720 <- read.csv("data/temp1500_time720.csv")
colnames(temp1500_time720) <- c("time", depth)
# temp 1500
# length at peak 480 seconds 

temp2000_time720.long <-  temp2000_time720 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp2000_time720.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none") + geom_hline(yintercept=107, linetype="dashed")
# top 1.1 cm exposed to lethal heat 



# find time for each temperature 450 and 700 at which the top cm is burned 

temp450_time900 <- read.csv("data/temp450_time90.csv")
colnames(temp450_time900) <- c("time", depth)
# temp 1500
# length at peak 480 seconds 

temp450_time900.long <-  temp450_time900 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp450_time900.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none") + geom_hline(yintercept=107, linetype="dashed")
# top 0.5 cm exposed to lethal heat 


temp450_time1200 <- read.csv("data/temp450_time1200.csv")
colnames(temp450_time1200) <- c("time", depth)
# temp 1500
# length at peak 480 seconds 

temp450_time1200.long <-  temp450_time1200 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp450_time1200.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none") + geom_hline(yintercept=107, linetype="dashed")
# top 0.6 cm exposed to lethal heat 

temp450_time1800 <- read.csv("data/temp450_time1800.csv")
colnames(temp450_time1800) <- c("time", depth)
# temp 1500
# length at peak 480 seconds 

temp450_time1800.long <-  temp450_time1800 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp450_time1800.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none") + geom_hline(yintercept=107, linetype="dashed")
# top 0.6 cm exposed to lethal heat...I am going to stop there and move onto the 700 


temp700_time1800 <- read.csv("data/temp700_time1800.csv")
colnames(temp700_time1800) <- c("time", depth)
# temp 1500
# length at peak 480 seconds 

temp700_time1800.long <-  temp700_time1800 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp700_time1800.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none") + geom_hline(yintercept=107, linetype="dashed")
# top 0.8 cm exposed 



# make models for actula experiment you are running 

temp450time10min <- read.csv("data/temp450time10min.csv")
colnames(temp450time10min) <- c("time", depth)


temp450time10min <-  temp450time10min %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp450time10min, aes(x = time, y = value, color = depth))+ geom_line()  + geom_hline(yintercept=107, linetype="dashed") + theme_light() + theme(legend.position = "none") + labs(title = "Burn at 450 C for 10 Minutes", ylab = "Temperature")

temp450time20min <- read.csv("data/temp450time20min.csv")
colnames(temp450time20min) <- c("time", depth)


temp450time20min <-  temp450time20min %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp450time20min, aes(x = time, y = value, color = depth))+ geom_line()  + geom_hline(yintercept=107, linetype="dashed") + theme_light() + theme(legend.position = "none") + labs(title = "Burn at 450 C for 20 Minutes", ylab = "Temperature")



temp650time10min <- read.csv("data/temp650time600.csv")
colnames(temp650time10min) <- c("time", depth)


temp650time10min <-  temp650time10min %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp650time10min, aes(x = time, y = value, color = depth))+ geom_line()  + geom_hline(yintercept=107, linetype="dashed") + theme_light() + theme(legend.position = "none") + labs(title = "Burn at 650 C for 10 Minutes", y = "Temperature")


temp650time20min <- read.csv("data/temp650time20min.csv")
colnames(temp650time20min) <- c("time", depth)


temp650time20min <-  temp650time20min %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp650time20min, aes(x = time, y = value, color = depth))+ geom_line()  + geom_hline(yintercept=107, linetype="dashed") + theme_light() + theme(legend.position = "none") + labs(title = "Burn at 650 C for 20 Minutes", y = "Temperature")


# change water content to 0.042 to match values from Sherly's Paper 

temp450_time480_wc0.042 <- read.csv("data/temp450_time480_vwc004_tfile.csv")
colnames(temp450_time480_wc0.042) <- c("time", depth)


temp450.long <-  temp450_time480_wc0.042 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp450.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none") + geom_hline(yintercept=107, linetype="dashed")
# top 0.5 cm exposed to lethal heat 

temp650_time480_wc0.042 <- read.csv("data/temp650_time480_vwc004_tfile.csv")
colnames(temp650_time480_wc0.042) <- c("time", depth)


temp650.long <-  temp650_time480_wc0.042 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp650.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none") + geom_hline(yintercept=107, linetype="dashed") + theme_bw()
# top 0.6 cm exposed to lethal heat 



temp480_time600_wc0.042 <- read.csv("data/temp480_time600_vwc004_tfile.csv")
depth <- as.numeric(depth)
colnames(temp480_time600_wc0.042) <- c("time", depth)
temp480_time600_wc0.042 <- temp480_time600_wc0.042[c(1:102)]

temp480.long <-  temp480_time600_wc0.042 %>% pivot_longer(-time, names_to = "depth") 
temp480.long$depth <- as.numeric(temp480.long$depth)

                                        
ggplot(temp480.long %>% filter(depth %in% c(0:10)), aes(x = time, y = value, color = as.factor(depth)))+ geom_line(size = 1) + geom_hline(yintercept=107, linetype="dashed") + theme_bw()+ labs(y = "Temp (C)", x = "Time (sec)") 

ggplot(temp480.long %>% filter(depth %in% c(4:10)), aes(x = time, y = value, color = as.factor(depth)))+ geom_line(size = 1) + theme_bw()+ labs(y = "Temp (C)", x = "Time (sec)") 

temp630_time600_wc0.042 <- read.csv("data/temp630_time600_vwc004_tfile.csv")
colnames(temp630_time600_wc0.042) <- c("time", depth)


temp630.long <-  temp630_time600_wc0.042 %>% pivot_longer(-time, names_to = "depth") 

ggplot(temp630.long, aes(x = time, y = value, color = depth))+ geom_line() + theme(legend.position = "none") + geom_hline(yintercept=107, linetype="dashed")


realData <- read.csv("data/based_on_real_bulkDensity.csv")
depth <- (seq(0,20, by=0.1))
colnames(realData) <- c("time", depth)

realData.long <-  realData %>% pivot_longer(-time, names_to = "depth") 

realData.long$depth <- as.numeric(realData.long$depth)


p450 <- ggplot(realData.long %>% filter(depth %in% c(0:10)), aes(x = time, y = value, color = as.factor(depth)))+ geom_line(size = 1)  + theme_bw()+ labs(y = "Temp (C)", x = "Time (sec)") + ggtitle("Max Temp 450")


realData650 <- read.csv("data/based_on_real_bulkDensity650.csv")
depth <- (seq(0,20, by=0.1))
colnames(realData650) <- c("time", depth)

realData650.long <-  realData650 %>% pivot_longer(-time, names_to = "depth") 

realData650.long$depth <- as.numeric(realData650.long$depth)


p600 <- ggplot(realData650.long %>% filter(depth %in% c(0:10)), aes(x = time, y = value, color = as.factor(depth)))+ geom_line(size = 1)  + theme_bw()+ labs(y = "Temp (C)", x = "Time (sec)") + ggtitle("Max Temp 600")
p600


library(patchwork)

pdf("figures/model graph.pdf")
p450/p600
dev.off()