#################### Exploratory Data Analysis ##################

###### Reading the Data ######

lead <- read.table("C:/Users/Dell/Dropbox/NCSU/ST 537/Project/lead.full.txt", header = F)
colnames(lead) = c("id", "ind.age", "sex", "week", "blood", "trt")
head(lead)

Y= lead$blood
week = lead$week
ind.age = lead$ind.age
sex= lead$sex

### Indicator variable for groups
C1= as.numeric(lead$trt==1)
C2= as.numeric(lead$trt==2)
C3= as.numeric(lead$trt==3)
lead=cbind(lead,C1,C2,C3)


### Groups
lead1 = lead[C1 == 1,1:7]
lead2 = lead[C2 == 1,1:7]
lead3 = lead[C3 == 1,1:7]

### Converting data to wide format
library(tidyr)
lead_wide<- spread(lead, week, blood)
head(lead_wide)


######################## Missing Values ###########################
###################### Number of missing values####################

sum(is.na(lead_wide))

################ Number of missing values across each week and their pattern ##############
library(naniar)

gg_miss_upset(lead_wide[,8:12])

vis_miss(lead_wide[,8:12])

######## Number of children with same number of  missing observations ######################

lead_wide$MissCount = rowSums(is.na(lead_wide))
head(lead_wide)
M=matrix(data=c(0,0,1,0,2,0,3,0,4,0,5,0), nrow=6, byrow = TRUE)
for (i in 0:5){
 M[i+1,2]= nrow(subset(lead_wide,MissCount == i))
  }
colnames(M)=c('No of missing values','No of Entries')
M


################## Box Plots ###############

boxplot(blood~trt,
        data=lead,
        main="Different boxplots for each treatment group",
        xlab="Treatment group",
        ylab="Blood lead level",
        col="orange",
        border="brown"
)


################## Histograms ###############


############################## Group 1 Profile plots ################################################

## Group 1, age<24 months and Male
dat1=lead[(C1 == 1 & ind.age == 0 & sex == 1),]
dat1_wide<- spread(dat1, week, blood)
dat1_wide
y1=colMeans(dat1_wide[,8:12],na.rm = TRUE)


## Group 1, age<24 months and Female
dat2=lead[(C1 == 1 & ind.age == 0 & sex == 0),]
dat2_wide<- spread(dat2, week, blood)
dat2_wide
y2=colMeans(dat2_wide[,8:12],na.rm = TRUE)

### Plotting

xdata= c(0,2,4,6,8)
plot(xdata, y1, type="b", col="blue", pch="M", lty=1, ylim=c(20,30), 
     xlab = "week", ylab = "lead", main = 'Estimated mean trends for \n 
     Group 1 (Placebo) & age<24 months')

lines(xdata, y2,type = "b", lty = 2, col = "red", pch="F")
legend(6,30,legend=c('Male','Female'),
       col=c("blue", "red"), lty=c(1,1))

####################################################################################################################

## Group 1, age>24 months and Male
dat3=lead[(C1 == 1 & ind.age == 1 & sex == 1),]
dat3_wide<- spread(dat3, week, blood)
dat3_wide
y3=colMeans(dat3_wide[,8:12],na.rm = TRUE)


## Group 1, age>24 months and Female
dat4=lead[(C1 == 1 & ind.age == 1 & sex == 0),]
dat4_wide<- spread(dat4, week, blood)
dat4_wide
y4=colMeans(dat4_wide[,8:12],na.rm = TRUE)

### Plotting

xdata= c(0,2,4,6,8)
plot(xdata, y3, type="b", col="blue", pch="M", lty=1, ylim=c(25,32), 
     xlab = "week", ylab = "lead", main = 'Estimated mean trends for 
     \n Group 1(Placebo) & age>24 months')

lines(xdata, y4,type = "b", lty = 2, col = "red", pch="F")
legend(6,31.5,legend=c('Male','Female'),
       col=c("blue", "red"), lty=c(1,1))

###################################### Group 2 Profile plots###########################################################

## Group 2, age<24 months and Male
dat5=lead[(C2 == 1 & ind.age == 0 & sex == 1),]
dat5_wide<- spread(dat5, week, blood)
dat5_wide
y5=colMeans(dat5_wide[,8:12],na.rm = TRUE)


## Group 2, age<24 months and Female
dat6=lead[(C2 == 1 & ind.age == 0 & sex == 0),]
dat6_wide<- spread(dat6, week, blood)
dat6_wide
y6=colMeans(dat6_wide[,8:12],na.rm = TRUE)

### Plotting

xdata= c(0,2,4,6,8)
plot(xdata, y5, type="b", col="blue", pch="M", lty=1, ylim=c(15,26), 
     xlab = "week", ylab = "lead", main = 'Estimated mean trends for 
     \n Group 2(Succimer-Low) & age<24 months')

lines(xdata, y6,type = "b", lty = 2, col = "red", pch="F")
legend(6,26,legend=c('Male','Female'),
       col=c("blue", "red"), lty=c(1,1))

####################################################################################################################

## Group 2, age>24 months and Male
dat7=lead[(C2 == 1 & ind.age == 1 & sex == 1),]
dat7_wide<- spread(dat7, week, blood)
dat7_wide
y7=colMeans(dat7_wide[,8:12],na.rm = TRUE)


## Group 2, age>24 months and Female
dat8=lead[(C2 == 1 & ind.age == 1 & sex == 0),]
dat8_wide<- spread(dat8, week, blood)
dat8_wide
y8=colMeans(dat8_wide[,8:12],na.rm = TRUE)

### Plotting

xdata= c(0,2,4,6,8)
plot(xdata, y7, type="b", col="blue", pch="M", lty=1, ylim=c(20,33), 
     xlab = "week", ylab = "lead", main = 'Estimated mean trends for 
     \n Group 2(Succimer-Low) & age>24 months')

lines(xdata, y8,type = "b", lty = 2, col = "red", pch="F")
legend(6,32.8,legend=c('Male','Female'),
       col=c("blue", "red"), lty=c(1,1))


#################################### Group 3 Profile plots #######################################################

## Group 3, age<24 months and Male
dat9=lead[(C3 == 1 & ind.age == 0 & sex == 1),]
dat9_wide<- spread(dat9, week, blood)
dat9_wide
y9=colMeans(dat9_wide[,8:12],na.rm = TRUE)


## Group 3, age<24 months and Female
dat10=lead[(C3 == 1 & ind.age == 0 & sex == 0),]
dat10_wide<- spread(dat10, week, blood)
dat10_wide
y10=colMeans(dat10_wide[,8:12],na.rm = TRUE)

### Plotting

xdata= c(0,2,4,6,8)
plot(xdata, y9, type="b", col="blue", pch="M", lty=1, ylim=c(15,26), 
     xlab = "week", ylab = "lead", main = 'Estimated mean trends for 
     \n Group 3(Succimer- High) & age<24 months')

lines(xdata, y10,type = "b", lty = 2, col = "red", pch="F")
legend(6,26,legend=c('Male','Female'),
       col=c("blue", "red"), lty=c(1,1))

####################################################################################################################

## Group 3, age>24 months and Male
dat11=lead[(C3 == 1 & ind.age == 1 & sex == 1),]
dat11_wide<- spread(dat11, week, blood)
dat11_wide
y11=colMeans(dat11_wide[,8:12],na.rm = TRUE)


## Group 3, age>24 months and Female
dat12=lead[(C3 == 1 & ind.age == 1 & sex == 0),]
dat12_wide<- spread(dat12, week, blood)
dat12_wide
y12=colMeans(dat12_wide[,8:12],na.rm = TRUE)

### Plotting

xdata= c(0,2,4,6,8)
plot(xdata, y11, type="b", col="blue", pch="M", lty=1, ylim=c(20,30), 
     xlab = "week", ylab = "lead", main = 'Estimated mean trends for
     \n Group 3(Succimer-High) & age>24 months')

lines(xdata, y12,type = "b", lty = 2, col = "red", pch="F")
legend(6,30,legend=c('Male','Female'),
       col=c("blue", "red"), lty=c(1,1))

################################## Profile plot for males with age > 24 ######################################################

plot(xdata, y4, type="b", col="blue", pch="o", lty=1, ylim=c(20,35), 
     xlab = "week", ylab = "lead", main = 'Mean trends for 
     \n males & age>24 months')

lines(xdata, y8,type = "b", lty = 2, col = "red", pch=19)
lines(xdata, y12,type = "b", lty = 2, col = "green", pch=19)

legend(5.5,35,legend=c('Placebo','Succimer-Low','Succimer-High'),
       col=c("green","red", "blue"), lty=c(1,1))
####################################################################################################################


### Plot given data

par(mfrow=c(1,3))
matplot(lead1$id,lead1$blood,  pch=16, col="green", xlab="ID", 
        ylab="lead content", main = 'Placebo group')
plot(lead2$id,lead2$blood,  pch=16, col="red",xlab="ID", 
     ylab="lead content",main = 'Succimer-Low group')
plot(lead3$id,lead3$blood,  pch=16, col="blue",xlab="ID", 
     ylab="lead content",main = 'Succimer- High group')

require(ggplot2)

## define base for the graphs and store in object 'p'
p <- ggplot(data = lead1, aes(x = week, y = blood, group = id))
## just plotting points (a.k.a., scatterplot)
p + geom_point()
## simple spaghetti plot
p + geom_line()

interaction.plot(lead1[1:60,]$week,
                 lead1[1:60,]$id, lead1[1:60,]$blood,
                 xlab="week", ylab="lead content", col=c(1:10), legend=F,pch = c(1:9, 0, letters)) 


interaction.plot(lead2[1:60,]$week,
                 lead2[1:60,]$id, lead2[1:60,]$blood,
                 xlab="week", ylab="lead content", col=c(1:10), legend=F,pch = c(1:9, 0, letters)) 


interaction.plot(lead3[1:60,]$week,
                 lead3[1:60,]$id, lead3[1:60,]$blood,
                 xlab="week", ylab="lead content", col=c(1:10), legend=F,pch = c(1:9, 0, letters)) 



############################ Comparison of three treatments ###############################


### Converting data to wide format
library(tidyr)
lead_wide<- spread(lead, week, blood)
lead_wide

## Number of rows for each treatment
table(lead_wide$trt)


y.p=colMeans(lead_wide[1:40,8:12],na.rm = TRUE)

y.l=colMeans(lead_wide[41:80,8:12],na.rm = TRUE)

y.h=colMeans(lead_wide[81:120,8:12],na.rm = TRUE)



### Plotting

xdata= c(0,2,4,6,8)
plot(xdata, y.p, type="b", col="green", pch="P", lty=1, ylim=c(20,30), 
     xlab = "week", ylab = "lead", main = 'Mean trends for 3 treatment groups')

lines(xdata, y.l,type = "b", lty = 2, col = "red", pch="L")
lines(xdata, y.h,type = "b", lty = 2, col = "blue", pch="H")
legend(5.5,30,legend=c('Placebo','Succimer-Low','Succimer-High'),
       col=c("green","red", "blue"), lty=c(1,1))













