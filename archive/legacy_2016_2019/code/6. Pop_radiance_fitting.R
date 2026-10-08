library(ggplot2)
library(tidyverse)
library(plyr)
library(broom)
library(viridis)
library(drc)
library(aomisc)
library(nlme)
library(lme4)
library(ggExtra)
library(broom)
library(data.table)
library(lubridate)
library(lmerTest)
library(broom.mixed)
library(Matrix)


######################### finds the fit per country of the 
## correlation between the av radiance in each NUTS3 area 
## and the population density

data <- read.csv("outputs/2.rad_dataframe.csv", header = T)
data <- data[order(data$ID),]
colnames(data)[2] <- '2019'
colnames(data)[3] <- '2018'
colnames(data)[4] <- '2017'
colnames(data)[5] <- '2016'


cd <- data.frame(ar = substr(data$ID, start = 1, stop = 2), NUTS_ID = data$ID, av_rad = rowMeans(data[,2:5]))
cd <- melt(data,id.vars="ID")
cd$ar <- substr(cd$ID, start = 1, stop = 2)

colnames(cd)[2] <- 'date'
colnames(cd)[3] <- 'av_rad'

pop <- read.table("NUTS3data/NUTS3populations.csv", header=T, sep="\t")
pop <- pop[,1:5]
area <- read.table("NUTS3data/NUTS3area.csv", header=T, sep="\t")
area <- area[,1:2]
pro <- melt(pop,id.vars="id")
are <- melt(area,id.vars="id")
#str(are)

pro$value <- sub(" p", "", pro$value)
pro$value <- as.numeric(pro$value)

#create dates
pro$date <- sub("^.*?X", "", pro$variable)
#pro$date <- as.Date(ISOdate(pro$date, 1, 1))
pro$nuts <- sub('.*\\,', '', pro$id)
pro$ar <- substr(pro$nuts, start = 1, stop = 2)

are$value <- as.numeric(are$value)
#create dates
are$date <- sub("^.*?X", "", are$variable)
#are$date <- as.Date(ISOdate(are$date, 1, 1))
are$nuts <- sub('.*\\,', '', are$id)
are$nuts_cat <- substr(are$nuts, start = 1, stop = 2)


tom <- data.frame(nuts = are$nuts, value = are$value)
pro <- merge(tom, pro, by = "nuts")

pro$value.x <- as.numeric(pro$value.x)
pro$density <- pro$value.y/pro$value.x


names(pro)[1] <- "NUTS_ID"
names(cd)[1] <- "NUTS_ID"
pd <- merge(cd, pro, by = c("NUTS_ID", "date"))
pd$NUTS_ID <- as.factor(pd$NUTS_ID)
pd$ar.x <- as.factor(pd$ar.x)
names(pd)[4] <- "country"


#analysis using lmer
m1 <- lmer(log10(av_rad)~ log10(density)+date+country+country:log10(density) + (1|NUTS_ID), na.omit(pd), REML=F)

m2 <- update(m1, ~.-country:log10(density))
anova(m1, m2, test="Chisq")

m3 <- update(m1, ~.-date)
anova(m1, m3, test="Chisq")




#plot all data and all years
ggplot(data= pd, aes(x = log10(density), y = log10(av_rad), color=country)) + geom_point() + scale_color_viridis_d() + 
  geom_smooth(method = lm, se = FALSE)


#fit data with lm controlling for for year and use the setDT commmand to do it for each country - in the data table are the slopes and intercepts
tn <- subset(pd, pd$date == '2019')

slopes <- setDT(tn)[, list(Slope = summary(lm(log10(av_rad) ~ log10(density)))$coeff[2], Intercept = summary(lm(log10(av_rad)~log10(density)))$coeff[1]), country]
new2019 <- arrange(slopes, Slope)

tn <- subset(pd, pd$date == '2018')

slopes <- setDT(tn)[, list(Slope = summary(lm(log10(av_rad) ~ log10(density)))$coeff[2], Intercept = summary(lm(log10(av_rad)~log10(density)))$coeff[1]), country]
new2018 <- arrange(slopes, Slope)

tn <- subset(pd, pd$date == '2017')
tn <- tn[complete.cases(tn), ]

slopes <- setDT(tn)[, list(Slope = summary(lm(log10(av_rad) ~ log10(density)))$coeff[2], Intercept = summary(lm(log10(av_rad)~log10(density)))$coeff[1]), country]
new2017 <- arrange(slopes, Slope)

tn <- subset(pd, pd$date == '2016')
tn <- tn[complete.cases(tn), ]
slopes <- setDT(tn)[, list(Slope = summary(lm(log10(av_rad) ~ log10(density)))$coeff[2], Intercept = summary(lm(log10(av_rad)~log10(density)))$coeff[1]), country]
new2016 <- arrange(slopes, Slope)

### merge alll these data together
data_CP <-  merge(new2019, new2018, by = 'country', all = TRUE)
colnames(data_CP)[2] <- 'slope2019'
colnames(data_CP)[3] <- 'intercept2019'
colnames(data_CP)[4] <- 'slope2018'
colnames(data_CP)[5] <- 'intercept2018'
data_CP <-  merge(data_CP, new2017, by = 'country', all = TRUE)
colnames(data_CP)[6] <- 'slope2017'
colnames(data_CP)[7] <- 'intercept2017'
data_CP <-  merge(data_CP, new2016, by = 'country', all = TRUE)
colnames(data_CP)[8] <- 'slope2016'
colnames(data_CP)[9] <- 'intercept2016'

data_CP <- data.frame(data_CP)

data_CP$avSlope <- rowMeans(subset(data_CP, select = c(slope2019, slope2018, slope2017, slope2016)),na.rm=TRUE)
data_CP$avIntercept <- rowMeans(subset(data_CP, select = c(intercept2019, intercept2018, intercept2017, intercept2016)),na.rm=TRUE)

### at this point data_CP has all the details of the fits for all the years and means in.

#vis <- melt(10^(subset(data_CP, select = c(intercept2019, intercept2018, intercept2017, intercept2016))))
#ggplot(vis, aes(x = variable, y = value))+geom_boxplot() + ylim(0,0.1)

### plot the fits for a few example countries and original data
pg <- subset(pd, pd$country == 'UK' | pd$country == 'DE' | 
               pd$country == 'NL' | pd$country == 'SE' |
               pd$country == 'FR' | pd$country == 'IT' |
               pd$country == 'NO' | pd$country == 'ES')
#dq <- subset(data_CP,data_CP$country == 'UK' | data_CP$country == 'DE' | data_CP$country == 'NL' | data_CP$country == 'SE')

pg$date <- as.Date(pg$date, format='%Y')

ggplot(data= pg, aes(x = log10(density), y = log10(av_rad), colour = as.factor(date))) + geom_point() + scale_color_viridis_d() + 
  facet_wrap(~ country, nrow = 2) +
  geom_smooth(method=lm, se=F)

str(pg)

#save the data
write.csv(data_CP, "outputs/6.rad_density_relationships.csv")

#################################################################################################

### just a few countries

pg <- subset(pd, pd$ar.x == 'LI' | pd$ar.x == 'MT' | pd$ar.x == 'LU' | pd$ar.x == 'SK')

ggplot(data= pg, aes(x = log10(density), y = log10(av_rad), color=ar.x)) + geom_point() + scale_color_viridis_d() + 
  geom_smooth(method = lm, se = FALSE)















##subset justa few countries
ggplot(data= pg, aes(x = log10(denmean), y = log10(radmean), color=ar, size = (gdpmean))) + geom_point() + scale_color_viridis_d()

#maximum model
m1 <- glm(log10(radmean) ~  log10(denmean)+ar, data = pd)
summary(m1)

m2 <- update(m1, ~.-log10(denmean))
anova(m1, m2, test = 'Chisq')
summary(m1)


m3 <- update(m1, ~.-ar)
anova(m1, m3, test = 'Chisq')
summary(m1)



#Germany
test <- subset(pg, pg$ar == 'DE')
newdata = data.frame(denmean = test$denmean, ar= rep('DE', 401), gdpmean = test$gdpmean)
gj <-predict(m1, newdata, type="response")
dt <- data.frame(ar = rep('DE', 401), denmean = test$denmean, gdpmean = test$gdpmean, pro_rad = gj, radmean = log10(test$radmean))

final <- dt
#UK
test <- subset(pg, pg$ar == 'UK')
newdata = data.frame(denmean = test$denmean, ar= rep('UK', 179), gdpmean = test$gdpmean)
gj <-predict(m1, newdata, type="response")
dt <- data.frame(ar = rep('UK', 179), denmean = test$denmean, gdpmean = test$gdpmean, pro_rad = gj, radmean = log10(test$radmean))

final <- rbind(final, dt)

#IT
test <- subset(pg, pg$ar == 'IT')
newdata = data.frame(denmean = test$denmean, ar= rep('IT', 110), gdpmean = test$gdpmean)
gj <-predict(m1, newdata, type="response")
dt <- data.frame(ar = rep('IT', 110), denmean = test$denmean, gdpmean = test$gdpmean, pro_rad = gj, radmean = log10(test$radmean))

final <- rbind(final, dt)


##plot real versus predicted for 3 countries, DE, UK nad IT
ggplot(data = final) + geom_point(aes(x = log10(denmean), y = radmean, color = ar, size = gdpmean)) +
  geom_line(aes(x = log10(denmean), y = pro_rad, color = ar, size = 8)) + scale_color_viridis_d() + scale_x_continuous(name="log Population Density [km^-2]") +
  scale_y_continuous(name="log Mean Radiance [nanoWatts cm^-2 sr-1] ")
d