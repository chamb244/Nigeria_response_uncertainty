# new spatial time-series extractions
# this file: Nigeria_stochastic_response_250813.R

#### setup ####
setwd("C:/Github/chamb244/EiA2030-ex-ante/Nigeria_response_uncertainty/")

library(terra)
library(randomForest)
library(geodata)

# Nigeria code NGA
adm0 <- gadm(country="NGA", level=0, path=getwd())
adm1 <- gadm(country="NGA", level=1, path=getwd())
adm2 <- gadm(country="NGA", level=2, path=getwd())

#### load data ####
# Bring in price prediction surfaces for 10 years
# this is maize local market price in Dec for 2015-2024, from WB source data
# generated at C:/DATA/prices/WB/NGA_RTFP_mkt_2007_2025-07-28/maize_price_surfaces_250813.R
mprices <- rast("C:/DATA/prices/WB/NGA_RTFP_mkt_2007_2025-07-28/mpkg_Dec_2015-2024.tif")
names(mprices)

# Bring in rain surfaces for 10 years
# these are from CHIRPS dekadal data for Africa
# generated at C:/DATA/prices/WB/NGA_RTFP_mkt_2007_2025-07-28/maize_price_surfaces_250813.R
rain <- rast("C:/DATA/Nigeria/stack/rain/Nigeria_rain_summaries.tif")
names(rain)

#compareGeom(mprices,rain)
mprices <- resample(mprices,rain, method="average")
compareGeom(mprices,rain)

plot(c(mprices[[1]],rain[[1]]))


# bring in response data
response_data_dir <- c("C:/Github/chamb244/EiA2030-ex-ante/Nigeria_response_uncertainty")
list.files(response_data_dir)
trials <- read.csv(paste(response_data_dir, "NGA_Jordan.csv", sep="/"))
trials_pts <- vect(trials, geom=c("longitude", "latitude"))


# bring in components of predictor stack
stack_data_dir <- c("C:/Github/chamb244/EiA2030-ex-ante/Nigeria_response_uncertainty")
soil  <- rast(paste(stack_data_dir, "Nigeria_soil_layers.tif", sep="/"))

# bring in N price (time-invariant for now)
npkg <- rast("C:/Github/chamb244/EiA2030-ex-ante/Nigeria_response_uncertainty/npkg.tif")
# convert to Naira & resample to same extent/resolution as maize prices
nprice <- resample(npkg, mprices)
nprice <- nprice*1535
names(nprice) <- "N price (NGN/kg)"
compareGeom(nprice,mprices)

#### fortify training data with year-appropriate rainfall measures ####
# set up predictor stack 
stack <- c(soil, rain)
names(stack)
# extract values to trials data
xtra <- extract(stack, trials_pts)
trials_pts <- cbind(trials_pts,xtra)
# check years in data
sort(unique(trials_pts$year)) # 1971-2017
dim(trials_pts) # 6750 obs
# keep only years for which we have rainfall (currently from 2001) 
trials_pts.recent <- trials_pts[trials_pts$year > 2000]
dim(trials_pts.recent) # we still have 5187 obs, so not bad; this is 77% of all trials

# fix year-appropriate rain
# define new variable with correct rainfall to use for each year
trials_pts.recent[,"rain.sum"] <- 0
trials_pts.recent[trials_pts.recent$year==2001 ,"rain.sum"] <- trials_pts.recent$rain.sum.2001[trials_pts.recent$year==2001]
trials_pts.recent[trials_pts.recent$year==2002 ,"rain.sum"] <- trials_pts.recent$rain.sum.2002[trials_pts.recent$year==2002]
trials_pts.recent[trials_pts.recent$year==2003 ,"rain.sum"] <- trials_pts.recent$rain.sum.2003[trials_pts.recent$year==2003]
trials_pts.recent[trials_pts.recent$year==2004 ,"rain.sum"] <- trials_pts.recent$rain.sum.2004[trials_pts.recent$year==2004]
trials_pts.recent[trials_pts.recent$year==2005 ,"rain.sum"] <- trials_pts.recent$rain.sum.2005[trials_pts.recent$year==2005]
trials_pts.recent[trials_pts.recent$year==2006 ,"rain.sum"] <- trials_pts.recent$rain.sum.2006[trials_pts.recent$year==2006]
trials_pts.recent[trials_pts.recent$year==2007 ,"rain.sum"] <- trials_pts.recent$rain.sum.2007[trials_pts.recent$year==2007]
trials_pts.recent[trials_pts.recent$year==2008 ,"rain.sum"] <- trials_pts.recent$rain.sum.2008[trials_pts.recent$year==2008]
trials_pts.recent[trials_pts.recent$year==2009 ,"rain.sum"] <- trials_pts.recent$rain.sum.2009[trials_pts.recent$year==2009]
trials_pts.recent[trials_pts.recent$year==2010 ,"rain.sum"] <- trials_pts.recent$rain.sum.2010[trials_pts.recent$year==2010]
trials_pts.recent[trials_pts.recent$year==2011 ,"rain.sum"] <- trials_pts.recent$rain.sum.2011[trials_pts.recent$year==2011]
trials_pts.recent[trials_pts.recent$year==2012 ,"rain.sum"] <- trials_pts.recent$rain.sum.2012[trials_pts.recent$year==2012]
trials_pts.recent[trials_pts.recent$year==2013 ,"rain.sum"] <- trials_pts.recent$rain.sum.2013[trials_pts.recent$year==2013]
trials_pts.recent[trials_pts.recent$year==2014 ,"rain.sum"] <- trials_pts.recent$rain.sum.2014[trials_pts.recent$year==2014]
trials_pts.recent[trials_pts.recent$year==2015 ,"rain.sum"] <- trials_pts.recent$rain.sum.2015[trials_pts.recent$year==2015]
trials_pts.recent[trials_pts.recent$year==2016 ,"rain.sum"] <- trials_pts.recent$rain.sum.2016[trials_pts.recent$year==2016]
trials_pts.recent[trials_pts.recent$year==2017 ,"rain.sum"] <- trials_pts.recent$rain.sum.2017[trials_pts.recent$year==2017]
trials_pts.recent[trials_pts.recent$year==2018 ,"rain.sum"] <- trials_pts.recent$rain.sum.2018[trials_pts.recent$year==2018]
trials_pts.recent[trials_pts.recent$year==2019 ,"rain.sum"] <- trials_pts.recent$rain.sum.2019[trials_pts.recent$year==2019]
trials_pts.recent[trials_pts.recent$year==2020 ,"rain.sum"] <- trials_pts.recent$rain.sum.2020[trials_pts.recent$year==2020]

trials_pts.recent[,"rain.avg"] <- 0
trials_pts.recent[trials_pts.recent$year==2001 ,"rain.avg"] <- trials_pts.recent$rain.avg.2001[trials_pts.recent$year==2001]
trials_pts.recent[trials_pts.recent$year==2002 ,"rain.avg"] <- trials_pts.recent$rain.avg.2002[trials_pts.recent$year==2002]
trials_pts.recent[trials_pts.recent$year==2003 ,"rain.avg"] <- trials_pts.recent$rain.avg.2003[trials_pts.recent$year==2003]
trials_pts.recent[trials_pts.recent$year==2004 ,"rain.avg"] <- trials_pts.recent$rain.avg.2004[trials_pts.recent$year==2004]
trials_pts.recent[trials_pts.recent$year==2005 ,"rain.avg"] <- trials_pts.recent$rain.avg.2005[trials_pts.recent$year==2005]
trials_pts.recent[trials_pts.recent$year==2006 ,"rain.avg"] <- trials_pts.recent$rain.avg.2006[trials_pts.recent$year==2006]
trials_pts.recent[trials_pts.recent$year==2007 ,"rain.avg"] <- trials_pts.recent$rain.avg.2007[trials_pts.recent$year==2007]
trials_pts.recent[trials_pts.recent$year==2008 ,"rain.avg"] <- trials_pts.recent$rain.avg.2008[trials_pts.recent$year==2008]
trials_pts.recent[trials_pts.recent$year==2009 ,"rain.avg"] <- trials_pts.recent$rain.avg.2009[trials_pts.recent$year==2009]
trials_pts.recent[trials_pts.recent$year==2010 ,"rain.avg"] <- trials_pts.recent$rain.avg.2010[trials_pts.recent$year==2010]
trials_pts.recent[trials_pts.recent$year==2011 ,"rain.avg"] <- trials_pts.recent$rain.avg.2011[trials_pts.recent$year==2011]
trials_pts.recent[trials_pts.recent$year==2012 ,"rain.avg"] <- trials_pts.recent$rain.avg.2012[trials_pts.recent$year==2012]
trials_pts.recent[trials_pts.recent$year==2013 ,"rain.avg"] <- trials_pts.recent$rain.avg.2013[trials_pts.recent$year==2013]
trials_pts.recent[trials_pts.recent$year==2014 ,"rain.avg"] <- trials_pts.recent$rain.avg.2014[trials_pts.recent$year==2014]
trials_pts.recent[trials_pts.recent$year==2015 ,"rain.avg"] <- trials_pts.recent$rain.avg.2015[trials_pts.recent$year==2015]
trials_pts.recent[trials_pts.recent$year==2016 ,"rain.avg"] <- trials_pts.recent$rain.avg.2016[trials_pts.recent$year==2016]
trials_pts.recent[trials_pts.recent$year==2017 ,"rain.avg"] <- trials_pts.recent$rain.avg.2017[trials_pts.recent$year==2017]
trials_pts.recent[trials_pts.recent$year==2018 ,"rain.avg"] <- trials_pts.recent$rain.avg.2018[trials_pts.recent$year==2018]
trials_pts.recent[trials_pts.recent$year==2019 ,"rain.avg"] <- trials_pts.recent$rain.avg.2019[trials_pts.recent$year==2019]
trials_pts.recent[trials_pts.recent$year==2020 ,"rain.avg"] <- trials_pts.recent$rain.avg.2020[trials_pts.recent$year==2020]

trials_pts.recent[,"rain.std"] <- 0
trials_pts.recent[trials_pts.recent$year==2001 ,"rain.std"] <- trials_pts.recent$rain.std.2001[trials_pts.recent$year==2001]
trials_pts.recent[trials_pts.recent$year==2002 ,"rain.std"] <- trials_pts.recent$rain.std.2002[trials_pts.recent$year==2002]
trials_pts.recent[trials_pts.recent$year==2003 ,"rain.std"] <- trials_pts.recent$rain.std.2003[trials_pts.recent$year==2003]
trials_pts.recent[trials_pts.recent$year==2004 ,"rain.std"] <- trials_pts.recent$rain.std.2004[trials_pts.recent$year==2004]
trials_pts.recent[trials_pts.recent$year==2005 ,"rain.std"] <- trials_pts.recent$rain.std.2005[trials_pts.recent$year==2005]
trials_pts.recent[trials_pts.recent$year==2006 ,"rain.std"] <- trials_pts.recent$rain.std.2006[trials_pts.recent$year==2006]
trials_pts.recent[trials_pts.recent$year==2007 ,"rain.std"] <- trials_pts.recent$rain.std.2007[trials_pts.recent$year==2007]
trials_pts.recent[trials_pts.recent$year==2008 ,"rain.std"] <- trials_pts.recent$rain.std.2008[trials_pts.recent$year==2008]
trials_pts.recent[trials_pts.recent$year==2009 ,"rain.std"] <- trials_pts.recent$rain.std.2009[trials_pts.recent$year==2009]
trials_pts.recent[trials_pts.recent$year==2010 ,"rain.std"] <- trials_pts.recent$rain.std.2010[trials_pts.recent$year==2010]
trials_pts.recent[trials_pts.recent$year==2011 ,"rain.std"] <- trials_pts.recent$rain.std.2011[trials_pts.recent$year==2011]
trials_pts.recent[trials_pts.recent$year==2012 ,"rain.std"] <- trials_pts.recent$rain.std.2012[trials_pts.recent$year==2012]
trials_pts.recent[trials_pts.recent$year==2013 ,"rain.std"] <- trials_pts.recent$rain.std.2013[trials_pts.recent$year==2013]
trials_pts.recent[trials_pts.recent$year==2014 ,"rain.std"] <- trials_pts.recent$rain.std.2014[trials_pts.recent$year==2014]
trials_pts.recent[trials_pts.recent$year==2015 ,"rain.std"] <- trials_pts.recent$rain.std.2015[trials_pts.recent$year==2015]
trials_pts.recent[trials_pts.recent$year==2016 ,"rain.std"] <- trials_pts.recent$rain.std.2016[trials_pts.recent$year==2016]
trials_pts.recent[trials_pts.recent$year==2017 ,"rain.std"] <- trials_pts.recent$rain.std.2017[trials_pts.recent$year==2017]
trials_pts.recent[trials_pts.recent$year==2018 ,"rain.std"] <- trials_pts.recent$rain.std.2018[trials_pts.recent$year==2018]
trials_pts.recent[trials_pts.recent$year==2019 ,"rain.std"] <- trials_pts.recent$rain.std.2019[trials_pts.recent$year==2019]
trials_pts.recent[trials_pts.recent$year==2020 ,"rain.std"] <- trials_pts.recent$rain.std.2020[trials_pts.recent$year==2020]

trials_pts.recent[,"rain.cv"] <- 0
trials_pts.recent[trials_pts.recent$year==2001 ,"rain.cv"] <- trials_pts.recent$rain.cv.2001[trials_pts.recent$year==2001]
trials_pts.recent[trials_pts.recent$year==2002 ,"rain.cv"] <- trials_pts.recent$rain.cv.2002[trials_pts.recent$year==2002]
trials_pts.recent[trials_pts.recent$year==2003 ,"rain.cv"] <- trials_pts.recent$rain.cv.2003[trials_pts.recent$year==2003]
trials_pts.recent[trials_pts.recent$year==2004 ,"rain.cv"] <- trials_pts.recent$rain.cv.2004[trials_pts.recent$year==2004]
trials_pts.recent[trials_pts.recent$year==2005 ,"rain.cv"] <- trials_pts.recent$rain.cv.2005[trials_pts.recent$year==2005]
trials_pts.recent[trials_pts.recent$year==2006 ,"rain.cv"] <- trials_pts.recent$rain.cv.2006[trials_pts.recent$year==2006]
trials_pts.recent[trials_pts.recent$year==2007 ,"rain.cv"] <- trials_pts.recent$rain.cv.2007[trials_pts.recent$year==2007]
trials_pts.recent[trials_pts.recent$year==2008 ,"rain.cv"] <- trials_pts.recent$rain.cv.2008[trials_pts.recent$year==2008]
trials_pts.recent[trials_pts.recent$year==2009 ,"rain.cv"] <- trials_pts.recent$rain.cv.2009[trials_pts.recent$year==2009]
trials_pts.recent[trials_pts.recent$year==2010 ,"rain.cv"] <- trials_pts.recent$rain.cv.2010[trials_pts.recent$year==2010]
trials_pts.recent[trials_pts.recent$year==2011 ,"rain.cv"] <- trials_pts.recent$rain.cv.2011[trials_pts.recent$year==2011]
trials_pts.recent[trials_pts.recent$year==2012 ,"rain.cv"] <- trials_pts.recent$rain.cv.2012[trials_pts.recent$year==2012]
trials_pts.recent[trials_pts.recent$year==2013 ,"rain.cv"] <- trials_pts.recent$rain.cv.2013[trials_pts.recent$year==2013]
trials_pts.recent[trials_pts.recent$year==2014 ,"rain.cv"] <- trials_pts.recent$rain.cv.2014[trials_pts.recent$year==2014]
trials_pts.recent[trials_pts.recent$year==2015 ,"rain.cv"] <- trials_pts.recent$rain.cv.2015[trials_pts.recent$year==2015]
trials_pts.recent[trials_pts.recent$year==2016 ,"rain.cv"] <- trials_pts.recent$rain.cv.2016[trials_pts.recent$year==2016]
trials_pts.recent[trials_pts.recent$year==2017 ,"rain.cv"] <- trials_pts.recent$rain.cv.2017[trials_pts.recent$year==2017]
trials_pts.recent[trials_pts.recent$year==2018 ,"rain.cv"] <- trials_pts.recent$rain.cv.2018[trials_pts.recent$year==2018]
trials_pts.recent[trials_pts.recent$year==2019 ,"rain.cv"] <- trials_pts.recent$rain.cv.2019[trials_pts.recent$year==2019]
trials_pts.recent[trials_pts.recent$year==2020 ,"rain.cv"] <- trials_pts.recent$rain.cv.2020[trials_pts.recent$year==2020]


#### estimate response model ####
mypredictors <- c("N_fertilizer", "P_fertilizer", "K_fertilizer", "OC", "pH", "sand", "clay", "rain.sum", "rain.cv")
mydepvar <- c("yield") 
myvars <- append(mydepvar, mypredictors)
names(trials_pts.recent)

df <- as.data.frame(trials_pts.recent)
df <- df[complete.cases(df[, myvars]), myvars]

# tune the forest
trf <- tuneRF(x=df[, -1], y=df[,1])
trf
mt <- trf[which.min(trf[,2]),1]
mt

# fit the model
crf <- randomForest(x=df[, -1], y=df[,1],
                    mtry=mt,importance = TRUE)
crf
plot(crf)
importance(crf)
varImpPlot(crf)


#### estimate variability in predicted value coming from joint rainfall & price uncertainty ####

set.seed(1492)

# list of rainfall seasonal totals 
tlist <- c("rain.sum.2015", "rain.sum.2016", "rain.sum.2017", "rain.sum.2018", "rain.sum.2019", "rain.sum.2020", "rain.sum.2021", "rain.sum.2022", "rain.sum.2023", "rain.sum.2024")
# list of rainfall seasonal dekadal CVs 
clist <- c("rain.cv.2015", "rain.cv.2016", "rain.cv.2017", "rain.cv.2018", "rain.cv.2019", "rain.cv.2020", "rain.cv.2021", "rain.cv.2022", "rain.cv.2023", "rain.cv.2024")
# list of predicted prices 
mlist <- c("maipkg_Dec_2015", "maipkg_Dec_2016", "maipkg_Dec_2017", "maipkg_Dec_2018", "maipkg_Dec_2019", "maipkg_Dec_2020", "maipkg_Dec_2021", "maipkg_Dec_2022", "maipkg_Dec_2023", "maipkg_Dec_2024")
# probabilities for each year
# defined such that oldest year is half as likely to be chosen as most recent year, and where all probabilities sum to 1
# see C:\DATA\Nigeria\EiA\rainfall_year_probabillities.xlsx
# this version for 20 years: plist <- c(0.033, 0.035, 0.037, 0.039, 0.040, 0.042, 0.044, 0.046, 0.047, 0.049, 0.051, 0.053, 0.054, 0.056, 0.058, 0.060, 0.061, 0.063, 0.065, 0.067)
# this version for 10 years:
plist <- c(0.067,0.074,0.081,0.089,0.096,0.104,0.111,0.119,0.126,0.133)



# predict using raster stack
# note: we must set the rainfall variables equal to a given year
rain.sum <- stack["rain.sum.2020"]
names(rain.sum) <- c("rain.sum")
rain.cv  <- stack["rain.cv.2020"]
names(rain.cv) <- c("rain.cv")

newstack <- c(stack, rain.sum, rain.cv)

pr <- predict(newstack, crf, const=data.frame(N_fertilizer=100, P_fertilizer=50, K_fertilizer=15), na.rm=TRUE)
names(pr) <- "predicted fertilized yield"
plot(pr) # this is first iteration

# control yield 
pr0 <- predict(newstack, crf, const=data.frame(N_fertilizer=0, P_fertilizer=0, K_fertilizer=0), na.rm=TRUE)
names(pr0) <- "predicted unfertilized yield"

# gains relative to control 
prdif <- pr - pr0
names(prdif) <- "predicted yield gain"

# profit (convert npkg to USD for comparison with maize prices, currently in Naira)
pi <- prdif*mprices["maipkg_Dec_2020"] - 100*nprice
names(pi) <- "predicted profit"

plot(c(pr,pr0,prdif,pi))


rm(pr.X)
rm(pr.0)
rm(pr.dif)
rm(pr.pro)

for (i in 1:10) {
  print(paste("Iteration:", i))  
  # predict using raster stack
  # note: we set the rainfall variables with random draws from last 10 years, as defined in the vectors above
  
  tchoice <-  sample(tlist, 1, replace=TRUE, prob=plist)
  cchoice <-  clist[which(tchoice==tlist)]
  mchoice <-  mlist[which(tchoice==tlist)]
  
  print(tchoice)
  print(cchoice)
  print(mchoice)
  
  rain.sum <- stack[tchoice]
  names(rain.sum) <- c("rain.sum")
  rain.cv  <- stack[cchoice]
  names(rain.cv) <- c("rain.cv")
  mprice  <- mprices[mchoice]
  names(mprice) <- c("mprice")
  
  # replace newstack for each iteration 
  newstack <- c(stack, rain.sum, rain.cv)
  # predict and send each prediction as new layer of output stack
  if (i==1) {
      pr.X <- predict(newstack, crf, const=data.frame(N_fertilizer=100, P_fertilizer=50, K_fertilizer=15), na.rm=TRUE)
      names(pr.X[[i]]) <- paste0("x",i)
  } else {
      add(pr.X) <- predict(newstack, crf, const=data.frame(N_fertilizer=100, P_fertilizer=50, K_fertilizer=15), na.rm=TRUE)
      names(pr.X[[i]]) <- paste0("x",i)
    }
  # control yield 
  if (i==1) {
    pr.0 <- predict(newstack, crf, const=data.frame(N_fertilizer=0, P_fertilizer=0, K_fertilizer=0), na.rm=TRUE)
    names(pr.0[[i]]) <- paste0("x",i)
  } else {
    add(pr.0) <- predict(newstack, crf, const=data.frame(N_fertilizer=0, P_fertilizer=0, K_fertilizer=0), na.rm=TRUE)
    names(pr.0[[i]]) <- paste0("x",i)
  }    
  # gains relative to control 
  if (i==1) {
    pr.dif <- pr.X[[i]] - pr.0[[i]]
    names(pr.dif[[i]]) <- paste0("x",i)
  } else {
    add(pr.dif) <- pr.X[[i]] - pr.0[[i]]
    names(pr.dif[[i]]) <- paste0("x",i)
  }    
  # profit (convert npkg to USD for comparison with maize prices, currently in Naira)
  if (i==1) {
    pr.pro <- pr.dif[[i]]*mprice - 100*nprice
    names(pr.pro[[i]]) <- paste0("x",i)
  } else {
    add(pr.pro) <- pr.dif[[i]]*mprice - 100*nprice
    names(pr.pro[[i]]) <- paste0("x",i)
  }        
}

summary(c(pr,pr0,prdif,pi))

#define threshold profit
threshold_profit <- 5000

rcl_matrix_continuous <- matrix(c(
  0, threshold_profit, 0,    # Values below X become 0
  threshold_profit, Inf, 1   # Values above X become 1
), ncol=3, byrow=TRUE)

junk <- pr.pro[[6]]/1535 # extract one layer and convert to USD
junk <- classify(junk, rcl_matrix_continuous)
plot(junk)
# profitable almost everywhere

##### 251104 ####

# trying now with spatial copula
library(copula)
