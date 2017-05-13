### Basic random effects estimator for survey biomass smoothing in ADMB
### read in your survey estimates and CV into temp and tempcv however you want...
### start year and end year are your predictions
#### yrs_srv are your survey observation years
#### Example is "roughly" rougheye biomass estimates in the GOA
## Set working directory to location of script and "re.exe"
### use this if your default path is My Documents
#setwd<-path.expand("~/SimpleRE")
### Otherwise set path to working directory
### get the data
temp<-c(5091.1, 43680.7, 44836.6, 61862.6, 45913.1, 39559.5, 43202.1, 47862.3, 59880.1, 50773.6, 50000, 44115.4, 27580.5, 35000)
tempcv<-rep(0.2,14)
styr <-1984
endyr <-2016
yrs_srv<-c(1984,1987,1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015)
### assign variables and write dat file
nobs<-length(yrs_srv)
srv_est<-temp
srv_cv<-tempcv
cat(styr,"\n",endyr,"\n",nobs,"\n",yrs_srv,"\n",srv_est,"\n",srv_cv,"\n",sep=" ",file="re.dat")
### run compiled ADMB model
system("re.exe")
## read the stuf back in 
srv_re<-scan(file="rwout.rep",nlines=1,skip=11)
srv_recv<-scan(file="rwout.rep",nlines=1,skip=21)
## make a data frame of results
ranests<-data.frame(rbind(srv_re,srv_recv))
rownames(ranests)<-c("Estimate","CV")
all_yrs<-seq(styr,endyr)
colnames(ranests)<-all_yrs
ranests
## write them to current direcotry
write.csv(ranests, "ranests.csv")
