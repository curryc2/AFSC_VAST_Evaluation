## Compute GOA survey biomass estimates from scratch

##### Generalize
## Compute GOA survey biomass estimates from scratch
#### Prep the data ##############
HAUL<-read.csv("HAUL.csv",header=T)
CATCH<-read.csv("CATCH.csv",header=T)
STRATUM<-read.csv("STRATUM2.csv",header=T)
BS<-read.csv("BIENNIAL_SURVEYS.csv",header=T)
BSGOA<-BS[BS$SURVEY_NAME=="GOA BIENNIAL SURVEY"|BS$SURVEY_NAME=="GOA TRIENNIAL SURVEY"|BS$SURVEY_NAME=="2009 Gulf of Alaska Bottom Trawl Survey"|BS$SURVEY_NAME=="2013 Gulf of Alaska Bottom Trawl Survey",]
CPUE<-read.csv("CPUE.csv",header=T)
catch<-CATCH[CATCH$CRUISEJOIN%in%BSGOA$CRUISEJOIN,]
haul<-HAUL[HAUL$CRUISEJOIN%in%BSGOA$CRUISEJOIN,]
cpue<-CPUE[CPUE$HAULJOIN%in%haul$HAULJOIN,]
catch$YEAR<-floor(catch$CRUISE/100)
haul$YEAR<-floor(haul$CRUISE/100)
dim(CATCH)
dim(catch)
dim(CPUE)
dim(cpue)
dim(HAUL)
dim(haul)
haul<-haul[!is.na(haul$STRATUM),]
haul<-haul[haul$PERFORMANCE>-0.01,]
haul<-haul[haul$HAUL_TYPE==3,]
haul<-haul[!is.na(haul$NET_WIDTH),]
# get rid of some more stuff, groundfish only
catch<-catch[catch$SPECIES_CODE>20000,]
catch<-catch[catch$SPECIES_CODE<40000,]
# calculate CPUE per catch records
goastrat<-STRATUM
#year<-c(2013) # Choose number of survey years to estimate from
year<-as.numeric(names(table(haul$YEAR)))
species<-c(30052,30051) #choose species group to choose from
species<-c(30060)
# Orock:
#species<-c(30100,30120,30170,30190,30200,30220,30240,30260,30270,30320,30340,30350,30370,30380,30400,30410,30420,30430,30470,30475,30490,30535,30550,30560,30600)
#species<-c(30576) #shortraker is 30576
### Get haul data and catch data in one object
catch2<-merge(catch,haul,by=c("CRUISEJOIN","HAULJOIN"),all.x=TRUE,all.y=TRUE) # outer join
### Get rid of hauls with no stratum association
cpue<-catch2[!is.na(catch2$STRATUM),]
#### Get rid of "unstatisfactory hauls"
## some changes for orocks
#cpue$ew<-(cpue$STRATUM+100)/100-floor((cpue$STRATUM+100)/100)
#cpue<-cpue[-which(cpue$SPECIES_CODE==30420&cpue$ew<0.4,),]
#dsr<-c(30120,30270,30320,30340,30370,30380,30410,30470)
#cpue<-cpue[-which(cpue$SPECIES_CODE%in%dsr&cpue$ew>0.39),]
############
######
cpue<-cpue[cpue$PERFORMANCE>-0.01,]
#cpue<-cpue[cpue$GEAR!=174,] weird gear type in 1999
### Extract years
temp<-cpue[cpue$YEAR.y%in%year,]
### Get an object with all hauljoins, but no duplicates
temp2<-temp[!duplicated(temp$HAULJOIN),]
### Calculated combined catch for multiple species
cpue2<-cpue[!cpue$SPECIES_CODE%in%species,]
cpue<-cpue[cpue$SPECIES_CODE%in%species,]
cpue2$WEIGHT<-0
cpue2$SPECIES_CODE<-99999
cpue2<-cpue2[!duplicated(cpue2$HAULJOIN),]
cpue3<-rbind(cpue,cpue2)
cpue3$SPECIES_CODE<-"POP"

h2<-haul[!duplicated(haul$HAULJOIN),]

### sum up catch across species for each haul
sumcatch<-as.matrix(tapply(cpue$WEIGHT,cpue$HAULJOIN,sum)) 
### sloppy matrix work
x<-as.numeric(rownames(sumcatch))
sumcatch<-cbind(x,as.numeric(sumcatch[,1]))
colnames(sumcatch)<-c("HAULJOIN","sumcatch")
### name all the species codes to one...
cpue$SPECIES_CODE<-species[1]
#### add new column onto table with combined catch
cpue<-merge(cpue,sumcatch,by="HAULJOIN")
### calculate CPUE (kg/km^2)
cpue$cpue<-cpue$sumcatch/(cpue$DISTANCE_FISHED*cpue$NET_WIDTH/1000)
#### again getting rid of extra records
cpue<-cpue[!duplicated(cpue$HAULJOIN),]
cpue<-cpue[cpue$SPECIES_CODE%in%species[1],]
cpue<-cpue[!is.na(cpue$SPECIES_CODE),]
#pc2013<-pc[pc$YEAR.y%in%year,]
cpue<-cpue[cpue$YEAR.y%in%year,]
temp2$cpue<-NA
temp2[match(cpue$HAULJOIN,temp2$HAULJOIN),]$cpue<-cpue$cpue
temp2[is.na(temp2$cpue),]$cpue<-0
temp2[is.infinite(temp2$cpue),]$cpue<-0
####


cstrat<-tapply(temp2$cpue,temp2$STRATUM,sum,na.rm=T)
vstrat<-tapply(temp2$cpue,temp2$STRATUM,var)
h<-haul[haul$YEAR%in%year,]
hstrat<-tapply(h$HAUL,h$STRATUM,length)
### build stratified esitmates table

xx<-data.frame(cbind(rownames(cstrat),as.numeric(cstrat),as.numeric(vstrat)))
yy<-data.frame(cbind(rownames(hstrat),as.numeric(hstrat)))
names(xx)<-c("STRATUM","SUMC","var")
names(yy)<-c("STRATUM","n")
zz<-merge(yy,xx,all.x=TRUE,all.y=TRUE)
zzz<-merge(zz,goastrat,all.x=TRUE)
zzz$bio<-as.numeric(as.character(zzz[,3]))/as.numeric(as.character(zzz[,2]))*as.numeric(as.character(zzz$AREA))
zzz$var2<-as.numeric(as.character(zzz$AREA))^2*(as.numeric(as.character(zzz[,4]))/as.numeric(as.character(zzz[,2])))

### Output results
print("Biomass")
sum(zzz$bio,na.rm=T)/1000
print("Variance")
sum(zzz$var2,na.rm=T)/1000000
print("SD")
sqrt(sum(zzz$var2,na.rm=T))/1000
print("CV")
#cv
sqrt(sum(zzz$var2,na.rm=T))/1000/(sum(zzz$bio,na.rm=T)/1000)

#write.csv(zzz,"c://sa//orox_strata.csv")

####Cindy's addition
library(plyr)
#for whole GOA
cstrat<-ddply(temp2,c("YEAR.x","STRATUM"),summarize,CPUE=sum(cpue),CPUEvar=var(cpue))
hstrat<-ddply(haul,c("YEAR","STRATUM"),summarize,n_sta=length(unique(HAULJOIN)))
biomvar<-merge(cstrat,hstrat,by.x=c("YEAR.x","STRATUM"),by.y=c("YEAR","STRATUM"),all.x=T)
colnames(biomvar)<-c("YEAR","STRATUM","CPUE","VAR","n_stations")
stratarea<-STRATUM[,c("STRATUM","AREA","INPFC_AREA")]
biomvar<-merge(biomvar,stratarea,by=c("STRATUM"),all.x=T)
biomvar$BIOMASS<-(biomvar$CPUE/biomvar$n_stations)*biomvar$AREA
biomvar$VAR2<-biomvar$AREA^2*(biomvar$VAR/biomvar$n_stations)
BiomassGOA<-ddply(biomvar,c("YEAR"),summarize,Biomass=sum(BIOMASS,na.rm=T)/1000,Variance=sum(VAR2,na.rm=T)/1000000,
                  SD=sqrt(sum(VAR2,na.rm=T))/1000,CV=SD/(sum(BIOMASS,na.rm=T)/1000))

cstrat<-ddply(temp4,c("YEAR.x","STRATUM"),summarize,CPUE=sum(cpue),CPUEvar=var(cpue))
hstrat<-ddply(haul,c("YEAR","STRATUM"),summarize,n_sta=length(unique(HAULJOIN)))
biomvar<-merge(cstrat,hstrat,by.x=c("YEAR.x","STRATUM"),by.y=c("YEAR","STRATUM"),all.x=T)
colnames(biomvar)<-c("YEAR","STRATUM","CPUE","VAR","n_stations")
stratarea<-STRATUM[,c("STRATUM","AREA","INPFC_AREA")]
biomvar<-merge(biomvar,stratarea,by=c("STRATUM"),all.x=T)
biomvar$BIOMASS<-(biomvar$CPUE/biomvar$n_stations)*biomvar$AREA
biomvar$VAR2<-biomvar$AREA^2*(biomvar$VAR/biomvar$n_stations)
BiomassGOA<-ddply(biomvar,c("YEAR"),summarize,Biomass=sum(BIOMASS,na.rm=T)/1000,Variance=sum(VAR2,na.rm=T)/1000000,
                  SD=sqrt(sum(VAR2,na.rm=T))/1000,CV=SD/(sum(BIOMASS,na.rm=T)/1000))



library(ggmap)
library(plyr)
library(ggplot2)
library(reshape2)
library(sp)
library(grid)
#library(geoshape)
testmap <- get_googlemap(c(lon=-155,lat=58) ,maptype='satellite',zoom=4, xlim=c(-135,-170), ylim=c(51,63)) # set up mapframe

species<-30060
species_name<-"POP"
c2<-catch2
y2<-(unique(c2$YEAR.x))
y2<-sort(y2[!is.na(y2)])
c2<-c2[c2$SPECIES_CODE==species,]
for(i in 1:length(y2)){
  c3<-c2[c2$YEAR.x==y2[i],]
  c3<-c3[!is.na(c3$YEAR.x),]
  c3<-c3[!is.na(c3$START_LONGITUDE),]
  
  ggmap(testmap,legend="topleft") + geom_point(data=c3,aes(x= START_LONGITUDE, y=START_LATITUDE,size=log(WEIGHT)),alpha=0.62,colour="red")+
    xlab("Longitude") +ylab("Latitude") +
    ggtitle(paste("Gulf of Alaska Survey Biomass of",species_name))+
    geom_text(mapping=aes(label = c(y2[i]), x = -140, y = 51.5),colour="white",size=10) +
    theme_grey(base_size=20) 
   ggsave(paste("GOABiomass_",y2[i],".png",sep=""), dpi=150, width=10, height=12)
}

