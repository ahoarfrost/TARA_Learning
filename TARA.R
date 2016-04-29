require(vegan)
require(scales)
require(dplyr)
require(reshape2)
require(ggplot2)
require(stringr)
require(geosphere)
require(ecodist)
require(gridExtra)
#load map-related R libraries
libs <- c("ggmap", "rgdal", "rgeos", "maptools", "dplyr", "tidyr", "tmap")
lapply(libs, require, character.only = TRUE) 

set.seed(1009)

samples <- read.csv("SampleDescriptionTable.csv",header=TRUE)
samples2 <- samples[,c(1,2,3,5,7,9,10,11,12,13,14,15,16,17)]
colnames(samples2) <- c("sample.label","INSDC.sample.id","INSDC.run.id.s","PANGAEA.sample.id","stn","latN","longE","sample.depth.m","depth.id","size.frac.lower","size.frac.upper","longhurst.biome","marine.region","pelagic.biome")
samples2$sample.label <- as.character(samples2$sample.label)
samples2$INSDC.sample.id <- as.character(samples2$INSDC.sample.id)
samples2$INSDC.run.id.s <- as.character(samples2$INSDC.run.id.s)
samples2$PANGAEA.sample.id <- as.character(samples2$PANGAEA.sample.id)
levels(samples2$depth.id) <- list(DCM="(DCM) deep chlorophyll maximum layer (ENVO:01000326)",DCM.OMZ="(DCM) deep chlorophyll maximum layer (ENVO:01000326) & marine oxygen minimum zone (ENVO:01000065)",MES="(MES) mesopelagic zone (ENVO:00000213)",MES.OMZ="(MES) mesopelagic zone (ENVO:00000213) & marine oxygen minimum zone (ENVO:01000065)",SRF="(SRF) surface water layer (ENVO:00002042)")
levels(samples2$longhurst.biome) <- list(Coastal="Coastal Biome",Polar="Polar Biome",Trades="Trades Biome",Westerlies="Westerlies Biome")
levels(samples2$marine.region) <- list(IO="(IO) Indian Ocean [MRGID:1904]",MS="(MS) Mediterranean Sea [MRGID:1905]",NAO="(NAO) North Atlantic Ocean [MRGID:1912]",NPO="(NPO) North Pacific Ocean [MRGID:1908]",RS="(RS) Red Sea [MRGID:4264]",SAO="(SAO) South Atlantic Ocean [MRGID:1914]",SO="(SO) Southern Ocean [MRGID:1907]",SPO="(SPO) South Pacific Ocean [MRGID:1910]")
levels(samples2$pelagic.biome) <- list(ANTA="(ANTA) Antarctic Province [MRGID:21502]",ARAB="(ARAB) Northwest Arabian Sea Upwelling Province [MRGID:21475]",BENG="(BENG) Benguela Current Coastal Province [MRGID:21470]",CAMR="(CAMR) Central American Coastal Province [MRGID:21494]",CARB="(CARB) Caribbean Province [MRGID:21466]",CHIL="(CHIL) Chile-Peru Current Coastal Province [MRGID:21495]",EAFR="(EAFR) Eastern Africa Coastal Province [MRGID:21473]",FKLD="(FKLD) Southwest Atlantic Shelves Province [MRGID:21469]",GFST="(GFST) Gulf Stream Province [MRGID:21454]",GUIA="(GUIA) Guianas Coastal Province [MRGID:21463]",ISSG="(ISSG) Indian South Subtropical Gyre Province [MRGID:21472]",MEDI="(MEDI) Mediterranean Sea, Black Sea Province [MRGID:21465]",MONS="(MONS) Indian Monsoon Gyres Province [MRGID:21471]",NASTE="(NAST-E) North Atlantic Subtropical Gyral Province [MRGID:21467]",NASTW="(NAST-W) North Atlantic Subtropical Gyral Province [MRGID:21455]",NPST="(NPST) North Pacific Subtropical and Polar Front Provinces [MRGID:21484]",PEOD="(PEOD) Pacific Equatorial Divergence Province [MRGID:21489]",PNEC="(PNEC) North Pacific Equatorial Countercurrent Province [MRGID:21488]",REDS="(REDS) Red Sea, Persian Gulf Province [MRGID:21474]",SATL="(SATL) South Atlantic Gyral Province [MRGID:21459]",SPSG="(SPSG) South Pacific Subtropical Gyre Province, North and South [MRGID:21486]")
write.csv(samples2,"SampleDescriptionTableMunged.csv",row.names=FALSE)

##extract run IDs to download raw runs for mapping
#make sample description table for prokaryotes only
proksamples <- samples2[samples2$size.frac.lower==0.22&(samples2$size.frac.upper==1.60|samples2$size.frac.upper==3.00),]
#remove four depths in-betweeners
proksamples <- proksamples[-c(100,102,104,138),]
#create a vector of runs
runs <- proksamples$INSDC.run.id.s
runs <- unlist(strsplit(runs,split="|",fixed=TRUE))
runs <- data.frame(runs=runs,runssix=substr(runs,start=1,stop=6))
write.csv(runs,file="RunsToDownload.csv",row.names=FALSE)

key <- read.table("SampleAssemblyKey.txt",header=TRUE)
new <- merge(samples2,key,by="INSDC.sample.id",incomparables=NA)
new2 <- new[,c(1,14,2,15,3:13)]
new2$wgs <- tolower(new2$wgs)
new2$wgs.url <- substr(new2$wgs,1,6)
write.csv(new2,"SampleDescriptions.csv",row.names=FALSE)

#subset new2 to be only prokaryote-enriched fractions
new3 <- new2[new[,"size.frac.lower"]==0.22&(new[,"size.frac.upper"]==1.60|new[,"size.frac.upper"]==3.00),]
#subset new3 so only surface, dcm, and meso depths (are four in-betweeners)
new4 <- new3[-c(102,103,105,139),]
write.csv(new4,"AssembliesToDownload.csv",row.names=FALSE)

#load KO annotations
KOannotations <- read.table("TARA_KO_annotations.txt",header=TRUE,row.names=1)
#only look at prokaryote-enriched sites
ProkSites <- new4$sample.label
#convert new4 sample labels to match KOannotations colnames format
ProkSites <- gsub(pattern="TARA_([0-9]+)_([A-Z]+)_([0-9.]+)-([0-9.]+)",replacement="TARA_\\1_\\2_\\3.\\4",ProkSites)
#note lost a couple sites in this. fix that later. 133 sites we're looking at now. 
KOannotationsProk <- KOannotations[,colnames(KOannotations)%in%ProkSites]
#save file
write.csv(KOannotationsProk,"KOannotations_ProkaryotesOnly.csv")


######VISUALIZATIONS######

melt <- melt(KOannotationsProk)
#ggplot(melt,aes(x=value))+geom_bar()
melt$log10 <- log10(melt$value) #base 10 log of values looks better
ggplot(melt,aes(x=log10))+geom_bar()

#make master table with metadata and counts both 
#transpose count table
ko <- read.csv("KOannotations_ProkaryotesOnly.csv",header=TRUE,row.names=1)
tko <- t(ko)
metadata <- read.csv("SampleDescriptionTableMunged.csv",header=TRUE,row.names=1)
#fix rownames to match ko (lose the dashes)
row.names(metadata) <- gsub(pattern="TARA_([0-9]+)_([A-Z]+)_([0-9.]+)-([0-9.]+)",replacement="TARA_\\1_\\2_\\3.\\4",row.names(metadata))
#make stations prettier
metadata$stn <- gsub(pattern="TARA_([0-9]+)",replacement="\\1",metadata$stn)
#subset metadata to only columns want to add to master, and rows also in tko
new.metadata <- metadata[,c("stn","latN","longE","sample.depth.m","depth.id","marine.region","longhurst.biome","pelagic.biome")]
new.metadata <- new.metadata[row.names(new.metadata)%in%row.names(tko),]
#add total gene abundance at each site to master
KOP <- ko[,order(colnames(ko))]
sums <- colSums(KOP)
new.metadata$totalGeneAbun <- sums
#reorder tko to be in same order new.metadata
tko <- tko[order(row.names(tko)),]
#cbind together
master <- cbind(new.metadata,tko)


####bray-curtis dissimilarities####
#calculate dissimilarities
bc <- metaMDS(tko,distance="bray")
png(file="figures/bc_nmds_general.png",width=800,height=650,pointsize=18)
#look at things in general
plot(bc)
#stress is only 0.105
text(1,1,labels=paste("stress=",round(bc$stress,3)))
dev.off()
#the fit looks pretty good; original dissimilarities are represented pretty well in the reduced dimensions
png(file="figures/stressplot_BCnmds.png")
stressplot(bc)
#stress is only 0.105
text(0.1,3,labels=paste("stress=",round(bc$stress,3)))
dev.off()

#plot by lat first
#colors want factor to be (in order of levels of factor)
latid <- function(lat) {
    if(abs(lat) > 66.5) {
        id="polar"
    } else if(abs(lat)>23.5 & abs(lat)<=66.5) {
        id="temperate"
    } else if(abs(lat)>=0 & abs(lat)<=23.5) {
        id="tropics"
    }
    return(id)
}
master$lat.id <- sapply(master$latN,latid)
master$lat.id <- factor(master$lat.id,levels=c("tropics","temperate","polar"))

cols.latid <- c("pink","palegreen2","cornflowerblue")
png(file="figures/BCDistances_LatitudeRegion.png",width=800,height=600,pointsize=18)
plot(bc,type="n",xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
points(bc,pch=16,col=alpha(cols.latid[master$lat.id],1))
legend("topright",legend=levels(master$lat.id),fill=cols.latid,border="white",bg="transparent",bty="n",cex=0.75)
#for (i in 1:length(levels(master$lat.id))) {
#    ordiellipse(bc,groups=master$lat.id,label=TRUE,col=cols.latid[i],show.groups=levels(master$lat.id)[i],draw="polygon",alpha=70)
#}
dev.off()
#are bc distances significantly different comparing across stations
#P=0.02 - surprises me, look like overlap a lot
adonis(tgh~lat.id,data=master)

#by marine region
png(file="figures/BCdistances_marineRegion.png",width=800,height=600,pointsize=18)
cols.region <- c("orange","pink","blue","darkgreen","red","cornflowerblue","black","green")
plot(bc,type="n",xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
points(bc,pch=16,col=alpha(cols.region[master$marine.region],0.7))
legend("topright",legend=levels(master$marine.region),fill=cols.region,border="white",bg="transparent",bty="n",cex=0.75)
#for (i in length(levels(master$marine.region))) {
#    ordiellipse(bc,groups=master$marine.region,label=TRUE,col=cols.region[i],show.groups=levels(master$marine.region)[i],draw="polygon",alpha=70)
#}
dev.off()
adonis(tko~marine.region,data=master)


#by depth id
png(file="figures/BCdistances_depth.png",width=800,height=600,pointsize=18)
cols.depthid <- c("tomato3","tomato1","slateblue4","slateblue3","goldenrod")
plot(bc,type="n",xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
points(bc,pch=16,col=alpha(cols.depthid[master$depth.id],0.7))
legend("topright",legend=levels(master$depth.id),fill=cols.depthid,border="white",bg="transparent",bty="n",cex=0.75)
#have to ignore DCM.OMZ for ellipses because only one point
#master$depth.id2 <- as.character(master$depth.id)
#master$depth.id2[80] <- "DCM"
#master$depth.id2 <- as.factor(master$depth.id2)
#cols.depthid2 <- c("springgreen3","dodgerblue2","dodgerblue1","pink")
#for (i in 1:length(levels(master$depth.id2))) {
#    ordiellipse(bc,groups=master$depth.id2,label=TRUE,col=cols.depthid2[i],show.groups=levels(master$depth.id2)[i],draw="polygon",alpha=70)
#}
dev.off()


#calculate shannon diversity for every station, append to master
H <- diversity(tko)
master$H <- H

###map visualizations###
#load world dataset that has continent polygons for context
data(World)
#create spatial points object containing data associated with spatial points
latlong <- as.matrix(master[,c("longE","latN")])
#subset data from master want to store in spatial dataset 
spdata <- master[,c("stn","sample.depth.m","depth.id","marine.region","longhurst.biome","pelagic.biome","lat.id","depth.id","H","totalGeneAbun")]
#create spatial object with coordinates and associated data; use CRS for global geography
spdata <- SpatialPointsDataFrame(coords=latlong,data=spdata,proj4string=CRS("+init=epsg:4326"))

#map Shannon Diversity, color points by H; this is for all depths
png("figures/map_H.png",width=1200,height=750,pointsize=24)
qtm(World) + tm_shape(spdata) + tm_bubbles(col="H",style="kmeans",alpha=0.5,size=0.2,border.col="grey20")
dev.off()

#same thing, color by TotalGeneAbun
png("figures/map_TotalGeneAbun.png",width=1200,height=750,pointsize=24)
qtm(World) + tm_shape(spdata) + tm_bubbles(col="totalGeneAbun",style="kmeans",alpha=0.5,size=0.2,border.col="grey20")
dev.off()

#looks sort of like negative relationship btwn diversity and totalGeneAbun, but at least overall not the case. 
#Maybe look at specific regions individually like IO
HvsGene <- lm(master$totalGeneAbun~master$H)
png("figures/HvsTotalGeneAbun.png",width=800,height=600,pointsize=18)
ggplot(master,aes(x=totalGeneAbun,y=H)) + geom_point(aes(color=marine.region)) + scale_color_manual(values=cols.region) + annotate("text",x=4e+07,y=5,label=paste("R2=",round(summary(HvsGene)$r.squared,3),"; P=",round(summary(HvsGene)$coef[2,4],2))) + geom_smooth(method="lm")
dev.off()

#map some things depth-specific
masterSRF <- master[master$depth.id=="SRF",]
masterDCM <- master[master$depth.id=="DCM"|master$depth.id=="DCM.OMZ",]
masterMES <- master[master$depth.id=="MES"|master$depth.id=="MES.OMZ",]

srf <- spdata$depth.id=="SRF"
dcm <- spdata$depth.id=="DCM"|spdata$depth.id=="DCM.OMZ"
meso <- spdata$depth.id=="MES"|spdata$depth.id=="MES.OMZ"

#H in surface only
png("figures/map_Surface_H.png",width=1200,height=750,pointsize=24)
qtm(World) + tm_shape(spdata[srf,]) + tm_bubbles(col="H",size=0.2, border.col="grey20",style="fixed",alpha=0.85,breaks=c(3.5,3.85,4.2,4.55,4.9,5.25)) + tm_layout(title="H in Surface")
dev.off()
#H in DCM
png("figures/map_DCM_H.png",width=1200,height=750,pointsize=24)
qtm(World) + tm_shape(spdata[dcm,]) + tm_bubbles(col="H",size=0.2, border.col="grey20",style="fixed",alpha=0.85,breaks=c(3.5,3.85,4.2,4.55,4.9,5.25)) + tm_layout(title="H in DCM")
dev.off()
#H in meso
png("figures/map_Meso_H.png",width=1200,height=750,pointsize=24)
qtm(World) + tm_shape(spdata[meso,]) + tm_bubbles(col="H",size=0.2, border.col="grey20",style="fixed",alpha=0.85,breaks=c(3.5,3.85,4.2,4.55,4.9,5.25)) + tm_layout(title="H in Meso")
dev.off()


#Read in module counts
mod <- read.table("TARA_KeggModule_annotations.txt",header=TRUE,row.names=1)
modref <- read.csv("KeggModule_Hierarchy.csv",fill=TRUE,strip.white=TRUE,header=FALSE)
modref$V1 <- str_trim(modref$V1,side="right")
modref$V3 <- str_trim(modref$V3,side="right")
#remove KO numbers
modref <- modref[modref$V1!="E",]

#which higher-level module each M number corresponds to
modcut <- read.csv("cutModules.csv")
modcut <- modcut[modcut$V1=="D",]
#sort modcut in order of mod
modcut <- modcut[order(modcut$Mnum),]
#remove M numbers from modcut not represented in mod
mod$Blevel <- "NA"
mod$Clevel <- "NA"
modcut$Blevel <- as.character(modcut$Blevel)
modcut$Clevel <- as.character(modcut$Clevel)
for (rownum in 2:nrow(mod)) {
    #find matching annotations
    row <- row.names(mod)[rownum]
    inserts <- modcut[modcut$Mnum==row,c("Blevel","Clevel")]
    if (nrow(inserts)==0) {
        next
    } else {
        #insert into mod
        mod[row,c("Blevel","Clevel")] <- inserts
    }
}
mod$Blevel <- as.factor(mod$Blevel)
mod$Clevel <- as.factor(mod$Clevel)

##relationship geographic distance to BC distance
#geographic distance matrix
geogdist <- distm(latlong)
distances <- as.vector(geogdist)
#make bc distance matrix in square symmetrical matrix form, then vectorize
bcdist <- vegdist(tko,method="bray")
bcmat <- as.matrix(bcdist,method="dist")
bcvect <- as.vector(bcmat)
#lm between distances and bcvect
geog <- cbind(distances,bcvect)
geog <- as.data.frame(geog)
#not great fit
fit.geog <- lm(bcvect~distances)
png("figures/BC_vs_GeographicDistance.png",width=600)
ggplot(geog,aes(x=distances,y=bcvect)) + geom_point(color="darkmagenta",line="black",alpha=0.2) + geom_smooth(method="lm") + theme(axis.text=element_text(size=20)) + ggtitle("Great-Circle Distance vs BC Dissimilarity") + annotate("text",x=1.5e+07,y=0.6,label=paste("R2=",round(summary(fit.geog)$r.squared,3)))
dev.off()
#do for surface, dcm, and meso
latlong.srf <- latlong[grep(pattern="SRF",x=row.names(latlong)),]
geogdist.srf <- distm(latlong.srf)
distances.srf <- as.vector(geogdist.srf)
#remove diagonals of matrices
distances.srf <- distances.srf[which(distances.srf!=0)]
tko.srf <- tko[grep(pattern="SRF",x=row.names(tko)),]
bcdist.srf <- vegdist(tko.srf,method="bray")
bcmat.srf <- as.matrix(bcdist.srf,method="dist")
bcvect.srf <- as.vector(bcmat.srf)
#remove diagonals of matrices
bcvect.srf <- bcvect.srf[which(bcvect.srf!=0)]
geog.srf <- cbind(distances.srf,bcvect.srf)
geog.srf <- as.data.frame(geog.srf)
fit.geog.srf <- lm(bcvect.srf~distances.srf)
srf <- ggplot(geog.srf,aes(x=distances.srf,y=bcvect.srf)) + geom_point(color="goldenrod",line="black",alpha=0.2) + geom_smooth(method="lm") + theme(axis.text=element_text(size=20)) + ggtitle("Great-Circle Distance vs BC Dissimilarity SRF") + annotate("text",x=1.5e+07,y=0.6,label=paste("R2=",round(summary(fit.geog.srf)$r.squared,4)))


latlong.dcm <- latlong[grep(pattern="DCM",x=row.names(latlong)),]
geogdist.dcm <- distm(latlong.dcm)
distances.dcm <- as.vector(geogdist.dcm)
distances.dcm <- distances.dcm[which(distances.dcm!=0)]
tko.dcm <- tko[grep(pattern="DCM",x=row.names(tko)),]
bcdist.dcm <- vegdist(tko.dcm,method="bray")
bcmat.dcm <- as.matrix(bcdist.dcm,method="dist")
bcvect.dcm <- as.vector(bcmat.dcm)
bcvect.dcm <- bcvect.dcm[which(bcvect.dcm!=0)]
geog.dcm <- cbind(distances.dcm,bcvect.dcm)
geog.dcm <- as.data.frame(geog.dcm)
fit.geog.dcm <- lm(bcvect.dcm~distances.dcm)
dcm <- ggplot(geog.dcm,aes(x=distances.dcm,y=bcvect.dcm)) + geom_point(color="tomato3",line="black",alpha=0.2) + geom_smooth(method="lm") + theme(axis.text=element_text(size=20)) + ggtitle("Great-Circle Distance vs BC Dissimilarity dcm") + annotate("text",x=1.5e+07,y=0.6,label=paste("R2=",round(summary(fit.geog.dcm)$r.squared,3)))

latlong.mes <- latlong[grep(pattern="MES",x=row.names(latlong)),]
geogdist.mes <- distm(latlong.mes)
distances.mes <- as.vector(geogdist.mes)
distances.mes <- distances.mes[which(distances.mes!=0)]
tko.mes <- tko[grep(pattern="MES",x=row.names(tko)),]
bcdist.mes <- vegdist(tko.mes,method="bray")
bcmat.mes <- as.matrix(bcdist.mes,method="dist")
bcvect.mes <- as.vector(bcmat.mes)
bcvect.mes <- bcvect.mes[which(bcvect.mes!=0)]
geog.mes <- cbind(distances.mes,bcvect.mes)
geog.mes <- as.data.frame(geog.mes)
fit.geog.mes <- lm(bcvect.mes~distances.mes)
mes <- ggplot(geog.mes,aes(x=distances.mes,y=bcvect.mes)) + geom_point(color="slateblue3",line="black",alpha=0.2) + geom_smooth(method="lm") + theme(axis.text=element_text(size=20)) + ggtitle("Great-Circle Distance vs BC Dissimilarity MES") + annotate("text",x=1.5e+07,y=0.6,label=paste("R2=",round(summary(fit.geog.mes)$r.squared,3)))

#total gene abundance by depth
png("figures/TotalGeneAbun_bydepth_violin.png",width=800,height=600)
ggplot(master,aes(x=depth.id,y=totalGeneAbun)) + geom_violin(aes(fill=depth.id),scale="width") + ggtitle("distribution total gene abundance by depth")
dev.off()

#H by depth violin plot
png("figures/H_bydepth_violin.png",width=800,height=600)
ggplot(master,aes(x=depth.id,y=H)) + geom_violin(aes(fill=depth.id),scale="width") + ggtitle("distribution total gene abundance by depth")
dev.off()

#
ggplot(mod,aes(x=Blevel)) + geom_violin(scale="width") + ggtitle("distribution Kegg modules")

