#download MODIS brdf to VI data for BAED 

# load libraries
library("MODIS")
library(rgdal)
library(raster)


################# section 1: download MODIS NBAR data  ############
MODISoptions(localArcPath="D:\\users\\Zhihua\\MODIS",
             outDirPath="D:\\users\\Zhihua\\MODIS",
             gdalPath='c:/OSGeo4W64/bin')

getProduct() # list available products

dates <- as.POSIXct( as.Date(c("1/1/2000","30/7/2015"),format = "%d/%m/%Y") )
dates2 <- transDate(dates[1],dates[2]) # Transform input dates from before
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "


#download reflectance
runGdal(product="MCD43A4",  #Nadir BRDF-Adjusted Reflectance, 1000 m reso
        begin=dates2$beginDOY,
        end = dates2$endDOY,
        tileH = 23,tileV = 2,
        #tileH = 23:24,tileV = 2:3,
        #SDSstring = "1", #only extract the first layers
        outProj=proj.geo,
        job = "NBAR_larch")

#download quality data
runGdal(product="MCD43A2",  #Nadir BRDF-Adjusted Reflectance, 1000 m reso
        begin=dates2$beginDOY,
        end = dates2$endDOY,
        tileH = 23,tileV = 2,
        #tileH = 23:24,tileV = 2:3,
        #SDSstring = "1", #only extract the first layers
        outProj=proj.geo,
        job = "NBAR_larch_qa")


################# section 2: extract fire parameters and determine disturbance time############



###############################################################################################
#### use bfast methods and all available landsat images to find fire disturbance ##
# need R3.1 to work

#download landsat from glovis
#1 download all the jpg and meta data to select scenes
setwd("D:\\dataset\\XinganTM\\p122024")
scene.meta = list.files(path = ".", pattern = "*.meta")

#read into meta files, and extract cloud cover info
scene.id = c()
cloud.cover = c()
for (i in 1:length(scene.meta)){
  d <- scan(scene.meta[i], what=character() )
  scene.id = c(scene.id, d[3])
  cloud.cover = c(cloud.cover, as.numeric(d[18]))

}

scene.df = data.frame(scene_id = scene.id, cloud_cover = cloud.cover)
scene1.df <- scene.df[which(as.numeric(substr(scene.df$scene_id,10,13)) >= 2000 &   #select year
                             as.numeric(substr(scene.df$scene_id,14,16)) >= 106 &   #select doy
                             as.numeric(substr(scene.df$scene_id,14,16)) <= 290 & 
                             as.numeric(scene.df$cloud_cover) < 70),]

scene1.lc = scene1.df$scene_id[substr(scene1.df$scene_id,1,3) == "LC8"]
scene1.lt = scene1.df$scene_id[substr(scene1.df$scene_id,1,3) == "LT5"]
scene1.le_slc_on = scene1.df$scene_id[substr(scene1.df$scene_id,1,3) == "LE7" & as.numeric(substr(scene1.df$scene_id,10,16)) < 2003152]
scene1.le_slc_off = scene1.df$scene_id[substr(scene1.df$scene_id,1,3) == "LE7" & as.numeric(substr(scene1.df$scene_id,10,16)) > 2003152]

scene1.lc <- factor(scene1.lc)
scene1.lt <- factor(scene1.lt)
scene1.le_slc_on <- factor(scene1.le_slc_on)
scene1.le_slc_off <- factor(scene1.le_slc_off)


sink('p122024_since2000_2.txt') # starting write out into a text files

cat("GloVis Scene List","\n")
cat("sensor=Landsat 8 OLI","\n")
cat("ee_dataset_name=LANDSAT_8","\n")
for (i in 1:length(scene1.lc)){
  cat(levels(scene1.lc)[i],"\n")
}

cat("sensor=L7 SLC-off (2003->)","\n")
cat("ee_dataset_name=LANDSAT_ETM_SLC_OFF","\n")
for (i in 1:length(scene1.le_slc_off)){
  cat(levels(scene1.le_slc_off)[i],"\n")
}

cat("sensor=L7 SLC-on (1999-2003)","\n")
cat("ee_dataset_name=LANDSAT_ETM","\n")
for (i in 1:length(scene1.le_slc_on)){
  cat(levels(scene1.le_slc_on)[i],"\n")
}

cat("sensor=Landsat 4-5 TM","\n")
cat("ee_dataset_name=LANDSAT_TM","\n")
for (i in 1:length(scene1.lt)){
  cat(levels(scene1.lt)[i],"\n")
}

sink() # Stop writing to the file



#using firefox download all app to download all the landsat images


#install packages
install.packages('devtools') 

library(devtools)

install_github('dutri001/bfastSpatial')
install_github('dutri001/bfastSpatial')
install_github('bendv/Rgrowth')

# load the package
library(bfastSpatial)
library(rgdal)
library(Rgrowth)

####################  Data pre-processing  #######################
# list files
Work.Dir <- "D:/users/Zhihua/Landsat/p122024"
setwd(Work.Dir)

testr1 = readOGR("D:\\users\\Zhihua\\Landsat\\p122024\\test_region", "testr3")

#list names
fn = list.files(path = "./images.zipped", pattern = "*.tar.gz")

#create folder for each images, and extract into folder
for (i in 1:length(fn)){
  dir.create(paste("./images.unzipped/", substr(fn[i], 1, 16), sep = ""))
  folder.name = substr(fn[i], 1, 16)
  untar(paste("./images.zipped/", fn[i], sep = ""),
              exdir = paste("./images.unzipped/", folder.name, sep = ""))
  
  print(paste("Finish Extracting", i," of ", length(fn), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}

# get year and day of image acquisitions
Year = substr(fn, 10, 13); DOY = substr(fn, 14, 16)
YearDOY = as.numeric(substr(fn, 10, 16))
YearDOY_sorted = sort(YearDOY)

#rearrange the order of the files
Index = c(); for (i in 1:length(YearDOY_sorted)){Index = c(Index, which(YearDOY_sorted[i]==YearDOY))}
folder.name = substr(fn, 1, 16)    #get folder names
folder.name = folder.name[Index]   #rearrange the order of the files

#save vegetation index file into a list files
evi.list = list()
ndmi.list = list()
scene.name = rep("tmp",length(folder.name)) 


for (i in 1:length(folder.name)){
  path1 = paste("./images.unzipped/", folder.name[i], sep = "")
  
  #cloud, 0 - clear data; 1 - water; 2-shadow, 3-snow; 4-cloud
  cloud.nm = list.files(path = path1, pattern = "*cfmask.tif")[1]
  cloudm = raster(paste(path1, "\\", cloud.nm, sep = ""))
  cloudm = crop(cloudm, testr1)
  
  scene.name[i] <- substr(cloud.nm, 1, 21)
  
  #list NDVI
  evi.nm = list.files(path = path1, pattern = "*evi.tif")[1]
  evi = raster(paste(path1,"\\", evi.nm, sep = ""))
  evi = crop(evi, testr1)
  evi[cloudm > 0] = NA
  evi.list[[i]] <- evi
  
  #list NDMI
  ndmi.nm = list.files(path = path1, pattern = "*ndmi.tif")[1]
  ndmi = raster(paste(path1,"\\", ndmi.nm, sep = ""))
  ndmi = crop(ndmi, testr1)
  ndmi[cloudm > 0] = NA
  ndmi.list[[i]] <- ndmi
  
  print(paste("Finish Listing files ", i," of ", length(folder.name), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}

evi.list = stack(evi.list)
ndmi.list = stack(ndmi.list)

names(evi.list) <- scene.name
names(ndmi.list) <- scene.name

writeRaster(evi.list,"D:\\users\\Zhihua\\Landsat\\p122024\\results\\evi.grd", overwrite=TRUE) #write a raster stack files
writeRaster(ndmi.list,"D:\\users\\Zhihua\\Landsat\\p122024\\results\\ndmi.grd", overwrite=TRUE) #write a raster stack files

############## data inspections ####################
#get scene info
evi = stack("D:\\users\\Zhihua\\Landsat\\p122024\\results\\evi.grd")

s <- getSceneinfo(names(evi))
s$year <- as.numeric(substr(s$date, 1, 4))
hist(s$year, breaks=c(2000:2015), main="p122r24: Scenes per Year", 
     xlab="year", ylab="# of scenes")

#Valid Observations

obs <- countObs(evi)
plot(obs)

# valid observations
obs <- countObs(evi, as.perc=TRUE)
summary(obs)
# % NA per pixel
percNA <- 100 - countObs(evi, as.perc=TRUE)
plot(percNA, main="percent NA per pixel")

#summary
meanVI <- summaryBrick(evi, fun=mean, na.rm=TRUE)
plot(meanVI)

#annual summary
annualMed <- annualSummary(evi, fun=median, na.rm=TRUE)
plot(annualMed)

############### Spatial BFASTMonitor  #######################
#need to plot a layer first to use interactive functions
plot(evi, 25)

bfm <- bfmPixel(evi, start=c(2010, 1), interactive=TRUE)

windows()
plot(bfm$bfm)

plot(evi, which(s$date == "2011-09-05"))
bfm4 <- bfmPixel(evi, start=c(2011, 1), monend=c(2012, 1), 
                 interactive=TRUE,
                 plot=TRUE)


bfm2 <- bfmPixel(evi, cell=targcell, start=c(2009, 1), 
                 formula=response~harmon, order=1, plot=TRUE)

bfm <- bfmSpatial(evi, start=c(2009, 1), order=1)

# extract change raster
change <- raster(bfm, 1)

months <- changeMonth(change)
# set up labels and colourmap for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# extract magn raster
magn <- raster(bfm, 2)
# make a version showing only breakpoing pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")

################use one year monitoring period   #############
bfm.list = list()

for (i in 2009:2014) {
  bfm1 <- bfmSpatial(evi, start=c(i, 1), monend=c(i+1, 1), order=1)
  bfm.list[[i-2008]] <- bfm1
  
  print(paste("Finish monitoring ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}


for(i in 1:length(bfm.list)){
  
  writeRaster(bfm.list[[i]],paste(".\\results\\bfm.", i+2008, ".grd", sep = ""), overwrite=TRUE) #write a raster stack files
  
 print(paste("Finish writing ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
              
}




bfm09 <- bfmSpatial(evi, start=c(2010, 1), monend=c(2011, 1), order=1)
# extract change
change09 <- raster(bfm09, 1)
# extract and rescale magnitude
magn09 <- raster(bfm09, 2) / 10000
# remove all non-breakpoints
magn09[is.na(change09)] <- NA

# extract and rescale magnitude and apply a -500 threshold
magn09thresh <- magn09
magn09thresh[magn09 > -0.05] <- NA
# compare
op <- par(mfrow=c(1, 2))
plot(magn09, main="magnitude")
plot(magn09thresh, main="magnitude < -0.05")


##detect post-disturbance regrowth
# 1 extract a value
evi.1 <- extract(evi , 1139808)
evi.1 <- zoo(as.numeric(evi.1), as.Date(s$date))
plot(evi.1, type = 'b', cex = 0.5)

bts <- bfastts(evi.1, dates = time(evi.1), type = 'irregular')
bfm <- bfastmonitor(bts, start = c(2011, 1), formula = response~harmon, order = 1, plot = TRUE)

reg <- tsreg(evi.1, change = bfm$breakpoint, h = 0.5, plot = TRUE)
reg2 <- tsreg(evi.1, change = bfm$breakpoint, startOffset = "floor", h = 0.5, plot = TRUE)
reg3 <- tsreg(evi.1, change = bfm$breakpoint, startOffset = "floor", h = 0.5, history='all', plot=TRUE)
reg4 <- tsreg(evi.1, change = bfm$breakpoint, startOffset = "floor", h = 0.5, history='all', s=0, plot=TRUE)


#12/9/2015 

# load the package
library(bfastSpatial)
library(rgdal)
library(raster)
library(Rgrowth)

#setwd directories
Work.Dir <- "D:/users/Zhihua/Landsat/p122024"
setwd(Work.Dir)

testr1 = readOGR("D:\\users\\Zhihua\\Landsat\\p122024\\test_region", "testr3")
fn = list.files(path = "./images.zipped", pattern = "*.tar.gz")
#only select the peak season images from Jun 1st[doy = 152] - Oct 1st[doy = 275]
fn = fn[which(as.numeric(substr(fn, 14, 16)) >= 152 & as.numeric(substr(fn, 14, 16)) <= 275)]
folder.name = substr(fn, 1, 16)

# get year and day of image acquisitions, 
Year = substr(fn, 10, 13); DOY = substr(fn, 14, 16)
YearDOY = as.numeric(substr(fn, 10, 16))
YearDOY_sorted = sort(YearDOY)

#rearrange the order of the files
Index = c(); for (i in 1:length(YearDOY_sorted)){Index = c(Index, which(YearDOY_sorted[i]==YearDOY))}
folder.name = substr(fn, 1, 16)    #get folder names
folder.name = folder.name[Index]   #rearrange the order of the files

#get annual best composite using Griffths et th 2013 (IEEE journal)
#step 1: get band 5 values, for stable (NDVI>0.8) for all clear observations

evi.list = list() #store real evi data
evi.list2 = list() #store if it mature forest: evi > 5000
b5.list = list()
cloud.list = list()
scene.name = rep("tmp",length(folder.name)) 

for (i in 1:length(folder.name)){
  path1 = paste("./images.unzipped/", folder.name[i], sep = "")
  
  #cloud, 0 - clear data; 1 - water; 2-shadow, 3-snow; 4-cloud
  cloud.nm = list.files(path = path1, pattern = "*cfmask.tif")[1]
  cloudm = raster(paste(path1, "\\", cloud.nm, sep = ""))
  cloudm = crop(cloudm, testr1)
  cloud.list[[i]] <- cloudm
  
  scene.name[i] <- substr(cloud.nm, 1, 21)
  
  #list EVI
  evi.nm = list.files(path = path1, pattern = "*evi.tif")[1]
  evi = raster(paste(path1,"\\", evi.nm, sep = ""))
  evi = crop(evi, testr1)
  evi[cloudm > 0] = NA
  evi2[] = evi > 2500
  
  evi.list[[i]] <- evi
  evi.list2[[i]] <- evi2
  
  #list B5
  b5.nm = list.files(path = path1, pattern = "*sr_band5.tif")[1]
  b5 = raster(paste(path1,"\\", b5.nm, sep = ""))
  b5 = crop(b5, testr1)
  b5[cloudm > 0] = NA
  b5.list[[i]] <- b5
  
  print(paste("Finish Listing files ", i," of ", length(folder.name), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}

evi.list = stack(evi.list)
b5.list = stack(b5.list)
cloud.list = stack(cloud.list)

names(evi.list) <- scene.name
names(b5.list) <- scene.name
names(cloud.list) <- scene.name

#select stable, mature forest
evi.list2 = stack(evi.list2)
evi.list2.sum = calc(evi.list2, sum, na.rm = TRUE)

#extract the stable, mature forest point and b5 value
Ex.pts = function(x, threshold){ #x is raster, threshold is how to extract data points
  require(raster)
  proj.geo = projection(x)
  x[x < threshold] = NA
  pts = rasterToPoints(x) #get raster coordinate, from left to right, top to bottom
  pts = data.frame(pts)
  pts <- SpatialPoints(coords = cbind(pts$x,pts$y),proj4string = CRS(proj.geo))
  pts.sp = SpatialPoints(coords = pts, proj4string = CRS(proj.geo))
  return(pts.sp)
}

pts.sp = Ex.pts(evi.list2.sum,  threshold = 60)
pts.b5.df = extract(b5.list, pts.sp)

#estimate gussian model, 
#the reflectable were relative constant during the growing season, 
#therefore, do not need this parameters
pts.b5.df.mean = apply(pts.b5.df, 2, mean, na.rm = TRUE)
plot(pts.b5.df.mean, type = "b")
pts.b5.df.mean = data.frame(b5 = pts.b5.df.mean, doy = as.numeric(substr(scene.name, 14,16)))
plot(pts.b5.df.mean$doy, pts.b5.df.mean$b5)
#however, we also want to give a weight to favor selection of Mid Aug ( DOY = 228)
Score_doy = function(x, mu = 214, sigma = 35){exp(-0.5*(((x-mu)/sigma)^2)/(sigma*sqrt(2*3.14159)))}
#Std = sd(as.numeric(DOY)) #Std = 35.09
x = 154:275
plot(x, Score_doy(x,mu = 193, sigma = 34.33))

#estimate the distance functions
Score_could = function(x, d_req= 50, d_min=50){1/(1+exp(-0.2*min(c(x,d_req)) - 0.5*(d_req - d_min)))}
Score_could(x = 20,d_req = 50,  d_min = 30)
x = 1:100     demfile = 
y = c()
for (i in 1:length(x)){y=c(y, Score_could(x[i],d_req = 100,  d_min = 100))}
plot(x, y)

#composite for each year
Year_uni = sort(as.numeric(unique(Year)))
for (i in Year_uni){
  
  score.list = list()
  b1.list = list()
  b2.list = list()
  b3.list = list()
  b4.list = list()
  b5.list = list()
  b7.list = list()
  
  folder.name1 = folder.name[which(as.numeric(substr(folder.name, 10, 13)) == i)]
  for(j in 1:length(folder.name1)){
    path1 = paste("./images.unzipped/", folder.name1[j], sep = "")
    
    #calculate cloud score
    #cloud, 0 - clear data; 1 - water; 2-shadow, 3-snow; 4-cloud
    cloud.nm = list.files(path = path1, pattern = "*cfmask.tif")[1]
    cloudm = raster(paste(path1, "\\", cloud.nm, sep = ""))
    cloudm = crop(cloudm, testr1)
    cloudm = cloudm == 2 | cloudm == 4
    cloudm.dist = gridDistance(cloudm, origin=1)
    cloudm.score = calc(cloudm.dist, Score_could)
    
    #calculate DOY score
    doy1 = as.numeric(substr(cloud.nm, 14,16))
    doy.score = Score_doy(doy1)
    total.score = cloudm.score + doy.score
    total.score[cloudm > 0] = NA
    score.list[[j]] <- total.score
    
    #list band reflectance
    if (substr(cloud.nm,1,3) == "LT5"| substr(cloud.nm,1,3) == "LE7"){ #for LT and LE
    #b1
    b1.nm = list.files(path = path1, pattern = "*sr_band1.tif")[1]
    b1 = raster(paste(path1,"\\", b1.nm, sep = ""))
    b1 = crop(b1, testr1)
    b1[cloudm > 0] = NA
    b1[b1 > 10000] = NA
    b1.list[[j]] <- b1
    
    #b2
    b1.nm = list.files(path = path1, pattern = "*sr_band2.tif")[1]
    b1 = raster(paste(path1,"\\", b1.nm, sep = ""))
    b1 = crop(b1, testr1)
    b1[cloudm > 0] = NA
    b1[b1 > 10000] = NA
    b2.list[[j]] <- b1
    #b3
    b1.nm = list.files(path = path1, pattern = "*sr_band3.tif")[1]
    b1 = raster(paste(path1,"\\", b1.nm, sep = ""))
    b1 = crop(b1, testr1)
    b1[cloudm > 0] = NA
    b1[b1 > 10000] = NA
    b3.list[[j]] <- b1
    #b4
    b1.nm = list.files(path = path1, pattern = "*sr_band4.tif")[1]
    b1 = raster(paste(path1,"\\", b1.nm, sep = ""))
    b1 = crop(b1, testr1)
    b1[cloudm > 0] = NA
    b1[b1 > 10000] = NA
    b4.list[[j]] <- b1
    
    #b5
    b1.nm = list.files(path = path1, pattern = "*sr_band5.tif")[1]
    b1 = raster(paste(path1,"\\", b1.nm, sep = ""))
    b1 = crop(b1, testr1)
    b1[cloudm > 0] = NA
    b1[b1 > 10000] = NA
    b5.list[[j]] <- b1
    
    #b7
    b1.nm = list.files(path = path1, pattern = "*sr_band7.tif")[1]
    b1 = raster(paste(path1,"\\", b1.nm, sep = ""))
    b1 = crop(b1, testr1)
    b1[cloudm > 0] = NA
    b1[b1 > 10000] = NA
    b7.list[[j]] <- b1
    } else { #for LC
      
      b1.nm = list.files(path = path1, pattern = "*sr_band2.tif")[1]
      b1 = raster(paste(path1,"\\", b1.nm, sep = ""))
      b1 = crop(b1, testr1)
      b1[cloudm > 0] = NA
      b1[b1 > 10000] = NA
      b1.list[[j]] <- b1
      
      #b2
      b1.nm = list.files(path = path1, pattern = "*sr_band3.tif")[1]
      b1 = raster(paste(path1,"\\", b1.nm, sep = ""))
      b1 = crop(b1, testr1)
      b1[cloudm > 0] = NA
      b1[b1 > 10000] = NA
      b2.list[[j]] <- b1
      #b3
      b1.nm = list.files(path = path1, pattern = "*sr_band4.tif")[1]
      b1 = raster(paste(path1,"\\", b1.nm, sep = ""))
      b1 = crop(b1, testr1)
      b1[cloudm > 0] = NA
      b1[b1 > 10000] = NA
      b3.list[[j]] <- b1
      #b4
      b1.nm = list.files(path = path1, pattern = "*sr_band5.tif")[1]
      b1 = raster(paste(path1,"\\", b1.nm, sep = ""))
      b1 = crop(b1, testr1)
      b1[cloudm > 0] = NA
      b1[b1 > 10000] = NA
      b4.list[[j]] <- b1
      
      #b5
      b1.nm = list.files(path = path1, pattern = "*sr_band6.tif")[1]
      b1 = raster(paste(path1,"\\", b1.nm, sep = ""))
      b1 = crop(b1, testr1)
      b1[cloudm > 0] = NA
      b1[b1 > 10000] = NA
      b5.list[[j]] <- b1
      
      #b7
      b1.nm = list.files(path = path1, pattern = "*sr_band7.tif")[1]
      b1 = raster(paste(path1,"\\", b1.nm, sep = ""))
      b1 = crop(b1, testr1)
      b1[cloudm > 0] = NA
      b1[b1 > 10000] = NA
      b7.list[[j]] <- b1
            
    } #end of list
    
    print(paste("Finish listing files for Year ", i,", ", j," of ", length(folder.name1), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
        
  }
  
  score.list = stack(score.list)
  b1.list = stack(b1.list)
  b2.list = stack(b2.list)
  b3.list = stack(b3.list)
  b4.list = stack(b4.list)
  b5.list = stack(b5.list)
  b7.list = stack(b7.list)
  
  #compositing dataset
  #find the location of the highest score
  #rm(list = c("evi.list","cloud.list","evi.list2"))
  
  score.list.max = calc(score.list, find.max)

  b1.list <- addLayer(b1.list, score.list.max)
  b1 = calc(b1.list, fun = Rps.value )
  
  b2.list <- addLayer(b2.list, score.list.max)
  b2 = calc(b2.list, fun = Rps.value)
  
  b3.list <- addLayer(b3.list, score.list.max)
  b3 = calc(b3.list, fun = Rps.value)
  
  b4.list <- addLayer(b4.list, score.list.max)
  b4 = calc(b4.list, fun = Rps.value)
  
  b5.list <- addLayer(b5.list, score.list.max)
  b5 = calc(b5.list, fun = Rps.value)
  
  b7.list <- addLayer(b7.list, score.list.max)
  b7 = calc(b7.list, fun = Rps.value)
  
  writeRaster(b1,paste(Work.Dir, "/results/PCB_sr_band1_", i, ".tif", sep = ""),format="GTiff", overwrite=TRUE) #write a raster stack files
  writeRaster(b2,paste(Work.Dir, "/results/PCB_sr_band2_", i, ".tif", sep = ""),format="GTiff", overwrite=TRUE) #write a raster stack files
  writeRaster(b3,paste(Work.Dir, "/results/PCB_sr_band3_", i, ".tif", sep = ""),format="GTiff", overwrite=TRUE) #write a raster stack files
  writeRaster(b4,paste(Work.Dir, "/results/PCB_sr_band4_", i, ".tif", sep = ""),format="GTiff", overwrite=TRUE) #write a raster stack files
  writeRaster(b5,paste(Work.Dir, "/results/PCB_sr_band5_", i, ".tif", sep = ""),format="GTiff", overwrite=TRUE) #write a raster stack files
  writeRaster(b7,paste(Work.Dir, "/results/PCB_sr_band7_", i, ".tif", sep = ""),format="GTiff", overwrite=TRUE) #write a raster stack files
   
  print(paste("Finish compositing Year ", i," at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

#functions to use
find.max = function(x){
  if(length(which(is.na(x))) == length(x)){
    return(NA)
  } else {
    return(which.max(x)[1])
  }
}

#assign pixel values used on best pixels
Rps.value = function(x){ 
  if (is.na(x[length(x)])) {
    return(NA)
    } else {
      return(x[x[length(x)]])
    }
}

####################  Data pre-processing  #######################
#using LandsatLinkr: http://landsatlinkr.jdbcode.com/index.html

#download Landsat data, and copy to right directors
Work.Dir <- "D:/users/Zhihua/Landsat/XinganImages"
setwd(Work.Dir)

scene_id = c("122023","122024","121023","121024","120023","120024", "123023","119024","121025","120025")

fn = list.files(path = "./oli", pattern = "*.tar.gz")

for (i in 1:length(scene_id)){

fn1 = fn[which(substr(fn, 4,9) == scene_id[i])]

file.copy(from = paste("./oli/", fn1, sep = ""),
            to = paste("./oli/wrs2/",scene_id[i], "/targz/", sep = ""))

}

#In R3.2, run LandsatLinkr
download.file("https://github.com/jdbcode/LandsatLinkr/releases/download/0.2.2-beta/LandsatLinkr_0.2.2.zip", "LandsatLinkr")
install.packages("LandsatLinkr", repos=NULL, type="binary")
install.packages(c("raster","SDMTools","doParallel","foreach","ggplot2","MASS","gridExtra","zoo","segmented","plyr","ecp","gdalUtils","rgdal","igraph","hexbin"))

#setup: http://landsatlinkr.jdbcode.com/guide.html#setup
#remove file not in JUN JUL AUG
proj.utm51:EPSG:32651: WGS 84 / UTM zone 51N

proj.utm51 = "+proj=utm +zone=51 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

library(raster)
library(rgdal)
library(LandsatLinkr)
rasterOptions(tmpdir="D:/users/Zhihua/Landsat/XinganImages/TempDir2") 
rasterOptions(tmpdir="D:/dataset/XinganTM/TempDir")

run_landsatlinkr()

PROJ string: 
+proj=utm +zone=51 +ellps=WGS84 +datum=WGS84 +units=m +no_defs
  
#change study area into a binary files: 1: study area; 0:outside study area
library(raster)
library(rgdal)

setwd("D:\\users\\Zhihua\\Landsat\\XinganImages\\boundry")
  
b = readOGR(dsn="D:\\users\\Zhihua\\Landsat\\XinganImages\\boundry",layer="dxal_bj_proj_polygon")
wrs2 = readOGR(dsn="D:\\users\\Zhihua\\Landsat\\XinganImages\\boundry",layer="StudyAreaPart")

scene_id = c("122023","122024","121023","121024","120023","120024", "123023","119024","121025","120025")

for(i in 1:length(scene_id)){
  
t = readGDAL(paste("mask", scene_id[i], sep = ""))
t = raster(t)
t[1,1] = 0
writeRaster(t,
            paste("D:\\users\\Zhihua\\Landsat\\XinganImages\\boundry\\mask", scene_id[i], ".tif", sep = ""),
            format="GTiff", overwrite=TRUE) 
print(paste("Finish compositing Year ", i," at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}


#get the extent of each landsat 


#merging data myself

library(raster)
library(rgdal)

Work.Dir <- "D:/users/Zhihua/Landsat/XinganImages"
setwd(Work.Dir)

output.dir <- "D:/users/Zhihua/Landsat/XinganImages/output/"
scene_id = c("122023","122024","121023","120023","120024", "121025","120025","123023","121024")
#scene_id = rev(scene_id)
scenor = c("tm", "oli")

rasterOptions(tmpdir="D:/users/Zhihua/Landsat/XinganImages/TempDir2")
proj.utm = "+proj=utm +zone=51 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#for (s1 in 1:length(scenor)){
s1 = 1
#get year info
for (i in 1:length(scene_id)){

tm.year = list.dirs(paste('./', scenor[s1],'/wrs2/', scene_id[i], '/images', sep = ""), recursive=FALSE)

#get files info within each year
for (j in 1:length(tm.year)){

tm.year.file = list.files(path = tm.year[j], pattern = "*ledaps.tif$")
basename = substr(tm.year.file, 1,16)

#remove could-continimated values
ledaps.list = list()
#tca.list = list()
#tc.list = list()
ext = c()
for (k in 1:length(tm.year.file)){
  cloud = raster(paste(tm.year[j], "/", basename[k],"_cloudmask.tif", sep = ""))
  ledaps = stack(paste(tm.year[j], "/", basename[k],"_ledaps.tif", sep = ""))
  ledaps[cloud != 1] = NA
  
  #tca = raster(paste(tm.year[j], "/", basename[k],"_tca.tif", sep = ""))
  #tca[cloud != 1] = NA
  
  #tc = stack(paste(tm.year[j], "/", basename[k],"_tc.tif", sep = ""))
  #tc[cloud != 1] = NA
  
  ledaps.list[[k]] <- ledaps
  #tca.list[[k]] <- tca
  #tc.list[[k]] <- tc
  
  ext1 <-  extent(cloud)
  
  ext = rbind(ext, c(ext1@xmin, ext1@xmax,ext1@ymin, ext1@ymax))
  
  print(paste("Listing:::: Science ", scenor[s1], '::', scene_id[i], " for year ", tm.year[j], " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

#cropping
ext = c(max(ext[,1]),min(ext[,2]),max(ext[,3]),min(ext[,4]))

ledaps.list2 = list()
#tca.list2 = list()
#tc.list2 = list()
for (k in 1:length(ledaps.list)){
  ledaps = crop(ledaps.list[[k]], ext)
  #tca = crop(tca.list[[k]], ext)
  #tc = crop(tc.list[[k]], ext)

  extent(ledaps) <- ext
  #extent(tca) <- ext
  #extent(tc) <- ext
  
  ledaps.list2[[k]] <- ledaps
  #tca.list2[[k]] <- tca
  #tc.list2[[k]] <- tc
  
  print(paste("Cropping:::: ", scenor[s1], '::', scene_id[i], " for year ", tm.year[j], " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

#compositing
b1.list = list()
b2.list = list()
b3.list = list()
b4.list = list()
b5.list = list()
b6.list = list()

#tc1.list = list()
#tc2.list = list()
#tc3.list = list()

for (k in 1:length(ledaps.list)){
  b1 = ledaps.list2[[k]][[1]]
  b1[b1 > 10000 | b1 < 0] = NA
  b1.list[[i]] <- b1

  b2 = ledaps.list2[[k]][[2]]
  b2[b2 > 10000 | b2 < 0] = NA
  b2.list[[i]] <- b2

  b3 = ledaps.list2[[k]][[3]]
  b3[b3 > 10000 | b3 < 0] = NA
  b3.list[[i]] <- b3
  
  b4 = ledaps.list2[[k]][[4]]
  b4[b4 > 10000 | b4 < 0] = NA
  b4.list[[i]] <- b4

  b5 = ledaps.list2[[k]][[5]]
  b5[b5 > 10000 | b5 < 0] = NA
  b5.list[[i]] <- b5
  
  b6 = ledaps.list2[[k]][[6]]
  b6[b6 > 10000 | b6 < 0] = NA
  b6.list[[i]] <- b1

  #tc1.list = tc.list2[[k]][[1]]
  #tc2.list = tc.list2[[k]][[2]]
  #tc3.list = tc.list2[[k]][[3]]
  
  print(paste("Compositing:::: ", scenor[s1], '::', scene_id[i], " for year ", tm.year[j], " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

b1.list = stack(b1.list)
b2.list = stack(b2.list)
b3.list = stack(b3.list)
b4.list = stack(b4.list)
b5.list = stack(b5.list)
b6.list = stack(b6.list)

#tc1.list = stack(tc1.list)
#tc2.list = stack(tc2.list)
#tc3.list = stack(tc3.list)
#tca.list2 = stack(tca.list2)

if (nlayers(b1.list) > 1) {

b1.list = calc(b1.list, mean, na.rm = TRUE)
b2.list = calc(b2.list, mean, na.rm = TRUE)
b3.list = calc(b3.list, mean, na.rm = TRUE)
b4.list = calc(b4.list, mean, na.rm = TRUE)
b5.list = calc(b5.list, mean, na.rm = TRUE)
b6.list = calc(b6.list, mean, na.rm = TRUE)

#tc1.list = calc(tc1.list, mean, na.rm = TRUE)
#tc2.list = calc(tc2.list, mean, na.rm = TRUE)
#tc3.list = calc(tc3.list, mean, na.rm = TRUE)
#tca.list2 = calc(tca.list2, mean, na.rm = TRUE)

#stacking:

}

b = stack(b1.list,b2.list,b3.list,b4.list,b5.list,b6.list)
#tc = stack(tc1.list,tc2.list,tc3.list)

writeRaster(b,
            paste(output.dir, scene_id[i], "/",scenor[s1], substr(basename[1], 4,13), "_reflectance_composite.grd", sep = ""), 
            overwrite=TRUE) #write a raster stack files
#writeRaster(tc,paste(output.dir, scene_id[i], "/LT", substr(basename[1], 4,13), "_tc_composite.grd", sep = ""), overwrite=TRUE) #write a raster stack files
#writeRaster(tca,paste(output.dir, scene_id[i], "/LT", substr(basename[1], 4,13), "_tca_composite.tif", sep = ""), format="GTiff",overwrite=TRUE) #write a raster stack files

removeTmpFiles(h=24) # remove temporal files 24 hours before

} # end of year

} # end of scene

#} # end of scenor

#### calculate indices
b2001 = stack("./output/122023/LT1220232001_reflectance_composite.grd")
ndvi2000 = (b2001[[4]]-b2001[[3]])/(b2001[[4]]+b2001[[3]])


###########################################################################
## mosaic the results from Landsatlinkr

setwd("D:/users/Zhihua/Landsat/XinganImages/output")
library(raster)
library(rgdal)
#library(LandsatLinkr)
rasterOptions(tmpdir="D:/users/Zhihua/Landsat/XinganImages/TempDir")

proj.utm = "+proj=utm +zone=51 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

demfiles = c("R:/users/Zhihua/Landsat/XinganImages/output/120023/tca/2015_p120023_tca_composite.bsq",
             "R:/users/Zhihua/Landsat/XinganImages/output/120024/tca/2015_p120024_tca_composite.bsq",
             "R:/users/Zhihua/Landsat/XinganImages/output/120025/tca/2015_p120025_tca_composite.bsq")

#replace the extreme value and 0 with NA 
scene_id = c("122023","122024","121023","120023","120024", "121025","120025","123023","121024")
spectral.indices = c("tca", "tcb","tcg","tcw")

  for (i in 2:length(scene_id)){
    
    i = 9
    
    for (k in 2:length(spectral.indices)){
      
    files = list.files(path = paste('./', scene_id[i], '/',spectral.indices[k], "/", sep = ""), 
                       pattern = "*composite.bsq$")
    
      for (j in 1:length(files)){
      
      file = raster(paste('./', scene_id[i], '/',spectral.indices[k], "/",files[j], sep = ""))
      
      file.v = as.vector(as.matrix(file))
      file.v[which(file.v == 0)] = NA
      file.v1 = quantile(file.v, probs = c(0.0025,0.9975), na.rm = TRUE)
      file[file == 0 | file < file.v1[1] | file > file.v1[2]] = NA
      
      basename = substr(files[j], 1, nchar(files[j])-4)
      outfile = paste('./', scene_id[i], '/',spectral.indices[k], "/",basename, ".tif", sep = "")
      writeRaster(file, outfile, format="GTiff",overwrite=TRUE)

      print(paste("Replacing value:::: ", scene_id[i], " for indice ", spectral.indices[k], " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
      
    } # end of j
    
    removeTmpFiles(h=24) # remove temporal files 24 hours before
    
    } # end of k
  } # end of i



#merging files from Year 2000

setwd("R:/users/Zhihua/Landsat/XinganImages/output")
library(raster)
library(rgdal)
#library(LandsatLinkr)
rasterOptions(tmpdir="R:/users/Zhihua/Landsat/XinganImages/TempDir2")

scene_id = c("122023","122024","121023","120023","120024", "121025","120025","123023","121024")
spectral.indices = c("tca", "tcb","tcg","tcw")

output.dir = "R:/users/Zhihua/Landsat/XinganImages/output/merge_files"

  
for (k in 1:length(spectral.indices))
  
  for (p in 2000:2015){
      
  # list files START from here
  demfiles = c()
  
  for (j in 1:length(scene_id)){
       
     demfile = paste("./", scene_id[j], "/",spectral.indices[k],"/", 
                     p, "_p", scene_id[j],"_",spectral.indices[k], "_composite.tif", sep = "")  
  
    if (file.exists(demfile)) {demfiles = c(demfiles, demfile)}
  
 } # end of j
 # list files END from here
 
 #mosaic files start from here
 len = length(demfiles)
 refimg = raster(demfiles[1])
 
 for(i in 1:len){
   if(i == 1){mergeit = "r1"} else {mergeit = paste(mergeit,",r",i, sep="")}
   #open the image as raster and aligns it
   dothis = paste("r",i,"=align(demfiles[",i,"], refimg)", sep="")
   eval(parse(text=dothis))
   if(i == len){mergeit = paste("big = mosaic(",mergeit,",fun=mean,na.rm=TRUE,tolerance=0.5)", sep="")}
 }
 
 #run the merge function
 print("calculating mosiac")
 if(len == 1){big = r1} else {eval(parse(text=mergeit))} #only run merge it if there are multiple files to merge
 
 #write out the DEM mosiac
 print("writing mosiac")
  
 outfile = file.path(output.dir,
                     paste(spectral.indices[k], "_",p, "_merge.tif", sep = ""))
 writeRaster(big, outfile, format="GTiff", overwrite=TRUE)
 
 print(paste("Finishing writing :::: Year ", p, " for indice ", spectral.indices[k], " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
} #end of spectral indices

removeTmpFiles(h=24) # remove temporal files 24 hours before

} #End of year
 
 
#doing some plot: in ROSA
setwd("D:/users/Zhihua/Landsat/XinganImages/output")
library(raster)
library(rgdal)
library(rasterVis)
library(colorRamps)

dxal.sp = readOGR(dsn="D:/users/Zhihua/Landsat/XinganImages/boundry",layer="dxal_bj_proj_polygon")
wrs2.sp = readOGR(dsn="D:/users/Zhihua/Landsat/XinganImages/boundry",layer="wrs2_descending_dxal")
wrs2.sp = spTransform(wrs2.sp, projection(dxal.sp))

spectral.indices = c("tca", "tcb","tcg","tcw")
spectral.indices.names = c("TC Angle", "TC Brightness","TC Greenness","TC Wetness")

for (j in 1:length(spectral.indices)){

  j = 3
tca.files = list()
for (i in 2000:2015){
  tca.file = raster(paste("./merge_files/",spectral.indices[j], "_",i, "_merge.tif", sep = ""))
  #tca.file = crop(tca.file, dxal.sp)
  tca.file = crop(tca.file, dxal.sp[3,])
  tca.files[[i - 1999]] <- tca.file
  
  print(paste("Finishing listing :::: Year ", i, " for ", spectral.indices[j], " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
    
}

tca.files = stack(tca.files)
#rasterVis plot
names(tca.files) <- paste("Year", 2000:2015, sep = "")

color_tc3 = rev(rainbow(99, start=0,end=1))
color_tc32 = rev(blue2green2red(99))

qu.val = quantile(as.vector(as.matrix(tca.files)), prob = c(0.01, 0.99), na.rm = T)

breaks_tc3 <- round(seq(min(minValue(tca.files)),max(maxValue(tca.files)),length.out=100),3)  
legendbrks2_tc3 <- round(seq(min(minValue(tca.files)),max(maxValue(tca.files)),length.out=10),0)

breaks_tc3 <- round(seq(min(qu.val[1]),max(qu.val[2]),length.out=100),3)  
legendbrks2_tc3 <- round(seq(min(qu.val[1]),max(qu.val[2]),length.out=10),0)

png(paste("annual_", spectral.indices[j], ".png", sep = ""),height = 6000, width = 6000, res = 300, units = "px")

levelplot(tca.files, 
          main=paste("Annual ", spectral.indices[j], " Map", sep = ""),
          at= breaks_tc3, margin=FALSE,
          maxpixels = nrow(tca.files)*ncol(tca.files),
          col.regions=color_tc32,
          colorkey= list(labels= list(labels= legendbrks2_tc3,at= legendbrks2_tc3, cex = 1.5), space = "bottom"),
          layout=c(4, 4))+
  latticeExtra::layer(sp.polygons(dxal.sp, col = "black", lwd = 2)) 
#+ latticeExtra::layer(sp.polygons(wrs2.sp, col = "black", lty = 2, lwd = 1))
  

dev.off()

} # end of j


#some statistics: annual numbers of LT, LE, and LC for each scenes

setwd("D:/users/Zhihua/Landsat/XinganImages/")

scenor = c("tm", "oli")
scene_id = c("122023","122024","121023","120023","120024", "121025","120025","123023","121024")

# for le7 and lt5
#for (i in 1:length(scenor)){
  lt5.tmp = data.frame(year = 2000:2015, num = 0)
  le7.tmp = data.frame(year = 2000:2015, num = 0)
  
  for (j in 1:length(scene_id)){
    
    years = list.dirs(paste('./', scenor[i], '/', 'wrs2/', scene_id[j], '/images', sep = ""), recursive=FALSE)
    
    lt5.year.files = c()
    le7.year.files = c()
    
    for (k in 1:length(years)){
      files = list.files(path = years[k], pattern = "*cloudmask.tif$")
      lt5.year.files1 = data.frame(year = as.numeric(substr(years[k], nchar(years[k])-3, nchar(years[k]))),                               
                                   num = length(which(substr(files, 1,3) == "LT5")))

      lt5.year.files = rbind(lt5.year.files,lt5.year.files1)
      
      le7.year.files1 = data.frame(year = as.numeric(substr(years[k], nchar(years[k])-3, nchar(years[k]))),                               
                                   num = length(which(substr(files, 1,3) == "LE7")))
      
      le7.year.files = rbind(le7.year.files,le7.year.files1)
    }
    
    lt5.tmp = merge(lt5.tmp, lt5.year.files, by.x = "year", by.y = "year", all.x = TRUE)
    le7.tmp = merge(le7.tmp, le7.year.files, by.x = "year", by.y = "year", all.x = TRUE)
    
  } #end of Scene
  

lt5.tmp[is.na(lt5.tmp)] = 0; lt5.tmp = lt5.tmp[,-2]; colnames(lt5.tmp) <- c("year", scene_id)
le7.tmp[is.na(le7.tmp)] = 0; le7.tmp = le7.tmp[,-2]; colnames(le7.tmp) <- c("year", scene_id)

write.csv(lt5.tmp, "./output/lt5.num.csv")
write.csv(le7.tmp, "./output/le7.num.csv")

#}

# for LC8
#for (i in 1:length(scenor)){
i = 2
tmp = data.frame(year = 2000:2015, num = 0)

for (j in 1:length(scene_id)){
  
  years = list.dirs(paste('./', scenor[i], '/', 'wrs2/', scene_id[j], '/images', sep = ""), recursive=FALSE)
  
  year.files = c()
  
  for (k in 1:length(years)){
    files = list.files(path = years[k], pattern = "*cloudmask.tif$")
    year.files1 = data.frame(year = as.numeric(substr(years[k], nchar(years[k])-3, nchar(years[k]))),                               
                                 num = length(files))
    
    year.files = rbind(year.files,year.files1)
  }
  
  tmp = merge(tmp, year.files, by.x = "year", by.y = "year", all.x = TRUE)
  
} #end of Scene

tmp[is.na(tmp)] = 0; tmp = tmp[,-2]; colnames(tmp) <- c("year", scene_id)
write.csv(tmp, "./output/lc8.num.csv")

Sce.num = lt5.tmp + le7.tmp + tmp
Sce.num[,1] = lt5.tmp[,1]

Sce.st = apply(Sce.num[,c(2:9)], 1, function(x){c(median(x), min(x), max(x))})   
Year.st = apply(Sce.num[,c(2:9)], 2, function(x){c(median(x), min(x), max(x))})   

write.csv(Sce.num, "./output/image.total.csv")
write.csv(Year.st, "./output/image.scene.sta.csv")
write.csv(Sce.st, "./output/image.year.sta.csv")





len = length(demfiles)
refimg = raster(demfiles[1])

for(i in 1:len){
  if(i == 1){mergeit = "r1"} else {mergeit = paste(mergeit,",r",i, sep="")}
  #open the image as raster and aligns it
  dothis = paste("r",i,"=align(demfiles[",i,"], refimg)", sep="")
  eval(parse(text=dothis))
  if(i == len){mergeit = paste("big = mosaic(",mergeit,",fun=mean,na.rm=TRUE,tolerance=0.5)", sep="")}
}

#run the merge function
print("calculating DEM mosiac")
if(len == 1){big = r1} else {eval(parse(text=mergeit))} #only run merge it if there are multiple files to merge

#write out the DEM mosiac
print("writing DEM mosiac")
outfile = file.path(dir,"dem_mosaic.tif")
writeRaster(big, outfile, format="GTiff", datatype = "INT2S",overwrite=T, bandorder = "BSQ",options=c("COMPRESS=NONE"))

######################################
#get work 
t = stack("D:/users/Zhihua/Landsat/XinganImages/output/120025/tca/p120025_tca_composite_stack.bsq")
t.df = extract(t, 30452030)
plot(as.vector(t.df)[c(4:18)])



#1/4/2016
#unzip files, and organize the
#setwd("D:/users/Zhihua/Landsat/XinganImages2/")

scene_id = c("122023","122024","121023","120023","120024", "121025","120025","123023","121024")

#create folders, if folder NOT exist, create it
mainDir <- "R:/users/Zhihua/Landsat/XinganImages2/"

for (i in 1:length(scene_id)){
  for (j in 1999:2015){
    
    subDir <- paste(mainDir, "Images.unzipped/", scene_id[i], "/", j, sep = "")
    #subDir <- paste(mainDir, "Images.unzipped/", scene_id[i], sep = "")
  
    if (!file.exists(subDir)){
      dir.create(file.path(subDir))
   }
  }
}

#copy file and unzip the file into the folder
mainDir <- "D:/users/Zhihua/Landsat/XinganImages2/"

setwd(mainDir)

files = list.files(path = "./Images.zipped", pattern = "*.tar.gz")

for (i in 4:length(scene_id)){
  for (j in 1999:2015){
    
  file = files[which(substr(files, 4,9) == scene_id[i] & substr(files, 10,13) == as.character(j))]
  
  if (length(file) >= 1){
    for (k in 1:length(file)){
      untar(paste("./Images.zipped/", file[k],sep = ""), 
            exdir = paste("./Images.unzipped/",scene_id[i], "/", j, sep = ""))
  
  } #end of k
  } #end of length
  print(paste("Finishing unzipping :::: Year ", j, " for scene ", scene_id[i], " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
  }
  
}

#replace bad value with na
library(raster)
library(rgdal)
library(rasterVis)
rasterOptions(tmpdir="D:/users/Zhihua/Landsat/XinganImages/TempDir2") 

scene_id = c("122023","122024","121023","120023","120024", "121025","120025","123023","121024")
scene_id = rev(scene_id)
#create folders, if folder NOT exist, create it
mainDir <- "D:/users/Zhihua/Landsat/XinganImages2/"
setwd(mainDir)
for (i in 1:length(scene_id)){
  for (j in 1999:2015){
    
    files = list.files(path = paste("./Images.unzipped/",scene_id[i], "/", j, sep = ""), pattern = "*.cfmask.tif$")
    basenames = substr(files, 1, 21)
    if (length(files) >= 1){
      for (k in 1:length(basenames)){
        cloud = raster(paste("./Images.unzipped/",scene_id[i], "/", j, "/",basenames[k], "_cfmask.tif",sep = ""))
        cloud[cloud > 0] = NA
        cloud[cloud == 0] = 1
        
        evi = raster(paste("./Images.unzipped/",scene_id[i], "/", j, "/",basenames[k], "_sr_evi.tif",sep = ""))
        evi = evi*cloud
        
        nbr = raster(paste("./Images.unzipped/",scene_id[i], "/", j, "/",basenames[k], "_sr_nbr.tif",sep = ""))
        nbr = nbr*cloud
        
        nbr2 = raster(paste("./Images.unzipped/",scene_id[i], "/", j, "/",basenames[k], "_sr_nbr2.tif",sep = ""))
        nbr2 = nbr2*cloud
        
        writeRaster(evi,paste("./Images.unzipped/",scene_id[i], "/", j, "/",basenames[k], "_sr_evi-2.tif",sep = ""), format="GTiff", overwrite=TRUE)
        writeRaster(nbr,paste("./Images.unzipped/",scene_id[i], "/", j, "/",basenames[k], "_sr_nbr-2.tif",sep = ""), format="GTiff", overwrite=TRUE)
        writeRaster(nbr2,paste("./Images.unzipped/",scene_id[i], "/", j, "/",basenames[k], "_sr_nbr2-2.tif",sep = ""), format="GTiff", overwrite=TRUE)
        
        
      }
      
    }
    removeTmpFiles(h=24) # remove temporal files 24 hours before
    print(paste("Finishing replacing :::: Year ", j, " for scene ", scene_id[i], " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
    
  }
}





library(raster)
library(rgdal)
library(gdalUtils)

make_usearea_file_bsq = function(infile, projref){
  tempout = paste(infile,"_temp.tif",sep="")
  template = r = raster(infile)
  r = as.matrix(r)
  nas = which(is.na(r) == T)
  ones = which(r != 0)
  zero = which(r == 0)
  r[nas,zero] = 0
  r[ones] = 1
  r = setValues(template,r)
  projection(r) = set_projection(projref)
  r = as(r, "SpatialGridDataFrame")
  writeGDAL(r, tempout, drivername = "GTiff", type = "Byte", options="INTERLEAVE=BAND")
  
  outfile = paste(substr(infile,1,(nchar(infile)-3)),"bsq",sep="")
  
  s_srs = projection(raster(tempout))
  t_srs = set_projection(projref)
  gdalwarp(srcfile=tempout, dstfile=outfile, s_srs=s_srs, t_srs=s_srs, order=1, tr=c(30,30), r="near", of="ENVI")
  unlink(tempout)
}

################################################################################

infile = "K:/test/composite/useareafile.tif"
projref = "K:/test/mss/wrs1/041029/images/1976/LM10410291976181_archv.tif"

################################################################################

infile = "R:/users/Zhihua/Landsat/XinganImages/output/122024/p122024_usearea_sub.tif"
projref = "R:/users/Zhihua/Landsat/XinganImages\tm/wrs2/122024/images/2000/LE71220242000162_cloudmask.tif"



make_usearea_file_bsq(infile, projref)

.reset_session 
@"R:\users\Zhihua\Landsat\XinganImages\LLR-LandTrendr\output\run_llr_lt_seg_and_fit_batchfile.pro"

mag = stack("R:/users/Zhihua/Landsat/XinganImages/LLR-LandTrendr/output/tcangle/LT_v2.00_tcangle_p122024_mags.bsq")
yrs = stack("R:/users/Zhihua/Landsat/XinganImages/LLR-LandTrendr/output/tcangle/LT_v2.00_tcangle_p122024_vertyrs.bsq")


#########################################################################################
# 2/12/2016 
# processing climate data for each fires
# 

setwd("D:\\users\\Zhihua\\Landsat")

library(raster)
library(rgdal)
library("ncdf")
library("maptools")

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

fire.sp = readOGR(dsn=".\\XinganImages\\boundry",layer="FirePatch")
b.sp = readOGR(dsn=".\\XinganImages\\boundry",layer="dxal_bj_proj_polygon")

fire.geo.sp = spTransform(fire.sp, CRS(proj.geo))
b.geo.sp = spTransform(b.sp, CRS(proj.geo))

ext = extent(b.geo.sp)
ext = c(ext@xmin-1, ext@xmax+1,ext@ymin-1, ext@ymax+1)
treecover = raster(".\\XinganImages\\boundry\\Hansen_GFC2015_treecover2000_60N_120E_dxal.tif")
treecover = crop(treecover, fire.geo.sp)

#treelossyr = raster(".\\XinganImages\\boundry\\Hansen_GFC2015_lossyear_60N_120E_dxal.tif")

#calculate climate anomaly for each fires

#############################################################################
#### for minumum temperature   #########
#### for spring (March April May), summary (Jun, Jul, Aug), winter (Dec, Jan, Feb)
tmn1 = brick(".\\Reanalysis2\\CRU\\cru_ts3.23.1981.1990.tmn.dat.nc")
tmn2 = brick(".\\Reanalysis2\\CRU\\cru_ts3.23.1991.2000.tmn.dat.nc")
tmn3 = brick(".\\Reanalysis2\\CRU\\cru_ts3.23.2001.2010.tmn.dat.nc")
tmn4 = brick(".\\Reanalysis2\\CRU\\cru_ts3.23.2011.2014.tmn.dat.nc")

tmn =stack(tmn1,tmn2,tmn3,tmn4)
tmn = crop(tmn, ext)

rm(list = c('tmn1','tmn2','tmn3','tmn4'))

#long-term means for spring, summary, and winter
lynames = names(tmn)

#spring
idx = which(as.numeric(substr(lynames, 2,5)) >=1981 & 
              as.numeric(substr(lynames, 2,5)) <=2010 & 
              (as.numeric(substr(lynames, 7,8)) == 3 | 
               as.numeric(substr(lynames, 7,8)) == 4 |
               as.numeric(substr(lynames, 7,8)) == 5)) #find the index for each year
  tmn1 <- subset(tmn, idx, drop=FALSE)
  tmn1 = calc(tmn1, mean)

#summer
idx = which(as.numeric(substr(lynames, 2,5)) >=1981 & 
              as.numeric(substr(lynames, 2,5)) <=2010 & 
              (as.numeric(substr(lynames, 7,8)) == 6 | 
                 as.numeric(substr(lynames, 7,8)) == 7 |
                 as.numeric(substr(lynames, 7,8)) == 8)) #find the index for each year
tmn2 <- subset(tmn, idx, drop=FALSE)
tmn2 = calc(tmn2, mean)

#winter
idx = which(as.numeric(substr(lynames, 2,5)) >=1981 & 
              as.numeric(substr(lynames, 2,5)) <=2010 & 
              (as.numeric(substr(lynames, 7,8)) == 12 | 
                 as.numeric(substr(lynames, 7,8)) == 1 |
                 as.numeric(substr(lynames, 7,8)) == 2)) #find the index for each year
tmn3 <- subset(tmn, idx, drop=FALSE)
tmn3 = calc(tmn3, mean)

tmn.mn = stack(tmn1, tmn2,tmn3)

# calculate anomaly
spring.ano = list()
summer.ano = list()
winter.ano = list()


for (i in 1998:2014){
    
  #spring
  idx = which(as.numeric(substr(lynames, 2,5)) == i & 
                  (as.numeric(substr(lynames, 7,8)) == 3 | 
                   as.numeric(substr(lynames, 7,8)) == 4 |
                   as.numeric(substr(lynames, 7,8)) == 5)) #find the index for each year
  tmn2 <- subset(tmn, idx, drop=FALSE)
  tmn2 = calc(tmn2, mean)
  tmn2 = tmn2-tmn.mn[[1]]
  spring.ano[[i-1997]] <- tmn2
  
  #summer
  idx = which(as.numeric(substr(lynames, 2,5)) == i & 
                (as.numeric(substr(lynames, 7,8)) == 6 | 
                   as.numeric(substr(lynames, 7,8)) == 7 |
                   as.numeric(substr(lynames, 7,8)) == 8)) #find the index for each year
  tmn2 <- subset(tmn, idx, drop=FALSE)
  tmn2 = calc(tmn2, mean)
  tmn2 = tmn2-tmn.mn[[2]]
  summer.ano[[i-1997]] <- tmn2
  
  #winter
  idx = c(which(as.numeric(substr(lynames, 2,5)) == i & 
                (as.numeric(substr(lynames, 7,8)) == 1 | 
                   as.numeric(substr(lynames, 7,8)) == 2)),
          which(as.numeric(substr(lynames, 2,5)) == i-1 & 
                  (as.numeric(substr(lynames, 7,8)) == 12)))
                
                
                #find the index for each year
  tmn2 <- subset(tmn, idx, drop=FALSE)
  tmn2 = calc(tmn2, mean)
  tmn2 = tmn2-tmn.mn[[3]]
  winter.ano[[i-1997]] <- tmn2
  
}

spring.ano = stack(spring.ano)
summer.ano = stack(summer.ano)
winter.ano = stack(winter.ano)

names(spring.ano) <- paste("Y", 1998:2014, sep = "")
names(summer.ano) <- paste("Y", 1998:2014, sep = "")
names(winter.ano) <- paste("Y", 1998:2014, sep = "")

writeRaster(spring.ano,paste(".\\Reanalysis2\\CRU\\tmin.spring.ano.grd"),overwrite=TRUE)
writeRaster(summer.ano,paste(".\\Reanalysis2\\CRU\\tmin.summer.ano.grd"),overwrite=TRUE)
writeRaster(winter.ano,paste(".\\Reanalysis2\\CRU\\tmin.winter.ano.grd"),overwrite=TRUE)

#############################################################################
#### for mean temperature   #########
#### for spring (March April May), summary (Jun, Jul, Aug), winter (Dec, Jan, Feb)
tmn1 = brick(".\\Reanalysis2\\CRU\\cru_ts3.23.1981.1990.tmp.dat.nc")
tmn2 = brick(".\\Reanalysis2\\CRU\\cru_ts3.23.1991.2000.tmp.dat.nc")
tmn3 = brick(".\\Reanalysis2\\CRU\\cru_ts3.23.2001.2010.tmp.dat.nc")
tmn4 = brick(".\\Reanalysis2\\CRU\\cru_ts3.23.2011.2014.tmp.dat.nc")


tmn =stack(tmn1,tmn2,tmn3,tmn4)
tmn = crop(tmn, ext)

rm(list = c('tmn1','tmn2','tmn3','tmn4'))

#long-term means for spring, summary, and winter
lynames = names(tmn)

#spring
idx = which(as.numeric(substr(lynames, 2,5)) >=1981 & 
              as.numeric(substr(lynames, 2,5)) <=2010 & 
              (as.numeric(substr(lynames, 7,8)) == 3 | 
                 as.numeric(substr(lynames, 7,8)) == 4 |
                 as.numeric(substr(lynames, 7,8)) == 5)) #find the index for each year
tmn1 <- subset(tmn, idx, drop=FALSE)
tmn1 = calc(tmn1, mean)

#summer
idx = which(as.numeric(substr(lynames, 2,5)) >=1981 & 
              as.numeric(substr(lynames, 2,5)) <=2010 & 
              (as.numeric(substr(lynames, 7,8)) == 6 | 
                 as.numeric(substr(lynames, 7,8)) == 7 |
                 as.numeric(substr(lynames, 7,8)) == 8)) #find the index for each year
tmn2 <- subset(tmn, idx, drop=FALSE)
tmn2 = calc(tmn2, mean)

#winter
idx = which(as.numeric(substr(lynames, 2,5)) >=1981 & 
              as.numeric(substr(lynames, 2,5)) <=2010 & 
              (as.numeric(substr(lynames, 7,8)) == 12 | 
                 as.numeric(substr(lynames, 7,8)) == 1 |
                 as.numeric(substr(lynames, 7,8)) == 2)) #find the index for each year
tmn3 <- subset(tmn, idx, drop=FALSE)
tmn3 = calc(tmn3, mean)

tmn.mn = stack(tmn1, tmn2,tmn3)

# calculate anomaly
spring.ano = list()
summer.ano = list()
winter.ano = list()


for (i in 1998:2014){
  
  #spring
  idx = which(as.numeric(substr(lynames, 2,5)) == i & 
                (as.numeric(substr(lynames, 7,8)) == 3 | 
                   as.numeric(substr(lynames, 7,8)) == 4 |
                   as.numeric(substr(lynames, 7,8)) == 5)) #find the index for each year
  tmn2 <- subset(tmn, idx, drop=FALSE)
  tmn2 = calc(tmn2, mean)
  tmn2 = tmn2-tmn.mn[[1]]
  spring.ano[[i-1997]] <- tmn2
  
  #summer
  idx = which(as.numeric(substr(lynames, 2,5)) == i & 
                (as.numeric(substr(lynames, 7,8)) == 6 | 
                   as.numeric(substr(lynames, 7,8)) == 7 |
                   as.numeric(substr(lynames, 7,8)) == 8)) #find the index for each year
  tmn2 <- subset(tmn, idx, drop=FALSE)
  tmn2 = calc(tmn2, mean)
  tmn2 = tmn2-tmn.mn[[2]]
  summer.ano[[i-1997]] <- tmn2
  
  #winter
  idx = c(which(as.numeric(substr(lynames, 2,5)) == i & 
                  (as.numeric(substr(lynames, 7,8)) == 1 | 
                     as.numeric(substr(lynames, 7,8)) == 2)),
          which(as.numeric(substr(lynames, 2,5)) == i-1 & 
                  (as.numeric(substr(lynames, 7,8)) == 12)))
  
  
  #find the index for each year
  tmn2 <- subset(tmn, idx, drop=FALSE)
  tmn2 = calc(tmn2, mean)
  tmn2 = tmn2-tmn.mn[[3]]
  winter.ano[[i-1997]] <- tmn2
  
}

spring.ano = stack(spring.ano)
summer.ano = stack(summer.ano)
winter.ano = stack(winter.ano)

names(spring.ano) <- paste("Y", 1998:2014, sep = "")
names(summer.ano) <- paste("Y", 1998:2014, sep = "")
names(winter.ano) <- paste("Y", 1998:2014, sep = "")

writeRaster(spring.ano,paste(".\\Reanalysis2\\CRU\\tmean.spring.ano.grd"),overwrite=TRUE)
writeRaster(summer.ano,paste(".\\Reanalysis2\\CRU\\tmean.summer.ano.grd"),overwrite=TRUE)
writeRaster(winter.ano,paste(".\\Reanalysis2\\CRU\\tmean.winter.ano.grd"),overwrite=TRUE)

#############################################################################
#### for mean precipitation   #########
#### for spring (March April May), summary (Jun, Jul, Aug), winter (Dec, Jan, Feb)
tmn1 = brick(".\\Reanalysis2\\CRU\\cru_ts3.23.1981.1990.pre.dat.nc")
tmn2 = brick(".\\Reanalysis2\\CRU\\cru_ts3.23.1991.2000.pre.dat.nc")
tmn3 = brick(".\\Reanalysis2\\CRU\\cru_ts3.23.2001.2010.pre.dat.nc")
tmn4 = brick(".\\Reanalysis2\\CRU\\cru_ts3.23.2011.2014.pre.dat.nc")


tmn =stack(tmn1,tmn2,tmn3,tmn4)
tmn = crop(tmn, ext)

rm(list = c('tmn1','tmn2','tmn3','tmn4'))

#long-term means for spring, summary, and winter
lynames = names(tmn)

#spring
idx = which(as.numeric(substr(lynames, 2,5)) >=1981 & 
              as.numeric(substr(lynames, 2,5)) <=2010 & 
              (as.numeric(substr(lynames, 7,8)) == 3 | 
                 as.numeric(substr(lynames, 7,8)) == 4 |
                 as.numeric(substr(lynames, 7,8)) == 5)) #find the index for each year
tmn1 <- subset(tmn, idx, drop=FALSE)
tmn1 = calc(tmn1, mean)

#summer
idx = which(as.numeric(substr(lynames, 2,5)) >=1981 & 
              as.numeric(substr(lynames, 2,5)) <=2010 & 
              (as.numeric(substr(lynames, 7,8)) == 6 | 
                 as.numeric(substr(lynames, 7,8)) == 7 |
                 as.numeric(substr(lynames, 7,8)) == 8)) #find the index for each year
tmn2 <- subset(tmn, idx, drop=FALSE)
tmn2 = calc(tmn2, mean)

#winter
idx = which(as.numeric(substr(lynames, 2,5)) >=1981 & 
              as.numeric(substr(lynames, 2,5)) <=2010 & 
              (as.numeric(substr(lynames, 7,8)) == 12 | 
                 as.numeric(substr(lynames, 7,8)) == 1 |
                 as.numeric(substr(lynames, 7,8)) == 2)) #find the index for each year
tmn3 <- subset(tmn, idx, drop=FALSE)
tmn3 = calc(tmn3, mean)

tmn.mn = stack(tmn1, tmn2,tmn3)

# calculate anomaly
spring.ano = list()
summer.ano = list()
winter.ano = list()


for (i in 1998:2014){
  
  #spring
  idx = which(as.numeric(substr(lynames, 2,5)) == i & 
                (as.numeric(substr(lynames, 7,8)) == 3 | 
                   as.numeric(substr(lynames, 7,8)) == 4 |
                   as.numeric(substr(lynames, 7,8)) == 5)) #find the index for each year
  tmn2 <- subset(tmn, idx, drop=FALSE)
  tmn2 = calc(tmn2, mean)
  tmn2 = tmn2-tmn.mn[[1]]
  spring.ano[[i-1997]] <- tmn2
  
  #summer
  idx = which(as.numeric(substr(lynames, 2,5)) == i & 
                (as.numeric(substr(lynames, 7,8)) == 6 | 
                   as.numeric(substr(lynames, 7,8)) == 7 |
                   as.numeric(substr(lynames, 7,8)) == 8)) #find the index for each year
  tmn2 <- subset(tmn, idx, drop=FALSE)
  tmn2 = calc(tmn2, mean)
  tmn2 = tmn2-tmn.mn[[2]]
  summer.ano[[i-1997]] <- tmn2
  
  #winter
  idx = c(which(as.numeric(substr(lynames, 2,5)) == i & 
                  (as.numeric(substr(lynames, 7,8)) == 1 | 
                     as.numeric(substr(lynames, 7,8)) == 2)),
          which(as.numeric(substr(lynames, 2,5)) == i-1 & 
                  (as.numeric(substr(lynames, 7,8)) == 12)))
  
  
  #find the index for each year
  tmn2 <- subset(tmn, idx, drop=FALSE)
  tmn2 = calc(tmn2, mean)
  tmn2 = tmn2-tmn.mn[[3]]
  winter.ano[[i-1997]] <- tmn2
  
}

spring.ano = stack(spring.ano)
summer.ano = stack(summer.ano)
winter.ano = stack(winter.ano)

names(spring.ano) <- paste("Y", 1998:2014, sep = "")
names(summer.ano) <- paste("Y", 1998:2014, sep = "")
names(winter.ano) <- paste("Y", 1998:2014, sep = "")

writeRaster(spring.ano,paste(".\\Reanalysis2\\CRU\\prec.spring.ano.grd"),overwrite=TRUE)
writeRaster(summer.ano,paste(".\\Reanalysis2\\CRU\\prec.summer.ano.grd"),overwrite=TRUE)
writeRaster(winter.ano,paste(".\\Reanalysis2\\CRU\\prec.winter.ano.grd"),overwrite=TRUE)


###########################################################################
#calculate relationship among forest density and VIs

setwd("D:\\users\\Zhihua\\Landsat")

library(raster)
library(rgdal)
library("ncdf")
library("maptools")

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
proj.utm = "+proj=utm +zone=51 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

fire.sp = readOGR(dsn=".\\XinganImages\\boundry",layer="FirePatch")
b.sp = readOGR(dsn=".\\XinganImages\\boundry",layer="dxal_bj_proj_polygon")

field.df = read.csv(".\\XinganImages2\\chuli0828.csv")
##change to spatial points
field1.df = field.df[which(field.df$X < 10000),]
coordinates(field1.df) <- ~X+Y
proj4string(field1.df) <- proj.geo
field1.df = spTransform(field1.df, CRS(proj.utm))

field2.df = field.df[-which(field.df$X < 10000),]
coordinates(field2.df) <- ~X+Y
proj4string(field2.df) <- proj.utm

field.sp = rbind(field1.df, field2.df)

plot(fire.sp, add = T)
plot(field.sp, pch = 1, col = "red", add = T)

table(field.sp@data$Year)

#in landsat 122024 [in the center] and 121024[in the left side of scene]
rasterOptions(tmpdir="D:/users/Zhihua/Landsat/XinganImages/TempDir2")

#get files info for 2011 and 2012
year.file.2011 = list.files(path = ".\\XinganImages2\\Images.unzipped\\122024\\2011", pattern = "*cfmask.tif$")
basename.2011 = substr(year.file.2011, 1,21)
 
year.file.2012 = list.files(path = ".\\XinganImages2\\Images.unzipped\\122024\\2012", pattern = "*cfmask.tif$")
basename.2012 = substr(year.file.2012, 1,21)

#calculate vegetation indices for 2011

    file.2011 = list()
  
    for (i in 1:length(year.file.2011)){
      cloud = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2011", "\\", basename.2011[i],"_cfmask.tif", sep = ""))
      cloud = crop(cloud, fire.sp[1,])
      
      b1 = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2011", "\\", basename.2011[i],"_sr_band1.tif", sep = ""))
      b1 = crop(b1, fire.sp[1,])
      
      b2 = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2011", "\\", basename.2011[i],"_sr_band2.tif", sep = ""))
      b2 = crop(b2, fire.sp[1,])
      
      b3 = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2011", "\\", basename.2011[i],"_sr_band3.tif", sep = ""))
      b3 = crop(b3, fire.sp[1,])
      
      b4 = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2011", "\\", basename.2011[i],"_sr_band4.tif", sep = ""))
      b4 = crop(b4, fire.sp[1,])
      
      b5 = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2011", "\\", basename.2011[i],"_sr_band5.tif", sep = ""))
      b5 = crop(b5, fire.sp[1,])
      
      b6 = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2011", "\\", basename.2011[i],"_sr_band7.tif", sep = ""))
      b6 = crop(b6, fire.sp[1,])
      
      stk = stack(b1,b2,b3,b4,b5,b6)
      stk = stk/10000
      stk[cloud != 0] = NA
      
      ndvi = (stk[[4]]-stk[[3]])/(stk[[4]]+stk[[3]])
      evi = 2.5*(stk[[4]]-stk[[3]])/(stk[[4]]+2.4*stk[[3]]+1)
      savi = 1.5*(stk[[4]]-stk[[3]])/(stk[[4]]+stk[[3]]+0.5)
      ndwi = (stk[[4]]-stk[[5]])/(stk[[4]]+stk[[5]])
      nbr = (stk[[4]]-stk[[6]])/(stk[[4]]+stk[[6]])
      nbr2 = (stk[[5]]-stk[[6]])/(stk[[5]]+stk[[6]])
      
      #calculate TC index
      if(substr(basename.2011[i],1,2) == "LE"){
        brightness = 0.3561*stk[[1]]+0.3972*stk[[2]]+0.3904*stk[[3]]+0.6966*stk[[4]]+0.2286*stk[[5]]+0.1596*stk[[6]]  #brightness
        greenness = -0.3344*stk[[1]]-0.3544*stk[[2]]-0.4556*stk[[3]]+0.6966*stk[[4]]-0.0242*stk[[5]]-0.2630*stk[[6]]  #greenness
        wetness = 0.2626*stk[[1]]+0.2141*stk[[2]]+0.0926*stk[[3]]+0.0656*stk[[4]]-0.7629*stk[[5]]-0.5388*stk[[6]]  #wetness
      } else if (substr(basename.2011[i],1,2) == "LT"){
        brightness = 0.3037*stk[[1]]+0.2793*stk[[2]]+0.4343*stk[[3]]+0.5585*stk[[4]]+0.5082*stk[[5]]+0.1863*stk[[6]]  #brightness
        greenness = -0.2848*stk[[1]]-0.2435*stk[[2]]-0.5436*stk[[3]]+0.7246*stk[[4]]+0.0840*stk[[5]]-0.18*stk[[6]]  #greenness
        wetness = 0.1509*stk[[1]]+0.1793*stk[[2]]+0.3299*stk[[3]]+0.3406*stk[[4]]-0.7112*stk[[5]]-0.4572*stk[[6]]  #wetness
      
      } 
      
      tca = greenness/brightness
      
      file.2011.st = stack(stk, ndvi, evi,savi, ndwi, nbr, nbr2, brightness, greenness, wetness, tca)
      names(file.2011.st) <- c("b1","b2","b3","b4","b5","b7","ndvi", "evi","savi", "ndwi", "nbr", "nbr2", "tcb", "tcg", "tcw", "tca")
      
      file.2011[[i]] <- file.2011.st
      
    }
    
#calculate vegetation indices for 2012

file.2012 = list()

for (i in 1:length(year.file.2012)){
  cloud = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2012", "\\", basename.2012[i],"_cfmask.tif", sep = ""))
  cloud = crop(cloud, fire.sp[1,])
  
  b1 = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2012", "\\", basename.2012[i],"_sr_band1.tif", sep = ""))
  b1 = crop(b1, fire.sp[1,])
  
  b2 = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2012", "\\", basename.2012[i],"_sr_band2.tif", sep = ""))
  b2 = crop(b2, fire.sp[1,])
  
  b3 = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2012", "\\", basename.2012[i],"_sr_band3.tif", sep = ""))
  b3 = crop(b3, fire.sp[1,])
  
  b4 = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2012", "\\", basename.2012[i],"_sr_band4.tif", sep = ""))
  b4 = crop(b4, fire.sp[1,])
  
  b5 = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2012", "\\", basename.2012[i],"_sr_band5.tif", sep = ""))
  b5 = crop(b5, fire.sp[1,])
  
  b6 = raster(paste(".\\XinganImages2\\Images.unzipped\\122024\\2012", "\\", basename.2012[i],"_sr_band7.tif", sep = ""))
  b6 = crop(b6, fire.sp[1,])
  
  stk = stack(b1,b2,b3,b4,b5,b6)
  stk = stk/10000
  stk[cloud != 0] = NA
  
  ndvi = (stk[[4]]-stk[[3]])/(stk[[4]]+stk[[3]])
  evi = 2.5*(stk[[4]]-stk[[3]])/(stk[[4]]+2.4*stk[[3]]+1)
  savi = 1.5*(stk[[4]]-stk[[3]])/(stk[[4]]+stk[[3]]+0.5)
  ndwi = (stk[[4]]-stk[[5]])/(stk[[4]]+stk[[5]])
  nbr = (stk[[4]]-stk[[6]])/(stk[[4]]+stk[[6]])
  nbr2 = (stk[[5]]-stk[[6]])/(stk[[5]]+stk[[6]])
  
  #calculate TC index
  if(substr(basename.2012[i],1,2) == "LE"){
    brightness = 0.3561*stk[[1]]+0.3972*stk[[2]]+0.3904*stk[[3]]+0.6966*stk[[4]]+0.2286*stk[[5]]+0.1596*stk[[6]]  #brightness
    greenness = -0.3344*stk[[1]]-0.3544*stk[[2]]-0.4556*stk[[3]]+0.6966*stk[[4]]-0.0242*stk[[5]]-0.2630*stk[[6]]  #greenness
    wetness = 0.2626*stk[[1]]+0.2141*stk[[2]]+0.0926*stk[[3]]+0.0656*stk[[4]]-0.7629*stk[[5]]-0.5388*stk[[6]]  #wetness
  } else if (substr(basename.2012[i],1,2) == "LT"){
    brightness = 0.3037*stk[[1]]+0.2793*stk[[2]]+0.4343*stk[[3]]+0.5585*stk[[4]]+0.5082*stk[[5]]+0.1863*stk[[6]]  #brightness
    greenness = -0.2848*stk[[1]]-0.2435*stk[[2]]-0.5436*stk[[3]]+0.7246*stk[[4]]+0.0840*stk[[5]]-0.18*stk[[6]]  #greenness
    wetness = 0.1509*stk[[1]]+0.1793*stk[[2]]+0.3299*stk[[3]]+0.3406*stk[[4]]-0.7112*stk[[5]]-0.4572*stk[[6]]  #wetness
    
  } 
  
  tca = greenness/brightness
  
  file.2012.st = stack(stk, ndvi, evi,savi, ndwi, nbr, nbr2, brightness, greenness, wetness, tca)
  names(file.2012.st) <- c("b1","b2","b3","b4","b5","b7","ndvi", "evi","savi", "ndwi", "nbr", "nbr2", "tcb", "tcg", "tcw", "tca")
  
  file.2012[[i]] <- file.2012.st
  
}


#extract values
le.2011.df = data.frame(extract(file.2011[[1]],field.sp))
lt.2011.df = data.frame(extract(file.2011[[2]],field.sp))
le.2012.df = data.frame(extract(file.2012[[1]],field.sp))

field.df2 = field.sp@data[,c("CoverUnderstory", "DensityTotal","UnderstoryANPP","TotalANPP","DBH","TreeANPP")]
field.df2$UnderstoryANPPper = field.df2$UnderstoryANPP/field.df2$TotalANPP
field.df2$aboveANPPper = 1-field.df2$UnderstoryANPPper

field.df2$DensityTotal = log(field.df2$DensityTotal)
field.df2$TotalANPP = log(field.df2$TotalANPP)
field.df2$TreeANPP = log(field.df2$TreeANPP)


library(psych)
png(file = ".\\analysis_results2\\LT2011-field.data.relation.png", width = 4000, height = 4000, units = "px", res = 300)

#pairs.panels(cbind(lt.2011.df[,c("b4","b5","ndvi", "evi", "ndwi", "nbr", "tcb", "tcg", "tcw", "tca")],field.df2))
#pairs.panels(cbind(lt.2011.df[,c("ndvi", "ndwi", "nbr", "tcg", "tcw", "tca")],field.df2[,c("DensityTotal", "TotalANPP", "aboveANPPper")]))

pairs.panels(cbind(lt.2011.df[,c("ndvi", "ndwi", "tcw", "tca")],field.df2[,c("DensityTotal", "TreeANPP")]),
             lm=TRUE,
             labels = c("NDVI", "NDWI", "TCW", "TCA", "Tree Density", "Tree ANPP"),
             cex.labels = 3,
             cex.axis = 3,
             #main = "Correlation Between VIs and Tree Density and ANPP",
             )

dev.off()

png(file = ".\\analysis_results2\\LE2011-field.data.relation.png", width = 4000, height = 4000, units = "px", res = 300)

#pairs.panels(cbind(le.2011.df[,c("b4","b5","ndvi", "evi", "ndwi", "nbr", "tcb", "tcg", "tcw", "tca")],field.df2))

pairs.panels(cbind(le.2011.df[,c("ndvi", "ndwi", "tcw", "tca")],field.df2[,c("DensityTotal", "TreeANPP")]),
             lm=TRUE,
             labels = c("NDVI", "NDWI", "TCW", "TCA", "Tree Density", "Tree ANPP"),
             cex.labels = 3,
             cex.axis = 3,
             #main = "Correlation Between VIs and Tree Density and ANPP",
)

dev.off()

png(file = ".\\analysis_results2\\LE2012-field.data.relation.png", width = 4000, height = 4000, units = "px", res = 300)

#pairs.panels(cbind(le.2012.df[,c("b4","b5","ndvi", "evi", "ndwi", "nbr", "tcb", "tcg", "tcw", "tca")],field.df2))

pairs.panels(cbind(le.2012.df[,c("ndvi", "ndwi", "tcw", "tca")],field.df2[,c("DensityTotal", "TreeANPP")]),
             lm=TRUE,
             labels = c("NDVI", "NDWI", "TCW", "TCA", "Tree Density", "Tree ANPP"),
             cex.labels = 3,
             cex.axis = 3,
             #main = "Correlation Between VIs and Tree Density and ANPP",
)


dev.off()

write.csv(cbind(le.2011.df[,c("b4","b5","ndvi", "evi", "ndwi", "nbr", "tcb", "tcg", "tcw", "tca")],field.df2),
          ".\\analysis_results2\\LE2011-field.data.relation.csv")

#conclusion: NDWI is a good indicator of forest recovery than any other indices


##############################################################################
## sampling vegetation indices for each fire
plot(gbm.ndwi,
#composite the VI images for each fire

fire.sp = readOGR(dsn=".\\XinganImages\\boundry",layer="FirePatch")

fire.sp@data$Id2 = 1:nrow(fire.sp)
fire.sp@data$area = sapply(slot(fire.sp, "polygons"), slot, "area")


for (i in 1:nrow(fire.sp)){
  
  #go to path/rows, and list the directory
  dirs = list.dirs(paste("./XinganImages2/Images.unzipped/", fire.sp[i,]$wrspr, sep = ""), recursive=FALSE)
  #dirs = substr(dirs, nchar(dirs)-3, nchar(dirs))
  #generate some data point
  #Rpts1 = spsample(fire.sp[i,], n=ceiling(fire.sp[i,]$area/(30*30*81)), type='regular')
  
  for(j in 1:length(dirs)){ #for each year
    
    year.file = list.files(path = dirs[j], pattern = "*cfmask.tif$")
    
     if (length(year.file) > 0) { #must have files under the folder
    
    basename = substr(year.file, 1,21)
    
    #change projection if they are not the same
    proj = raster(paste(dirs[j], "\\", basename[1],"_cfmask.tif", sep = ""))
    if (projection(fire.sp) != projection(proj)) {fire.sp = spTransform(fire.sp, projection(proj))}  
        
    #calculate vegetation indices for 2011
    
    files = list()
    
    for (k in 1:length(year.file)){
      
      if(substr(basename[k],1,2) == "LE" | substr(basename[k],1,2) == "LT"){
      
      cloud = raster(paste(dirs[j], "\\", basename[k],"_cfmask.tif", sep = ""))
      cloud = crop(cloud, fire.sp[i,])
      
      b1 = raster(paste(dirs[j], "\\", basename[k],"_sr_band1.tif", sep = ""))
      b1 = crop(b1, fire.sp[i,])
      
      b2 = raster(paste(dirs[j], "\\", basename[k],"_sr_band2.tif", sep = ""))
      b2 = crop(b2, fire.sp[i,])
      
      b3 = raster(paste(dirs[j], "\\", basename[k],"_sr_band3.tif", sep = ""))
      b3 = crop(b3, fire.sp[i,])
      
      b4 = raster(paste(dirs[j], "\\", basename[k],"_sr_band4.tif", sep = ""))
      b4 = crop(b4, fire.sp[i,])
      
      b5 = raster(paste(dirs[j], "\\", basename[k],"_sr_band5.tif", sep = ""))
      b5 = crop(b5, fire.sp[i,])
      
      b6 = raster(paste(dirs[j], "\\", basename[k],"_sr_band7.tif", sep = ""))
      b6 = crop(b6, fire.sp[i,])
      } else if (substr(basename[k],1,2) == "LC"){
      
        cloud = raster(paste(dirs[j], "\\", basename[k],"_cfmask.tif", sep = ""))
        cloud = crop(cloud, fire.sp[i,])
        
        b1 = raster(paste(dirs[j], "\\", basename[k],"_sr_band2.tif", sep = ""))
        b1 = crop(b1, fire.sp[i,])
        
        b2 = raster(paste(dirs[j], "\\", basename[k],"_sr_band3.tif", sep = ""))
        b2 = crop(b2, fire.sp[i,])
        
        b3 = raster(paste(dirs[j], "\\", basename[k],"_sr_band4.tif", sep = ""))
        b3 = crop(b3, fire.sp[i,])
        
        b4 = raster(paste(dirs[j], "\\", basename[k],"_sr_band5.tif", sep = ""))
        b4 = crop(b4, fire.sp[i,])
        
        b5 = raster(paste(dirs[j], "\\", basename[k],"_sr_band6.tif", sep = ""))
        b5 = crop(b5, fire.sp[i,])
        
        b6 = raster(paste(dirs[j], "\\", basename[k],"_sr_band7.tif", sep = ""))
        b6 = crop(b6, fire.sp[i,])
        
      }
      
      stk = stack(b1,b2,b3,b4,b5,b6)
      stk[cloud != 0] = NA
      
      stk = stk/10000
      
      ndvi = (stk[[4]]-stk[[3]])/(stk[[4]]+stk[[3]])
      evi = 2.5*(stk[[4]]-stk[[3]])/(stk[[4]]+2.4*stk[[3]]+1)
      savi = 1.5*(stk[[4]]-stk[[3]])/(stk[[4]]+stk[[3]]+0.5)
      ndwi = (stk[[4]]-stk[[5]])/(stk[[4]]+stk[[5]])
      nbr = (stk[[4]]-stk[[6]])/(stk[[4]]+stk[[6]])
      nbr2 = (stk[[5]]-stk[[6]])/(stk[[5]]+stk[[6]])
      
      #calculate TC index
      if(substr(basename[k],1,2) == "LE"){
        brightness = 0.3561*stk[[1]]+0.3972*stk[[2]]+0.3904*stk[[3]]+0.6966*stk[[4]]+0.2286*stk[[5]]+0.1596*stk[[6]]  #brightness
        greenness = -0.3344*stk[[1]]-0.3544*stk[[2]]-0.4556*stk[[3]]+0.6966*stk[[4]]-0.0242*stk[[5]]-0.2630*stk[[6]]  #greenness
        wetness = 0.2626*stk[[1]]+0.2141*stk[[2]]+0.0926*stk[[3]]+0.0656*stk[[4]]-0.7629*stk[[5]]-0.5388*stk[[6]]  #wetness
      } else if (substr(basename[k],1,2) == "LT"){
        brightness = 0.3037*stk[[1]]+0.2793*stk[[2]]+0.4343*stk[[3]]+0.5585*stk[[4]]+0.5082*stk[[5]]+0.1863*stk[[6]]  #brightness
        greenness = -0.2848*stk[[1]]-0.2435*stk[[2]]-0.5436*stk[[3]]+0.7246*stk[[4]]+0.0840*stk[[5]]-0.18*stk[[6]]  #greenness
        wetness = 0.1509*stk[[1]]+0.1793*stk[[2]]+0.3299*stk[[3]]+0.3406*stk[[4]]-0.7112*stk[[5]]-0.4572*stk[[6]]  #wetness
        
      } else if (substr(basename[k],1,2) == "LC"){
      brightness = 0.3029*stk[[1]]+0.2786*stk[[2]]+0.4733*stk[[3]]+0.5599*stk[[4]]+0.508*stk[[5]]+0.1872*stk[[6]]  #brightness
      greenness = -0.2941*stk[[1]]-0.243*stk[[2]]-0.5424*stk[[3]]+0.7276*stk[[4]]+0.0713*stk[[5]]-0.1608*stk[[6]]  #greenness
      wetness = 0.1511*stk[[1]]+0.1973*stk[[2]]+0.3283*stk[[3]]+0.3407*stk[[4]]-0.7117*stk[[5]]-0.4559*stk[[6]]  #wetness
    }
      
      tca = greenness/brightness
      
      file.2011.st = stack(stk, ndvi, evi,savi, ndwi, nbr, nbr2, brightness, greenness, wetness, tca)
      names(file.2011.st) <- c("b1","b2","b3","b4","b5","b7","ndvi", "evi","savi", "ndwi", "nbr", "nbr2", "tcb", "tcg", "tcw", "tca")
      
      files[[k]] <- file.2011.st
      
    }

    #composite the vegetation indices
    if (length(files)==1){
        writeRaster(files[[1]],
                    paste("./analysis_results/fire_", i,"_", substr(dirs[j], nchar(dirs[j])-3, nchar(dirs[j])),".grd", sep = ""), 
                    overwrite=TRUE) #write a raster stack files
        
    } else {
      b1 = stack(files[[1]][[1]])
      b2 = stack(files[[1]][[2]])
      b3 = stack(files[[1]][[3]])
      b4 = stack(files[[1]][[4]])
      b5 = stack(files[[1]][[5]])
      b6 = stack(files[[1]][[6]])
      ndvi = stack(files[[1]][[7]])
      evi = stack(files[[1]][[8]])
      savi = stack(files[[1]][[9]])
      ndwi = stack(files[[1]][[10]])
      nbr = stack(files[[1]][[11]])
      nbr2 = stack(files[[1]][[12]])
      tcb = stack(files[[1]][[13]])
      tcg = stack(files[[1]][[14]])
      tcw = stack(files[[1]][[15]])
      tca = stack(files[[1]][[16]])
      
      for (k1 in 2:length(files)){
        b1 = addLayer(b1, files[[k1]][[1]])
        b2 = addLayer(b2, files[[k1]][[2]])
        b3 = addLayer(b3, files[[k1]][[3]])
        b4 = addLayer(b4, files[[k1]][[4]])
        b5 = addLayer(b5, files[[k1]][[5]])
        b6 = addLayer(b6, files[[k1]][[6]])
        ndvi = addLayer(ndvi, files[[k1]][[7]])
        evi = addLayer(evi, files[[k1]][[8]])
        savi = addLayer(savi, files[[k1]][[9]])
        ndwi = addLayer(ndwi, files[[k1]][[10]])
        nbr = addLayer(nbr, files[[k1]][[11]])
        nbr2 = addLayer(nbr2, files[[k1]][[12]])
        tcb = addLayer(tcb, files[[k1]][[13]])
        tcg = addLayer(tcg, files[[k1]][[14]])
        tcw = addLayer(tcw, files[[k1]][[15]])
        tca = addLayer(tca, files[[k1]][[16]])
        
      }
      
      b1 = calc(b1, mean, na.rm =TRUE)
      b2 = calc(b2, mean, na.rm =TRUE)
      b3 = calc(b3, mean, na.rm =TRUE)
      b4 = calc(b4, mean, na.rm =TRUE)
      b5 = calc(b5, mean, na.rm =TRUE)
      b6 = calc(b6, mean, na.rm =TRUE)
      ndvi = calc(ndvi, mean, na.rm =TRUE)
      evi = calc(evi, mean, na.rm =TRUE)
      savi = calc(savi, mean, na.rm =TRUE)
      ndwi = calc(ndwi, mean, na.rm =TRUE)
      nbr = calc(nbr, mean, na.rm =TRUE)
      nbr2 = calc(nbr2, mean, na.rm =TRUE)
      tcb = calc(tcb, mean, na.rm =TRUE)
      tcg = calc(tcg, mean, na.rm =TRUE)
      tcw = calc(tcw, mean, na.rm =TRUE)
      tca = calc(tca, mean, na.rm =TRUE)
      
      file = stack(b1,b2,b3,b4,b5,b6, ndvi, evi,savi, ndwi, nbr, nbr2, brightness, greenness, wetness, tca)
      names(file) <- c("b1","b2","b3","b4","b5","b7","ndvi", "evi","savi", "ndwi", "nbr", "nbr2", "tcb", "tcg", "tcw", "tca")
      
      writeRaster(file,
                  paste("./analysis_results/fire_", i,"_", substr(dirs[j], nchar(dirs[j])-3, nchar(dirs[j])),".grd", sep = ""), 
                  overwrite=TRUE) #write a raster stack files
    }
    
    print(paste("Finishing for :::: Year ", dirs[j], " for fire ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
    
     } # END of if (length(year.file)>0)
    
  } #end of J
  
} # end of i
  

###############################################################################################################
# calculate vegetation indice in the peak growing season
# DO NOT composite, only select a image most close to peak growing season
# preference on sensors: 7 > 5 > 8

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
proj.utm = "+proj=utm +zone=51 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj.utm52 = "+proj=utm +zone=52 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

fire.sp = readOGR(dsn=".\\XinganImages\\boundry",layer="FirePatch")

fire.sp@data$Id2 = 1:nrow(fire.sp)
fire.sp@data$area = sapply(slot(fire.sp, "polygons"), slot, "area")

fire.sp.utm52 = spTransform(fire.sp, proj.utm52)
fire.sp.utm51 = fire.sp

for (i in 1:nrow(fire.sp)){
  
  #go to path/rows, and list the directory
  dirs = list.dirs(paste("./XinganImages2/Images.unzipped/", fire.sp[i,]$wrspr, sep = ""), recursive=FALSE)

  for(j in 1:length(dirs)){ #for each year
    
    year.file = list.files(path = dirs[j], pattern = "*cfmask.tif$")
    
    if (length(year.file) > 0) { #must have files under the folder
      
      basename = substr(year.file, 1,21)
      basename.doy = as.numeric(substr(basename, 14, 16))
      basename.sensor = substr(basename, 1, 2)
            
      #change projection if they are not the same
      proj = raster(paste(dirs[j], "\\", basename[1],"_cfmask.tif", sep = ""))
      proj = projection(proj)
      if (proj == proj.utm) {fire.sp = fire.sp.utm51} else {fire.sp = fire.sp.utm52}  
      
      #select a image to calcuate the vegetation indices,
      #using a weighted methods to calculate the preference of images
      basename.doy = as.numeric(substr(basename, 14, 16))
      basename.sensor = substr(basename, 1, 2)
      
      w.df = data.frame(doy = basename.doy, sensor = basename.sensor)
      w.df$doy2 = 1
      w.df$doy2[w.df$doy >= 200 & w.df$doy <= 245] = 4
      w.df$doy2[w.df$doy < 200 | w.df$doy > 245] = 2
      
      w.df$sensor2 = 1
      w.df$sensor2[w.df$sensor == "LE"] = 3
      w.df$sensor2[w.df$sensor == "LT"] = 2
      
      w.df$weight = w.df$doy2 + w.df$sensor2
      
      k = which.max(w.df$weight)[1]
      
      files = list()
      
      #for (k in 1:length(year.file)){
        
        if(substr(basename[k],1,2) == "LE" | substr(basename[k],1,2) == "LT"){
          
          cloud = raster(paste(dirs[j], "\\", basename[k],"_cfmask.tif", sep = ""))
          cloud = crop(cloud, fire.sp[i,])
          
          b1 = raster(paste(dirs[j], "\\", basename[k],"_sr_band1.tif", sep = ""))
          b1 = crop(b1, fire.sp[i,])
          
          b2 = raster(paste(dirs[j], "\\", basename[k],"_sr_band2.tif", sep = ""))
          b2 = crop(b2, fire.sp[i,])
          
          b3 = raster(paste(dirs[j], "\\", basename[k],"_sr_band3.tif", sep = ""))
          b3 = crop(b3, fire.sp[i,])
          
          b4 = raster(paste(dirs[j], "\\", basename[k],"_sr_band4.tif", sep = ""))
          b4 = crop(b4, fire.sp[i,])
          
          b5 = raster(paste(dirs[j], "\\", basename[k],"_sr_band5.tif", sep = ""))
          b5 = crop(b5, fire.sp[i,])
          
          b6 = raster(paste(dirs[j], "\\", basename[k],"_sr_band7.tif", sep = ""))
          b6 = crop(b6, fire.sp[i,])
        } else if (substr(basename[k],1,2) == "LC"){
          
          cloud = raster(paste(dirs[j], "\\", basename[k],"_cfmask.tif", sep = ""))
          cloud = crop(cloud, fire.sp[i,])
          
          b1 = raster(paste(dirs[j], "\\", basename[k],"_sr_band2.tif", sep = ""))
          b1 = crop(b1, fire.sp[i,])
          
          b2 = raster(paste(dirs[j], "\\", basename[k],"_sr_band3.tif", sep = ""))
          b2 = crop(b2, fire.sp[i,])
          
          b3 = raster(paste(dirs[j], "\\", basename[k],"_sr_band4.tif", sep = ""))
          b3 = crop(b3, fire.sp[i,])
          
          b4 = raster(paste(dirs[j], "\\", basename[k],"_sr_band5.tif", sep = ""))
          b4 = crop(b4, fire.sp[i,])
          
          b5 = raster(paste(dirs[j], "\\", basename[k],"_sr_band6.tif", sep = ""))
          b5 = crop(b5, fire.sp[i,])
          
          b6 = raster(paste(dirs[j], "\\", basename[k],"_sr_band7.tif", sep = ""))
          b6 = crop(b6, fire.sp[i,])
          
        }
        
        stk = stack(b1,b2,b3,b4,b5,b6)
        stk[cloud != 0] = NA
        
        stk = stk/10000
        
        ndvi = (stk[[4]]-stk[[3]])/(stk[[4]]+stk[[3]])
        evi = 2.5*(stk[[4]]-stk[[3]])/(stk[[4]]+2.4*stk[[3]]+1)
        savi = 1.5*(stk[[4]]-stk[[3]])/(stk[[4]]+stk[[3]]+0.5)
        ndwi = (stk[[4]]-stk[[5]])/(stk[[4]]+stk[[5]])
        nbr = (stk[[4]]-stk[[6]])/(stk[[4]]+stk[[6]])
        nbr2 = (stk[[5]]-stk[[6]])/(stk[[5]]+stk[[6]])
        
        #calculate TC index
        if(substr(basename[k],1,2) == "LE"){
          brightness = 0.3561*stk[[1]]+0.3972*stk[[2]]+0.3904*stk[[3]]+0.6966*stk[[4]]+0.2286*stk[[5]]+0.1596*stk[[6]]  #brightness
          greenness = -0.3344*stk[[1]]-0.3544*stk[[2]]-0.4556*stk[[3]]+0.6966*stk[[4]]-0.0242*stk[[5]]-0.2630*stk[[6]]  #greenness
          wetness = 0.2626*stk[[1]]+0.2141*stk[[2]]+0.0926*stk[[3]]+0.0656*stk[[4]]-0.7629*stk[[5]]-0.5388*stk[[6]]  #wetness
        } else if (substr(basename[k],1,2) == "LT"){
          brightness = 0.3037*stk[[1]]+0.2793*stk[[2]]+0.4343*stk[[3]]+0.5585*stk[[4]]+0.5082*stk[[5]]+0.1863*stk[[6]]  #brightness
          greenness = -0.2848*stk[[1]]-0.2435*stk[[2]]-0.5436*stk[[3]]+0.7246*stk[[4]]+0.0840*stk[[5]]-0.18*stk[[6]]  #greenness
          wetness = 0.1509*stk[[1]]+0.1793*stk[[2]]+0.3299*stk[[3]]+0.3406*stk[[4]]-0.7112*stk[[5]]-0.4572*stk[[6]]  #wetness
          
        } else if (substr(basename[k],1,2) == "LC"){
          brightness = 0.3029*stk[[1]]+0.2786*stk[[2]]+0.4733*stk[[3]]+0.5599*stk[[4]]+0.508*stk[[5]]+0.1872*stk[[6]]  #brightness
          greenness = -0.2941*stk[[1]]-0.243*stk[[2]]-0.5424*stk[[3]]+0.7276*stk[[4]]+0.0713*stk[[5]]-0.1608*stk[[6]]  #greenness
          wetness = 0.1511*stk[[1]]+0.1973*stk[[2]]+0.3283*stk[[3]]+0.3407*stk[[4]]-0.7117*stk[[5]]-0.4559*stk[[6]]  #wetness
        }
        
        tca = greenness/brightness
        
        file.2011.st = stack(stk, ndvi, evi,savi, ndwi, nbr, nbr2, brightness, greenness, wetness, tca)
        names(file.2011.st) <- c("b1","b2","b3","b4","b5","b7","ndvi", "evi","savi", "ndwi", "nbr", "nbr2", "tcb", "tcg", "tcw", "tca")
        
        files[[1]] <- file.2011.st
        
      #} end of k
      
      #write the vegetation indiecs
      
        writeRaster(files[[1]],
                    paste("./analysis_results/fire_", i,"_", substr(dirs[j], nchar(dirs[j])-3, nchar(dirs[j])),
                          substr(basename[k], 1, 2),substr(basename[k], 14, 16),".grd", sep = ""), 
                    overwrite=TRUE) #write a raster stack files
        
            
      print(paste("Finishing for :::: Year ", dirs[j], " for fire ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
      
    } # END of if (length(year.file)>0)
    
  } #end of J
  
} # end of i

# END OF only using 1 image for each year

#plot image info for each fire

data.img = c()

for (i in 1:nrow(fire.sp)){

  files = list.files(path = "./analysis_results", pattern = paste("fire_", i, "_*", sep = ""))
  files = files[which(nchar(files) > 20)]
  
  data.img1 = cbind(year = substr(files, 9,12), sensor = substr(files, 13,14), doy = substr(files, 15,17))
  
  data.img = rbind(data.img, data.img1)
}

#collapse the duplicate images
data.img = unique(data.img)
data.img = data.frame(data.img)

data.img$year = as.numeric(levels(data.img$year))[data.img$year]
data.img$doy = as.numeric(levels(data.img$doy))[data.img$doy]

data.img1 = data.img[which(data.img$sensor == "LE"), ]
data.img2 = data.img[which(data.img$sensor == "LT"), ]

png(file = ".\\analysis_results2\\images.png", width = 3000, height = 2000, units = "px", res = 300)

plot(data.img1$year, data.img1$doy, pch = 0, col = "red", 
     xlim = c(1999,2015),ylim = c(185,max(data.img$doy)), cex = 2,
     xlab = "", ylab = "", 
     cex.lab = 2, cex.axis = 1.5)

points(data.img2$year, data.img2$doy, pch = 6, col = "blue", 
     xlim = c(1999,2015),ylim = c(185,max(data.img$doy)),cex = 2)

legend("bottomleft",
       #0.8,0.85, 
       #title = "Landsat sensors",
       legend = c("ETM+", "TM"), pch = c(0,6),col = c("red","blue"), cex = 2, , bty = "black")

dev.off()

########################################################################################

#######################################################################################
# sample variable for fire points
 
setwd("D:\\users\\Zhihua\\Landsat")

library(raster)
library(rgdal)
library("ncdf")
library("maptools")

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
proj.utm = "+proj=utm +zone=51 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

fire.sp = readOGR(dsn=".\\XinganImages\\boundry",layer="FirePatch")
fire.sp@data$Id2 = 1:nrow(fire.sp)
fire.sp@data$area = sapply(slot(fire.sp, "polygons"), slot, "area")

b.sp = readOGR(dsn=".\\XinganImages\\boundry",layer="dxal_bj_proj_polygon")

fire.geo.sp = spTransform(fire.sp, CRS(proj.geo))
b.geo.sp = spTransform(b.sp, CRS(proj.geo))

#read into tree cover
treecover = raster(".\\XinganImages\\boundry\\Hansen_GFC2015_treecover2000_60N_120E_dxal.tif")
treecover = crop(treecover, fire.geo.sp)

#read into biophysical factors
dem = raster(".\\DEM\\srtm_mosaic_utm.tif")
slp = raster(".\\DEM\\slp.tif")
rad = readGDAL(".\\DEM\\raditation")
rad = raster(rad)

# read into climate variable
tmean.spring.ano = stack(".\\Reanalysis2\\CRU\\tmean.spring.ano.grd")
tmean.summer.ano = stack(".\\Reanalysis2\\CRU\\tmean.summer.ano.grd")
tmean.winter.ano = stack(".\\Reanalysis2\\CRU\\tmean.winter.ano.grd")

prec.spring.ano = stack(".\\Reanalysis2\\CRU\\prec.spring.ano.grd")
prec.summer.ano = stack(".\\Reanalysis2\\CRU\\prec.summer.ano.grd")
prec.winter.ano = stack(".\\Reanalysis2\\CRU\\prec.winter.ano.grd")

data.df = c()

for (i in 1:nrow(fire.sp)){
  
  Year = fire.sp[i,]$Year
  
  if (Year <= 2010){ #only do fire between 2000-2010
  
  #generate some data point
  Rpts1 = spsample(fire.sp[i,], n=ceiling(fire.sp[i,]$area/(30*30*81)), type='regular')
  
  Rpts1.geo = spTransform(Rpts1, CRS(proj.geo))
  
  #sample tree cover
  tc.df = extract(treecover, Rpts1.geo)
  
  #sample biophysical factors
  dem.df = extract(dem, Rpts1)
  slp.df = extract(slp, Rpts1)
  rad.df = extract(rad, Rpts1)
  
  #calculate fire characteristics
  # disturbance to fire boundry
  dem.grd = crop(dem, extent(fire.sp[i,]) + c(-90,90,-90,90))
  fire.grd = rasterize(fire.sp[i,], dem.grd)
  dem.grd[is.na(fire.grd)] = 1
  dem.grd[fire.grd == 1] = 0
  d2b = gridDistance(dem.grd, origin=1)
  d2b.df = extract(d2b, Rpts1)
  
  # calculate shape complexity 
  # use shape index: Shape index corrects for the size problem of the perimeter-area ratio index (see
  # previous description) by adjusting for a square standard and, as a result, is the
  # simplest and perhaps most straightforward measure of shape complexity. see FRAGSTATS HELP McGarigal et al., (2012)
  require(SDMTools)
     fire.ID = ConnCompLabel(fire.grd)  ##get unique ID
     fire.shpidx = data.frame(PatchStat(fire.ID,cellsize=30,latlon=FALSE))$shape.index
      
  ####################################################
  ##sample climate
  
  cli.names = names(tmean.spring.ano)
  cli.idx = which(substr(cli.names,2,5) == as.character(Year))
  
  #the 0 year climate
  tmean.spring0.df = extract(tmean.spring.ano[[cli.idx]], Rpts1.geo)
  tmean.summer0.df = extract(tmean.summer.ano[[cli.idx]], Rpts1.geo)
  tmean.winter0.df = extract(tmean.winter.ano[[cli.idx]], Rpts1.geo)
  
  prec.spring0.df = extract(prec.spring.ano[[cli.idx]], Rpts1.geo)
  prec.summer0.df = extract(prec.summer.ano[[cli.idx]], Rpts1.geo)
  prec.winter0.df = extract(prec.winter.ano[[cli.idx]], Rpts1.geo)
  
  #the 1 year climate
  tmean.spring1.df = extract(tmean.spring.ano[[cli.idx+1]], Rpts1.geo)
  tmean.summer1.df = extract(tmean.summer.ano[[cli.idx+1]], Rpts1.geo)
  tmean.winter1.df = extract(tmean.winter.ano[[cli.idx+1]], Rpts1.geo)
  
  prec.spring1.df = extract(prec.spring.ano[[cli.idx+1]], Rpts1.geo)
  prec.summer1.df = extract(prec.summer.ano[[cli.idx+1]], Rpts1.geo)
  prec.winter1.df = extract(prec.winter.ano[[cli.idx+1]], Rpts1.geo)
  
  #the 0-4 year average
  tmean.spring.04avg.df = extract(calc(tmean.spring.ano[[c(cli.idx: (cli.idx+4))]], mean, na.rm = TRUE),Rpts1.geo)
  tmean.summer.04avg.df = extract(calc(tmean.summer.ano[[c(cli.idx: (cli.idx+4))]], mean, na.rm = TRUE),Rpts1.geo)
  tmean.winter.04avg.df = extract(calc(tmean.winter.ano[[c(cli.idx: (cli.idx+4))]], mean, na.rm = TRUE),Rpts1.geo)
  
  prec.spring.04avg.df = extract(calc(prec.spring.ano[[c(cli.idx: (cli.idx+4))]], mean, na.rm = TRUE),Rpts1.geo)
  prec.summer.04avg.df = extract(calc(prec.summer.ano[[c(cli.idx: (cli.idx+4))]], mean, na.rm = TRUE),Rpts1.geo)
  prec.winter.04avg.df = extract(calc(prec.winter.ano[[c(cli.idx: (cli.idx+4))]], mean, na.rm = TRUE),Rpts1.geo)
  
  climate.df = data.frame(tmean.spring0 = tmean.spring0.df, tmean.summer0 = tmean.summer0.df,tmean.winter0=tmean.winter0.df,
                          prec.spring0 = prec.spring0.df, prec.summer0 = prec.summer0.df,prec.winter0=prec.winter0.df,
                          tmean.spring1 = tmean.spring1.df, tmean.summer1 = tmean.summer1.df,tmean.winter1=tmean.winter1.df,
                          prec.spring1 = prec.spring1.df, prec.summer1 = prec.summer1.df,prec.winter1=prec.winter1.df,
                          tmean.spring04avg = tmean.spring.04avg.df, tmean.summer04avg = tmean.summer.04avg.df,tmean.winter04avg=tmean.winter.04avg.df,
                          prec.spring04avg = prec.spring.04avg.df, prec.summer04avg = prec.summer.04avg.df,prec.winter04avg=prec.winter.04avg.df)
  
  ####################################################
  ##sample vegetation indices
  # 0 year VI
  vi00 = stack(paste("./analysis_results/fire_", i,"_", Year,".grd", sep = "")) # composite

  filename = paste("fire_", i,"_", Year, sep = "")
  filename = list.files(path = "./analysis_results/", pattern = filename)
  vi00 = stack(paste("./analysis_results/", filename[3], sep = "")) # single images
  
  # -1 year VI
  vi01 = stack(paste("./analysis_results/fire_", i,"_", Year-1,".grd", sep = ""))
  
  filename = paste("fire_", i,"_", Year-1, sep = "")
  filename = list.files(path = "./analysis_results/", pattern = filename)
  vi01 = stack(paste("./analysis_results/", filename[3], sep = "")) # single images
    
  # 1 year VI
  vi10 = stack(paste("./analysis_results/fire_", i,"_", Year+1,".grd", sep = ""))   
  
  filename = paste("fire_", i,"_", Year+1, sep = "")
  filename = list.files(path = "./analysis_results/", pattern = filename)
  vi10 = stack(paste("./analysis_results/", filename[3], sep = "")) # single images
  
  #calculate the fire severity heterogeneity within a 5 by 5 window
  dnbr.grd = vi01$nbr - vi10$nbr
  dnbr.grd.sd <- focal(dnbr.grd, w=matrix(1, nc=5, nr=5), na.rm = TRUE, fun=sd)
  
  dnbrvar.df = extract(dnbr.grd.sd, Rpts1)
  # 5 year VI
  vi50 = stack(paste("./analysis_results/fire_", i,"_", Year+5,".grd", sep = ""))   
  
  filename = paste("fire_", i,"_", Year+5, sep = "")
  filename = list.files(path = "./analysis_results/", pattern = filename)
  vi50 = stack(paste("./analysis_results/", filename[3], sep = "")) # single images
  
  #extract data
  vi00.df = extract(vi00, Rpts1)
  colnames(vi00.df) <- paste(colnames(vi00.df), "curyr", sep = "")
  
  vi01.df = extract(vi01, Rpts1)
  colnames(vi01.df) <- paste(colnames(vi01.df), "preyr", sep = "")
  
  vi10.df = extract(vi10, Rpts1)
  colnames(vi10.df) <- paste(colnames(vi10.df), "oneyr", sep = "")
  
  vi50.df = extract(vi50, Rpts1)
  colnames(vi50.df) <- paste(colnames(vi50.df), "fiveyr", sep = "")
  
  #combine data together
 data.df1 = data.frame(tc = tc.df, dem = dem.df, slp = slp.df, rad = rad.df, year = fire.sp[i,]$Year,
                       d2b = d2b.df, shpidx = fire.shpidx,dnbrvar = dnbrvar.df,
                       climate.df, vi00.df, vi01.df, vi10.df, vi50.df)
  
 data.df = rbind(data.df, data.df1)
 
  } #end of if
 print(paste("Finishing for :::: Fire ", i, " of ", nrow(fire.sp), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
 
} #end of fire i
 
write.csv(data.df, "./analysis_results2/data.csv")
 
 
#do brt analysis

library(gbm)
library(randomForest)
library(dismo)

Variable.list1year = c("tc", "dem", "slp", "rad", "d2b","shpidx","dnbrvar",
                  "tmean.spring0", "tmean.summer0","tmean.winter0","prec.spring0", "prec.summer0","prec.winter0",
                  "tmean.spring1", "tmean.summer1","tmean.winter1","prec.spring1", "prec.summer1","prec.winter1",
                  "tmean.summer04avg","prec.summer04avg",
                  "ndvipreyr","ndwipreyr", "nbrpreyr","ndvicuryr","ndwicuryr","nbrcuryr","ndvioneyr","ndwioneyr","nbroneyr",
                  "ndvifiveyr","ndwifiveyr","nbrfiveyr")

data1 = data.df[ , which(names(data.df) %in% Variable.list1year)]
data1 = data1[complete.cases(data1),]
#data1 = data1[which(data1$tc >= 30), ]
#dNBR=prefireNBR???postfireNBR

data1$dnbr= data1$nbrpreyr - data1$nbrcuryr
data1$rdnbr= data1$dnbr/sqrt(abs(data1$nbrpreyr))

Variable.list.excu = c("tc", "ndvipreyr","ndwipreyr", "nbrpreyr","ndwicuryr","nbrcuryr","ndvioneyr","ndwioneyr","nbroneyr",
                       "ndvifiveyr","nbrfiveyr")
data1 = data1[ , -which(names(data1) %in% Variable.list.excu)]

#fitting gbm models for ndvi
#reorder the column
Variable.order = c("dem", "slp", "rad", 
                   "tmean.spring0", "tmean.summer0","tmean.winter0","prec.spring0", "prec.summer0","prec.winter0",
                   "tmean.spring1", "tmean.summer1","tmean.winter1","prec.spring1", "prec.summer1","prec.winter1",
                   "tmean.summer04avg","prec.summer04avg",
                   "d2b","shpidx",
                   "dnbr", "rdnbr", "dnbrvar","ndvicuryr",
                   "ndwifiveyr")
data1 = data1[Variable.order]

#data1.df1 = data1[sample(nrow(data1), 10000), ] 
#png(file = ".\\XinganImages2\\LT2011-field.data.relation.png", width = 4000, height = 4000, units = "px", res = 300)

pairs.panels(data1[,c("dnbr","rdnbr","ndvicuryr","ndwifiveyr")])
pairs.panels(data1[,c("tmean.spring0", "tmean.summer0","tmean.winter0","prec.spring0", "prec.summer0","prec.winter0",
                      "tmean.spring1", "tmean.summer1","tmean.winter1","prec.spring1", "prec.summer1","prec.winter1")])

pairs.panels(data1[,c( "tmean.summer0","prec.summer0",
                      "tmean.spring1", "tmean.summer1","tmean.winter1","prec.spring1", "prec.summer1","prec.winter1")])

#dev.off()
#remove significantly correlated variables: r > 0.6 
data1 = data1[ , -which(names(data1) %in% c("dnbr","ndvicuryr"))]
data1 = data1[ , -which(names(data1) %in% c("tmean.winter0", "prec.winter0","tmean.spring0","prec.spring0","tmean.spring1", "prec.spring1"))]

#remove some outliers
data1 = data1[-which(data1$rdnbr > 1.9 | data1$rdnbr < 0), ]
data1 = data1[-which(data1$dem > 1300), ]
data1 = data1[-which(data1$dnbrvar > 0.4), ]
data1 = data1[-which(data1$tmean.summer1 > 1 | data1$tmean.summer1 < -0.5), ]
data1 = data1[-which(data1$d2b > 10000), ]

gbm.ndwi=gbm.step(data = data1, gbm.x = 1:(ncol(data1)-1),
                 gbm.y = ncol(data1), 
                 family = "gaussian",
                 tree.complexity = 3, 
                 n.trees = 20, learning.rate = 0.1, bag.fraction = 0.5)

best.iter=gbm.perf(gbm.ndwi,method="OOB")
print(best.iter)

r <- data1$ndwifiveyr - predict.gbm(gbm.ndwi, data1,n.trees = gbm.ndwi$gbm.call$best.trees, type="response")
1-var(r)/var(data1$ndwifiveyr) 

#variation explained: 64.8%

summary(gbm.ndwi,n.trees=best.iter,cBars=10,las=1,cex.axis=1.5,cex.lab=1.6,mar=c(5,8,4,2)+0.1,axes=TRUE,axisnames=TRUE, plotit=TRUE)


plot(gbm.ndwi,
      i.var = 14,
      n.trees = gbm.ndwi$n.trees,
      continuous.resolution = 100,
      return.grid = FALSE,
      type = "link",
      cex.lab = 2,cex.axis = 2,
     lwd = 3,
     xlab = "")




###########################################################################
# plot partial plot within the barplot (only for variable > 3%)
###########################################################################

b = summary(gbm.ndwi)

Variable.list = colnames(data1)[-16]
b.col = c(rep("lightgreen", 3), #Topography 
                     rep("lightcoral", 8), #Climate
                     rep("skyblue", 2), #Fire Size
                     rep("yellow3", 2) # Fier Severity
                     )

b.col = cbind(Variable.list, b.col)

b.col = merge(b, b.col, by.x = "var", by.y = "Variable.list")
rownames(b.col) = b.col$var
b.col = b.col[rownames(b), ]
b.col = as.character(b.col[,3])

Var.mn = c("RdNBR", "DEM", "Shape Index","RdNBR Variability", "Distance to Perimeter", "Precipation Year+1")

png(file = ".\\analysis_results2\\VarImp.png", width = 2000, height = 2000, units = "px", res = 300)

barplot(height=b$rel.inf[1:6],horiz=TRUE,xlab="",
        #xlim = c(, 15),
        #ylim = c(-4.9, 15),
        cex.axis=1.5,
        cex.lab=1.5,
        col = b.col)
for(i in 1:6){text(5, 0.6+1.2*(i-1), labels=Var.mn[i],cex=1.5, srt=0)}

legend("topright",
        #0.8,0.85, 
       legend = c("Topography", "Climate", "Fire Size","Fire Severity" ), 
       fill = c("lightgreen", "lightcoral", "skyblue", "yellow3"), cex = 2, , bty = "n")


dev.off()


Var.mn1 = c("RdNBR", "DEM", "SHP_IDX","RdNBR_Var", "D2B", "PAno_Summer + 1")

png(file = ".\\analysis_results2\\VarImp2.png", width = 2000, height = 2000, units = "px", res = 300)

barplot(height=b$rel.inf[1:6],horiz=TRUE,xlab="",
        #xlim = c(, 15),
        #ylim = c(-4.9, 15),
        cex.axis=1.5,
        cex.lab=1.5,
        col = c("lightcoral", "yellow3", "lightcoral","lightcoral","lightcoral","skyblue"))
#for(i in 1:6){text(5, 0.6+1.2*(i-1), labels=Var.mn1[i],cex=1.5, srt=0)}

text(5, 0.6+1.2*(1-1), labels=Var.mn1[1],cex=1.5, srt=0)
text(4, 0.6+1.2*(2-1), labels=Var.mn1[2],cex=1.5, srt=0)
text(5, 0.6+1.2*(3-1), labels=Var.mn1[3],cex=1.5, srt=0)
text(5.5, 0.6+1.2*(4-1), labels=Var.mn1[4],cex=1.5, srt=0)
text(3, 0.6+1.2*(5-1), labels=Var.mn1[5],cex=1.5, srt=0)
text(8, 0.6+1.2*(6-1), labels=Var.mn1[6],cex=1.5, srt=0)


legend("topright",
       #0.8,0.85, 
       legend = c("Fire", "Topography", "Climate"), 
       fill = c("lightcoral", "yellow3", "skyblue"), cex = 2, , bty = "n")


dev.off()


#partial plot
png(file = ".\\analysis_results2\\Partplot.brt.png", width = 3000, height = 2000, units = "px", res = 300)

gbm.plot(gbm.ndwi, n.plots=6,
         smooth=TRUE,
         write.title = F, 
         show.contrib=T, 
         y.label="",
         plot.layout=c(2, 3),
         cex.lab = 2,cex.axis = 2,
         lwd = 3,
         common.scale=FALSE)

dev.off()

print(pretty.gbm.tree(gbm.ndwi,1))

#produce interactions
gbm.ndwi.int <- gbm.interactions(gbm.ndwi)
gbm.ndwi.int$interactions
gbm.ndwi.int$rank.list

> gbm.ndwi.int$rank.list
var1.index var1.names var2.index var2.names int.size
1          14      rdnbr          1        dem     0.08
2          12        d2b          1        dem     0.08
3          15    dnbrvar         14      rdnbr     0.07
4          15    dnbrvar          1        dem     0.05
5          15    dnbrvar         12        d2b     0.04
6          14      rdnbr         12        d2b     0.03
7          15    dnbrvar          3        rad     0.02
8          15    dnbrvar          2        slp     0.02
9          14      rdnbr         13     shpidx     0.02
10         14      rdnbr          2        slp     0.02
11         13     shpidx          1        dem     0.02

png(file = ".\\analysis_results2\\Partplot.brt.interaction.png", width = 2400, height = 2400, units = "px", res = 300)

par(mfrow=c(2,2),mar=c(0,0,0,0)+0.5,oma=c(0,0,0,0))

gbm.perspec(gbm.ndwi, 14, 1, z.range=c(-0.1,0.4),x.lab = "RdNBR", y.lab = "DEM", z.lab = "")
gbm.perspec(gbm.ndwi, 12, 1, z.range=c(-0.1,0.4),x.lab = "D2B", y.lab = "DEM", z.lab = "")
gbm.perspec(gbm.ndwi, 15, 14, z.range=c(-0.1,0.4),x.lab = "RdNBR_var", y.lab = "RdNBR", z.lab = "")
gbm.perspec(gbm.ndwi, 15, 1, z.range=c(-0.1,0.4),x.lab = "RdNBR_var", y.lab = "DEM", z.lab = "")

dev.off()

#do random forest analysis
library(randomForest)
rf <- randomForest(ndwifiveyr~., data=data1, na.action=na.omit)
rf
summary(rf)
names(rf)
rf$err.rate
rf$importance
rf$ntree
importance(rf) # importance of each predictor 

png(file = ".\\analysis_results2\\VarImp.rf.png", width = 2000, height = 2000, units = "px", res = 300)

varImpPlot(rf)

dev.off()


png(file = ".\\analysis_results2\\Partplot.rf.png", width = 3000, height = 2000, units = "px", res = 300)

par(mfrow=c(2,3),mar=c(4,3,1,0),oma=c(0,3,0,0))

partialPlot(rf, data1,  x.var = "rdnbr", xlab = "RdNBR",cex.lab = 2,cex.axis = 2,lwd = 3,main = "")
partialPlot(rf, data1,  x.var = "dem", xlab = "DEM",cex.lab = 2,cex.axis = 2,lwd = 3,main = "")
partialPlot(rf, data1,  x.var = "dnbrvar", xlab = "RdNBR Variability",cex.lab = 2,cex.axis = 2,lwd = 3,main = "")
partialPlot(rf, data1,  x.var = "shpidx", xlab = "Shape Index",cex.lab = 2,cex.axis = 2,lwd = 3,main = "")
partialPlot(rf, data1,  x.var = "d2b",xlab = "Distance to Boundry",cex.lab = 2,cex.axis = 2,lwd = 3,main = "")
partialPlot(rf, data1,  x.var = "prec.summer1",xlab = "Precipation Year+1",cex.lab = 2,cex.axis = 2,lwd = 3,main = "")

dev.off()

#plot the first tree
to.dendrogram <- function(dfrep,rownum=1,height.increment=0.1){
  
  if(dfrep[rownum,'status'] == -1){
    rval <- list()
    
    attr(rval,"members") <- 1
    attr(rval,"height") <- 0.0
    attr(rval,"label") <- dfrep[rownum,'prediction']
    attr(rval,"leaf") <- TRUE
    
  }else{##note the change "to.dendrogram" and not "to.dendogram"
    left <- to.dendrogram(dfrep,dfrep[rownum,'left daughter'],height.increment)
    right <- to.dendrogram(dfrep,dfrep[rownum,'right daughter'],height.increment)
    rval <- list(left,right)
    
    attr(rval,"members") <- attr(left,"members") + attr(right,"members")
    attr(rval,"height") <- max(attr(left,"height"),attr(right,"height")) + height.increment
    attr(rval,"leaf") <- FALSE
    attr(rval,"edgetext") <- dfrep[rownum,'split var']
  }
  
  class(rval) <- "dendrogram"
  
  return(rval)
}

tree <- getTree(rf,1,labelVar=TRUE)

d <- to.dendrogram(tree,rownum=50)
str(d)
plot(d,center=TRUE,leaflab='none',edgePar=list(t.cex=1,p.col=NA,p.lty=0))

#Plot a sample tree
cf = cforest(ndwifiveyr~., data=data1, controls=cforest_control(mtry=2, mincriterion=0))
plot(ct, main="Conditional Inference Tree for Kyphosis")



# Regression Tree
library(rpart)

# grow tree
rt <- rpart(ndwifiveyr~., 
            data=data1, 
            method="anova", #"class" for classification, "anova" for a regression tree
            )

printcp(rt) # display the results
plotcp(rt) # visualize cross-validation results
summary(rt) # detailed summary of splits
# create additional plots
par(mfrow=c(1,2)) # two plots on one page
rsq.rpart(rt) # visualize cross-validation results

# plot tree
plot(rt, uniform=TRUE,
     main="Regression Tree for Recovery")
text(rt, use.n=TRUE, all=TRUE, cex=.8)

# prune the tree
rt$cptable[which.min(rt$cptable[,"xerror"]),"CP"]
rt$cptable

prt<- prune(rt, cp=0.03) # from cptable   

# plot the pruned tree
plot(prt, uniform=TRUE,
     main="Pruned Regression Tree for Mileage")
text(prt, use.n=TRUE, all=TRUE, cex=.8)

###############################################################################
#post-fire recovery at fire levels
###############################################################################







###########################################################
#examining post-fire VIs trajectory

setwd("D:\\users\\Zhihua\\Landsat")

library(raster)
library(rgdal)
library("ncdf")
library("maptools")

proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
proj.utm = "+proj=utm +zone=51 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

fire.sp = readOGR(dsn=".\\XinganImages\\boundry",layer="FirePatch")
fire.sp@data$Id2 = 1:nrow(fire.sp)
fire.sp@data$area = sapply(slot(fire.sp, "polygons"), slot, "area")

data.vi.mean.df = list()
data.vi.sd.df = list()

for (i in 1:nrow(fire.sp)){
  
  vi.mn = c()
  vi.sd = c()
  
  Year = fire.sp[i,]$Year
  
    #generate some data point
    Rpts1 = spsample(fire.sp[i,], n=ceiling(fire.sp[i,]$area/(30*30*81)), type='regular')
    
    if (Year >= 2000) {
  
    for (j in (Year-1):2015) {
    # 0 year VI
    
    filename = paste("fire_", i,"_", j, sep = "")
    filename = list.files(path = "./analysis_results/", pattern = filename)
    vi00 = stack(paste("./analysis_results/", filename[3], sep = "")) # single images
    
    #extract data
    vi00.df = extract(vi00, Rpts1)
    
    vi00.df.mn = apply(vi00.df, 2, mean, na.rm = TRUE)
    vi00.df.sd = apply(vi00.df, 2, sd, na.rm = TRUE)
    
    vi.mn = rbind(vi.mn,vi00.df.mn)
    vi.sd = rbind(vi.sd,vi00.df.sd)
    
    } #end of j
    
    if (nrow(vi.mn) < 17) {
      na.matrix = matrix(NA, nrow = 17 - nrow(vi.mn), ncol = ncol(vi.mn))
      vi.mn = rbind(vi.mn, na.matrix)
      vi.sd = rbind(vi.sd, na.matrix)
      
    }
    
    } #end of if
  print(paste("Finishing for :::: Fire ", i, " of ", nrow(fire.sp), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
  data.vi.mean.df[[i]] <- data.frame(vi.mn)
  data.vi.sd.df[[i]] <- data.frame(vi.sd)
  
} #end of fire i


#calculate means of all fires
vi.mn = as.data.frame(vi.mn)

plot(data.vi.mean.df[[1]]$ndvi, type = "b")

z<-matrix(0,nrow(data.vi.mean.df[[1]]),ncol(data.vi.mean.df[[1]]))

for(l in 1:nrow(data.vi.mean.df[[1]])){
  for(m in 1:ncol(data.vi.mean.df[[1]])){
    z[l,m]<-mean(unlist(lapply(data.vi.mean.df, `[`, i =l, j = m)), na.rm = TRUE)
  }
}

z = data.frame(z)
colnames(z) <- colnames(data.vi.mean.df[[1]]) 
plot(z$ndvi, type = "b")
plot(z$ndwi, type = "b")
z = z[, c("ndvi", "ndwi","tcw","tca")]
z = data.frame(z, year = -1:15)

z.sd<-matrix(0,nrow(data.vi.mean.df[[1]]),ncol(data.vi.mean.df[[1]]))

for(l in 1:nrow(data.vi.sd.df[[1]])){
  for(m in 1:ncol(data.vi.sd.df[[1]])){
    z.sd[l,m]<-mean(unlist(lapply(data.vi.sd.df, `[`, i =l, j = m)), na.rm = TRUE)
  }
}

z.sd = data.frame(z.sd)
colnames(z.sd) <- colnames(data.vi.mean.df[[1]]) 
plot(z.sd$ndvi, type = "b")

z.sd = z.sd[, c("ndvi", "ndwi","tcw","tca")]
z.sd = data.frame(z.sd, year = -1:15)

library("reshape2")
library("ggplot2")

z.sd.long = melt(z.sd, id.vars=c("year"))
colnames(z.sd.long) <- c("year", "vi","sd")

z.long = melt(z, id.vars=c("year"))
colnames(z.long) <- c("year", "vi","mean")

long.df = data.frame(z.sd.long, mean = z.long$mean)

levels(long.df$vi) <- c("NDVI","NDWI", "TCW","TCA")

ggplot(long.df, aes(x=year, y=mean, colour=vi)) + 
  facet_wrap( ~ vi, scales = "free",ncol=2) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, lwd = 1) +
  geom_line(lwd = 2) +
  geom_point()+
  xlab("Year") + ylab("VI value") + 
  #theme(legend.position="none")+
  theme(legend.position=c(0.95,0.65))+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title=element_blank()) +
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  theme(strip.text.x = element_text(size=18))+
  theme(strip.text.y = element_text(size=18)) 

ggsave(".\\analysis_results2\\vi_recovery.png", width = 12, height = 9, units = "in")


save.image(".\\analysis_results2\\workimage1.RData")
save.image(".\\analysis_results2\\workimage2.RData")


#4/21/2016
#downland modis burned area data to validate the results
#product == MCD45A1, version = 5.1
#MCD45A1.051
  
  
# load libraries
library(MODIS)
library(rgdal)
library(raster)

wkdir <- "D:/users/Zhihua/Landsat/MODIS_validation"
setwd(wkdir)

#set spatial extent
dxal <- readOGR(dsn = "D:/users/Zhihua/Landsat/XinganImages/boundry", layer = "dxal_bj_proj_polygon")
proj.utm = projection(dxal)
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

####################### downlaod MODIS burned area data ##########################
#set MODIS path, see MODIS package for details
MODISoptions(localArcPath="D:/users/Zhihua/Landsat/MODIS_validation",
             outDirPath="D:/users/Zhihua/Landsat/MODIS_validation",
             gdalPath='c:/OSGeo4W64/bin')

#set date period for data
dates <- as.POSIXct( as.Date(c("1/1/2000","1/11/2015"),format = "%d/%m/%Y") )
dates2 <- transDate(dates[1],dates[2]) # Transform input dates from before

# getProduct() # list available MODIS products
#download MOD13Q1 data
runGdal(product="MCD45A1",  #VI/combined/Tile/500m/monthly
        begin=dates2$beginDOY, #start date
        end = dates2$endDOY,#end date
        tileH = 25,tileV = 3, #tile
        SDSstring = "11", #extract the first 3 layers
        extent = dxal, #crop to extent
        outProj=proj.utm, #reproject to UTM
        job = "MCD45A1") #download folder


################## extracting individual fires #############################
# 6/11/2016 at U of M

# load libraries
library(MODIS)
library(rgdal)
library(raster)

wkdir <- "F:/Rosa/Landsat/MODIS_validation"
setwd(wkdir)

#set spatial extent
dxal <- readOGR(dsn = "F:/Rosa/Landsat/XinganImages/boundry", layer = "dxal_bj_proj_polygon")
proj.utm = projection(dxal)
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

ba_date = preStack(path=".\\MCD45A1", pattern="*.burndate.tif$")
ba_quality = preStack(path=".\\MCD45A1", pattern="*.ba_qa.tif$")

#1. produce yearly burned area data
YearDOY = substr(ba_date, 20, 26)

ba.list = list()

for (yr in 2000:2015){
Year2010 = which(substr(YearDOY, 1, 4)==as.character(yr))

ba_date2010 = stack(ba_date[Year2010])
ba_quality2010 = stack(ba_quality[Year2010])

# remove bad data based on data quality
for (i in 1:nlayers(ba_date2010)) {
  
  ba_date2010[[i]][ba_quality2010[[i]] != 1] = 0
  
}

# to assess whether there are multiple burn or not within one year
Non.zero.length = function(x){length(which(x > 0))}
date2010_freq = calc(ba_date2010, Non.zero.length)
#freq(ba_date2010_1)

# pixels with multple date usually have the same date, but may contain some have different date, 
# use the first date to assign the date information

#get the location of the non-zero values
date2010_freq2 = date2010_freq > 0
date2010_freq2[date2010_freq2==0] = NA

Ex.pts.all.nonNA = function(x){
  proj.geo = projection(x)
  pts = rasterToPoints(x) #get raster coordinate, from left to right, top to bottom
  pts = data.frame(pts)
  pts <- SpatialPoints(coords = cbind(pts$x,pts$y),proj4string = CRS(proj.geo))
  pts.sp = SpatialPoints(coords = pts, proj4string = CRS(proj.geo))
  return(pts.sp)
}

#
date2010.sp = Ex.pts.all.nonNA(date2010_freq2)

#extract burn date values
date2010.date.df = extract(ba_date2010, date2010.sp)
#get the burn date for each pixels
date2010.date.df2 = apply(date2010.date.df, 1, function(x){x[which(x>1)[1]]})

#convert from point to raster
date2010.sp2 = data.frame(rasterToPoints(date2010_freq2)[,c(1,2)], date = date2010.date.df2)
coordinates(date2010.sp2) <- ~x+y
proj4string(date2010.sp2) = projection(date2010_freq)
date2010.sp2 = rasterize(date2010.sp2,date2010_freq2, field = "date")

date2010.sp2[date2010.sp2 > 273] = NA

ba.list[[yr - 1999]] = date2010.sp2

print(paste("Finish Extracting Year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}

ba.list = stack(ba.list)
ba.list = crop(ba.list, dxal[3,])

#plot
names(ba.list) <- paste("Year", 2000:2015, sep = "")
ba.list2 = ba.list

#produce one raster, and get the burn date from the rasterStack
ba.grd = ba.list2[[1]]; ba.grd[] = 0
for (i in 1:nlayers(ba.list2)){
ba.grd[!is.na(ba.list2[[i]])] = ba.list2[[i]][!is.na(ba.list2[[i]])] + (i+1999)*1000
}

writeRaster(ba.grd,"F:/Rosa/Landsat/MODIS_validation/results/ba00-15.tif",format="GTiff", overwrite=TRUE) #write a raster stack files

#plot to see the figure
plot(date2010.sp2)
plot(dxal, add = T)

#########START DO NOT USE: BLOCK 1 ############ 
################## extracting individual fires #############################
########### batch processing start from here
Year = 2001:2015

Year_burndate = list()
for (year in 2001:2015) {
Year2010 = which(as.numeric(substr(YearDOY, 1, 4))==year)

ba_date2010 = stack(ba_date[Year2010])
ba_quality2010 = stack(ba_quality[Year2010])

}

# remove bad data based on data quality
for (i in 1:nlayers(ba_date2010)) {
  
#  ba_date2010[[i]][ba_quality2010[[i]] < 1 | ba_quality2010[[i]] > 2] = 0 # for quality 1 + 2
  ba_date2010[[i]][ba_quality2010[[i]] != 1] = 0 # only quality 1 
  
  
}

# to assess whether there are multiple burn or not within one year
Non.zero.length = function(x){length(which(x > 0))}
date2010_freq = calc(ba_date2010, Non.zero.length)
freq(ba_date2010_1)

# pixels with multple date usually have the same date, but may contain some have different date, 
# use the first date to assign the date information

#get the location of the non-zero values
date2010_freq2 = date2010_freq > 0
date2010_freq2[date2010_freq2==0] = NA
date2010.sp = Ex.pts.all.nonNA(date2010_freq2)

#extract burn date values
date2010.date.df = extract(ba_date2010, date2010.sp)
#get the burn date for each pixels
date2010.date.df2 = apply(date2010.date.df, 1, function(x){x[which(x>1)[1]]})

#convert from point to raster
date2010.sp2 = data.frame(rasterToPoints(date2010_freq2)[,c(1,2)], date = date2010.date.df2)
coordinates(date2010.sp2) <- ~x+y
proj4string(date2010.sp2) = projection(date2010_freq)
date2010.sp2 = rasterize(date2010.sp2,date2010_freq2, field = "date")

date2010.sp2[date2010.sp2 > 273] = NA

Year_burndate[[year - 2000]] <-date2010.sp2

print(paste("Finishing for :::: Year ", year, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}

Year_burndate = stack(Year_burndate)
names(Year_burndate) <- paste("Y", 2001:2015, sep = "")

#for Huzhong
Year_burndate_hz = crop(Year_burndate, dxal[3,])


plot(Year_burndate_hz[[3]])
plot(dxal, add = T)

#########END DO NOT USE: BLOCK 1 ############ 

fire.sp = readOGR(dsn="F:/Elias/XinganTM/GIS_data","FirePatch") # manual digitized fire polygon

# PLOTING
#for test area 1
#read into test areas
require("rgeos")
r1 = readOGR(dsn=".\\BEAD_ploys",layer="testr1")

fire.sp1 = crop(fire.sp, r1)
ba.grd1 = crop(ba.grd, r1)

r1_firepoly1 = readOGR(dsn=".\\BEAD_ploys",layer="testr1_4_r4_ndvi")
r1_firepoly2 = readOGR(dsn=".\\BEAD_ploys",layer="testr1_4_r4_ndvi_90_2")

#for test area 3
#read into test areas
r3 = readOGR(dsn=".\\BEAD_ploys",layer="testr3")
fire.sp3 = crop(fire.sp, r3)
ba.grd3 = crop(ba.grd, r3)

r3_firepoly1 = readOGR(dsn=".\\BEAD_ploys",layer="testr3_1_r4_ndvi_90_2")
r3_firepoly2 = readOGR(dsn=".\\BEAD_ploys",layer="testr3_4_r4_ndvi_90_2")


#MODIS burned area data did not detect any burned area in this area


#############downland modis fire hotspot in this area to validate the results
#download data from http://modis-fire.umd.edu/pages/ActiveFire.php?target=GetData
# ftp://fuoco.geog.umd.edu/

#decompress the zip files
install.packages("R.utils")

library(R.utils)

# file.names = list.files(path = "./MCD14MLV5", pattern = "*.gz$")
file.names = list.files(path = "./MCD14MLV5", pattern = "*.asc$")

fire.df = c()

for (i in 1:length(file.names)){
#gunzip(paste("./MCD14MLV5/", file.names[i], sep = ""))

#file.tmp = substr(file.names[i], 1, nchar(file.names[i])-3)
file.tmp = data.frame(read.table(paste("./MCD14MLV5/", file.names[i], sep = ""), header = TRUE))

file.tmp = file.tmp[which(file.tmp$lon > 121 & file.tmp$lon < 127 & file.tmp$lat > 50 & file.tmp$lat < 53.5),]

fire.df = rbind(fire.df, file.tmp)

print(paste("Finishing for extracting ", i, " of ", length(file.names), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}

plot(density(fire.df$FRP), xlim = c(0,120))
plot(density(fire.df$conf), xlim = c(0,120))

#change into spatial database
library(sp)
fire.sp <- fire.df[which(fire.df$conf > 50),]
coordinates(fire.sp) <- ~lon+lat
projection(fire.sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

#dxal <- readOGR(dsn = "D:/users/Zhihua/Landsat/XinganImages/boundry", layer = "dxal_bj_proj_polygon")
proj.utm = projection(dxal)
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
dxal.geo = spTransform(dxal, CRS(proj.geo))

fire.sp.hz = fire.sp[dxal.geo[3,], ]

fire.sp.hz = spTransform(fire.sp.hz, CRS(proj.utm))

#crop into test area
fire.sp.hz1 = fire.sp.hz[r1, ]
fire.sp.hz3 = fire.sp.hz[r3, ]

######## plot final ######
#par(mfrow=c(2,1),mar=c(0,0,0,0))

png("region1-2.png",height = 3000, width = 3000, res = 300, units = "px")

plot(ba.grd1, legend=FALSE, axes=FALSE, box=FALSE)
plot(r1_firepoly2[!is.na(r1_firepoly2$dist_time),], border = "black", lwd = 2, add = TRUE)
plot(fire.sp1[1,], border = "black", lwd = 2, add = TRUE)
#plot(fire.sp1, border = "red", lwd = 2, add = TRUE)
plot(fire.sp.hz1, add = TRUE, pch = 20, col = "blue")

scalebar(10000, xy=c(505000, 5742500), type='bar', divs=4,below = "Meter")
#text(x=510000, y=5767000, "Region 1", cex = 1.5)
dev.off()

png("region2-2.png",height = 3000, width = 4500, res = 300, units = "px")

plot(ba.grd3, legend=FALSE, axes=FALSE, box=FALSE)
plot(r3_firepoly2[!is.na(r3_firepoly2$dist_time),], border = "black", lwd = 2, add = T)
#plot(fire.sp3, border = "red", lwd = 2, add = T)
plot(fire.sp.hz3, add = TRUE, pch = 20, col = "blue")

scalebar(10000, xy=c(498000, 5680000), type='bar', divs=4,below = "Meter")
#text(x=500000, y=5700000, "Region 2", cex = 1.5)
dev.off()


#read into test areas
r1 = readOGR(dsn=".\\BEAD_ploys",layer="testr1")

r1_firepoly1 = readOGR(dsn=".\\BEAD_ploys",layer="testr1_4_r4_ndvi")
r1_firepoly2 = readOGR(dsn=".\\BEAD_ploys",layer="testr1_4_r4_ndvi_90_2")

r1_firepoly1 = r1_firepoly1[which(r1_firepoly1@data$dist_time > 0),]
r1_firepoly2 = r1_firepoly2[which(r1_firepoly2@data$dist_time > 0),]

r1_val = readOGR(dsn=".\\BEAD_ploys",layer="testr1_validation")

#for rest area 3
#read into test areas
r3 = readOGR(dsn=".\\BEAD_ploys",layer="testr3")

r3_firepoly1 = readOGR(dsn=".\\BEAD_ploys",layer="testr3_1_r4_ndvi_90_2")
r3_firepoly2 = readOGR(dsn=".\\BEAD_ploys",layer="testr3_4_r4_ndvi_90_2")

r3_firepoly1 = r3_firepoly1[which(r3_firepoly1@data$dist_time > 0),]
r3_firepoly2 = r3_firepoly2[which(r3_firepoly2@data$dist_time > 0),]

r3_val = readOGR(dsn=".\\BEAD_ploys",layer="testr3_validation")

fire.sp.hz.r1 = fire.sp.hz[r1,]
fire.sp.hz.r3 = fire.sp.hz[r3,]

#plot to see the results
plot(r3)
plot(fire.sp.hz.r3, add = T)
plot(r3_firepoly2, add =T, border = "blue")
plot(r3_firepoly1, add =T, border = "red")
plot(r3_val, add =T, border = "green")

plot(r1)
plot(fire.sp.hz.r1, add = T)
plot(r1_firepoly2, add =T, border = "blue")
plot(r1_firepoly1, add =T, border = "red")

plot(r1_val, add =T, border = "green")


########################################################
#6/11/2016 at U of Montana
#
setwd("F:/Elias/XinganTM/GIS_data")

library(rgdal)
library(raster)
library(maptools)

########################################################################################################
############### begin of reading into validation fire  and bead-extracted fire patches #################
#read tested regions
testr1 = readOGR(dsn=".",layer="testr1")
testr3 = readOGR(dsn=".",layer="testr3")

#read validation fire patches 
fire.sp <- readOGR(".", "FirePatch")


