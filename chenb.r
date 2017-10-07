#download MODIS brdf to VI data for BAED 

# load libraries
library(rgdal)
library(raster)

####################  Data pre-processing  #######################
# list files
Work.Dir <- "D:/users/Zhihua/Landsat/p122024"
setwd(Work.Dir)

# read a polygon shapefile
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
  
 
#数据裁切 crop
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
