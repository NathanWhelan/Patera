library(ENMTools)
library(raster)
library(rgbif)
library(maptools)
library(rgdal)
print("libraries loaded")

######This script will take a long time to run. I suggest splitting up the ID tests into seperate scripts so multiple ID tests can run at once if you want to speed things up.


setwd("~/Patera/Niche-modeling/thinning/") #Change to for your file system.
rm(list=ls())
unixtools::set.tempdir("/home/labgroup/Patera/Niche-modeling/thinning/tmp4") ##May want to use non-standard tmp folder. Some steps create massive temporary folders that may be bigger than you system's tmp folder can handle.
#wclim<-getData("worldclim", var="bio", res=0.5, lon = -83, lat =35, path = "~/Patera/Niche-modeling/thinning/") #Step for downloading bioclim data. I did this and then pulled it into arcGIS Pro, so it's really pre-processing step.

###Data Wrangling####
##get wclim data from arcgis output or data uploaded to FigShare
wclim_high_res1<-raster("soils/Soils/bio1_13_resize2.tif")
wclim_high_res2<-raster("soils/Soils/bio2_13_2.tif")
wclim_high_res3<-raster("soils/Soils/bio3_13_resize2.tif")
wclim_high_res4<-raster("soils/Soils/bio4_13_reszie2.tif")
wclim_high_res5<-raster("soils/Soils/bio5_13_resize2.tif")
wclim_high_res6<-raster("soils/Soils/bio6_13_revised2.tif")
wclim_high_res7<-raster("soils/Soils/bio7_13_resize2.tif")
wclim_high_res8<-raster("soils/Soils/bio8_13_resize2.tif")
wclim_high_res9<-raster("soils/Soils/bio9_13_resize2.tif")
wclim_high_res10<-raster("soils/Soils/bio10_13_resize2.tif")
wclim_high_res11<-raster("soils/Soils/bio11_13_resize2.tif")
wclim_high_res12<-raster("soils/Soils/bio12_13_resize2.tif")
wclim_high_res13<-raster("soils/Soils/bio13_13_1.tif")
wclim_high_res14<-raster("soils/Soils/bio14_13_resize2.tif")
wclim_high_res15<-raster("soils/Soils/bio15_13_1.tif")
wclim_high_res16<-raster("soils/Soils/bio16_13_resize2.tif")
wclim_high_res17<-raster("soils/Soils/bio17_13_resize2.tif")
wclim_high_res18<-raster("soils/Soils/bio18_13_resize2.tif")
wclim_high_res19<-raster("soils/Soils/bio19_13_resize2.tif")
# 
# 
# #get USGS data (elevation, aspect, slope)
aspect<-raster(x="soils/Soils/aspect_usgs1_Band_resize2.tif")
elevation<-raster(x="soils/Soils/usgselevation_resize2.tif")
slope<-raster(x="soils/Soils/Slope_usgsel1_resize2.tif")

# ##USA Soils Data
erodability<-raster(x="soils/Soils/USA SSURGO  Erodibil_resize2.tif")
albedo<-raster(x="soils/Soils/USA_Soils_Albedo_resize2.tif")
erodability[is.na(erodability)]<-0
albedo[is.na(albedo)]<-0

##Landfire
canopy_cover<-raster(x="soils/Soils/LC20_CC_200_resize2.tif")
canopy_base_height<-raster(x="soils/Soils/LC20_CBH_200_resize2.tif")
vegetation_height<-raster(x="soils/Soils/LC16_EVH_200_resize2.tif")
vegetation_cover<-raster(x="soils/Soils/LC16_EVC_200_resize2.tif")

print("tifs read")
print("start stacking")

###Correlation Testing###  
##Done prior to trimming down to the varibales placed in all_cont_variables

all_cont_variables_correlation_test<-stack(wclim_high_res1,
                          wclim_high_res2,
                          wclim_high_res3,
                          wclim_high_res4,
                          wclim_high_res5,
                          wclim_high_res6
                          wclim_high_res7,
                          wclim_high_res8,
                          wclim_high_res9,
                          wclim_high_res10,
                          wclim_high_res11,
                          wclim_high_res12,
                          wclim_high_res13,
                          wclim_high_res14,
                          wclim_high_res15,
                          wclim_high_res16,
                          wclim_high_res17,
                          wclim_high_res18,
                          wclim_high_res19,
                          albedo,
                          aspect,
                          elevation,
                          slope,
                          erodability,
                          canopy_cover,
                          canopy_base_height,
                          vegetation_cover)

#raster.cor.matrix(all_cont_variables_correlation_test)
#raster.cor.plot(all_cont_variables_correlation_test)




all_cont_variables<-stack(wclim_high_res1,
                          wclim_high_res2,
                          wclim_high_res3,
                          wclim_high_res4,
                          wclim_high_res7,
                          wclim_high_res8,
                          wclim_high_res9,
                          wclim_high_res12,
                          wclim_high_res15,
                          albedo,
                          aspect,
                          elevation,
                          slope,
                          erodability,
                          canopy_cover,
                          canopy_base_height,
                          vegetation_cover)

bioclimate_variables<-stack(wclim_high_res1,
                          wclim_high_res2,
                          wclim_high_res3,
                          wclim_high_res4,
                          wclim_high_res7,
                          wclim_high_res8,
                          wclim_high_res9,
                          wclim_high_res12,
                          wclim_high_res15)

non_climate<-stack(albedo,
                   aspect,
                   elevation,
                   slope,
                   erodability,
                   canopy_cover,
                   canopy_base_height,
                   vegetation_cover)

                   ##Did not end up using categorical variables
#non_climate_cat<-stack(albedo,
#                       aspect,
#                       elevation,
#                       slope,
#                       erodability,
#                       canopy_cover,
#                       canopy_base_height,
#                       vegetation_cover,
#                       drainage_categorical,
#                       erosion_categorical,
#                       hydrologic_categorical)

#all_variables<-stack(wclim_high_res1,
#                          wclim_high_res2,
#                          wclim_high_res3,
#                          wclim_high_res4,
#                          wclim_high_res7,
#                          wclim_high_res8,
#                         wclim_high_res9,
#                          wclim_high_res12,
#                          wclim_high_res15,
#                          albedo,
#                          aspect,
#                          elevation,
#                          slope,
#                          erodability,
#                          canopy_cover,
#                          canopy_base_height,
#                          vegetation_cover,
#                          drainage_categorical,
#                          erosion_categorical,
#                          hydrologic_categorical)


print("stacking done")








##CSV should be "species,longitude,latitude"
<-read.csv("~/Patera/Niche-modeling/thinning/clarki-nantahala_1KM-thinned_rarefied_points.csv")

##Format for ENMTools
colnames(x)[1:3]<-c("species","Longitude","Latitude")

# ##Write each species to a separate file, with only "Longitude" and "Latitude columns
write.csv(x[which(x$species == names(table(x$species))[1]),2:3], "clarki.csv", row.names = FALSE)
write.csv(x[which(x$species == names(table(x$species))[2]),2:3], "nantahala.csv", row.names = FALSE)

#Import the two species occurrences as enmtools.species objects
sp1.path <- 'clarki.csv'
sp2.path <- 'nantahala.csv'

#Import the two species occurrences as enmtools.species objects
sp1 <- enmtools.species(species.name = "sp1", presence.points = read.csv(sp1.path))
sp2 <- enmtools.species(species.name = "sp2", presence.points = read.csv(sp2.path))

####Plotting Environmental Niche Models###########

##Only bioclimate glm####
clarki.glm.climate.only<-enmtools.glm(sp1,env=bioclimate_variables,test.prop=0.2)
pdf("clarki-glm-climate-only.pdf")
clarki.glm.climate.only
dev.off()

nantahala.glm.climate.only<-enmtools.glm(sp2,env=bioclimate_variables,test.prop=0.2)
pdf("nantahala-glm-climate-only.pdf")
nantahala.glm.climate.only
dev.off()



####Includes all variables GLM
clarki.glm.climate<-enmtools.glm(sp1,env=all_cont_variables,test.prop=0.2)
pdf("clarki-glm-climate.pdf")
clarki.glm.climate
dev.off()
 
nantahala.glm.climate<-enmtools.glm(sp2,env=all_cont_variables,test.prop=0.2)
pdf("nantahala-glm-climate.pdf")
nantahala.glm.climate
dev.off()
 
### Non-bioclimatic variables GLM
clarki.glm.no_worldclim<-enmtools.glm(sp1,env=non_climate,test.prop=0.2)
pdf("clarki-glm_no-worldclim.pdf")
clarki.glm.no_worldclim
dev.off()


nantahala.glm.no_worldclim<-enmtools.glm(sp2,env=non_climate,test.prop=0.2)
pdf("nantahala-glm-no-worldclim.pdf")
nantahala.glm.no_worldclim
dev.off()


####Maxent with all variables
 
clarki.mx.climate<-enmtools.maxent(sp1,env=all_cont_variables,test.prop=0.2)
pdf("clarki-mx-climate.pdf")
plot(clarki.mx.climate)
dev.off()
 
nantahala.mx.climate<-enmtools.maxent(sp2,env =all_cont_variables,test.prop=0.2)
pdf("nantahala-mx-climate.pdf")
plot(nantahala.mx.climate)
dev.off()
 
#####Maxent without bioclimatic variables
clarki.mx.no_worldclim<-enmtools.maxent(sp1,env=non_climate,test.prop=0.2)
pdf("clarki-mx-no-worldclim")
plot(clarki.mx.no_worldclim)
dev.off()
 
nantahala.mx.no_worldclim<-enmtools.maxent(sp2,env=non_climate,test.prop=0.2)
pdf("nantahala-mx-no-worldclim.pdf")
plot(nantahala.mx.no_worldclim)
dev.off()

#######GLM and MAXENT ID Test using only bioclimating data###########
dir.create("bioclimate_glm")
 
# ##Set the number of replicate runs of the identity.test() function

reps <- 100 ##should normally be 100

##Set the number of background points to be used

back <- 10000 ##should normally be 10000

##Perform the identity test

idtest_glm <- identity.test(
  sp1,
  sp2,
  bioclimate_variables,
  type = 'glm',
  f = NULL,
  nreps = reps,
  nback = back,
  low.memory = FALSE,
  bg.source="env",
  rep.dir = "bioclimate_glm",
  verbose = FALSE,
  clamp = TRUE,
)

###Write the raw results to a file
write.csv(idtest_glm$reps.overlap, "bioclimate_glm/reps_overlap.csv")

##Split the permuted results from the empirical results
permute.reps <- idtest_glm$reps.overlap[2:nrow(idtest_mx$reps.overlap),]

##Set the percentile you want to use from your permuted results to obtain the critical values
critical.percent <- 95
critical.number <- round((critical.percent/100)*(reps))
d.critical <- as.vector(permute.reps[order(permute.reps[,"D"], decreasing = TRUE),"D"])[critical.number]
i.critical <- as.vector(permute.reps[order(permute.reps[,"I"], decreasing = TRUE),"I"])[critical.number]
rank.cor.critical <- as.vector(permute.reps[order(permute.reps[,"rank.cor"], decreasing = TRUE),"rank.cor"])[critical.number]

##Create a summary results object
results.summary <- rbind(

  idtest_glm$reps.overlap[1,1:3],
  c(d.critical, i.critical, rank.cor.critical)

)

##Create the proper row names for the summary results object
rownames(results.summary) <- c("empirical.value", "permuted.critical.value")

##Write the summary results object to a file
write.csv(results.summary, "bioclimate_glm/results_summary.csv")



dir.create("bioclimate_mx")
 
##Perform the identity test

idtest_mx <- identity.test(
  sp1,
  sp2,
  bioclimate_variables,
  type = 'mx',
  f = NULL,
  nreps = reps,
  nback = back,
  low.memory = FALSE,
  bg.source="env",
  rep.dir = "bioclimate_mx",
  verbose = FALSE,
  clamp = TRUE,
)

# ##Write the raw results to a file
write.csv(idtest_mx$reps.overlap, "bioclimate_mx/reps_overlap.csv")

##Split the permuted results from the empirical results
permute.reps <- idtest_mx$reps.overlap[2:nrow(idtest_mx$reps.overlap),]

##Set the percentile you want to use from your permuted results to obtain the critical values
critical.percent <- 95
critical.number <- round((critical.percent/100)*(reps))
d.critical <- as.vector(permute.reps[order(permute.reps[,"D"], decreasing = TRUE),"D"])[critical.number]
i.critical <- as.vector(permute.reps[order(permute.reps[,"I"], decreasing = TRUE),"I"])[critical.number]
rank.cor.critical <- as.vector(permute.reps[order(permute.reps[,"rank.cor"], decreasing = TRUE),"rank.cor"])[critical.number]

##Create a summary results object
results.summary <- rbind(

  idtest_mx$reps.overlap[1,1:3],
  c(d.critical, i.critical, rank.cor.critical)

)

##Create the proper row names for the summary results object
rownames(results.summary) <- c("empirical.value", "permuted.critical.value")

##Write the summary results object to a file
write.csv(results.summary, "bioclimate_mx/results_summary.csv")
 



############################GLM ID Test with Worldlcim and Climate#############################
##GLM ID TEST with all variables worldclim
dir.create("soil_climate_glm2")
idtest_glm <- identity.test(
  sp1,
  sp2,
  all_cont_variables,
  type = 'glm',
  nreps = reps,
  nback = back,
  low.memory = FALSE,
  bg.source="env",
  rep.dir = "soil_climate_glm2",
  verbose = FALSE,
  clamp = TRUE,
)
# 
# 
# 
# ##Write the raw results to a file
write.csv(idtest_glm$reps.overlap, "soil_climate_glm2/reps_overlap.csv")

##Split the permuted results from the empirical results
permute.reps <- idtest_glm$reps.overlap[2:nrow(idtest_glm$reps.overlap),]
 
##Set the percentile you want to use from your permuted results to obtain the critical values
critical.percent <- 95
critical.number <- round((critical.percent/100)*(reps))
d.critical <- as.vector(permute.reps[order(permute.reps[,"D"], decreasing = TRUE),"D"])[critical.number]
i.critical <- as.vector(permute.reps[order(permute.reps[,"I"], decreasing = TRUE),"I"])[critical.number]
rank.cor.critical <- as.vector(permute.reps[order(permute.reps[,"rank.cor"], decreasing = TRUE),"rank.cor"])[critical.number]
 
# ##Create a summary results object
results.summary <- rbind(
 idtest_glm$reps.overlap[1,1:3],
 c(d.critical, i.critical, rank.cor.critical)
)
 
###Create the proper row names for the summary results object
 
rownames(results.summary) <- c("empirical.value", "permuted.critical.value")

##Write the summary results object to a file 
write.csv(results.summary, "soil_climate_glm2/results_summary.csv")


###GLM ID Test with non-bioclimatic variables
dir.create("soil_glm")
idtest_soil_glm <- identity.test(
  sp1,
  sp2,
  non_climate,
  type = 'glm',
  nreps = reps,
  nback = back,
  low.memory = FALSE,
  bg.source="env",
  rep.dir = "soil_glm",
  verbose = FALSE,
  clamp = TRUE,
)
####Write the raw results to a file
write.csv(idtest_soil_glm$reps.overlap, "soil_glm/reps_overlap.csv")
##Split the permuted results from the empirical results
permute.reps <- idtest_soil_glm$reps.overlap[2:nrow(idtest_soil_glm$reps.overlap),]

##Set the percentile you want to use from your permuted results to obtain the critical values
 
critical.percent <- 95
critical.number <- round((critical.percent/100)*(reps))
d.critical <- as.vector(permute.reps[order(permute.reps[,"D"], decreasing = TRUE),"D"])[critical.number]
i.critical <- as.vector(permute.reps[order(permute.reps[,"I"], decreasing = TRUE),"I"])[critical.number]
rank.cor.critical <- as.vector(permute.reps[order(permute.reps[,"rank.cor"], decreasing = TRUE),"rank.cor"])[critical.number]
 
##Create a summary results object
 
results.summary <- rbind(
   
  idtest_soil_glm$reps.overlap[1,1:3],
  c(d.critical, i.critical, rank.cor.critical)
)
 
##Create the proper row names for the summary results object
rownames(results.summary) <- c("empirical.value", "permuted.critical.value")
##Write the summary results object to a file 
write.csv(results.summary, "soil_glm/results_summary.csv")
 
########################Maxent ID Test with Worldlcim and Climate#############
##Maxent without worldclim
dir.create("soil_maxent")
idtest_soil_maxent <- identity.test(
  sp1,
  sp2,
  non_climate,
  type = 'mx',
  f = NULL,
  nreps = reps,
  nback = back,
  low.memory = FALSE,
  bg.source="env",
  rep.dir = "soil_maxent",
  verbose = FALSE,
  clamp = TRUE,
)
 
 
 
##Write the raw results to a file
write.csv(idtest_soil_maxent$reps.overlap, "soil_maxent/reps_overlap.csv")
##Split the permuted results from the empirical results
permute.reps <- idtest_soil_maxent$reps.overlap[2:nrow(idtest_soil_maxent$reps.overlap),]
 
###Set the percentile you want to use from your permuted results to obtain the critical values
critical.percent <- 95
critical.number <- round((critical.percent/100)*(reps))
d.critical <- as.vector(permute.reps[order(permute.reps[,"D"], decreasing = TRUE),"D"])[critical.number]
i.critical <- as.vector(permute.reps[order(permute.reps[,"I"], decreasing = TRUE),"I"])[critical.number]
rank.cor.critical <- as.vector(permute.reps[order(permute.reps[,"rank.cor"], decreasing = TRUE),"rank.cor"])[critical.number]
 
##Create a summary results object
results.summary <- rbind(
  idtest_soil_maxent$reps.overlap[1,1:3],
  c(d.critical, i.critical, rank.cor.critical)
   
)
# 
# ##Create the proper row names for the summary results object
rownames(results.summary) <- c("empirical.value", "permuted.critical.value")
##Write the summary results object to a file 
write.csv(results.summary, "soil_maxent/results_summary.csv")




######Maxent ID Test with all continuous variables###
dir.create("soil_climate2")
 
# ##Set the number of replicate runs of the identity.test() function

reps <- 100 ##should normally be 100

##Set the number of background points to be used

back <- 10000 ##should normally be 10000

##Perform the identity test

idtest_mx <- identity.test(
  sp1,
  sp2,
  all_cont_variables,
  type = 'mx',
  f = NULL,
  nreps = reps,
  nback = back,
  low.memory = FALSE,
  bg.source="env",
  rep.dir = "soil_climate2",
  verbose = FALSE,
  clamp = TRUE,
)


##Write the raw results to a file
write.csv(idtest_mx$reps.overlap, "soil_climate2/reps_overlap.csv")

##Split the permuted results from the empirical results
permute.reps <- idtest_mx$reps.overlap[2:nrow(idtest_mx$reps.overlap),]

##Set the percentile you want to use from your permuted results to obtain the critical values
critical.percent <- 95
critical.number <- round((critical.percent/100)*(reps))
d.critical <- as.vector(permute.reps[order(permute.reps[,"D"], decreasing = TRUE),"D"])[critical.number]
i.critical <- as.vector(permute.reps[order(permute.reps[,"I"], decreasing = TRUE),"I"])[critical.number]
rank.cor.critical <- as.vector(permute.reps[order(permute.reps[,"rank.cor"], decreasing = TRUE),"rank.cor"])[critical.number]

##Create a summary results object
results.summary <- rbind(

  idtest_mx$reps.overlap[1,1:3],
  c(d.critical, i.critical, rank.cor.critical)

)

##Create the proper row names for the summary results object
rownames(results.summary) <- c("empirical.value", "permuted.critical.value")

##Write the summary results object to a file
write.csv(results.summary, "soil_climate2/results_summary.csv")