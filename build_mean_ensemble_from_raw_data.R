# This script build the mean across the model ensemble for heat, drought and flood based 
rm(list=ls())


packages <- c("raster", "sf", "terra", "ncdf4")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


setwd(dirname(rstudioapi::getSourceEditorContext()$path))


##############################################
#defining all functions 
##################################################

compute_weights_drought<-function(model_rasters,observations,D_const,model_names){
  all_distance<-matrix(NA,nrow=length(model_rasters),ncol=length(model_names))
  for (y in 1:length(model_rasters)){
    predicted<-model_rasters[[y]]
    observation<-observations[[y]]
    all_points<-values((observation-predicted)^2)
    obs_number<-sapply(1:ncol(all_points),FUN=function(x) length((all_points[,x])[which(!is.na(all_points[,x]))]) )
    annual_distance<-sqrt(colSums(all_points,na.rm=T)/obs_number)
    all_distance[y,]<-annual_distance
    
  }
  each_distance<-colMeans(all_distance)
  names(each_distance)<-model_names
  the_weights<-exp(-(each_distance^2/D_const^2))
  the_weights_nrm_avrg<-the_weights/sum(the_weights)
  
  each_distance2<-apply(all_distance, 2, FUN = min)
  names(each_distance2)<-model_names
  the_weights<-exp(-(each_distance2^2/D_const^2))
  the_weights_nrm_min<-the_weights/sum(the_weights)
  return (data.frame(each_distance,each_distance2,the_weights_nrm_avrg,the_weights_nrm_min))
}



formatRast <- function(rst,the_shape) {
  r_rotated <- raster::rotate(rst)
  r_croped <- raster::crop(r_rotated, the_shape)
  rst <- raster::mask(r_croped, the_shape)
  rst
}


#input the word shape file containing country boundaries
input_fold<-'./Input/'
output_fold<-'./Output2/'
dir.create(paste0(output_fold,'Intermediate/'),recursive=T)
world_shape <- st_read(paste0(input_fold,'world_shape2.geojson'))


###drought#########################
#for drought, weighted average is computed. Weights are computed based on calibration from observation data

#read the observations from the baseline period and store them in a list
base_SPEI<-raster::brick(paste0(input_fold,'./SPEI/SPEI12_observation.nc'))
all_the_years<-format(getZ(base_SPEI),'%Y')
year_ind<-match(as.character(1995:2014),all_the_years)
year_ind<-year_ind+11
base_SPEI<-base_SPEI[[year_ind]]

#read the generated SPEI for the baseline period
all_f_hist<-list.files(paste0(input_fold,'SPEI/hist/'),pattern='.nc$')
model_names<-sapply(all_f_hist,FUN=function(x) substring(x,17,nchar(x)-14),USE.NAMES = FALSE )

example_brick<-raster::brick(paste0(input_fold,'SPEI/hist/spei12_december_',model_names[1],'_historical.nc'))
all_the_years<-format(getZ(example_brick),'%Y')
year_ind<-match(as.character(1995:2014),all_the_years)

modeled_per_year<-list()
for (y in year_ind){
  print (y)
  model_stack<-raster::stack()
  for (f in model_names){
    model_stack<-addLayer(model_stack,raster::brick(paste0(input_fold,'SPEI/hist/spei12_december_',f,'_historical.nc'))[[y]])
  }
  model_stack<-formatRast(model_stack,world_shape)
  
  modeled_per_year[[which(year_ind==y)]]<-model_stack
}

#saveRDS(modeled_per_year,paste0(output_fold,'/Intermediate/baseline_drought_per_model_per_year.Rds'))
#modeled_per_year<-readRDS(paste0(output_fold,'/Intermediate/baseline_drought_per_model_per_year.Rds'))

#aggregate the observations to the model resolution
base_SPEI2<-projectRaster(base_SPEI,modeled_per_year[[1]][[1]])
#compute weight based on validation from the observed SPEI
dist_and_weights<-compute_weights_drought(modeled_per_year,base_SPEI2,0.2,model_names)
model_weights<-dist_and_weights$the_weights_nrm_avrg

all_years<-raster::stack()
for (y in year_ind){
  yearly_rast_stack<-raster::stack()
  for (n in (1:length(model_names))){
    each_model_rast<-raster::brick(paste0(input_fold,'SPEI/hist/spei12_december_',model_names[n],'_historical.nc'))[[y]]
    each_model_rast<-model_weights[n]*each_model_rast
    yearly_rast_stack<-addLayer(yearly_rast_stack,each_model_rast)
  }
  the_rast <- raster::calc(yearly_rast_stack,fun=function(x) if(all(is.na(x))) NA_integer_ else sum(x, na.rm=TRUE))
  all_years<-addLayer(all_years,the_rast)
}
#saveRDS(all_years,paste0(output_fold,'Intermediate/hist_ensemble_drought.RData'))
writeRaster(all_years, paste0(output_fold,'Intermediate/hist_ensemble_drought.tif'), format = "GTiff", overwrite = TRUE)

#future projection (SSP585)
example_brick<-raster::brick(paste0(input_fold,'SPEI/ssp585/spei12_december_',model_names[1],'_ssp585.nc'))
all_the_years<-format(getZ(example_brick),'%Y')
year_ind<-match(as.character(2020:2040),all_the_years)
all_years<-raster::stack()
for (y in year_ind){
  yearly_rast_stack<-raster::stack()
  for (n in (1:length(model_names))){
    each_model_rast<-raster::brick(paste0(input_fold,'SPEI/ssp585/spei12_december_',model_names[n],'_ssp585.nc'))[[y]]
    each_model_rast<-model_weights[n]*each_model_rast
    yearly_rast_stack<-addLayer(yearly_rast_stack,each_model_rast)
  }
  the_rast <- raster::calc(yearly_rast_stack,fun=function(x) if(all(is.na(x))) NA_integer_ else sum(x, na.rm=TRUE))
  all_years<-addLayer(all_years,the_rast)
}
#saveRDS(all_years,paste0(output_fold,'Intermediate/ssp_ensemble_drought.RData'))
writeRaster(all_years, paste0(output_fold,'Intermediate/ssp_ensemble_drought.tif'), format = "GTiff", overwrite = TRUE)


#future projection (SSP126)
all_f_ssp1<-list.files(paste0(input_fold,'SPEI/ssp126/'),pattern='.nc4$')
model_names_ssp1<-sapply(all_f_ssp1,FUN=function(x) substring(x,8,nchar(x)-11),USE.NAMES = FALSE )

example_brick<-raster::brick(paste0(input_fold,'SPEI/ssp126/spei12_',model_names_ssp1[1],'_ssp126.nc4'))
all_the_years<-format(getZ(example_brick),'%Y')
year_ind<-match(as.character(2020:2040),all_the_years)+11
modeled_per_year<-list()
for (y in year_ind){
  print (y)
  model_stack<-raster::stack()
  for (f in model_names_ssp1){
    model_stack<-addLayer(model_stack,raster::brick(paste0(input_fold,'SPEI/ssp126/spei12_',f,'_ssp126.nc4'))[[y]])
  }
  model_stack<-formatRast(model_stack,world_shape)
  modeled_per_year[[which(year_ind==y)]]<-model_stack
}

mean_base <- raster::calc(base_SPEI,fun=function(x) if(all(is.na(x))) NA_integer_ else mean(x, na.rm=TRUE))
mean_base2<-projectRaster(mean_base,modeled_per_year[[1]][[1]])
all_years<-raster::stack()
for (y in (1:length(modeled_per_year))){
  dist_and_weights2<-compute_weights_drought(modeled_per_year[y],mean_base2,0.2,model_names_ssp1)
  model_weights_year<-dist_and_weights2$the_weights_nrm_avrg
  yearly_rast_stack<-raster::stack()
  for (n in (1:length(model_names_ssp1))){
    each_model_rast<-raster::brick(paste0(input_fold,'SPEI/ssp126/spei12_',model_names_ssp1[n],'_ssp126.nc4'))[[(year_ind[y])]]
    each_model_rast<-model_weights_year[n]*each_model_rast
    yearly_rast_stack<-addLayer(yearly_rast_stack,each_model_rast)
  }
  the_rast <- raster::calc(yearly_rast_stack,fun=function(x) if(all(is.na(x))) NA_integer_ else sum(x, na.rm=TRUE))
  all_years<-addLayer(all_years,the_rast)
}
#saveRDS(all_years,paste0(output_fold,'Intermediate/ssp126_ensemble_drought.RData'))
writeRaster(all_years, paste0(output_fold,'Intermediate/ssp126_ensemble_drought.tif'), format = "GTiff", overwrite = TRUE)



######################################################
#heat
####################################
#baseline
all_f<-list.files(paste0(input_fold,'HI/hist/'),pattern='HI_da_n_day')
#find the model with the highest resolution, then aggregate the rest to it
all_res<-c()
for (b in (1:length(all_f))){
  our_rast<-raster::brick(paste0(input_fold,'HI/hist/',all_f[b]))
  all_res<-c(all_res,min(res(our_rast)))
}
in_min<-which.min(all_res)
the_brick<-raster::brick(paste0(input_fold,'HI/hist/',all_f[in_min])) #the one with the highest resolution
all_bricks<-list()
for (b in (1:length(all_f))){
  our_rast<-raster::brick(paste0(input_fold,'HI/hist/',all_f[b]))
  br<-projectRaster(our_rast,the_brick[[1]])
  all_bricks<-append(all_bricks,br)
}
#build the mean ensemble and make it to the same resolution as the SPEI
#drought_ref<-readRDS(paste0(output_fold,'Intermediate/hist_ensemble_drought.RData'))
drought_ref<-raster(paste0(output_fold,'Intermediate/hist_ensemble_drought.tif'))
all_years<-raster::stack()
for (y in 1:nlayers(the_brick)){
  yearly_rast_stack<-raster::stack()
  for (n in (1:length(all_f))){
    yearly_rast_stack<-addLayer(yearly_rast_stack,all_bricks[[n]][[y]])
  }
  the_rast <- stackApply(yearly_rast_stack, indices =  rep(1,nlayers(yearly_rast_stack)), fun = "mean", na.rm = T)
  the_rast<-projectRaster(the_rast,drought_ref[[1]])
  all_years<-addLayer(all_years,the_rast)
}
#saveRDS(all_years,paste0(output_fold,'Intermediate/hist_ensemble_heat.RData'))
writeRaster(all_years, paste0(output_fold,'Intermediate/hist_ensemble_heat.tif'), format = "GTiff", overwrite = TRUE)


#heat future (ssp585)
future_heat<-raster::brick(paste0(input_fold,'HI/ssp585/ENS_mean_out.nc'),varname='da_n_day')
future_heat<-future_heat[[6:26]]
future_heat<-projectRaster(future_heat,drought_ref[[1]])
future_heat[future_heat<0]<-0
#saveRDS(future_heat,paste0(output_fold,'Intermediate/ssp_ensemble_heat.RData'))
writeRaster(future_heat, paste0(output_fold,'Intermediate/ssp_ensemble_heat.tif'), format = "GTiff", overwrite = TRUE)


#heat future (ssp126)
all_f<-list.files(paste0(input_fold,'HI/ssp126/'))
#find the model with the highest resolution, then aggregate the rest to it
all_res<-c()
for (b in (1:length(all_f))){
  our_rast<-terra::rast(paste0(input_fold,'HI/ssp126/',all_f[b]),subds = "da_n_day")
  all_res<-c(all_res,min(res(our_rast)))
}
in_min<-which.min(all_res)
the_brick<-raster::brick(paste0(input_fold,'HI/ssp126/',all_f[in_min]),varname='da_n_day') #the one with the highest resolution
all_bricks<-list()
for (b in (1:length(all_f))){
  
  if (all_f[b]=='HI_thres_n_day_FGOALS-g3_ssp126_r1i1p1f1_2020-2040.nc'){
    #build points data since the points are not equally spaced for this specific model
    our_rast<-rast(paste0(input_fold,'HI/ssp126/',all_f[b]),subds = "da_n_day")
    rr<-nc_open(paste0(input_fold,'HI/ssp126/',all_f[b]))
    the_values<-ncvar_get(rr, varid = "da_n_day")
    the_list<-list(NA,21)
    for (t in(1:21)){
      lon_matrix<-replicate(length(rr$dim$lat$vals), rr$dim$lon$vals)-180
      lat_matrix<-t(replicate(length(rr$dim$lon$vals),rr$dim$lat$vals))
      dat<-data.frame('lat'=as.vector(lat_matrix),'lon'=as.vector(lon_matrix),'value'=as.vector(the_values[,,t]))
      the_list[[t]]<-dat
    }
    nc_close(rr)
    #e <- ext(min(rr$dim$lon$vals), max(rr$dim$lon$vals),  min(rr$dim$lat$vals), max(rr$dim$lat$vals)) # xmin, xmax, ymin, ymax
    e<-ext(-180, 180, -90, 90)
    res_x <- max(diff(rr$dim$lon$vals))
    res_y <- max(diff(rr$dim$lat$vals))
    r_template <- rast(e, res = c(res_x, res_y))
    our_spatrast<-rast()
    for (t in(1:21)){
      points_sf <- vect(the_list[[t]], geom = c("lon", "lat"))
      rasterized_data <- rasterize(points_sf, r_template, field = "value", fun = mean)
      our_spatrast<-c(our_spatrast,rasterized_data)
    }
    our_rast<-brick(our_spatrast)
  }else{
    our_rast<-brick(paste0(input_fold,'HI/ssp126/',all_f[b]),varname = "da_n_day")
  }
  br<-projectRaster(our_rast,the_brick[[1]])
  all_bricks<-append(all_bricks,br)
}
#build the mean ensemble and make it to the same resolution as the SPEI
drought_ref<-raster(paste0(output_fold,'Intermediate/hist_ensemble_drought.tif'))
all_years<-raster::stack()
for (y in 1:nlayers(the_brick)){
  yearly_rast_stack<-raster::stack()
  for (n in (1:length(all_f))){
    yearly_rast_stack<-addLayer(yearly_rast_stack,all_bricks[[n]][[y]])
  }
  the_rast <- stackApply(yearly_rast_stack, indices =  rep(1,nlayers(yearly_rast_stack)), fun = "mean", na.rm = T)
  the_rast<-projectRaster(the_rast,drought_ref[[1]])
  all_years<-addLayer(all_years,the_rast)
}
#saveRDS(all_years,paste0(output_fold,'Intermediate/ssp126_ensemble_heat.RData'))
writeRaster(all_years, paste0(output_fold,'Intermediate/ssp126_ensemble_heat.tif'), format = "GTiff", overwrite = TRUE)


#################################
##flood
######################################
all_years_period<-list(1995:2014,2020:2040,2020:2040)
folders=c('hist','ssp585','ssp126')
names<-c('hist','ssp','ssp126')
for (period in 1:3){
  all_f<-list.files(paste0(input_fold,'./SDII/',folders[period],'/'),pattern='.nc$')
  the_brick<-raster::brick(paste0(input_fold,'./SDII/',folders[period],'/',all_f[1]))
  all_years<-format(getZ(the_brick),'%Y')
  year_ind<-match(as.character(all_years_period[[period]]),all_years)
  all_years<-raster::stack()
  for (y in year_ind){
    yearly_rast_stack<-raster::stack()
    for (n in (1:length(all_f))){
      yearly_rast_stack<-addLayer(yearly_rast_stack,raster::brick(paste0(input_fold,'./SDII/',folders[period],'/',all_f[n]))[[y]])
    }
    the_rast <- stackApply(yearly_rast_stack, indices =  rep(1,nlayers(yearly_rast_stack)), fun = "mean", na.rm = T)
    all_years<-addLayer(all_years,the_rast)
  }
  #saveRDS(all_years,paste0(output_fold,'Intermediate/',names[period],'_ensemble_flood.RData'))
  writeRaster(all_years, paste0(output_fold,'Intermediate/',names[period],'_ensemble_flood.tif'), format = "GTiff", overwrite = TRUE)

}



