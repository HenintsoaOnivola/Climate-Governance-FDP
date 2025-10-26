# This script has two parts (1) building the composite hazard index based on the individual hazards indices
# (2) investigate on the individual hazard that is mostly driving the high compounding effects and how these hazards are distributed over regions 
rm(list=ls())

packages <- c("DescTools", "raster", "tibble", "ggplot2",'sf','scales','RColorBrewer','exactextractr','geomtextpath','ggpubr')

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

##############################################
#defining all functions and plot formatting
theme_MAP <-   theme_bw() +
  theme(title = element_text(size = 15),
        plot.title=element_text(hjust=0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(size=10,face='bold'), 
        
  )

overlayRast <- function(rst) {
  tbl <- raster::rasterToPoints(rst, spatial = FALSE)
  tbl <- as_tibble(tbl)
  names(tbl) <- c('x', 'y', 'value')
  tbl
}

Visu_raster<-function (output_name,the_min,the_max,title,bar_title,color_palette,all_ticks,all_labels,the_tbl,shape_delimiter){
  all_val<-na.omit(sort(unique(the_tbl$value)))
  my_cols<-color_palette
  if (all(all_val >= -0.01)){
    res<-seq(0,1,1/length(my_cols))
  }else
  {
    my_cols<- c('#767676','#FFFFFF',my_cols)
    positive_val=all_val[all_val>0]
    val_after_zero=positive_val[1]
    scaling_colormap<-c(the_min,0,seq(val_after_zero,the_max,(the_max-val_after_zero)/(length(my_cols)-3)))
    res<-scales::rescale(scaling_colormap)
  }
  tiff(output_name, units="in", width=10, height=6, res=500,compression = 'lzw')
  print(ggplot() +
          geom_tile(data = the_tbl, aes(x = x, y = y, fill = value)) +
          geom_sf(data = shape_delimiter, fill = NA,lwd=0.5) +
          scale_fill_gradientn(colors = my_cols,name=bar_title, limits=c(the_min,the_max),values=res,breaks=all_ticks,labels=all_labels)+   
          labs(title = title) +
          theme_MAP)
  dev.off()
}

Visu_raster_discrete <-function (output_name,title,bar_title,color_palette,the_lim,all_ticks,all_labels,the_tbl,shape_delimiter){
  tiff(output_name, units="in", width=10, height=6, res=500,compression = 'lzw')
  print(ggplot() +
          geom_tile(data = the_tbl, aes(x = x, y = y, fill = factor(value))) +
          geom_sf(data = shape_delimiter, fill = NA,lwd=0.5) +
          scale_fill_manual(values=color_palette, name=bar_title,limits=the_lim,breaks=all_ticks,labels=all_labels) +   
          labs(title = title) +
          theme_MAP)
  dev.off()
}


input_fold<-'./Input/'
output_fold<-'./Output/'
world_shape <- st_read(paste0(input_fold,'world_shape2.geojson'))

ff<-''
##reading individual indices from files
heat_hist_index<-raster(paste0(output_fold,"heat_hist_nrm.tif"))
heat_ssp_index<-raster(paste0(output_fold,"heat_ssp",ff,"_nrm.tif"))

drought_hist_index<-raster(paste0(output_fold,"drought_hist_nrm.tif"))
drought_ssp_index<-raster(paste0(output_fold,"drought_ssp",ff,"_nrm.tif"))

flood_hist_index<-raster(paste0(output_fold,"flood_hist_nrm.tif"))
flood_ssp_index<-raster(paste0(output_fold,"flood_ssp",ff,"_nrm.tif"))


#########################part 1: building the composite hazard index ############
#################################################################################
#build the non-normalized composite index and get the overall min and max of the composite index
ter_hist_ensemble_nrm<-readRDS(paste0(output_fold,'Intermediate/ter_hist_ensemble_nrm.RData'))
ter_ssp_ensemble_nrm<-readRDS(paste0(output_fold,'Intermediate/ter_ssp_ensemble_nrm.RData'))

#baseline
ter_hist_ensemble_nrm<-c(ter_hist_ensemble_nrm,raster::stack())
ter_hist_ensemble_nrm[[4]]<-raster::stack()
all_min_hist<-c()
all_max_hist<-c()
for (y in 1:nlayers(ter_hist_ensemble_nrm[[1]])){
  heat_raster<-ter_hist_ensemble_nrm[[1]][[y]]
  drought_raster<-ter_hist_ensemble_nrm[[2]][[y]]
  flood_raster<-ter_hist_ensemble_nrm[[3]][[y]]
  composite_raster<- log(sqrt((heat_raster^2)+(drought_raster^2)+(flood_raster^2)))
  
  all_min_hist<-c(all_min_hist,min(composite_raster[composite_raster!=-Inf],na.rm=T))
  all_max_hist<-c(all_max_hist,maxValue(composite_raster))
  ter_hist_ensemble_nrm[[4]]<-addLayer(ter_hist_ensemble_nrm[[4]],composite_raster)
}
ter_hist_ensemble_nrm<-ter_hist_ensemble_nrm[1:4]
#future
ter_ssp_ensemble_nrm<-c(ter_ssp_ensemble_nrm,raster::stack())
ter_ssp_ensemble_nrm[[4]]<-raster::stack()
all_min_ssp<-c()
all_max_ssp<-c()
for (y in 1:nlayers(ter_ssp_ensemble_nrm[[1]])){
  heat_raster<-ter_ssp_ensemble_nrm[[1]][[y]]
  drought_raster<-ter_ssp_ensemble_nrm[[2]][[y]]
  flood_raster<-ter_ssp_ensemble_nrm[[3]][[y]]
  composite_raster<- log(sqrt((heat_raster^2)+(drought_raster^2)+(flood_raster^2)))
  all_min_ssp<-c(all_min_ssp,min(composite_raster[composite_raster!=-Inf],na.rm=T))
  all_max_ssp<-c(all_max_ssp,maxValue(composite_raster))
  ter_ssp_ensemble_nrm[[4]]<-addLayer(ter_ssp_ensemble_nrm[[4]],composite_raster)
}
ter_ssp_ensemble_nrm<-ter_ssp_ensemble_nrm[1:4]

all_min<-c(all_min_hist,all_min_ssp)
all_min<-all_min[all_min!=-Inf]
all_max<-c(all_max_hist,all_max_ssp)
mi<-min(all_min)
ma<-max(all_max)
hist_br<-ter_hist_ensemble_nrm[[4]]
hist_br[hist_br==-Inf]<-mi
ter_hist_ensemble_nrm[[4]]<-hist_br
ssp_br<-ter_ssp_ensemble_nrm[[4]]
ssp_br[ssp_br==-Inf]<-mi
ter_ssp_ensemble_nrm[[4]]<-ssp_br

#store into a list, and plot the composite index per year, after normalization
#baseline
the_years<-1995:2014
l_plot<-list(1:12,13:20)
for (ll in 1:length(l_plot)){
  tiff(paste0('./Plots/time_series/composite_hist_',ll,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    ter_hist_ensemble_nrm[[4]][[y]]<-(ter_hist_ensemble_nrm[[4]][[y]]-mi)/(ma-mi)
    plot(ter_hist_ensemble_nrm[[4]][[y]],col=terrain.colors(255),zlim=c(0,1))
    title(the_years[y])
    writeRaster(ter_hist_ensemble_nrm[[4]][[y]], filename = (paste0(output_fold,"Yearly_hazard/composite_hist_nrm_",the_years[y],".tif")), format = "GTiff",overwrite=TRUE)
  }
  dev.off()
}
#future
the_years<-2020:2040
l_plot<-list(1:12,13:21)
for (ll in 1:length(l_plot)){
  tiff(paste0('./Plots/time_series/composite_ssp_',ll,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    ter_ssp_ensemble_nrm[[4]][[y]]<-(ter_ssp_ensemble_nrm[[4]][[y]]-mi)/(ma-mi)
    plot(ter_ssp_ensemble_nrm[[4]][[y]],col=terrain.colors(255),zlim=c(0,1))
    title(the_years[y])
    writeRaster(ter_ssp_ensemble_nrm[[4]][[y]], filename = (paste0(output_fold,"Yearly_hazard/composite_ssp_nrm_",the_years[y],".tif")), format = "GTiff",overwrite=TRUE)
  }
  dev.off()
}
saveRDS(ter_hist_ensemble_nrm,paste0(output_fold,'Intermediate/ter_hist_ensemble_nrm_with_composite.RData'))
saveRDS(ter_ssp_ensemble_nrm,paste0(output_fold,'/Intermediate/ter_ssp_ensemble_nrm_with_composite.RData'))

#######
#plotting yearly side by side with individual hazards
the_years<-1995:2014
l_plot<-seq(1,nlayers(ter_hist_ensemble_nrm[[1]]),2)
for (each_page in l_plot){
  tiff(paste0('./Plots/time_series/all_hist_',each_page,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,2),mar=c(2,3,2,3)+0.1)
  for (y in each_page:(each_page+1)){
    plot(ter_hist_ensemble_nrm[[1]][[y]],col=terrain.colors(255),zlim=c(0,1))
    title(paste0('heat-',the_years[y]))
    plot(ter_hist_ensemble_nrm[[2]][[y]],col=terrain.colors(255),zlim=c(0,1))
    title(paste0('drought-',the_years[y]))
    plot(ter_hist_ensemble_nrm[[3]][[y]],col=terrain.colors(255),zlim=c(0,1))
    title(paste0('flood-',the_years[y]))
    plot(ter_hist_ensemble_nrm[[4]][[y]],col=terrain.colors(255),zlim=c(0,1))
    title(paste0('compound-',the_years[y]))
  }
  dev.off()
}
the_years<-2020:2040
l_plot<-seq(1,nlayers(ter_ssp_ensemble_nrm[[1]]),2)
for (each_page in l_plot){
  tiff(paste0('./Plots/time_series/all_ssp_',each_page,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,2),mar=c(2,3,2,3)+0.1)
  for (y in each_page:(each_page+1)){
    if (y<=l_plot[length(l_plot)]){
      plot(ter_ssp_ensemble_nrm[[1]][[y]],col=terrain.colors(255),zlim=c(0,1))
      title(paste0('heat-',the_years[y]))
      plot(ter_ssp_ensemble_nrm[[2]][[y]],col=terrain.colors(255),zlim=c(0,1))
      title(paste0('drought-',the_years[y]))
      plot(ter_ssp_ensemble_nrm[[3]][[y]],col=terrain.colors(255),zlim=c(0,1))
      title(paste0('flood-',the_years[y]))
      plot(ter_ssp_ensemble_nrm[[4]][[y]],col=terrain.colors(255),zlim=c(0,1))
      title(paste0('compound-',the_years[y]))
    }
    
  }
  dev.off()
}

##compute the composite hazard that represents each period (baseline and future) and plot
composite_hist_index<-stackApply(ter_hist_ensemble_nrm[[4]], indices =  rep(1,nlayers(ter_hist_ensemble_nrm[[4]])), fun = "mean", na.rm = T)
composite_ssp_index<-stackApply(ter_ssp_ensemble_nrm[[4]], indices =  rep(1,nlayers(ter_ssp_ensemble_nrm[[4]])), fun = "mean", na.rm = T)

writeRaster(composite_hist_index, filename = (paste0(output_fold,"composite_hist_nrm.tif")), format = "GTiff",overwrite=TRUE)
writeRaster(composite_ssp_index, filename = (paste0(output_fold,"composite_ssp_nrm.tif")), format = "GTiff",overwrite=TRUE)


#plotting the spatial distribution of compound hazard
composite_change_index<-composite_ssp_index-composite_hist_index
composite_hist_tbl<- overlayRast(composite_hist_index)
composite_ssp_tbl <- overlayRast(composite_ssp_index)
composite_change_tbl <- overlayRast(composite_change_index)

the_min<-min(c(minValue(composite_hist_index),minValue(composite_ssp_index)))
the_max<-max(c(maxValue(composite_hist_index),maxValue(composite_ssp_index)))
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c(the_ticks[1],rep('',length(the_ticks)-2),the_ticks[length(the_ticks)])
col_map<-rev(brewer.pal(11,'Spectral'))
Visu_raster("./plots/composite_hist.tiff",color_bar_min,color_bar_max,"composite hazard - baseline",'composite index',col_map,the_ticks,the_labels,composite_hist_tbl,world_shape)
Visu_raster("./plots/composite_ssp.tiff",color_bar_min,color_bar_max,"composite hazard - future",'composite index',col_map,the_ticks,the_labels,composite_ssp_tbl,world_shape)
Visu_raster("./plots/composite_change.tiff",minValue(composite_change_index),maxValue(composite_change_index),"composite hazard - change",'composite change',col_map,waiver(),waiver(),composite_change_tbl,world_shape)

# build the 5 different severity categories and plot the resulting maps
brk  <- seq(0,1,0.2)
composite_hist_class <- raster::cut(composite_hist_index, breaks=brk) 
composite_hist_class_tbl<-overlayRast(composite_hist_class)

composite_ssp_class <- raster::cut(composite_ssp_index, breaks=brk) 
composite_ssp_class_tbl<-overlayRast(composite_ssp_class)

spect_pal<-rev(brewer.pal(5,'Spectral'))
my_pals=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])
Visu_raster_discrete("./plots/composite_hist_class.tiff","composite hazard - baseline",'hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),composite_hist_class_tbl,world_shape)
Visu_raster_discrete("./plots/composite_ssp_class.tiff","composite hazard - future",'hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),composite_ssp_class_tbl,world_shape)

#saving the composite hazard class into rasters
writeRaster(composite_hist_class, filename = paste0(output_fold,"composite_hist_class.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(composite_ssp_class, filename = paste0(output_fold,"composite_ssp",ff,"_class.tif"), format = "GTiff",overwrite=TRUE)


#######################  part 2: investigate on the individual hazard that is mostly driving the high compounding effects 
#          and how these hazards are distributed over regions 
######################

#plot the individual hazards under the high, severe and extreme compound hazard 

comp_hist_from_high<-composite_hist_class
comp_hist_from_high[comp_hist_from_high<=2]<-NA
comp_ssp_from_high<-composite_ssp_class
comp_ssp_from_high[comp_ssp_from_high<=2]<-NA

the_point<-0.4

high_ind_hist <- comp_hist_from_high
high_ind_hist[,]<-0
high_ind_hist[heat_hist_index>=the_point & drought_hist_index<the_point & flood_hist_index<the_point & !is.na(comp_hist_from_high)]<-1
high_ind_hist[heat_hist_index<the_point & drought_hist_index>=the_point & flood_hist_index<the_point & !is.na(comp_hist_from_high)]<-2
high_ind_hist[heat_hist_index<the_point & drought_hist_index<the_point & flood_hist_index>=the_point & !is.na(comp_hist_from_high)]<-3
high_ind_hist[heat_hist_index>=the_point & drought_hist_index>=the_point & flood_hist_index<the_point & !is.na(comp_hist_from_high)]<-4
high_ind_hist[heat_hist_index>=the_point & drought_hist_index<the_point & flood_hist_index>=the_point & !is.na(comp_hist_from_high)]<-5
high_ind_hist[heat_hist_index<the_point & drought_hist_index>=the_point & flood_hist_index>=the_point & !is.na(comp_hist_from_high)]<-6
high_ind_hist[heat_hist_index>=the_point & drought_hist_index>=the_point & flood_hist_index>=the_point & !is.na(comp_hist_from_high)]<-7
high_ind_hist[heat_hist_index>=the_point & drought_hist_index>=the_point & flood_hist_index>=the_point & !is.na(comp_hist_from_high)]<-7
high_ind_hist[heat_hist_index<the_point & drought_hist_index<the_point & flood_hist_index<the_point & !is.na(comp_hist_from_high)  ]<-7


high_ind_hist_tbl<-overlayRast(high_ind_hist)


high_ind_ssp <- comp_ssp_from_high
high_ind_ssp[,]<-0
high_ind_ssp[heat_ssp_index>=the_point & drought_ssp_index<the_point & flood_ssp_index<the_point & !is.na(comp_ssp_from_high)]<-1
high_ind_ssp[heat_ssp_index<the_point & drought_ssp_index>=the_point & flood_ssp_index<the_point & !is.na(comp_ssp_from_high)]<-2
high_ind_ssp[heat_ssp_index<the_point & drought_ssp_index<the_point & flood_ssp_index>=the_point & !is.na(comp_ssp_from_high)]<-3
high_ind_ssp[heat_ssp_index>=the_point & drought_ssp_index>=the_point & flood_ssp_index<the_point & !is.na(comp_ssp_from_high)]<-4
high_ind_ssp[heat_ssp_index>=the_point & drought_ssp_index<the_point & flood_ssp_index>=the_point & !is.na(comp_ssp_from_high)]<-5
high_ind_ssp[heat_ssp_index<the_point & drought_ssp_index>=the_point & flood_ssp_index>=the_point & !is.na(comp_ssp_from_high)]<-6
high_ind_ssp[heat_ssp_index>=the_point & drought_ssp_index>=the_point & flood_ssp_index>=the_point & !is.na(comp_ssp_from_high)]<-7
high_ind_ssp[heat_ssp_index<the_point & drought_ssp_index<the_point & flood_ssp_index<the_point & !is.na(comp_ssp_from_high)  ]<-7


high_ind_ssp_tbl<-overlayRast(high_ind_ssp)

my_pals=c('0'='white','1'='orangered1','2'='goldenrod3','3'='steelblue3','4'='saddlebrown','5'='olivedrab1','6'='purple','7'='grey')

Visu_raster_discrete("./Plots/high_to_extreme_compound_hist.tiff","high to extreme hazard - baseline",'hazard type',my_pals,c('0','1','2','3','4','5','6','7'),c('0','1','2','3','4','5','6','7'),c('','heat','drought','flood','heat and drougth','heat and flood','drought and flood','all'),high_ind_hist_tbl,world_shape)
Visu_raster_discrete("./Plots/high_to_extreme_compound_ssp.tiff","high to extreme hazard - future",'hazard type',my_pals,c('0','1','2','3','4','5','6','7'),c('0','1','2','3','4','5','6','7'),c('','heat','drought','flood','heat and drougth','heat and flood','drought and flood','all'),high_ind_ssp_tbl,world_shape)

#plot the distribution of underlying driving hazard per region
##baseline
regional_files<-c('AMERICAs_15m_UNHCR.gpkg',"ASIA_15m_UNHCR.gpkg","East_and_horn_of_Africa_shapefile.gpkg",
                  "Europe_new.gpkg","MENA_polygon_2.gpkg","Southern_Africa_Ash.gpkg","West_and_central_new.gpkg")

all_region<-c('Americas','Asia','East and\nHorn of Africa','Europe','MENA','Southern\nAfrica','West and\nCentral Africa')
all_types=c('heat','drought','flood','heat and drougth','heat and flood','drought and flood','all')

percentage=c()
type=c()
region=c()
for (count_region in 1:length(regional_files)){
  polygon_rst <- st_read(paste0(input_fold,"spatial_data/",regional_files[count_region]))
  r <- raster::crop(high_ind_hist, polygon_rst)
  rst_ini <- raster::mask(r, polygon_rst)
  rst_all_values<-rst_ini
  rst_all_values[rst_ini>=0]=1
  rst_all_high<-rst_ini
  rst_all_high[rst_ini>=1]=1
  all_type_high<-cellStats(rst_all_high, stat='sum', na.rm=TRUE)/cellStats(rst_all_values, stat='sum', na.rm=TRUE)
  for (type_ind in (1:length(all_types))){
    type_rast<-rst_ini
    type_rast[rst_ini!=type_ind]=0
    type_rast[type_rast==type_ind]=1
    proportion_type<-cellStats(type_rast, stat='sum', na.rm=TRUE)/cellStats(rst_all_values, stat='sum', na.rm=TRUE)
    
    percentage<-c(percentage,proportion_type*100)
    type<-c(type,all_types[type_ind])
    region<-c(region,all_region[count_region])
  }
}
all_colors=c('heat'='orangered1','drought'='goldenrod3','flood'='steelblue3','heat and drougth'='saddlebrown','heat and flood'='olivedrab1','drought and flood'='purple','all'='grey')
the_df<-data.frame(percentage,type,region)
the_df$region<-factor(the_df$region)
the_df2<-the_df[the_df$percentage!=0,]
tiff('./Plots/from_high_per_region_hist.tiff', units="in", width=11, height=6, res=500,compression = 'lzw')
ggplot(the_df,aes(x=region,y=percentage,fill=type))+
  geom_bar(stat='identity')+
  labs(y='area under high, severe and extreme (%)')+
  scale_fill_manual(values=all_colors,breaks=unique(the_df2$type),labels=all_types[match(unique(the_df2$type),all_types)])+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  ylim(0,100)+
  ggtitle('baseline')+
  theme(axis.text=element_text(size=16),axis.title = element_text(size = 20),plot.title = element_text(hjust = 0.1,size=25),legend.text = element_text(size = 15),legend.title = element_text(size = 18))
dev.off()

#future
percentage=c()
type=c()
region=c()
for (count_region in 1:length(regional_files)){
  polygon_rst <- st_read(paste0(input_fold,"spatial_data/",regional_files[count_region]))
  r <- raster::crop(high_ind_ssp, polygon_rst)
  rst_ini <- raster::mask(r, polygon_rst)
  rst_all_values<-rst_ini
  rst_all_values[rst_ini>=0]=1
  rst_all_high<-rst_ini
  rst_all_high[rst_ini>=1]=1
  all_type_high<-cellStats(rst_all_high, stat='sum', na.rm=TRUE)/cellStats(rst_all_values, stat='sum', na.rm=TRUE)
  for (type_ind in (1:length(all_types))){
    type_rast<-rst_ini
    type_rast[rst_ini!=type_ind]=0
    type_rast[type_rast==type_ind]=1
    proportion_type<-cellStats(type_rast, stat='sum', na.rm=TRUE)/cellStats(rst_all_values, stat='sum', na.rm=TRUE)
    
    percentage<-c(percentage,proportion_type*100)
    type<-c(type,all_types[type_ind])
    region<-c(region,all_region[count_region])
  }
}
the_df<-data.frame(percentage,type,region)
the_df$region<-factor(the_df$region)
the_df2<-the_df[the_df$percentage!=0,]

tiff('./Plots/from_high_per_region_ssp.tiff', units="in", width=11, height=6, res=500,compression = 'lzw')
ggplot(the_df,aes(x=region,y=percentage,fill=type))+
  geom_bar(stat='identity')+
  labs(y='area under high, severe and extreme (%)')+
  scale_fill_manual(values=all_colors,breaks=unique(the_df2$type),labels=all_types[match(unique(the_df2$type),all_types)])+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  ylim(0,100)+
  ggtitle('future')+
  theme(axis.text=element_text(size=16),axis.title = element_text(size = 20),plot.title = element_text(hjust = 0.1,size=25),legend.text = element_text(size = 15),legend.title = element_text(size = 18))

dev.off()


#######
###side visualization of the composite index per country
#average the composite hazard per country 
sp_data_hist<-exact_extract(composite_hist_index, world_shape,fun='mean') 
world_shape$hist_average<-sp_data_hist
the_class<-cut(sp_data_hist,breaks=seq(0,1,0.2),labels=1:5)
the_class<-as.numeric(levels(the_class))[the_class] 
world_shape$hist_class<-the_class

sp_data_ssp<-exact_extract(composite_ssp_index, world_shape,fun='mean') 
world_shape$ssp_average<-sp_data_ssp
the_class<-cut(sp_data_ssp,breaks=seq(0,1,0.2),labels=1:5)
the_class<-as.numeric(levels(the_class))[the_class] 
world_shape$ssp_class<-the_class

spect_pal<-rev(brewer.pal(5,'Spectral'))
tiff("./Plots/composite_hist_class_avrg.tiff", units="in", width=10, height=6, res=500,compression = 'lzw')
plot(world_shape['hist_class'],pal=spect_pal,breaks=0:5)
dev.off()
tiff("./Plots/composite_ssp_class_avrg.tiff", units="in", width=10, height=6, res=500,compression = 'lzw')
plot(world_shape['ssp_class'],pal=spect_pal,breaks=0:5)
dev.off()

#Euclidiean norm vs average
EuclideanNorm <- function(x) {
  if (length(which(is.na(x)))==length(x)){
    val<-NA
  }else{
    val<-sqrt(sum(x^2,na.rm=T))
  }
  val
}
par(mfrow=c(1,2))
x<-seq(0.01,1,0.01)
y<-x
rr<-length(x)
z<-matrix(NA,length(x),length(y))
z_raster<-matrix(NA,length(x),length(y))
for (xx in x){
  for (yy in y){
    z[which(x==xx),which(y==yy)]<-EuclideanNorm(c(xx,yy))
    z_raster[rr+1-which(x==xx),which(y==yy)]<-EuclideanNorm(c(xx,yy))
  }
}
z_tri<-matrix(NA,length(x),length(y))
z_tri_raster<-matrix(NA,length(x),length(y))
for (xx in x){
  for (yy in y){
    z_tri[which(x==xx),which(y==yy)]<-(z[which(x==xx),which(y==yy)] - min(z))/(max(z)-min(z))
    z_tri_raster[which(x==xx),which(y==yy)]<-(z_raster[which(x==xx),which(y==yy)] - min(z))/(max(z)-min(z))
    
  }
}
z_tri_frame<-as_tibble(raster::rasterToPoints(raster( z_tri_raster), spatial = FALSE))
names(z_tri_frame) <- c('x', 'y', 'value')
p1 <- ggplot(z_tri_frame, aes(x, y, z = value))+
  geom_raster(aes(fill = value))+
  geom_contour(color='black')+
  geom_textcontour(aes(z = value),hjust = 0.5) +
  scale_fill_gradient(low = "yellow", high = "red")+
  theme(
    panel.background = element_blank(),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    axis.line = element_line(colour = "black"),
    text = element_text(size = 20)
  )

z_mean<-matrix(NA,length(x),length(y))
z_mean_raster<-matrix(NA,length(x),length(y))
for (xx in x){
  for (yy in y){
    z_mean[which(x==xx),which(y==yy)]<-mean(c(xx,yy))
    z_mean_raster[rr+1-which(x==xx),which(y==yy)]<-mean(c(xx,yy))
  }
}
z_mean_frame<-as_tibble(raster::rasterToPoints(raster( z_mean_raster), spatial = FALSE))
names(z_mean_frame) <- c('x', 'y', 'value')
p2 <- ggplot(z_mean_frame, aes(x, y, z = value))+
  geom_raster(aes(fill = value))+
  geom_contour(color='black')+
  geom_textcontour(aes(z = value),hjust = 0.5) +
  scale_fill_gradient(low = "yellow", high = "red")+
  theme(
    panel.background = element_blank(),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(colour = "black"),
    text = element_text(size = 20)
  )
tiff("./Plots/mean_vs_euclidean.tiff", units="in", width=11, height=5, res=500,compression = 'lzw')
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="right")
dev.off()


