# This script generates results under SSP126 for both individual climate hazards and composite hazard
rm(list=ls())


packages <- c("raster", "tidyverse", 'tibble',"tibble", "ggplot2",'sf','RColorBrewer','colorspace','DescTools',
              'cowplot','exactextractr','readxl','ggrepel','ggnewscale')

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


setwd(dirname(rstudioapi::getSourceEditorContext()$path))


##############################################
#defining all functions and plot formatting
##################################################


overlayRast <- function(rst) {
  tbl <- raster::rasterToPoints(rst, spatial = FALSE)
  tbl <- as_tibble(tbl)
  names(tbl) <- c('x', 'y', 'value')
  tbl
}


formatRast <- function(rst,shape_delimiter) {
  r_rotated <- raster::rotate(rst)
  r_croped <- raster::crop(r_rotated, shape_delimiter)
  rst <- raster::mask(r_croped, shape_delimiter)
  rst
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

Visu_raster_discrete_with_points<-function (output_name,title,bar_title,color_palette,the_lim,all_ticks,all_labels,the_tbl,the_points,shape_delimiter){
  tiff(output_name, units="in", width=10, height=6, res=500,compression = 'lzw')
  
  print(
    ggplot() +
      geom_tile(data = the_tbl, aes(x = x, y = y, fill = factor(value))) +
      geom_sf(data = shape_delimiter, fill = NA,lwd=0.5) +
      scale_fill_manual(values=color_palette, name=bar_title,limits=the_lim,breaks=all_ticks,labels=all_labels) +
      new_scale_fill()+  
      geom_sf(data=the_points,size=0.5,aes(fill=new_col))+
      scale_fill_manual(name='',values=c('black'),labels='Location of forcibly displaced people')+
      guides(fill = guide_legend(override.aes = list(size = 2))) +
      labs(title = title) +
      theme_MAP)
  dev.off()
}


#input and pre-process individual hazards
#####################################################
input_fold<-'./Input/'
output_fold<-'./Output/'
world_shape <- st_read(paste0(input_fold,'world_shape2.geojson'))


rast_list=vector("list", length = 6)
rast_list[[1]]<-readRDS(paste0(output_fold,'./Intermediate/hist_ensemble_heat.RData'))
rast_list[[2]]<-readRDS(paste0(output_fold,'./Intermediate/hist_ensemble_drought.RData'))
rast_list[[3]]<-readRDS(paste0(output_fold,'./Intermediate/hist_ensemble_flood.RData'))
rast_list[[4]]<-readRDS(paste0(output_fold,'./Intermediate/ssp126_ensemble_heat.RData'))
rast_list[[5]]<-readRDS(paste0(output_fold,'./Intermediate/ssp126_ensemble_drought.RData'))
rast_list[[6]]<-readRDS(paste0(output_fold,'./Intermediate/ssp126_ensemble_flood.RData'))

#consider only terrestrial data (crop and mask marine data)
ter_ensemble<-list()
for (haz in 1:6){
  ter_ensemble[[haz]]<-formatRast(rast_list[[haz]],world_shape)
}
ter_hist_ensemble<-ter_ensemble[1:3]
ter_ssp_ensemble<-ter_ensemble[4:6]

#remove values in greenland
gr<-world_shape[which(world_shape$gis_name=='Greenland (DNK)'),]
for (h in 1:length(ter_hist_ensemble)){
  ter_hist_ensemble[[h]]<-mask(ter_hist_ensemble[[h]],gr,inverse=TRUE)
  ter_ssp_ensemble[[h]]<-mask(ter_ssp_ensemble[[h]],gr,inverse=TRUE)
}


##transforming raw indices to be comparable (normalization and distribution shift)
#######################################################################
##heat
#############
extreme_heat<-69.75972 #value defined from the standard analysis (baseline), same approach applies to drought and flood
dist_heat<-as.vector(values(ter_hist_ensemble[[1]]))
dist_heat_ssp<-as.vector(values(ter_ssp_ensemble[[1]]))
#normalization
min_heat=-0.4394801  #values defined from the standard analysis (baseline and ssp585) , same approach applies to drought and flood
max_heat=148.6373448
dist_heat_nrm<-(dist_heat-min_heat)/(max_heat-min_heat)
extreme_heat_nrm<-(extreme_heat-min_heat)/(max_heat-min_heat)
dist_heat_ssp_nrm<-(dist_heat_ssp-min_heat)/(max_heat-min_heat)
dist_heat_ssp_nrm[which(dist_heat_ssp_nrm>1)]<-1
#shifting
dist_heat_sft<-dist_heat_nrm+(0.8-extreme_heat_nrm)
dist_heat_sft[which(dist_heat_sft>1)]<-1
dist_heat_sft[which(dist_heat_sft<0)]<-0
dist_heat_sft[which(dist_heat<1)]<-0.1 
dist_heat_ssp_sft<-dist_heat_ssp_nrm+(0.8-extreme_heat_nrm)
dist_heat_ssp_sft[which(dist_heat_ssp_sft>1)]<-1
dist_heat_ssp_sft[which(dist_heat_ssp_sft<0)]<-0
dist_heat_ssp_sft[which(dist_heat_ssp<1)]<-0.1 
#####################
##drought
#####################
dist_drought<-as.vector(values(ter_hist_ensemble[[2]]))
dist_drought[dist_drought>0]<-0 #only taking the negative ones
dist_drought<--dist_drought
dist_drought_ssp<-as.vector(values(ter_ssp_ensemble[[2]]))
dist_drought_ssp[dist_drought_ssp>0]<-0 #only taking the negative ones
dist_drought_ssp<--dist_drought_ssp
extreme_drought=1.94869 
#normalization
min_drought=0
max_drought=4.366562
dist_drought_nrm<-(dist_drought-min_drought)/(max_drought-min_drought)
extreme_drought_nrm<-(extreme_drought-min_drought)/(max_drought-min_drought) 
dist_drought_ssp_nrm<-(dist_drought_ssp-min_drought)/(max_drought-min_drought)
dist_drought_ssp_nrm[dist_drought_ssp_nrm>1]<-1
#shifting
dist_drought_sft<-dist_drought_nrm+(0.8-extreme_drought_nrm)
dist_drought_sft[which(dist_drought_sft>1)]<-1
dist_drought_sft[which(dist_drought_sft<0)]<-0
dist_drought_sft[which(dist_drought==0)]=0.1
dist_drought_ssp_sft<-dist_drought_ssp_nrm+(0.8-extreme_drought_nrm)
dist_drought_ssp_sft[which(dist_drought_ssp_sft>1)]<-1
dist_drought_ssp_sft[which(dist_drought_ssp_sft<0)]<-0
dist_drought_ssp_sft[which(dist_drought_ssp==0)]=0.1
########################
#flood
############################
dist_flood<-as.vector(values(ter_hist_ensemble[[3]]))
dist_flood_ssp<-as.vector(values(ter_ssp_ensemble[[3]]))
extreme_flood=14.78567
#normalization
min_flood=1.044164
max_flood=19.173838
dist_flood_nrm<-(dist_flood-min_flood)/(max_flood-min_flood)
extreme_flood_nrm<-(extreme_flood-min_flood)/(max_flood-min_flood) 
dist_flood_ssp_nrm<-(dist_flood_ssp-min_flood)/(max_flood-min_flood)
dist_flood_ssp_nrm[dist_flood_ssp_nrm>1]<-1
#shifting
dist_flood_sft<-dist_flood_nrm+(0.8-extreme_flood_nrm)
dist_flood_sft[which(dist_flood_sft>1)]<-1
dist_flood_sft[which(dist_flood_sft<0)]<-0
dist_flood_ssp_sft<-dist_flood_ssp_nrm+(0.8-extreme_flood_nrm)
dist_flood_ssp_sft[which(dist_flood_ssp_sft>1)]<-1
dist_flood_ssp_sft[which(dist_flood_ssp_sft<0)]<-0

#plotting raw hazard distribution
tiff(paste0('./Plots/raw_hazard_distribution_ssp126.tiff'),units="in", width=7, height=4, res=300,compression = 'lzw')
par(mfrow=c(2,3))

hist(dist_heat[dist_heat>=1],xlab='number of days with HI>41°C',ylab='frequency',breaks = "Scott",main='',xlim=c(min_heat,200),ylim=c(0,50000))
title(main = "heat",adj = 0)
abline(v=extreme_heat,col='red',lwd=1.5,lty=2)
hist(dist_drought[dist_drought>0],xlab='-SPEI',ylab='frequency',breaks=75,main='',xlim=c(min_drought,max(dist_drought_ssp,na.rm=T)))
title(main = "drought",adj = 0)
abline(v=extreme_drought,col='red',lwd=1.5,lty=2)
mtext("Baseline", side = 3, line = 3, outer =FALSE)
hist(dist_flood,xlab='SDII',ylab='frequency',breaks = "Scott",main='',xlim=c(min_flood,max_flood))
title(main = "flood",adj = 0)
abline(v=extreme_flood,col='red',lwd=1.5,lty=2)
legend('topright',inset=c(-0.1, -0.5), legend="threshold for extreme",
       col=c("red"), lty=2, cex=0.8,xpd=TRUE)


hist(dist_heat_ssp[dist_heat_ssp>=1],xlab='number of days with HI>41°C',ylab='frequency',breaks = "Scott",main='',xlim=c(min_heat,200),ylim=c(0,50000))
title(main = "heat",adj = 0)
abline(v=extreme_heat,col='red',lwd=1.5,lty=2)
hist(dist_drought_ssp[dist_drought_ssp>0],xlab='-SPEI',ylab='frequency',breaks=75,main='',xlim=c(min_drought,max(dist_drought_ssp,na.rm=T)))
title(main = "drought",adj = 0)
abline(v=extreme_drought,col='red',lwd=1.5,lty=2)
mtext("Future", side = 3, line = 3, outer =FALSE)
hist(dist_flood_ssp,xlab='SDII',ylab='frequency',breaks = "Scott",main='',xlim=c(min_flood,max_flood))
title(main = "flood",adj = 0)
abline(v=extreme_flood,col='red',lwd=1.5,lty=2)
dev.off()

#plotting hazard distribution after normalizing and shifting
tiff(paste0('./Plots/normalized_shifted_hazard_distribution_ssp126.tiff'),units="in", width=7, height=4, res=300,compression = 'lzw')
par(mfrow=c(2,3))

hist(dist_heat_sft,xlab='heat index',ylab='frequency',breaks = "Scott",main='',xlim=c(0,1))
title(main = "heat",adj = 0)
abline(v=seq(0.2,0.8,0.2),col='blue',lwd=1.5,lty=2)
hist(dist_drought_sft,xlab='drought index',ylab='frequency',breaks=75,main='',xlim=c(0,1))
title(main = "drought",adj = 0)
abline(v=seq(0.2,0.8,0.2),col='blue',lwd=1.5,lty=2)
mtext("Baseline", side = 3, line = 3, outer =FALSE)
hist(dist_flood_sft,xlab='flood index',ylab='frequency',breaks = "Scott",main='',xlim=c(0,1))
title(main = "flood",adj = 0)
abline(v=seq(0.2,0.8,0.2),col='blue',lwd=1.5,lty=2)
legend('topright',inset=c(-0.1, -0.5), legend=c("class threshold"),
       col='blue', lty=2, cex=0.8,xpd=TRUE)

hist(dist_heat_ssp_sft,xlab='heat index',ylab='frequency',breaks = "Scott",main='',xlim=c(0,1))
title(main = "heat",adj = 0)
abline(v=seq(0.2,0.8,0.2),col='blue',lwd=1.5,lty=2)
hist(dist_drought_ssp_sft,xlab='drought index',ylab='frequency',breaks='Scott',main='',xlim=c(0,1))
title(main = "drought",adj = 0)
abline(v=seq(0.2,0.8,0.2),col='blue',lwd=1.5,lty=2)
mtext("Future", side = 3, line = 3, outer =FALSE)
hist(dist_flood_ssp_sft,xlab='flood index',ylab='frequency',breaks = "Scott",main='',xlim=c(0,1))
title(main = "flood",adj = 0)
abline(v=seq(0.2,0.8,0.2),col='blue',lwd=1.5,lty=2)
dev.off()


###########################
#store the final indices in raster bricks
all_years<-nlayers(ter_hist_ensemble[[1]])
yearly_index=seq(0,length(dist_heat_sft),ncell(ter_hist_ensemble[[1]][[1]]))
yearly_index=cbind(c((yearly_index[1:length(yearly_index)-1])+1),yearly_index[2:length(yearly_index)])

ter_hist_ensemble_nrm<-list(raster::stack(),raster::stack(),raster::stack())
#hist
l_plot<-list(1:12,13:20)
the_years<-1995:2014
for (ll in 1:length(l_plot)){
  for (y in l_plot[[ll]]){
    heat_val<-dist_heat_sft[yearly_index[y,1]:yearly_index[y,2]]
    heat_raster<-ter_hist_ensemble[[1]][[1]]
    values(heat_raster)<-heat_val
    ter_hist_ensemble_nrm[[1]]<-addLayer(ter_hist_ensemble_nrm[[1]],heat_raster)
  }

  for (y in l_plot[[ll]]){
    drought_val<-dist_drought_sft[yearly_index[y,1]:yearly_index[y,2]]
    drought_raster<-ter_hist_ensemble[[2]][[1]]
    values(drought_raster)<-drought_val
    ter_hist_ensemble_nrm[[2]]<-addLayer(ter_hist_ensemble_nrm[[2]],drought_raster)
  }

  for (y in l_plot[[ll]]){
    flood_val<-dist_flood_sft[yearly_index[y,1]:yearly_index[y,2]]
    flood_raster<-ter_hist_ensemble[[2]][[1]]
    values(flood_raster)<-flood_val
    ter_hist_ensemble_nrm[[3]]<-addLayer(ter_hist_ensemble_nrm[[3]],flood_raster)
  }

}
#ssp
all_years<-nlayers(ter_ssp_ensemble[[1]])
yearly_index=seq(0,length(dist_heat_ssp_sft),ncell(ter_ssp_ensemble[[1]][[1]]))
yearly_index=cbind(c((yearly_index[1:length(yearly_index)-1])+1),yearly_index[2:length(yearly_index)])

the_years<-2020:2040
ter_ssp_ensemble_nrm<-list(raster::stack(),raster::stack(),raster::stack())
l_plot<-list(1:12,13:21)
for (ll in 1:length(l_plot)){
  tiff(paste0('./Plots/time_series/heat_ssp126_',ll,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    heat_val<-dist_heat_ssp_sft[yearly_index[y,1]:yearly_index[y,2]]
    heat_raster<-ter_ssp_ensemble[[1]][[1]]
    values(heat_raster)<-heat_val
    ter_ssp_ensemble_nrm[[1]]<-addLayer(ter_ssp_ensemble_nrm[[1]],heat_raster)
    plot(heat_raster,col=terrain.colors(255),zlim=c(0,1))
    title(the_years[y])
    writeRaster(heat_raster, filename = (paste0(output_fold,"Yearly_hazard/heat_ssp126_nrm_",the_years[y],".tif")), format = "GTiff",overwrite=TRUE)
  }
  dev.off()
  
  tiff(paste0('./Plots/time_series/drought_ssp126_',ll,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    drought_val<-dist_drought_ssp_sft[yearly_index[y,1]:yearly_index[y,2]]
    drought_raster<-ter_ssp_ensemble[[2]][[1]]
    values(drought_raster)<-drought_val
    ter_ssp_ensemble_nrm[[2]]<-addLayer(ter_ssp_ensemble_nrm[[2]],drought_raster)
    plot(drought_raster,col=terrain.colors(255),zlim=c(0,1))
    title(the_years[y])
    writeRaster(heat_raster, filename = (paste0(output_fold,"Yearly_hazard/drought_ssp126_nrm_",the_years[y],".tif")), format = "GTiff",overwrite=TRUE)
  }
  dev.off()
  
  tiff(paste0('./Plots/time_series/flood_ssp126_',ll,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    flood_val<-dist_flood_ssp_sft[yearly_index[y,1]:yearly_index[y,2]]
    flood_raster<-ter_ssp_ensemble[[3]][[1]]
    values(flood_raster)<-flood_val
    ter_ssp_ensemble_nrm[[3]]<-addLayer(ter_ssp_ensemble_nrm[[3]],flood_raster)
    plot(flood_raster,col=terrain.colors(255),zlim=c(0,1))
    title(the_years[y])
    writeRaster(heat_raster, filename = (paste0(output_fold,"Yearly_hazard/drought_ssp126_nrm_",the_years[y],".tif")), format = "GTiff",overwrite=TRUE)
  }
  dev.off()
  
}

#####################################
##compute the individual hazard that represents each period (baseline and future) and plot
########################################
avr_hist<-raster::stack()
avr_ssp<-raster::stack()
for (h in 1:3){
  the_m<-stackApply(ter_hist_ensemble_nrm[[h]], indices =  rep(1,nlayers(ter_hist_ensemble_nrm[[h]])), fun = "mean", na.rm = T)
  the_m_ssp<-stackApply(ter_ssp_ensemble_nrm[[h]], indices =  rep(1,nlayers(ter_ssp_ensemble_nrm[[h]])), fun = "mean", na.rm = T)
  avr_hist<-addLayer(avr_hist,the_m)
  avr_ssp<-addLayer(avr_ssp,the_m_ssp)
}
heat_hist_index<-avr_hist[[1]]
drought_hist_index<-avr_hist[[2]]
flood_hist_index<-avr_hist[[3]]
heat_ssp_index<-avr_ssp[[1]]
drought_ssp_index<-avr_ssp[[2]]
flood_ssp_index<-avr_ssp[[3]]


writeRaster(heat_ssp_index, filename = (paste0(output_fold,"heat_ssp126_nrm.tif")), format = "GTiff",overwrite=TRUE)
writeRaster(drought_ssp_index, filename = (paste0(output_fold,"drought_ssp126_nrm.tif")), format = "GTiff",overwrite=TRUE)
writeRaster(flood_ssp_index, filename = (paste0(output_fold,"flood_ssp126_nrm.tif")), format = "GTiff",overwrite=TRUE)

##visualizing individual climate hazard indices (current, future, change)
heat_ssp5_index<-raster(paste0(output_fold,"heat_ssp_nrm.tif"))
drought_ssp5_index<-raster(paste0(output_fold,"drought_ssp_nrm.tif"))
flood_ssp5_index<-raster(paste0(output_fold,"flood_ssp_nrm.tif"))

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
        legend.title = element_text(size=15,face='bold'),
        legend.text = element_text(size=13)
        
  )

##visualizing heat hazard
heat_ssp_tbl <- overlayRast(heat_ssp_index)
heat_change_index<-heat_ssp_index-heat_hist_index
heat_change_tbl<-overlayRast(heat_change_index)
the_min<-min(c(minValue(heat_hist_index),minValue(heat_ssp_index),minValue(heat_ssp5_index)))
the_max<-max(c(maxValue(heat_hist_index),maxValue(heat_ssp_index),maxValue(heat_ssp5_index)))
col_map<-rev(sequential_hcl(27, palette = "Reds3"))
col_map<-col_map[2:length(col_map)]
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c('low',rep('',length(the_ticks)-2),'high')
Visu_raster(paste0("./Plots/heat_ssp126_index.tiff"),color_bar_min,color_bar_max,"heat hazard - future",'heat hazard',col_map,the_ticks,the_labels,heat_ssp_tbl,world_shape)
col_map<-rev(brewer.pal(11,'Spectral')) 
Visu_raster(paste0("./Plots/heat_change_index126.tiff"),minValue(heat_change_index),maxValue(heat_change_index),"heat hazard - change",'heat change',col_map,waiver(),waiver(),heat_change_tbl,world_shape)

##visualizing drought hazard
drought_ssp_tbl <- overlayRast(drought_ssp_index)
drought_change_index<-drought_ssp_index-drought_hist_index
drought_change_tbl<-overlayRast(drought_change_index)
the_min<-min(c(minValue(drought_hist_index),minValue(drought_ssp_index),minValue(drought_ssp5_index)))
the_max<-max(c(maxValue(drought_hist_index),maxValue(drought_ssp_index),maxValue(drought_ssp5_index)))
col_map<-rev(sequential_hcl(40, palette = "YlOrBr"))
col_map<-c(rep(col_map[1],3),rep(col_map[2],4),rep(col_map[3],5),col_map[4:length(col_map)])
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c('low',rep('',length(the_ticks)-2),'high')
Visu_raster(paste0("./Plots/drought_ssp126_index.tiff"),color_bar_min,color_bar_max,"drought hazard - future",'drought hazard',col_map,the_ticks,the_labels,drought_ssp_tbl,world_shape)
col_map<-rev(brewer.pal(11,'Spectral')) 
Visu_raster(paste0("./Plots/drought_change_index126.tiff"),minValue(drought_change_index),maxValue(drought_change_index),"drought hazard - change",'drought change',col_map,waiver(),waiver(),drought_change_tbl,world_shape)

##visualizing flood hazard
flood_ssp_tbl <- overlayRast(flood_ssp_index)
flood_change_index<-flood_ssp_index-flood_hist_index
flood_change_tbl<-overlayRast(flood_change_index)
the_min<-min(c(minValue(flood_hist_index),minValue(flood_ssp_index),minValue(flood_ssp5_index)))
the_max<-max(c(maxValue(flood_hist_index),maxValue(flood_ssp_index),maxValue(flood_ssp5_index)))
col_map<-rev(sequential_hcl(40, palette = "Blues3"))
col_map<-c(rep(col_map[1],3),rep(col_map[2],4),rep(col_map[3],5),col_map[4:length(col_map)])
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c('low',rep('',length(the_ticks)-2),'high')
Visu_raster(paste0("./Plots/flood_ssp126_index.tiff"),color_bar_min,color_bar_max,"flood hazard - future",'flood hazard',col_map,the_ticks,the_labels,flood_ssp_tbl,world_shape)
col_map<-rev(brewer.pal(11,'Spectral')) 
Visu_raster(paste0("./Plots/flood_change_index126.tiff"),minValue(flood_change_index),maxValue(flood_change_index),"flood hazard - change",'flood change',col_map,waiver(),waiver(),flood_change_tbl,world_shape)

##building the composite hazard index ############
#################################################################################
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


mi=-1.72013
ma=0.3983822
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
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    ter_hist_ensemble_nrm[[4]][[y]]<-(ter_hist_ensemble_nrm[[4]][[y]]-mi)/(ma-mi)
  }
}
#future
the_years<-2020:2040
l_plot<-list(1:12,13:21)
for (ll in 1:length(l_plot)){
  tiff(paste0('./Plots/time_series/composite_ssp126_',ll,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    ter_ssp_ensemble_nrm[[4]][[y]]<-(ter_ssp_ensemble_nrm[[4]][[y]]-mi)/(ma-mi)
    plot(ter_ssp_ensemble_nrm[[4]][[y]],col=terrain.colors(255),zlim=c(0,1))
    title(the_years[y])
    writeRaster(ter_ssp_ensemble_nrm[[4]][[y]], filename = (paste0(output_fold,"Yearly_hazard/composite_ssp126_nrm_",the_years[y],".tif")), format = "GTiff",overwrite=TRUE)
  }
  dev.off()
}
saveRDS(ter_ssp_ensemble_nrm,paste0(output_fold,'Intermediate/ter_ssp126_ensemble_nrm_with_composite.RData'))

#######
#plotting yearly side by side with individual hazards
the_years<-2020:2040
l_plot<-seq(1,nlayers(ter_ssp_ensemble_nrm[[1]]),2)
for (each_page in l_plot){
  tiff(paste0('./Plots/time_series/all_ssp126_',each_page,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
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
composite_hist_index<-raster(paste0(output_fold,"composite_hist_nrm.tif"))
composite_ssp5_index<-raster(paste0(output_fold,"composite_ssp_nrm.tif"))

composite_ssp_index<-stackApply(ter_ssp_ensemble_nrm[[4]], indices =  rep(1,nlayers(ter_ssp_ensemble_nrm[[4]])), fun = "mean", na.rm = T)
writeRaster(composite_ssp_index, filename = (paste0(output_fold,"composite_ssp126_nrm.tif")), format = "GTiff",overwrite=TRUE)

#mapping the spatial distribution of compound hazard
composite_change_index<-composite_ssp_index-composite_hist_index
composite_ssp_tbl <- overlayRast(composite_ssp_index)
composite_change_tbl <- overlayRast(composite_change_index)
the_min<-min(c(minValue(composite_hist_index),minValue(composite_ssp_index),minValue(composite_ssp5_index)))
the_max<-max(c(maxValue(composite_hist_index),maxValue(composite_ssp_index),maxValue(composite_ssp5_index)))
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c(the_ticks[1],rep('',length(the_ticks)-2),the_ticks[length(the_ticks)])
col_map<-rev(brewer.pal(11,'Spectral'))
Visu_raster("./plots/composite_ssp126.tiff",color_bar_min,color_bar_max,"composite hazard - future",'composite index',col_map,the_ticks,the_labels,composite_ssp_tbl,world_shape)
Visu_raster("./plots/composite_change126.tiff",minValue(composite_change_index),maxValue(composite_change_index),"composite hazard - change",'composite change',col_map,waiver(),waiver(),composite_change_tbl,world_shape)

# build the 5 different severity categories and plot the resulting maps
brk  <- seq(0,1,0.2)
composite_ssp_class <- raster::cut(composite_ssp_index, breaks=brk) 
composite_ssp_class_tbl<-overlayRast(composite_ssp_class)

spect_pal<-rev(brewer.pal(5,'Spectral'))
my_pals=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])
Visu_raster_discrete("./plots/composite_ssp126_class.tiff","composite hazard - future",'hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),composite_ssp_class_tbl,world_shape)

#saving the composite hazard class into rasters
writeRaster(composite_ssp_class, filename = paste0(output_fold,"composite_ssp126_class.tif"), format = "GTiff",overwrite=TRUE)


###########################################
####FDP exposure
#########################################
country_data<-readRDS(paste0(output_fold,'Intermediate/country_hazard_values_full.Rds'))
new_world_shape<-world_shape
new_world_shape<-cbind(new_world_shape,country_data[,10:ncol(country_data)])

rast_base<-composite_ssp_index
values(rast_base)<-NA
Bgov_rast<-rasterize(new_world_shape,rast_base,field='misgov_ind')

#build the new risk index of hazard and misgovernance
risk_ssp<-(composite_ssp_index+Bgov_rast)/2
#categorize the new risk maps
brk  <- seq(0,1,0.2)
risk_ssp_class <- raster::cut(risk_ssp, breaks=brk) 
#saving data of the new risk 
writeRaster(risk_ssp, filename = paste0(output_fold,"composite_gov_ssp126.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(risk_ssp_class, filename = paste0(output_fold,"composite_gov_ssp126_class.tif"), format = "GTiff",overwrite=TRUE)

##add climate ssp126 hazards into the world shape
haz_ssp126<-c()
haz_ssp126_class<-c()
hazgov_ssp126<-c()
hazgov_ssp126_class<-c()
region<-c()
for (c in 1:nrow(new_world_shape)){
  haz_ssp126[c]<-exact_extract(composite_ssp_index,world_shape[c,],'mean')
  ssp126_cl<-cut(haz_ssp126[c],breaks=seq(0,1,0.2),label=1:5)
  ssp126_cl<-as.numeric(levels(ssp126_cl))[ssp126_cl]
  haz_ssp126_class[c]<-ssp126_cl
  
  hazgov_ssp126[c]<-exact_extract(risk_ssp,world_shape[c,],'mean')
  ssp126_cl<-cut(hazgov_ssp126[c],breaks=seq(0,1,0.2),label=1:5)
  ssp126_cl<-as.numeric(levels(ssp126_cl))[ssp126_cl]
  hazgov_ssp126_class[c]<-ssp126_cl
}
new_world_shape$haz_ssp126<-haz_ssp126
new_world_shape$haz_ssp126_class<-haz_ssp126_class
new_world_shape$hazgov_ssp126<-hazgov_ssp126
new_world_shape$hazgov_ssp126_class<-hazgov_ssp126_class


#add ssp126 hazard intensity and class, risk intensity and class for each country within the PoC data 
#########################
PoC_data<-readRDS(paste0(output_fold,'Intermediate/FDP_population_data_frame_full.Rds'))

PoC_data<-PoC_data[PoC_data$country_code!='UKN',]
the_hazards<-matrix(NA,nrow=nrow(PoC_data),ncol=4)
PoC_dict<-as.data.frame(read_excel(paste0(input_fold,"PoC_country_region_dictionary.xlsx")))
for (c in (1:nrow(PoC_data))){
  the_ind<-which(new_world_shape$iso3==PoC_data$country_code[c])
  if (length(the_ind)!=1){
    
    if (PoC_data$country_code[c]=='PRT'){
      the_ind<-the_ind[2]
      the_hazards[c,]<-as.numeric(new_world_shape[the_ind,c('haz_ssp126', 'haz_ssp126_class', 'hazgov_ssp126', 'hazgov_ssp126_class')])[1:4]
    }
    if (PoC_data$country_code[c]=='GBR'){
      the_ind<-the_ind[1]
      the_hazards[c,]<-as.numeric(new_world_shape[the_ind,c('haz_ssp126', 'haz_ssp126_class', 'hazgov_ssp126', 'hazgov_ssp126_class')])[1:4]
    }
    if (PoC_data$country_code[c]=='AB9'){
      the_ind<-which(new_world_shape$gis_name=='Abyei')
      the_hazards[c,]<-as.numeric(new_world_shape[the_ind,c('haz_ssp126', 'haz_ssp126_class', 'hazgov_ssp126', 'hazgov_ssp126_class')])[1:4]
    }
    
    if (PoC_data$country_code[c]=='PSE'){
      the_ind<-the_ind[2]
      the_hazards[c,]<-as.numeric(new_world_shape[the_ind,c('haz_ssp126', 'haz_ssp126_class', 'hazgov_ssp126', 'hazgov_ssp126_class')])[1:4]
    }
    if (!(PoC_data$country_code[c] %in% c('PRT','GBR','AB9','PSE'))){
      proxy_iso<-PoC_dict[(PoC_dict$country)==(PoC_data$country_name)[c],3]
      if (!is.na(proxy_iso)){
        the_hazards[c,]<-as.numeric(new_world_shape[which(new_world_shape$iso3==proxy_iso),c('haz_ssp126', 'haz_ssp126_class', 'hazgov_ssp126', 'hazgov_ssp126_class')])[1:4]
      }else{
        the_hazards[c,]<-rep(NA,4)
      }
    }
    
    
  }else{
    the_hazards[c,]<-as.numeric(new_world_shape[the_ind,c('haz_ssp126', 'haz_ssp126_class', 'hazgov_ssp126', 'hazgov_ssp126_class')])[1:4]
  }
}
PoC_data<-cbind(PoC_data,as.data.frame(the_hazards))
names(PoC_data)[(ncol(PoC_data)-3):ncol(PoC_data)]<-c('haz_ssp126', 'haz_class_ssp126', 'risk_ssp126', 'risk_class_ssp126')


#Plot FDP numbers per hazard class, overly each point with distribution over regions
annotation_custom_list <- function(pie.names,p_coords){
  result <- vector("list", length(pie.names) + 1)
  for(i in seq_along(pie.names)){
    pie <- pie.names[i]
    
    result[[i]] <- annotation_custom(
      grob = pie_list[[pie]],
      xmin = p_coords$xmin[p_coords$pie == pie],
      xmax = p_coords$xmax[p_coords$pie == pie],
      ymin = p_coords$ymin[p_coords$pie == pie],
      ymax = p_coords$ymax[p_coords$pie == pie])
  }
  
  # add a blank geom layer to ensure the resulting ggplot's
  # scales extend sufficiently to show each pie
  result[[length(result)]] <- geom_blank(
    data = p_coords %>% filter(pie %in% pie.names),
    aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
  )
  return(result)
}
#build a combined shapefile with all regions
regional_files<-c('AMERICAs_15m_UNHCR.gpkg',"ASIA_15m_UNHCR.gpkg","East_and_horn_of_Africa_shapefile.gpkg",
                  "Europe_new.gpkg","MENA_polygon_2.gpkg","Southern_Africa_Ash.gpkg","West_and_central_new.gpkg")
all_region<-c('Americas','Asia','East and Horn of Africa','Europe','MENA','Southern Africa','West and Central Africa')
all_region_short<-c('Americas','Asia','EHAGL-Africa','Europe','MENA','S-Africa','WC-Africa')
all_region_short_bis<-c('Americas','Asia','EHAGL\nAfrica','Europe','MENA','S-Africa','WC\nAfrica')

for (count_region in 1:length(all_region)){
  polygon_rst <- st_read(paste0(input_fold,"spatial_data/",regional_files[count_region]))
  the_region<-polygon_rst
  the_region$region_name<-rep(all_region[count_region],nrow(the_region))
  the_region$region_name_short<-rep(all_region_short[count_region],nrow(the_region))
  
  if (count_region==1)
  {combined_region<-the_region}
  if (count_region>1){
    combined_region<-rbind(combined_region,the_region)
  }
}

region_col<-hcl.colors(length(all_region),palette='Dynamic')
names(region_col)<-all_region_short
all_type<-c('haz','risk')
type_names<-c('','gov_')
all_names<-c('low','moderate','high','severe','extreme') 
for (t in (1:length(all_type))){
  class_count_ar<-c()
  pie_list<-list()
  the_names<-c()
  count_hazard<-1
  for (cla in 1:5){
    all_PoC_class<-eval(parse(text=paste0('PoC_data[which(PoC_data$',all_type[t],'_class_ssp126==cla),]')))
    PoC_class_count<-sum(all_PoC_class$PoC_number,na.rm=TRUE)
    class_count_ar<-c(class_count_ar,PoC_class_count)
    
    # #for each hazard class, plot a pie chart showing region
    reg_array<-c()
    for (the_reg in all_region){
      reg_array<-c(reg_array,sum(all_PoC_class$PoC_number[all_PoC_class$region==the_reg]))
    }
    the_df<-data.frame(all_region_short,reg_array)
    if (sum(reg_array)!=0){
      the_lab=paste0(round(reg_array / sum(reg_array) * 100, 1), "%")
    }else
    {
      the_lab=paste0(all_region_short, "0%")
    }
    the_df<-the_df[the_df$reg_array!=0 & the_lab!='0%',]
    the_df$fraction<-the_df$reg_array/sum(the_df$reg_array)
    the_lab<-the_lab[reg_array!=0 & the_lab!='0%']

    if (nrow(the_df)!=0){
      
      the_df2 <- the_df %>%
        mutate(csum = cumsum(reg_array),
               y_pos = csum-(reg_array/2))
      the_df2$the_lab<-the_lab
      
      the_df$ymax<-cumsum(the_df$reg_array)
      the_df$ymin<-c(0,head(the_df$ymax,n=-1))
      
      p<-ggplot() +
        geom_rect(the_df,mapping=aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=as.factor(all_region_short)),color='black',lwd=0.2) +
        coord_polar(theta="y") + 
        xlim(c(3, 4))+
        scale_fill_manual(values = region_col,name='region')+
        geom_label_repel(data=the_df2,aes(x=3.5,y = y_pos, label = the_lab,fill=as.factor(all_region_short)),
                         size = 8, show.legend = FALSE,segment.alpha = 0,segment.colour="black",force=10,force_pull =1)+
        theme_void()+
        theme(legend.position = "none")
      
      pie_list[[count_hazard]]<-ggplotGrob(p)
      rm(p)
      
      the_names<-c(the_names,all_names[[cla]])
      count_hazard<-count_hazard+1
      
    }
    
  }
  names(pie_list)<-the_names
  
  pie_coords <- data.frame(
    pie = names(pie_list),
    center.x = match(the_names,all_names),
    center.y = (class_count_ar/1e6)[match(the_names,all_names)]+22,
    radius.x = 2,
    radius.y = 25
  )
  pie_coords <- pie_coords %>%
    mutate(xmin = center.x - radius.x,
           xmax = center.x + radius.x,
           ymin = center.y - radius.y,
           ymax = center.y + radius.y)
  
  the_bar <- data.frame(
    name=factor(c('low','moderate','high','severe','extreme'),levels = c('low','moderate','high','severe','extreme') ),
    value=class_count_ar/1e6,on_x<-1:5)

  df_col<-data.frame('xmin'=1,'ymin'=rep(-10,length(region_col)), 'xmax'=2,'ymax'=rep(-20,length(region_col)),'region'=as.factor(names(region_col)))
  
  tiff(paste0('./Plots/with_FDP/FDP_number_',type_names[t],'per_class_ssp126.tiff'),units="in", width=12, height=8, res=500,compression = 'lzw')
  the_bar_plot<-ggplot()+
    geom_bar(data=the_bar, aes(x=on_x, y=value),stat = "identity",fill='grey59')+
    annotation_custom_list(names(pie_list),pie_coords)+
    xlab("Hazard class")+
    ylab('FDP number (million)')+
    geom_rect(data=df_col,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=region))+
    scale_fill_manual(name = "Region",
                      values = region_col)+
    scale_y_continuous(limits = c(0, 145),expand=c(0,0))+
    scale_x_continuous(breaks = 1:5, labels = c('low','moderate','high','severe','extreme'),limits = c(0.5, 5.5),expand=c(0,0))+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text = element_text(size = 25),
          axis.title = element_text(size = 30),
          axis.line = element_line(color = "black"),
          legend.position = c(0.1,0.8),
          legend.text = element_text(size =15 ), # Adjust legend label text size
          legend.title = element_text(size = 15, face = "bold"))
  
  print(the_bar_plot)
  dev.off()
  
}

#Plot number of FDP per regions, overly each point with distribution over hazard classes
spect_pal<-rev(brewer.pal(5,'Spectral'))
my_pals=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])
all_type<-c('haz','risk')
type_names<-c('','gov_')
for (t in (1:length(all_type))){
  region_count_ar<-c()
  pie_list<-list()
  the_names<-c()
  count_region<-1
  for (the_reg in 1:length(all_region)){
    reg<-all_region[the_reg]
    region_count_ar<-c(region_count_ar,sum(PoC_data$PoC_number[PoC_data$region==reg]))
    
    ##for each region, plot a pie chart showing hazard class
    the_reg_data<-PoC_data[PoC_data$region==reg,]
    
    class_array<-c()
    hazard_class<-1:5
    for (cla in hazard_class){
      region_ind<-eval(parse(text=paste0('the_reg_data[which(the_reg_data$',all_type[t],'_class_ssp126==cla),]')))
      region_ind<-sum(region_ind$PoC_number)
      class_array<-c(class_array,region_ind)
      
    }
    class_array<-round(class_array / sum(class_array) * 100, 1)
    class_array[which.max(class_array)]<- round(100-sum(class_array[-which.max(class_array)]),1)
    hazard_class<-hazard_class[class_array!=0]
    class_array<-class_array[class_array!=0]
    
    
    the_df<-data.frame(hazard_class,class_array)
    the_lab=paste0(as.character(the_df$class_array), "%")
    the_df$fraction<-(the_df$class_array)/sum(the_df$class_array)
    
    if (nrow(the_df)!=0){
      
      the_df$ymax<-cumsum(the_df$class_array)
      the_df$ymin<-c(0,head(the_df$ymax,n=-1))
      
      the_df2 <- the_df %>%
        mutate(y_pos = ymax-(class_array/2))
      the_df2$the_lab<-the_lab
      
      
      p <- ggplot()+
        geom_rect(the_df,mapping=aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=as.factor(hazard_class)),color='black',lwd=0.2) +
        coord_polar(theta="y") + 
        xlim(c(3, 4))+
        scale_fill_manual(values=my_pals, name='hazard class')+
        geom_label_repel(data = the_df2,
                         aes(x=3.5,y = y_pos, label = the_lab,fill=as.factor(hazard_class)),
                         size = 10, show.legend = FALSE,segment.alpha = 0,segment.colour="black",force=20,force_pull =5) +
        theme_void()+
        theme(legend.position = "none")
      
      
      pie_list[[count_region]]<-ggplotGrob(p)
      rm(p)
      
      the_names<-c(the_names,all_region[the_reg])
      count_region<-count_region+1
      
    }
    
  }
  names(pie_list)<-the_names
  
  pie_coords <- data.frame(
    pie = names(pie_list),
    center.x = match(the_names,all_region),
    center.y = (region_count_ar/1e6)[match(the_names,all_region)]-3.2,
    radius.x = 0.5,
    radius.y = 5
  )
  pie_coords <- pie_coords %>%
    mutate(xmin = center.x - radius.x,
           xmax = center.x + radius.x,
           ymin = center.y - radius.y,
           ymax = center.y + radius.y)
  
  the_bar <- data.frame(
    name=factor(all_region_short_bis,levels = all_region_short_bis ),
    value=region_count_ar/1e6)
  
  
  tiff(paste0('./Plots/with_FDP/FDP_number_',type_names[t],'per_region_ssp126.tiff'),units="in", width=8, height=17, res=500,compression = 'lzw')
  the_bar_plot<-ggplot()+
    geom_bar(data=the_bar, aes(x=name, y=value),stat = "identity",fill='grey59')+
    annotation_custom_list(names(pie_list),pie_coords)+
    coord_flip() +
    ylab("FDP number (million)")+
    xlab('')+
    xlim(all_region_short_bis)+
    scale_y_continuous(limits = c(0,28), expand = c(0, 0))+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text.x= element_text(size = 25),
          axis.title = element_text(size = 35),
          axis.line = element_line(color = "black"),
          axis.text.y = element_text(angle = 90, hjust = 0.5,size = 35),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank())
  
  
  print(the_bar_plot)
  dev.off()
  
}

#overly FDP locations on the new risk map 
PoC_points<-st_read(paste0(input_fold,'wrl_prp_p_unhcr_PoC.geojson'))
PoC_points$new_col<-c('PoC')
PoC_points$new_col<-factor(PoC_points$new_col)
#visualize the climate hazard class maps, as well as the new risk class maps with location of FDP
spect_pal<-rev(brewer.pal(5,'Spectral'))
my_pals=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])

composite_ssp_class_tbl<-overlayRast(composite_ssp_class)
Visu_raster_discrete_with_points(paste0("./Plots/with_FDP/composite_ssp126_class_with_FDP.tiff"),"",'Hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),composite_ssp_class_tbl,PoC_points,world_shape)

risk_ssp_class_tbl<-overlayRast(risk_ssp_class)
Visu_raster_discrete_with_points(paste0("./Plots/with_FDP/climate_gov_ssp126_with_FDP.tiff"),"",'Hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),risk_ssp_class_tbl,PoC_points,world_shape)

#######
###side visualization of the composite index per country

tiff("./Plots/composite_ssp126_class_avrg.tiff", units="in", width=10, height=6, res=500,compression = 'lzw')
plot(new_world_shape['haz_ssp126_class'],pal=spect_pal,breaks=0:5)
dev.off()

tiff("./Plots/composite_gov_ssp126_class_avrg.tiff", units="in", width=10, height=6, res=500,compression = 'lzw')
plot(new_world_shape['hazgov_ssp126_class'],pal=spect_pal,breaks=0:5)
dev.off()

###########################################
#look at the exposure of FDP population to each individual hazard and their specific class
#########################################################

#add individual climate hazard intensity and class for countries in the world shape file
all_hazards<-matrix(NA,nrow=nrow(world_shape),ncol=6)
hazards<-c('heat','drought','flood')
count_column=1
for (h in (1:length(hazards))){
  haz_ssp<-c()
  haz_ssp_class<-c()
  for (c in 1:nrow(new_world_shape)){
    haz_ssp[c]<-eval(parse(text = paste0('exact_extract(',hazards[h],'_ssp_index,new_world_shape[c,],"mean")'))) 
    ssp_cl<-cut(haz_ssp[c],breaks=seq(0,1,0.2),label=1:5)
    ssp_cl<-as.numeric(levels(ssp_cl))[ssp_cl]
    haz_ssp_class[c]<-ssp_cl
  }
  all_hazards[,count_column:(count_column+1)]<-cbind(haz_ssp,haz_ssp_class)
  count_column<-count_column+2
  
}
all_hazards<-as.data.frame(all_hazards)
names(all_hazards)<-c('heat_ssp126','heat_ssp126_class','drought_ssp126','drought_ssp126_class',
                      'flood_ssp126','flood_ssp126_class')
new_world_shape<-cbind(new_world_shape,all_hazards)
saveRDS(st_drop_geometry(new_world_shape),paste0('./Output/Intermediate/country_hazard_values_full_with126.Rds'))

#add individual hazard class and value for each country within the FDP data for ssp126
new_world_data<-st_drop_geometry(new_world_shape)
the_hazards<-matrix(NA,nrow=nrow(PoC_data),ncol=6)
PoC_dict<-as.data.frame(read_excel(paste0(input_fold,"PoC_country_region_dictionary.xlsx")))
for (c in (1:nrow(PoC_data))){
  the_ind<-which(new_world_data$iso3==PoC_data$country_code[c])
  if (length(the_ind)!=1){
    
    if (PoC_data$country_code[c]=='PRT'){
      the_ind<-the_ind[2]
      the_hazards[c,]<-as.numeric(new_world_data[the_ind,(ncol(new_world_data)-5):ncol(new_world_data)])
    }
    if (PoC_data$country_code[c]=='GBR'){
      the_ind<-the_ind[1]
      the_hazards[c,]<-as.numeric(new_world_data[the_ind,(ncol(new_world_data)-5):ncol(new_world_data)])
    }
    if (PoC_data$country_code[c]=='AB9'){
      the_ind<-which(new_world_data$gis_name=='Abyei')
      the_hazards[c,]<-as.numeric(new_world_data[the_ind,(ncol(new_world_data)-5):ncol(new_world_data)])
    }
    
    if (PoC_data$country_code[c]=='PSE'){
      the_ind<-the_ind[2]
      the_hazards[c,]<-as.numeric(new_world_data[the_ind,(ncol(new_world_data)-5):ncol(new_world_data)])
    }
    if (!(PoC_data$country_code[c] %in% c('PRT','GBR','AB9','PSE'))){
      proxy_iso<-PoC_dict[(PoC_dict$country)==(PoC_data$country_name)[c],3]
      if (!is.na(proxy_iso)){
        the_hazards[c,]<-as.numeric(new_world_data[the_ind,(ncol(new_world_data)-5):ncol(new_world_data)])
      }else{
        the_hazards[c,]<-rep(NA,6)
      }
    }
  }else{
    the_hazards[c,]<-as.numeric(new_world_data[the_ind,(ncol(new_world_data)-5):ncol(new_world_data)])
  }
}
the_hazards<-as.data.frame(the_hazards)
names(the_hazards)<-c('heat_ssp126','heat_ssp126_class','drought_ssp126','drought_ssp126_class',
                      'flood_ssp126','flood_ssp126_class')
PoC_data<-cbind(PoC_data,the_hazards)
saveRDS(PoC_data,paste0(output_fold,'Intermediate/FDP_population_data_frame_full_with126.Rds'))


#Plot FDSP numbers per hazard class
exposure_ssp<-data.frame(hazard.class=rep(1:5,4),PoC.number=NA,hazard=c(rep('heat',5),rep('drought',5),rep('flood',5),rep('compound',5)))
for (h in c('heat','drought','flood','compound')){
  for (cl in 1:5){
    exposure_ssp[(exposure_ssp$hazard.class==cl) & exposure_ssp$hazard=='heat',2]<-PoC_data$PoC_number[which((PoC_data$heat_ssp126_class)==cl)]%>%sum
    exposure_ssp[(exposure_ssp$hazard.class==cl) & exposure_ssp$hazard=='drought',2]<-PoC_data$PoC_number[which((PoC_data$drought_ssp126_class)==cl)]%>%sum
    exposure_ssp[(exposure_ssp$hazard.class==cl) & exposure_ssp$hazard=='flood',2]<-PoC_data$PoC_number[which((PoC_data$flood_ssp126_class)==cl)]%>%sum
    exposure_ssp[(exposure_ssp$hazard.class==cl) & exposure_ssp$hazard=='compound',2]<-PoC_data$PoC_number[which((PoC_data$haz_class_ssp126)==cl)]%>%sum
  }
}

spect_pal<-rev(brewer.pal(5,'Spectral'))
all_colors=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])
all_class=names(all_colors)

exposure_ssp$hazard<-factor(exposure_ssp$hazard)
exposure_ssp$hazard.class<-factor(exposure_ssp$hazard.class)
the_df2<-exposure_ssp[exposure_ssp$PoC.number!=0,]
the_df2$PoC.number<-the_df2$PoC.number/1e6

tiff(paste0('./Plots/with_FDP/exposure_per_class_ssp126.tiff'), units="in", width=8, height=4, res=500,compression = 'lzw')
ggplot(the_df2,aes(x=hazard,y=PoC.number,fill=hazard.class))+
  geom_bar(stat='identity',color = "black",linewidth=0.2,position = position_stack(reverse = TRUE))+
  labs(y='FDP number (millions)')+
  scale_fill_manual(values=all_colors,name='hazard class')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  theme(axis.text=element_text(size=13),axis.title = element_text(size = 15),legend.text = element_text(size = 13),legend.title = element_text(size = 15))+
  coord_flip()
dev.off()

#plot bars wrt class
ssp_df<-exposure_ssp[exposure_ssp$hazard!='compound',]
ssp_df$PoC.number<-ssp_df$PoC.number/1e6
ssp_df$hazard.class<-factor(ssp_df$hazard.class)
ssp_df$hazard<-factor(ssp_df$hazard,levels=c('heat','drought','flood'))

my_pals=c('heat'='orangered1','drought'='goldenrod3','flood'='steelblue3')

p2<-ggplot() + geom_bar(data = ssp_df, aes(x = hazard.class, y = PoC.number, fill = hazard), position = "dodge", stat = "identity",width=0.8)+
  labs(x='hazard class', y='FDSP number (million)')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank())+
  theme(axis.text=element_text(size=12),axis.title = element_text(size = 12),legend.text = element_text(size = 9),legend.title = element_text(size = 12),legend.key.size = unit(.5, "cm"),
        legend.position='inside',legend.position.inside = c(0.9, 0.8),axis.text.y = element_text(angle=45,vjust = 0.5, hjust=0.5))+
  scale_x_discrete(labels=c('low','moderate','high','severe','extreme'),breaks=c("1","2","3",'4','5'))+
  scale_fill_manual( values = my_pals)
p2<-p2+coord_flip()
tiff(paste0('./Plots/with_FDP/bar_exposure_per_class_ssp126.tiff'),units="in", width=7, height=4, res=500,compression = 'lzw')
p2
dev.off() 

