# This script transforms raw climate data into climate hazard indices
# and create visualization of the individual climate hazards distribution
rm(list=ls())

packages <- c("raster", "tidyverse", 'tibble','ggplot2','sf','RColorBrewer','colorspace','DescTools')

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


#input and pre-process individual hazards
#####################################################
input_fold<-'./Input/'
output_fold<-'./Output/'
dir.create(paste0(output_fold,'Yearly_hazard/'),recursive=T)
dir.create('Plots',recursive=T)
dir.create('Plots/time_series',recursive=T)
world_shape <- st_read(paste0(input_fold,'world_shape2.geojson'))


rast_list=vector("list", length = 6)
rast_list[[1]]<-readRDS(paste0(output_fold,'./Intermediate/hist_ensemble_heat.RData'))
rast_list[[2]]<-readRDS(paste0(output_fold,'./Intermediate/hist_ensemble_drought.RData'))
rast_list[[3]]<-readRDS(paste0(output_fold,'./Intermediate/hist_ensemble_flood.RData'))
rast_list[[4]]<-readRDS(paste0(output_fold,'./Intermediate/ssp_ensemble_heat.RData'))
rast_list[[5]]<-readRDS(paste0(output_fold,'./Intermediate/ssp_ensemble_drought.RData'))
rast_list[[6]]<-readRDS(paste0(output_fold,'./Intermediate/ssp_ensemble_flood.RData'))

#consider only terrestrial data (crop and mask marine data)
ter_ensemble<-list()
for (haz in 1:6){
  ter_ensemble[[haz]]<-formatRast(rast_list[[haz]],world_shape)
}
ter_hist_ensemble<-ter_ensemble[1:3]
ter_ssp_ensemble<-ter_ensemble[4:6]

saveRDS(ter_hist_ensemble,paste0(output_fold,'Intermediate/ter_hist_ensemble.RData'))
saveRDS(ter_ssp_ensemble,paste0(output_fold,'Intermediate/ter_ssp_ensemble.RData'))

#remove values in greenland
gr<-world_shape[which(world_shape$gis_name=='Greenland (DNK)'),]
for (h in 1:length(ter_hist_ensemble)){
  ter_hist_ensemble[[h]]<-mask(ter_hist_ensemble[[h]],gr,inverse=TRUE)
  ter_ssp_ensemble[[h]]<-mask(ter_ssp_ensemble[[h]],gr,inverse=TRUE)
}


##transforming raw indices to be comparable (normalization and distribution shift)
#######################################################################
extreme_thresh<-0.99 #considered as 'extremes' are those values above the upper 99% of the distribution

##heat
#############

dist_heat<-as.vector(values(ter_hist_ensemble[[1]]))
dist_heat_ssp<-as.vector(values(ter_ssp_ensemble[[1]]))
considered_ref<-c(dist_heat,as.vector(ter_ssp_ensemble[[1]][[1:5]]))
extreme_heat<-quantile(considered_ref, extreme_thresh , na.rm = TRUE)
extreme_heat
#normalization
min_heat<-min(c(min(dist_heat,na.rm=T),min(dist_heat_ssp,na.rm=T)))
max_heat<-max(c(max(dist_heat,na.rm=T),max(dist_heat_ssp,na.rm=T)))
c(min_heat,max_heat)
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
dist_drought_ssp_ref<-as.vector(values(ter_ssp_ensemble[[2]][[1:5]]))
dist_drought_ssp_ref[dist_drought_ssp_ref>0]<-0 #only taking the negative ones
dist_drought_ssp_ref<--dist_drought_ssp_ref
extreme_drought<-quantile(c(dist_drought,dist_drought_ssp_ref), extreme_thresh , na.rm = TRUE)
extreme_drought
#normalization
min_drought<-min(c(min(dist_drought,na.rm=T),min(dist_drought_ssp,na.rm=T)))
max_drought<-max(c(max(dist_drought,na.rm=T),max(dist_drought_ssp,na.rm=T)))
c(min_drought,max_drought)
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
considered_ref<-c(dist_flood,as.vector(ter_ssp_ensemble[[3]][[1:5]]))
extreme_flood<-quantile(considered_ref[considered_ref>=10], extreme_thresh , na.rm = TRUE)
extreme_flood
#normalization
min_flood<-min(c(min(dist_flood,na.rm=T),min(dist_flood_ssp,na.rm=T)))
max_flood<-max(c(max(dist_flood,na.rm=T),max(dist_flood_ssp,na.rm=T)))
c(min_flood,max_flood)
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
tiff(paste0('./Plots/raw_hazard_distribution.tiff'),units="in", width=7, height=4, res=300,compression = 'lzw')
par(mfrow=c(2,3))

hist(dist_heat[dist_heat>=1],xlab='number of days with HI>41°C',ylab='frequency',breaks = "Scott",main='',xlim=c(min_heat,200),ylim=c(0,50000))
title(main = "heat",adj = 0)
abline(v=extreme_heat,col='red',lwd=1.5,lty=2)
abline(v=mean(dist_heat[dist_heat>=1],na.rm=T),col='blue',lwd=1.5,lty=2)

hist(dist_drought[dist_drought>0],xlab='-SPEI',ylab='frequency',breaks=75,main='',xlim=c(min_drought,max(dist_drought_ssp,na.rm=T)))
title(main = "drought",adj = 0)
abline(v=extreme_drought,col='red',lwd=1.5,lty=2)
abline(v=mean(dist_drought[dist_drought>0],na.rm=T),col='blue',lwd=1.5,lty=2)

mtext("Baseline", side = 3, line = 3, outer =FALSE)
hist(dist_flood,xlab='SDII',ylab='frequency',breaks = "Scott",main='',xlim=c(min_flood,max_flood))
title(main = "flood",adj = 0)
abline(v=extreme_flood,col='red',lwd=1.5,lty=2)
abline(v=mean(dist_flood,na.rm=T),col='blue',lwd=1.5,lty=2)
legend('topright',inset=c(-0.1, -0.5), legend=c("threshold for extreme",'average'),
       col=c("red",'blue'), lty=2, cex=0.8,xpd=TRUE)


hist(dist_heat_ssp[dist_heat_ssp>=1],xlab='number of days with HI>41°C',ylab='frequency',breaks = "Scott",main='',xlim=c(min_heat,200),ylim=c(0,50000))
title(main = "heat",adj = 0)
abline(v=extreme_heat,col='red',lwd=1.5,lty=2)
abline(v=mean(dist_heat_ssp[dist_heat_ssp>=1],na.rm=T),col='blue',lwd=1.5,lty=2)

hist(dist_drought_ssp[dist_drought_ssp>0],xlab='-SPEI',ylab='frequency',breaks=75,main='',xlim=c(min_drought,max(dist_drought_ssp,na.rm=T)))
title(main = "drought",adj = 0)
abline(v=extreme_drought,col='red',lwd=1.5,lty=2)
abline(v=mean(dist_drought_ssp[dist_drought_ssp>0],na.rm=T),col='blue',lwd=1.5,lty=2)

mtext("Future", side = 3, line = 3, outer =FALSE)
hist(dist_flood_ssp,xlab='SDII',ylab='frequency',breaks = "Scott",main='',xlim=c(min_flood,max_flood))
title(main = "flood",adj = 0)
abline(v=extreme_flood,col='red',lwd=1.5,lty=2)
abline(v=mean(dist_flood_ssp,na.rm=T),col='blue',lwd=1.5,lty=2)

dev.off()

#plotting hazard distribution after normalizing and shifting
spect_col<-rev(brewer.pal(5,'Spectral'))

tiff(paste0('./Plots/normalized_shifted_hazard_distribution.tiff'),units="in", width=7, height=5, res=300,compression = 'lzw')
par(mfrow=c(2,3))

hist(dist_heat_sft,xlab='heat index',ylab='frequency',breaks = "Scott",main='',xlim=c(0,1))
title(main = "heat",adj = 0)
abline(v=seq(0,0.8,0.2),col=spect_col,lwd=1.7,lty=2)
hist(dist_drought_sft,xlab='drought index',ylab='frequency',breaks=75,main='',xlim=c(0,1))
title(main = "drought",adj = 0)
abline(v=seq(0,0.8,0.2),col=spect_col,lwd=1.7,lty=2)
mtext("Baseline", side = 3, line = 3, outer =FALSE)
hist(dist_flood_sft,xlab='flood index',ylab='frequency',breaks = "Scott",main='',xlim=c(0,1))
title(main = "flood",adj = 0)
abline(v=seq(0,0.8,0.2),col=spect_col,lwd=1.7,lty=2)
legend('topright',inset=c(-0.15, -0.4), legend=c("low class","moderate class ",'high class ','severe class ','extreme class '),
       title='lower threshold',title.font=2,col=spect_col, lty=2, cex=1,xpd=TRUE)

hist(dist_heat_ssp_sft,xlab='heat index',ylab='frequency',breaks = "Scott",main='',xlim=c(0,1))
title(main = "heat",adj = 0)
abline(v=seq(0,0.8,0.2),col=spect_col,lwd=1.7,lty=2)
hist(dist_drought_ssp_sft,xlab='drought index',ylab='frequency',breaks='Scott',main='',xlim=c(0,1))
title(main = "drought",adj = 0)
abline(v=seq(0,0.8,0.2),col=spect_col,lwd=1.7,lty=2)
mtext("Future", side = 3, line = 3, outer =FALSE)
hist(dist_flood_ssp_sft,xlab='flood index',ylab='frequency',breaks = "Scott",main='',xlim=c(0,1))
title(main = "flood",adj = 0)
abline(v=seq(0,0.8,0.2),col=spect_col,lwd=1.7,lty=2)
dev.off()

#plotting distribution frequency of the hazards
theme_frequency<- theme(title = element_text(size = 15),
                        plot.title=element_text(hjust=0.5),
                        axis.text=element_text(size=10),
                        axis.line = element_line(colour = "black"),
                        plot.background = element_blank(),
                        panel.background=element_rect(fill = "white"),
                        panel.border = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.title = element_text(size=10,face='bold'),
                        legend.key = element_rect(fill = "white"))
the_x<-dist_heat_nrm
the_x[which(dist_heat<1)]<-NA
the_y<-dist_drought_nrm
the_y[which(dist_drought==0)]<-NA
the_z<-dist_flood_nrm
df_freq_hist<- data.frame(x = the_x, y = the_y,z=the_z) %>%
  gather(key, value)
tiff('./Plots/hazard_distribution_hist.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
ggplot(df_freq_hist, aes(value, y=after_stat(density),colour = key)) +
  stat_density(geom="line",position="identity") +
  geom_vline(xintercept = extreme_heat_nrm, linetype="dashed", 
             color = "red")+
  geom_vline(xintercept = extreme_drought_nrm, linetype="dashed", 
             color = "yellow")+
  geom_vline(xintercept = extreme_flood_nrm, linetype="dashed", 
             color = "blue")+
  theme_frequency+
  scale_color_manual(name='',values = c(x = "red", y = "yellow",z='blue'),labels=c('heat','drought','flood'))+
  labs(title='hazard distribution - baseline',x='hazard index value',y='density')
dev.off()

the_x<-dist_heat_ssp_nrm
the_x[which(dist_heat_ssp<1)]<-NA
the_y<-dist_drought_ssp_nrm
the_y[the_y==0]<-NA
the_z<-dist_flood_ssp_nrm
df_freq_ssp<- data.frame(x = the_x, y = the_y,z=the_z) %>%
  gather(key, value)
tiff('./Plots/hazard_distribution_ssp.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
ggplot(df_freq_ssp, aes(value,y=after_stat(density),colour = key)) +
  stat_density(geom="line",position="identity") +
  geom_vline(xintercept = extreme_heat_nrm, linetype="dashed", 
             color = "red")+
  geom_vline(xintercept = extreme_drought_nrm, linetype="dashed", 
             color = "yellow")+
  geom_vline(xintercept = extreme_flood_nrm, linetype="dashed", 
             color = "blue")+
  theme_frequency+
  scale_color_manual(name='',values = c(x = "red", y = "yellow",z='blue'),labels=c('heat','drought','flood'))+
  labs(title='hazard distribution - future',x='hazard index value',y='density')
dev.off()

#plot the new shifted distributions
the_x<-dist_heat_sft
the_x[which(dist_heat<1)]<-NA
the_y<-dist_drought_sft
the_y[which(dist_drought==0)]<-NA
the_z<-dist_flood_sft
df_freq_hist_sft<- data.frame(x = the_x, y = the_y,z=the_z) %>%
  gather(key, value)
tiff('./Plots/hazard_distribution_shifted_hist.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
ggplot(df_freq_hist_sft, aes(value,y=after_stat(density),colour = key)) +
  stat_density(geom="line",position="identity") +
  geom_vline(xintercept = 0.8, linetype="dashed", 
             color = "black")+
  geom_vline(xintercept = 0.8, linetype="dashed", 
             color = "black")+
  geom_vline(xintercept = 0.8, linetype="dashed", 
             color = "black")+
  theme_frequency+
  scale_color_manual(name='',values = c(x = "red", y = "yellow",z='blue'),labels=c('heat','drought','flood'))+
  labs(title='hazard distribution - baseline',x='hazard index value',y='density')
dev.off()

the_x<-dist_heat_ssp_sft
the_x[which(dist_heat_ssp<1)]<-NA
the_y<-dist_drought_ssp_sft
the_y[which(dist_drought_ssp==0)]<-NA
the_z<-dist_flood_ssp_sft
df_freq_ssp_sft<- data.frame(x = the_x, y = the_y,z=the_z) %>%
  gather(key, value)
tiff('./Plots/hazard_distribution_shifted_ssp.tiff', units="in", width=6, height=4, res=500,compression = 'lzw')
ggplot(df_freq_ssp_sft, aes(value,y=after_stat(density),colour = key)) +
  stat_density(geom="line",position="identity") +
  geom_vline(xintercept = 0.8, linetype="dashed", 
             color = "black")+
  geom_vline(xintercept = 0.8, linetype="dashed", 
             color = "black")+
  geom_vline(xintercept = 0.8, linetype="dashed", 
             color = "black")+
  theme_frequency+
  scale_color_manual(name='',values = c(x = "red", y = "yellow",z='blue'),labels=c('heat','drought','flood'))+
  labs(title='hazard distribution - future',x='hazard index value',y='density')
dev.off()



###########################
#plotting over years and store the final indices in raster bricks
all_years<-nlayers(ter_hist_ensemble[[1]])
yearly_index=seq(0,length(dist_heat_sft),ncell(ter_hist_ensemble[[1]][[1]]))
yearly_index=cbind(c((yearly_index[1:length(yearly_index)-1])+1),yearly_index[2:length(yearly_index)])

ter_hist_ensemble_nrm<-list(raster::stack(),raster::stack(),raster::stack())
#hist
l_plot<-list(1:12,13:20)
the_years<-1995:2014
for (ll in 1:length(l_plot)){
  tiff(paste0('./Plots/time_series/heat_hist_',ll,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    heat_val<-dist_heat_sft[yearly_index[y,1]:yearly_index[y,2]]
    heat_raster<-ter_hist_ensemble[[1]][[1]]
    values(heat_raster)<-heat_val
    ter_hist_ensemble_nrm[[1]]<-addLayer(ter_hist_ensemble_nrm[[1]],heat_raster)
    plot(heat_raster,col=terrain.colors(255),zlim=c(0,1))
    title(the_years[y])
    writeRaster(heat_raster, filename = (paste0(output_fold,"Yearly_hazard/heat_hist_nrm_",the_years[y],".tif")), format = "GTiff",overwrite=TRUE)
    
  }
  dev.off()
  
  tiff(paste0('./Plots/time_series/drought_hist_',ll,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    drought_val<-dist_drought_sft[yearly_index[y,1]:yearly_index[y,2]]
    drought_raster<-ter_hist_ensemble[[2]][[1]]
    values(drought_raster)<-drought_val
    ter_hist_ensemble_nrm[[2]]<-addLayer(ter_hist_ensemble_nrm[[2]],drought_raster)
    plot(drought_raster,col=terrain.colors(255),zlim=c(0,1))
    title(the_years[y])
    writeRaster(heat_raster, filename = (paste0(output_fold,"Yearly_hazard/drought_hist_nrm_",the_years[y],".tif")), format = "GTiff",overwrite=TRUE)
    
  }
  dev.off()
  
  tiff(paste0('./Plots/time_series/flood_hist_',ll,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    flood_val<-dist_flood_sft[yearly_index[y,1]:yearly_index[y,2]]
    flood_raster<-ter_hist_ensemble[[2]][[1]]
    values(flood_raster)<-flood_val
    ter_hist_ensemble_nrm[[3]]<-addLayer(ter_hist_ensemble_nrm[[3]],flood_raster)
    plot(flood_raster,col=terrain.colors(255),zlim=c(0,1))
    title(the_years[y])
    writeRaster(heat_raster, filename = (paste0(output_fold,"Yearly_hazard/flood_hist_nrm_",the_years[y],".tif")), format = "GTiff",overwrite=TRUE)
    
  }
  dev.off()
  
}
#ssp
all_years<-nlayers(ter_ssp_ensemble[[1]])
yearly_index=seq(0,length(dist_heat_ssp_sft),ncell(ter_ssp_ensemble[[1]][[1]]))
yearly_index=cbind(c((yearly_index[1:length(yearly_index)-1])+1),yearly_index[2:length(yearly_index)])

the_years<-2020:2040
ter_ssp_ensemble_nrm<-list(raster::stack(),raster::stack(),raster::stack())
l_plot<-list(1:12,13:21)
for (ll in 1:length(l_plot)){
  tiff(paste0('./Plots/time_series/heat_ssp_',ll,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    heat_val<-dist_heat_ssp_sft[yearly_index[y,1]:yearly_index[y,2]]
    heat_raster<-ter_ssp_ensemble[[1]][[1]]
    values(heat_raster)<-heat_val
    ter_ssp_ensemble_nrm[[1]]<-addLayer(ter_ssp_ensemble_nrm[[1]],heat_raster)
    plot(heat_raster,col=terrain.colors(255),zlim=c(0,1))
    title(the_years[y])
    writeRaster(heat_raster, filename = (paste0(output_fold,"Yearly_hazard/heat_ssp_nrm_",the_years[y],".tif")), format = "GTiff",overwrite=TRUE)
  }
  dev.off()
  
  tiff(paste0('./Plots/time_series/drought_ssp_',ll,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    drought_val<-dist_drought_ssp_sft[yearly_index[y,1]:yearly_index[y,2]]
    drought_raster<-ter_ssp_ensemble[[2]][[1]]
    values(drought_raster)<-drought_val
    ter_ssp_ensemble_nrm[[2]]<-addLayer(ter_ssp_ensemble_nrm[[2]],drought_raster)
    plot(drought_raster,col=terrain.colors(255),zlim=c(0,1))
    title(the_years[y])
    writeRaster(heat_raster, filename = (paste0(output_fold,"Yearly_hazard/drought_ssp_nrm_",the_years[y],".tif")), format = "GTiff",overwrite=TRUE)
  }
  dev.off()
  
  tiff(paste0('./Plots/time_series/flood_ssp_',ll,'.tiff'),units="in", width=8, height=8, res=300,compression = 'lzw')
  par(mfrow=c(4,3),mar=c(2,3,2,3)+0.1)
  for (y in l_plot[[ll]]){
    flood_val<-dist_flood_ssp_sft[yearly_index[y,1]:yearly_index[y,2]]
    flood_raster<-ter_ssp_ensemble[[3]][[1]]
    values(flood_raster)<-flood_val
    ter_ssp_ensemble_nrm[[3]]<-addLayer(ter_ssp_ensemble_nrm[[3]],flood_raster)
    plot(flood_raster,col=terrain.colors(255),zlim=c(0,1))
    title(the_years[y])
    writeRaster(heat_raster, filename = (paste0(output_fold,"Yearly_hazard/drought_ssp_nrm_",the_years[y],".tif")), format = "GTiff",overwrite=TRUE)
  }
  dev.off()
  
}
saveRDS(ter_hist_ensemble_nrm,paste0(output_fold,"Intermediate/ter_hist_ensemble_nrm.RData"))
saveRDS(ter_ssp_ensemble_nrm,paste0(output_fold,'Intermediate/ter_ssp_ensemble_nrm.RData'))


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

writeRaster(heat_hist_index, filename = (paste0(output_fold,"heat_hist_nrm.tif")), format = "GTiff",overwrite=TRUE)
writeRaster(drought_hist_index, filename = (paste0(output_fold,"drought_hist_nrm.tif")), format = "GTiff",overwrite=TRUE)
writeRaster(flood_hist_index, filename = (paste0(output_fold,"flood_hist_nrm.tif")), format = "GTiff",overwrite=TRUE)
writeRaster(heat_ssp_index, filename = (paste0(output_fold,"heat_ssp_nrm.tif")), format = "GTiff",overwrite=TRUE)
writeRaster(drought_ssp_index, filename = (paste0(output_fold,"drought_ssp_nrm.tif")), format = "GTiff",overwrite=TRUE)
writeRaster(flood_ssp_index, filename = (paste0(output_fold,"flood_ssp_nrm.tif")), format = "GTiff",overwrite=TRUE)

heat_hist_index<-raster(paste0(output_fold,"heat_hist_nrm.tif"))
heat_ssp_index<-raster(paste0(output_fold,"heat_ssp_nrm.tif"))
drought_hist_index<-raster(paste0(output_fold,"drought_hist_nrm.tif"))
drought_ssp_index<-raster(paste0(output_fold,"drought_ssp_nrm.tif"))
flood_hist_index<-raster(paste0(output_fold,"flood_hist_nrm.tif"))
flood_ssp_index<-raster(paste0(output_fold,"flood_ssp_nrm.tif"))


future=list('','126')
ff<-future[1]
##visualizing individual climate hazard indices (current, future, change)
##visualizing heat hazard
heat_hist_tbl <- overlayRast(heat_hist_index)
heat_ssp_tbl <- overlayRast(heat_ssp_index)
heat_change_index<-heat_ssp_index-heat_hist_index
heat_change_tbl<-overlayRast(heat_change_index)
the_min<-min(c(minValue(heat_hist_index),minValue(heat_ssp_index)))
the_max<-max(c(maxValue(heat_hist_index),maxValue(heat_ssp_index)))
col_map<-rev(sequential_hcl(27, palette = "Reds3"))
col_map<-col_map[2:length(col_map)]
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c('low',rep('',length(the_ticks)-2),'extreme')
Visu_raster("./Plots/heat_hist_index.tiff",color_bar_min,color_bar_max,"heat hazard - baseline",'heat hazard',col_map,the_ticks,the_labels,heat_hist_tbl,world_shape)
Visu_raster(paste0("./Plots/heat_ssp",ff,"_index.tiff"),color_bar_min,color_bar_max,"heat hazard - future",'heat hazard',col_map,the_ticks,the_labels,heat_ssp_tbl,world_shape)
col_map<-rev(brewer.pal(11,'Spectral')) 
Visu_raster(paste0("./Plots/heat_change_index",ff,".tiff"),minValue(heat_change_index),maxValue(heat_change_index),"heat hazard - change",'heat change',col_map,waiver(),waiver(),heat_change_tbl,world_shape)

##visualizing drought hazard
drought_hist_tbl <- overlayRast(drought_hist_index)
drought_ssp_tbl <- overlayRast(drought_ssp_index)
drought_change_index<-drought_ssp_index-drought_hist_index
drought_change_tbl<-overlayRast(drought_change_index)
the_min<-min(c(minValue(drought_hist_index),minValue(drought_ssp_index)))
the_max<-max(c(maxValue(drought_hist_index),maxValue(drought_ssp_index)))
col_map<-rev(sequential_hcl(40, palette = "YlOrBr"))
col_map<-c(rep(col_map[1],3),rep(col_map[2],4),rep(col_map[3],5),col_map[4:length(col_map)])
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c('low',rep('',length(the_ticks)-2),'extreme')
Visu_raster("./Plots/drought_hist_index.tiff",color_bar_min,color_bar_max,"drought hazard - baseline",'drought hazard',col_map,the_ticks,the_labels,drought_hist_tbl,world_shape)
Visu_raster(paste0("./Plots/drought_ssp",ff,"_index.tiff"),color_bar_min,color_bar_max,"drought hazard - future",'drought hazard',col_map,the_ticks,the_labels,drought_ssp_tbl,world_shape)
col_map<-rev(brewer.pal(11,'Spectral')) 
Visu_raster(paste0("./Plots/drought_change_index",ff,".tiff"),minValue(drought_change_index),maxValue(drought_change_index),"drought hazard - change",'drought change',col_map,waiver(),waiver(),drought_change_tbl,world_shape)

##visualizing flood hazard
flood_hist_tbl <- overlayRast(flood_hist_index)
flood_ssp_tbl <- overlayRast(flood_ssp_index)
flood_change_index<-flood_ssp_index-flood_hist_index
flood_change_tbl<-overlayRast(flood_change_index)
the_min<-min(c(minValue(flood_hist_index),minValue(flood_ssp_index)))
the_max<-max(c(maxValue(flood_hist_index),maxValue(flood_ssp_index)))
col_map<-rev(sequential_hcl(40, palette = "Blues3"))
col_map<-c(rep(col_map[1],3),rep(col_map[2],4),rep(col_map[3],5),col_map[4:length(col_map)])
color_bar_min<-RoundTo(the_min, multiple = 0.1, FUN = floor)
color_bar_max<-RoundTo(the_max, multiple = 0.1, FUN = ceiling)
the_ticks<-seq(color_bar_min,color_bar_max,0.1)
the_labels<-c('low',rep('',length(the_ticks)-2),'extreme')
Visu_raster("./Plots/flood_hist_index.tiff",color_bar_min,color_bar_max,"flood hazard - baseline",'flood hazard',col_map,the_ticks,the_labels,flood_hist_tbl,world_shape)
Visu_raster(paste0("./Plots/flood_ssp",ff,"_index.tiff"),color_bar_min,color_bar_max,"flood hazard - future",'flood hazard',col_map,the_ticks,the_labels,flood_ssp_tbl,world_shape)
col_map<-rev(brewer.pal(11,'Spectral')) 
Visu_raster(paste0("./Plots/flood_change_index",ff,".tiff"),minValue(flood_change_index),maxValue(flood_change_index),"flood hazard - change",'flood change',col_map,waiver(),waiver(),flood_change_tbl,world_shape)





