#This script looks at the exposure of FDP population to each individual hazard and their specific class

rm(list = ls())

rm(list=ls())

packages <- c("raster", "ggplot2", 'sf',"exactextractr", "RColorBrewer",'readxl','cowplot')

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


setwd(dirname(rstudioapi::getSourceEditorContext()$path))


input_fold<-'./Input/'
output_fold<-'./Output/'

world_shape <- st_read(paste0(input_fold,"world_shape2.geojson"))

ff<-''

#read raster layers of individual hazards
###########################################
heat_hist_index<-raster(paste0(output_fold,"heat_hist_nrm.tif"))
heat_ssp_index<-raster(paste0(output_fold,"heat_ssp_nrm.tif"))
drought_hist_index<-raster(paste0(output_fold,"drought_hist_nrm.tif"))
drought_ssp_index<-raster(paste0(output_fold,"drought_ssp_nrm.tif"))
flood_hist_index<-raster(paste0(output_fold,"flood_hist_nrm.tif"))
flood_ssp_index<-raster(paste0(output_fold,"flood_ssp_nrm.tif"))


#add individual climate hazard intensity and class for countries in the world shape file
all_hazards<-matrix(NA,nrow=nrow(world_shape),ncol=12)
hazards<-c('heat','drought','flood')
count_column=1
for (h in (1:length(hazards))){
  haz_hist<-c()
  haz_ssp<-c()
  haz_hist_class<-c()
  haz_ssp_class<-c()
  for (c in 1:nrow(world_shape)){
    haz_hist[c]<-eval(parse(text = paste0('exact_extract(',hazards[h],'_hist_index,world_shape[c,],"mean")'))) 
    haz_ssp[c]<-eval(parse(text = paste0('exact_extract(',hazards[h],'_ssp_index,world_shape[c,],"mean")'))) 
    hist_cl<-cut(haz_hist[c],breaks=seq(0,1,0.2),label=1:5)
    hist_cl<-as.numeric(levels(hist_cl))[hist_cl]
    haz_hist_class[c]<-hist_cl
    ssp_cl<-cut(haz_ssp[c],breaks=seq(0,1,0.2),label=1:5)
    ssp_cl<-as.numeric(levels(ssp_cl))[ssp_cl]
    haz_ssp_class[c]<-ssp_cl
  }
  all_hazards[,count_column:(count_column+3)]<-cbind(haz_hist,haz_ssp,haz_hist_class,haz_ssp_class)
  count_column<-count_column+4
  
}
world_data<-readRDS(paste0(output_fold,'Intermediate/country_hazard_values',ff,'.Rds'))
all_hazards<-as.data.frame(all_hazards)
names(all_hazards)<-c('heat_hist','heat_ssp','heat_hist_class','heat_ssp_class','drought_hist','drought_ssp','drought_hist_class','drought_ssp_class',
                      'flood_hist','flood_ssp','flood_hist_class','flood_ssp_class')
world_data_full<-cbind(world_data,all_hazards)
saveRDS(world_data_full,paste0(output_fold,'Intermediate/country_hazard_values_full.Rds'))

#read FDP data and add hazard class for each country within the FDP data 
PoC_data<-readRDS(paste0(output_fold,'Intermediate/FDP_population_data_frame',ff,'.Rds'))

the_hazards<-matrix(NA,nrow=nrow(PoC_data),ncol=12)
PoC_dict<-as.data.frame(read_excel(paste0(input_fold,"PoC_country_region_dictionary.xlsx")))
for (c in (1:nrow(PoC_data))){
  the_ind<-which(world_data_full$iso3==PoC_data$country_code[c])
  if (length(the_ind)!=1){
    
    if (PoC_data$country_code[c]=='PRT'){
      the_ind<-the_ind[2]
      the_hazards[c,]<-as.numeric(world_data_full[the_ind,(ncol(world_data_full)-11):ncol(world_data_full)])
    }
    if (PoC_data$country_code[c]=='GBR'){
      the_ind<-the_ind[1]
      the_hazards[c,]<-as.numeric(world_data_full[the_ind,(ncol(world_data_full)-11):ncol(world_data_full)])
    }
    if (PoC_data$country_code[c]=='AB9'){
      the_ind<-which(world_data_full$gis_name=='Abyei')
      the_hazards[c,]<-as.numeric(world_data_full[the_ind,(ncol(world_data_full)-11):ncol(world_data_full)])
    }
    
    if (PoC_data$country_code[c]=='PSE'){
      the_ind<-the_ind[2]
      the_hazards[c,]<-as.numeric(world_data_full[the_ind,(ncol(world_data_full)-11):ncol(world_data_full)])
    }
    if (!(PoC_data$country_code[c] %in% c('PRT','GBR','AB9','PSE'))){
      proxy_iso<-PoC_dict[(PoC_dict$country)==(PoC_data$country_name)[c],3]
      if (!is.na(proxy_iso)){
        the_hazards[c,]<-as.numeric(world_data_full[the_ind,(ncol(world_data_full)-11):ncol(world_data_full)])
      }else{
        the_hazards[c,]<-rep(NA,12)
      }
    }
  }else{
    the_hazards[c,]<-as.numeric(world_data_full[the_ind,(ncol(world_data_full)-11):ncol(world_data_full)])
  }
}
the_hazards<-as.data.frame(the_hazards)
names(the_hazards)<-c('heat_hist','heat_ssp','heat_hist_class','heat_ssp_class','drought_hist','drought_ssp','drought_hist_class','drought_ssp_class',
                      'flood_hist','flood_ssp','flood_hist_class','flood_ssp_class')
PoC_data<-cbind(PoC_data,the_hazards)
saveRDS(PoC_data,paste0(output_fold,'Intermediate/FDP_population_data_frame_full',ff,'.Rds'))


#Plot FDSP numbers per hazard class
exposure_hist<-data.frame(hazard.class=rep(1:5,4),PoC.number=NA,hazard=c(rep('heat',5),rep('drought',5),rep('flood',5),rep('compound',5)))
exposure_ssp<-data.frame(hazard.class=rep(1:5,4),PoC.number=NA,hazard=c(rep('heat',5),rep('drought',5),rep('flood',5),rep('compound',5)))
for (h in c('heat','drought','flood','compound')){
  for (cl in 1:5){
    exposure_hist[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='heat',2]<-PoC_data$PoC_number[which((PoC_data$heat_hist_class)==cl)]%>%sum
    exposure_hist[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='drought',2]<-PoC_data$PoC_number[which((PoC_data$drought_hist_class)==cl)]%>%sum
    exposure_hist[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='flood',2]<-PoC_data$PoC_number[which((PoC_data$flood_hist_class)==cl)]%>%sum
    exposure_hist[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='compound',2]<-PoC_data$PoC_number[which((PoC_data$haz_class_hist)==cl)]%>%sum
    
    exposure_ssp[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='heat',2]<-PoC_data$PoC_number[which((PoC_data$heat_ssp_class)==cl)]%>%sum
    exposure_ssp[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='drought',2]<-PoC_data$PoC_number[which((PoC_data$drought_ssp_class)==cl)]%>%sum
    exposure_ssp[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='flood',2]<-PoC_data$PoC_number[which((PoC_data$flood_ssp_class)==cl)]%>%sum
    exposure_ssp[(exposure_hist$hazard.class==cl) & exposure_hist$hazard=='compound',2]<-PoC_data$PoC_number[which((PoC_data$haz_class_ssp)==cl)]%>%sum
  }
}

spect_pal<-rev(brewer.pal(5,'Spectral'))
all_colors=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])
all_class=names(all_colors)

exposure_hist$hazard<-factor(exposure_hist$hazard)
exposure_hist$hazard.class<-factor(exposure_hist$hazard.class)
the_df<-exposure_hist[exposure_hist$PoC.number!=0,]
the_df$PoC.number<-the_df$PoC.number/1e6

tiff('./Plots/with_FDP/exposure_per_class.tiff', units="in", width=8, height=4, res=500,compression = 'lzw')
ggplot(the_df,aes(x=hazard,y=PoC.number,fill=hazard.class))+
  geom_bar(stat='identity',color = "black",linewidth=0.2,position = position_stack(reverse = TRUE))+
  labs(y='FDP number (millions)')+
  scale_fill_manual(values=all_colors,name='hazard class')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  theme(axis.text=element_text(size=13),axis.title = element_text(size = 15),legend.text = element_text(size = 13),legend.title = element_text(size = 15))+
  coord_flip()
dev.off()

exposure_ssp$hazard<-factor(exposure_ssp$hazard)
exposure_ssp$hazard.class<-factor(exposure_ssp$hazard.class)
the_df2<-exposure_ssp[exposure_ssp$PoC.number!=0,]
the_df2$PoC.number<-the_df2$PoC.number/1e6

tiff(paste0('./Plots/with_FDP/exposure_per_class_ssp',ff,'.tiff'), units="in", width=8, height=4, res=500,compression = 'lzw')
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
hist_df<-exposure_hist[exposure_hist$hazard!='compound',]
hist_df$PoC.number<-hist_df$PoC.number/1e6
hist_df$hazard.class<-factor(hist_df$hazard.class)
hist_df$hazard<-factor(hist_df$hazard,levels=c('heat','drought','flood'))

ssp_df<-exposure_ssp[exposure_ssp$hazard!='compound',]
ssp_df$PoC.number<-ssp_df$PoC.number/1e6
ssp_df$hazard.class<-factor(ssp_df$hazard.class)
ssp_df$hazard<-factor(ssp_df$hazard,levels=c('heat','drought','flood'))

my_pals=c('heat'='orangered1','drought'='goldenrod3','flood'='steelblue3')

p1<-ggplot() + geom_bar(data = hist_df, aes(x = hazard.class, y = PoC.number, fill = hazard), position = "dodge", stat = "identity",width=0.8)+
  labs(x='hazard class', y='FDP number (million)')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank())+
  theme(axis.text=element_text(size=12),axis.title = element_text(size = 12),legend.text = element_text(size = 9),legend.title = element_text(size = 12),legend.key.size = unit(.5, "cm"),
        legend.position='inside',legend.position.inside = c(0.9, 0.8),axis.text.y = element_text(angle=45,vjust = 0.5, hjust=0.5))+
  scale_x_discrete(labels=c('low','moderate','high','severe','extreme'),breaks=c("1","2","3",'4','5'))+
  scale_fill_manual( values = my_pals)

p1<-p1+coord_flip()
tiff('./Plots/with_FDP/bar_exposure_per_class_hist.tiff',units="in", width=7, height=4, res=500,compression = 'lzw')
p1
dev.off()   


p2<-ggplot() + geom_bar(data = ssp_df, aes(x = hazard.class, y = PoC.number, fill = hazard), position = "dodge", stat = "identity",width=0.8)+
  labs(x='hazard class', y='FDSP number (million)')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank())+
  theme(axis.text=element_text(size=12),axis.title = element_text(size = 12),legend.text = element_text(size = 9),legend.title = element_text(size = 12),legend.key.size = unit(.5, "cm"),
        legend.position='inside',legend.position.inside = c(0.9, 0.8),axis.text.y = element_text(angle=45,vjust = 0.5, hjust=0.5))+
  scale_x_discrete(labels=c('low','moderate','high','severe','extreme'),breaks=c("1","2","3",'4','5'))+
  scale_fill_manual( values = my_pals)
p2<-p2+coord_flip()
tiff(paste0('./Plots/with_FDP/bar_exposure_per_class_ssp',ff,'.tiff'),units="in", width=7, height=4, res=500,compression = 'lzw')
p2
dev.off() 

tiff('./Plots/with_FDP/bar_exposure_per_class.tiff',units="in", width=7, height=4, res=500,compression = 'lzw')
prow<-plot_grid(p1+theme(legend.position="none"),p2+theme(legend.position="none"),labels=c('A','B'),align = 'vh',hjust = -1,nrow = 1)
legend <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(prow, legend, rel_widths = c(2.7, .4))
dev.off()


