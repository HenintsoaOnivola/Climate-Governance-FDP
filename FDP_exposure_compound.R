#This script contains the statistical analyses with respect to the exposure of FDP to each class of climate hazards
#with and without governance taken into account
rm(list=ls())

packages <- c("DescTools", "raster", 'tidyverse',"tibble", "ggplot2",'sf','scales','RColorBrewer','exactextractr','ggrepel',
              'spatstat','maptools','readxl','readr','ggnewscale','cartography')

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
        legend.title = element_text(size=13,face='bold'),
        legend.position = 'bottom',
        legend.box="vertical", 
        legend.margin=margin(),
        legend.key.size = unit(1.5,"line"),
        legend.text=element_text(size=12)
  )


overlayRast <- function(rst) {
  
  tbl <- raster::rasterToPoints(rst, spatial = FALSE)
  tbl <- as_tibble(tbl)
  names(tbl) <- c('x', 'y', 'value')
  tbl
}


Visu_raster_discrete_with_points<-function (output_name,title,bar_title,color_palette,the_lim,all_ticks,all_labels,the_tbl,the_points,shape_delimiter){
  tiff(output_name, units="in", width=10, height=6, res=500,compression = 'lzw')
  
  print(
    ggplot() +
      geom_tile(data = the_tbl, aes(x = x, y = y, fill = factor(value))) +
      geom_sf(data = shape_delimiter, fill = NA,lwd=0.5) +
      scale_fill_manual(values=color_palette, name=bar_title,limits=the_lim,breaks=all_ticks,labels=all_labels) +
      new_scale_fill()+  
      geom_sf(data=the_points,size=0.4,aes(fill=new_col))+
      scale_fill_manual(name='',values=c('black'),labels='Location of forcibly displaced people')+
      guides(fill = guide_legend(override.aes = list(size = 2))) +
      labs(title = title) +
      theme_MAP)
  dev.off()
}

input_fold<-'./Input/'
output_fold<-'./Output/'
dir.create('Plots/with_FDP',recursive=T)
world_shape <- st_read("./Input/world_shape2.geojson")
ff<-''


#read the previously built composite hazard indices from rasters
composite_hist_index<-raster(paste0(output_fold,"composite_hist_nrm.tif"))
composite_ssp_index<-raster(paste0(output_fold,'composite_ssp',ff,"_nrm.tif"))

composite_hist_class<-raster(paste0(output_fold,"composite_hist_class.tif"))
composite_ssp_class<-raster(paste0(output_fold,"composite_ssp",ff,"_class.tif"))


#read the governance index per country
BTI_2022<-read_excel(paste0(input_fold,"BTI_2006-2022_Scores.xlsx"),sheet='BTI 2022')
country<-as.data.frame(BTI_2022[1:137,1])
names(country)<-'country'
BTI_ind_2022<-as.data.frame(BTI_2022[,54])
names(BTI_ind_2022)<-'GI_2022'
BTI_2020<-read_excel(paste0(input_fold,"BTI_2006-2022_Scores.xlsx"),sheet='BTI 2020')
BTI_ind_2020<-as.data.frame(BTI_2020[,54])
names(BTI_ind_2020)<-'GI_2020'
BTI_2024<-read_excel(paste0(input_fold,"BTI_2024_Scores.xlsx"))
BTI_ind_2024<-as.data.frame(BTI_2024[,52])
names(BTI_ind_2024)<-'GI_2024'

BTI<-cbind(BTI_ind_2024,BTI_ind_2022,BTI_ind_2020)
BTI[BTI=='-']<-NA
BTI<-apply(BTI, 2, FUN=as.numeric)
BTI<-apply(BTI,1,FUN=mean,na.rm = TRUE)
gov_data<-data.frame(country,BTI)

#scale the governance index to account for mis-governance and create new categories regularly
neg_gov=-gov_data$BTI
mis_gov<-(neg_gov-min(neg_gov))/(max(neg_gov)-min(neg_gov))
gov_data<-cbind(gov_data,misgovernance=mis_gov)

gov_ind<-NA
misgov_ind<-NA
world_shape<-cbind(world_shape,gov_ind,misgov_ind)
alt_country_name<-c()
c_code<-c()
for (each_c in (gov_data$country)){
  the_ind<-grep(each_c,world_shape$gis_name)
  if (length(the_ind)==0){
    var=FALSE
    if (each_c=='Turkey'){
      the_ind<-grep('Türkiye',world_shape$gis_name)
      alt_country_name<-c(alt_country_name,'Türkiye')
      c_code<-c(c_code,(world_shape$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='Vietnam'){
      the_ind<-grep('Viet Nam',world_shape$gis_name)
      alt_country_name<-c(alt_country_name,'Viet Nam')
      c_code<-c(c_code,(world_shape$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='South Korea'){
      the_ind<-grep('Republic of Korea',world_shape$gis_name)
      alt_country_name<-c(alt_country_name,'Republic of Korea')
      c_code<-c(c_code,(world_shape$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='North Korea'){
      the_ind<-grep('Democratic People\'s Rep. of Korea',world_shape$gis_name)
      alt_country_name<-c(alt_country_name,'Democratic People\'s Rep. of Korea')
      c_code<-c(c_code,(world_shape$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='Congo, DR'){
      the_ind<-grep('Democratic Republic of the Congo',world_shape$gis_name)
      alt_country_name<-c(alt_country_name,'Democratic Republic of the Congo')
      c_code<-c(c_code,(world_shape$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='Congo, Rep.'){
      the_ind<-grep('Republic of the Congo',world_shape$gis_name)
      alt_country_name<-c(alt_country_name,'Republic of the Congo')
      c_code<-c(c_code,(world_shape$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='Laos'){
      the_ind<-grep('Lao People\'s Democratic Republic',world_shape$gis_name)
      alt_country_name<-c(alt_country_name,'Lao People\'s Democratic Republic')
      c_code<-c(c_code,(world_shape$iso3[the_ind])[1])
      var=TRUE
    }
    if (each_c=='Czech Republic'){
      the_ind<-grep('Czechia',world_shape$gis_name)
      alt_country_name<-c(alt_country_name,'Czechia')
      c_code<-c(c_code,(world_shape$iso3[the_ind])[1])
      var=TRUE
    }
    if (var==FALSE){
      alt_country_name<-c(alt_country_name,each_c)
      c_code<-c(c_code,NA)
    }
    
  }
  else{
    alt_country_name<-c(alt_country_name,world_shape$gis_name[the_ind[1]])
    c_code<-c(c_code,(world_shape$iso3[the_ind])[1])
  }
  
  world_shape$gov_ind[the_ind]<-gov_data$BTI[gov_data$country==each_c]
  world_shape$misgov_ind[the_ind]<-gov_data$misgovernance[gov_data$country==each_c]
}
gov_data<-cbind(gov_data,alt_country_name,c_code)
gov_data<-data.frame(gov_data$country,gov_data$alt_country_name,gov_data$c_code, gov_data$BTI,gov_data$misgovernance)
names(gov_data)<-c('country','country_name_UNHCR','country_code','governance','misgovernance')

rast_base<-composite_hist_index
values(rast_base)<-NA
Bgov_rast<-rasterize(world_shape,rast_base,field='misgov_ind')

#############################################
#build the new risk index of hazard and misgovernance (resuming with the BTI index only)
risk_hist<-(composite_hist_index+Bgov_rast)/2
risk_ssp<-(composite_ssp_index+Bgov_rast)/2
#categorize the new risk maps
brk  <- seq(0,1,0.2)
risk_hist_class <- raster::cut(risk_hist, breaks=brk) 
risk_ssp_class <- raster::cut(risk_ssp, breaks=brk) 

#saving data of the new risk 
writeRaster(risk_hist, filename = paste0(output_fold,"composite_gov_hist.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(risk_ssp, filename = paste0(output_fold,"composite_gov_ssp",ff,".tif"), format = "GTiff",overwrite=TRUE)
writeRaster(risk_hist_class, filename = paste0(output_fold,"composite_gov_hist_class.tif"), format = "GTiff",overwrite=TRUE)
writeRaster(risk_ssp_class, filename = paste0(output_fold,"composite_gov_ssp",ff,"_class.tif"), format = "GTiff",overwrite=TRUE)

#overly FDP locations on the new risk map 
PoC_points<-st_read(paste0(input_fold,'wrl_prp_p_unhcr_PoC.geojson'))
PoC_points<-PoC_points[PoC_points$pop_type %in% c("Refugee",'IDP','Asylum-seeker',"Others of concern",'Unknown'),]
PoC_points$new_col<-c('PoC')
PoC_points$new_col<-factor(PoC_points$new_col)
#visualize the climate hazard class maps, as well as the new risk class maps with location of FDP
spect_pal<-rev(brewer.pal(5,'Spectral'))
my_pals=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])

composite_hist_class_tbl<-overlayRast(composite_hist_class)
Visu_raster_discrete_with_points("./Plots/with_FDP/composite_hist_class_with_FDP.tiff","",'Hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),composite_hist_class_tbl,PoC_points,world_shape)

composite_ssp_class_tbl<-overlayRast(composite_ssp_class)
Visu_raster_discrete_with_points(paste0("./Plots/with_FDP/composite_ssp",ff,"_class_with_FDP.tiff"),"",'Hazard class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),composite_ssp_class_tbl,PoC_points,world_shape)

risk_hist_class_tbl<-overlayRast(risk_hist_class)
Visu_raster_discrete_with_points("./Plots/with_FDP/climate_gov_hist_with_FDP.tiff","",'Risk class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),risk_hist_class_tbl,PoC_points,world_shape)

risk_ssp_class_tbl<-overlayRast(risk_ssp_class)
Visu_raster_discrete_with_points(paste0("./Plots/with_FDP/climate_gov_ssp",ff,"_with_FDP.tiff"),"",'Risk class',my_pals,c('1','2','3','4','5'),c('1','2','3','4','5'),c('low','moderate','high','severe','extreme'),risk_ssp_class_tbl,PoC_points,world_shape)


#read the FDP number per country
#extract UNHCR data
PoC<-read_csv(paste0(input_fold,"persons_of_concern.csv"))
PoC_data<-data.frame(PoC)
PoC_data<-PoC_data[1:(nrow(PoC_data)-1),]
PoC_data<-PoC_data[,c(1,2,5:7,9)]
names(PoC_data)[1:2]<-c('country_name','country_code')
PoC_data[PoC_data=='-']<-NA
#extract IDP data from IDMC
IDP=as.data.frame(read_excel(paste0(input_fold,"IDMC_GIDD_Internal_Displacement_Disaggregated.xlsx")))
IDP<-IDP[IDP$'Figure category'=='IDPs',]
IDP<-IDP[which(IDP$Year==2024),c('ISO3','Country','Total figures')]
names(IDP)<-c('country_code','country_name','IDP_number')
IDP<-aggregate(IDP_number~(country_code+country_name),data=IDP,FUN=function(x) (sum(x,na.rm=TRUE)),na.action=na.pass)
IDMC_IDP<-rep(NA,nrow(PoC_data))
new_values<-c()
new_codes<-c()
new_names<-c()
for (n in 1:nrow(IDP)){
  the_ind<-which(PoC_data$country_code==IDP$country_code[n])
  if (length(the_ind)==0){
    new_values<-c(new_values,IDP$IDP_number[n])
    new_codes<-c(new_codes,IDP$country_code[n])
    new_names<-c(new_names,IDP$country_name[n])
  }else{
    IDMC_IDP[the_ind]<-IDP$IDP_number[n]
  }
}
PoC_data$IDMC_IDPs<-IDMC_IDP
the_append<-data.frame(new_names,new_codes,rep(0,length(new_names)),rep(0,length(new_names)),
                       rep(0,length(new_names)), rep(0,length(new_names)),new_values)
names(the_append)<-names(PoC_data)
PoC_data<-rbind(PoC_data,the_append)
PoC_data<-PoC_data[PoC_data$country_code!='PSE',]#remove Palestinian IDP data as this will be replaced by UNRWA data
PoC_data<-subset(PoC_data, select = -c(IDPs.of.concern.to.UNHCR) )
PoC_number<-rowSums(PoC_data[,3:ncol(PoC_data)],na.rm=TRUE)
PoC_data<-data.frame(country_name=PoC_data[,1],country_code=PoC_data[,2],PoC_number)
#get UNRWA data for Palestinian refugees
PoC_p<-as.data.frame(read_csv(paste0(input_fold,"unrwa.csv")))
p_data<-PoC_p[PoC_p$`Country of Asylum`=='State of Palestine', c(2,4,6)]
names(p_data)<-names(PoC_data)
PoC_data<-rbind(PoC_data,p_data)
for (n in 2:nrow(PoC_p)){
  ind<-which(PoC_data$country_code==PoC_p$`Country of Asylum ISO`[n])
  PoC_data[ind,3]<-PoC_data[ind,3]+PoC_p$Total[n]
}

#####
##insert PoC numbers into the world shape file
#note that the PoC numbers into the world shape file are at country level, so duplicated for territories with the same iso3 code
################
PoC_num<-NA
world_shape<-cbind(world_shape,PoC_num)
for (country_ind in (1:length(PoC_data$country_code))){
  the_ind<-which(world_shape$iso3==(PoC_data$country_code[country_ind]))
  world_shape$PoC_num[the_ind]<-PoC_data[country_ind,3]
  if (length(the_ind)==0){
    the_country<-strsplit((PoC_data$country_name)[country_ind],' ')[[1]][1]
    the_ind<-which(grepl(the_country,world_shape$gis_name))
    world_shape$PoC_num[the_ind]<-PoC_data[country_ind,3]
  }
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

##add climate (including governance-corrected) variables as well as regions into the world shape for each country
haz_hist<-c()
haz_ssp<-c()
haz_hist_class<-c()
haz_ssp_class<-c()
hazgov_hist<-c()
hazgov_ssp<-c()
hazgov_hist_class<-c()
hazgov_ssp_class<-c()
region<-c()
for (c in 1:nrow(world_shape)){
  haz_hist[c]<-exact_extract(composite_hist_index,world_shape[c,],'mean')
  haz_ssp[c]<-exact_extract(composite_ssp_index,world_shape[c,],'mean')
  hist_cl<-cut(haz_hist[c],breaks=seq(0,1,0.2),label=1:5)
  hist_cl<-as.numeric(levels(hist_cl))[hist_cl]
  haz_hist_class[c]<-hist_cl
  ssp_cl<-cut(haz_ssp[c],breaks=seq(0,1,0.2),label=1:5)
  ssp_cl<-as.numeric(levels(ssp_cl))[ssp_cl]
  haz_ssp_class[c]<-ssp_cl
  
  hazgov_hist[c]<-exact_extract(risk_hist,world_shape[c,],'mean')
  hazgov_ssp[c]<-exact_extract(risk_ssp,world_shape[c,],'mean')
  hist_cl<-cut(hazgov_hist[c],breaks=seq(0,1,0.2),label=1:5)
  hist_cl<-as.numeric(levels(hist_cl))[hist_cl]
  hazgov_hist_class[c]<-hist_cl
  ssp_cl<-cut(hazgov_ssp[c],breaks=seq(0,1,0.2),label=1:5)
  ssp_cl<-as.numeric(levels(ssp_cl))[ssp_cl]
  hazgov_ssp_class[c]<-ssp_cl
  
  region_ind<-which((combined_region$iso3)==world_shape$iso3[c])
  if (length(region_ind)!=0){
    region[c]<-(combined_region$region_name)[region_ind]
  }else{
    region[c]<-''
  }
}
world_shape$haz_hist<-haz_hist
world_shape$haz_ssp<-haz_ssp
world_shape$haz_hist_class<-haz_hist_class
world_shape$haz_ssp_class<-haz_ssp_class
world_shape$hazgov_hist<-hazgov_hist
world_shape$hazgov_ssp<-hazgov_ssp
world_shape$hazgov_hist_class<-hazgov_hist_class
world_shape$hazgov_ssp_class<-hazgov_ssp_class
world_shape$region<-region

world_shape$hazgov_hist[is.na(world_shape$misgov_ind)]<-NA
world_shape$hazgov_ssp[is.na(world_shape$misgov_ind)]<-NA
world_shape$hazgov_hist_class[is.na(world_shape$misgov_ind)]<-NA
world_shape$hazgov_ssp_class[is.na(world_shape$misgov_ind)]<-NA

#add the region, hazard intensity and class, risk intensity and class for each country within the PoC data 
#########################
PoC_data<-PoC_data[PoC_data$country_code!='UKN',]
the_hazards<-matrix(NA,nrow=nrow(PoC_data),ncol=8)
the_region<-c()
PoC_dict<-as.data.frame(read_excel(paste0(input_fold,"PoC_country_region_dictionary.xlsx")))
for (c in (1:nrow(PoC_data))){
  the_ind<-which(world_shape$iso3==PoC_data$country_code[c])
  if (length(the_ind)!=1){
    
    if (PoC_data$country_code[c]=='PRT'){
      the_ind<-the_ind[2]
      the_hazards[c,]<-as.numeric(world_shape[the_ind,(ncol(world_shape)-8):(ncol(world_shape)-1)])[1:8]
      the_region<-c(the_region,world_shape$region[the_ind])
    }
    if (PoC_data$country_code[c]=='GBR'){
      the_ind<-the_ind[1]
      the_hazards[c,]<-as.numeric(world_shape[the_ind,(ncol(world_shape)-8):(ncol(world_shape)-1)])[1:8]
      the_region<-c(the_region,world_shape$region[the_ind])
    }
    if (PoC_data$country_code[c]=='AB9'){
      the_ind<-which(world_shape$gis_name=='Abyei')
      the_hazards[c,]<-as.numeric(world_shape[the_ind,(ncol(world_shape)-8):(ncol(world_shape)-1)])[1:8]
      the_region<-c(the_region,world_shape$region[the_ind])
    }
    
    if (PoC_data$country_code[c]=='PSE'){
      the_ind<-the_ind[2]
      the_hazards[c,]<-as.numeric(world_shape[the_ind,(ncol(world_shape)-8):(ncol(world_shape)-1)])[1:8]
      the_region<-c(the_region,world_shape$region[the_ind])
    }
    if (!(PoC_data$country_code[c] %in% c('PRT','GBR','AB9','PSE'))){
      the_region<-c(the_region,PoC_dict[(PoC_dict$country)==(PoC_data$country_name)[c],2])
      proxy_iso<-PoC_dict[(PoC_dict$country)==(PoC_data$country_name)[c],3]
      if (!is.na(proxy_iso)){
        the_hazards[c,]<-as.numeric(world_shape[which(world_shape$iso3==proxy_iso),(ncol(world_shape)-8):(ncol(world_shape)-1)])[1:8]
      }else{
        the_hazards[c,]<-rep(NA,8)
      }
    }
    
    
  }else{
    the_hazards[c,]<-as.numeric(world_shape[the_ind,(ncol(world_shape)-8):(ncol(world_shape)-1)])[1:8]
    the_region<-c(the_region,world_shape$region[the_ind])
  }
}
PoC_data$region<-the_region
PoC_data<-cbind(PoC_data,as.data.frame(the_hazards))
names(PoC_data)[5:12]<-c('haz_hist','haz_ssp','haz_class_hist','haz_class_ssp',
                         'risk_value_hist','risk_value_ssp','risk_class_hist','risk_class_ssp')

saveRDS(PoC_data,paste0(output_fold,'Intermediate/FDP_population_data_frame',ff,'.Rds'))
saveRDS(st_drop_geometry(world_shape),paste0(output_fold,'Intermediate/country_hazard_values',ff,'.Rds'))

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

region_col<-hcl.colors(length(all_region),palette='Dynamic')
names(region_col)<-all_region_short
temporal<-c('hist','ssp')
all_names<-c('low','moderate','high','severe','extreme') 


for (t in (1:length(temporal))){
  class_count_ar<-c()
  pie_list<-list()
  the_names<-c()
  count_hazard<-1
  for (cla in 1:5){
    
    all_PoC_class<-eval(parse(text=paste0('PoC_data[which(PoC_data$haz_class_',temporal[t],'==cla),]')))
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
      if (t==1){
        the_pulling=5
        the_pushing=30
      }
      if (t==2){
        the_pulling=30
        the_pushing=15
      }
      
      
      p<-ggplot() +
        geom_rect(the_df,mapping=aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=as.factor(all_region_short)),color='black',lwd=0.2) +
        coord_polar(theta="y") + 
        xlim(c(3, 4))+
        scale_fill_manual(values = region_col,name='region')+
        geom_label_repel(data=the_df2,aes(x=3.5,y = y_pos, label = the_lab,fill=as.factor(all_region_short)),
                         size = 8, show.legend = FALSE,segment.alpha = 0,segment.colour="black",force=the_pushing,force_pull =the_pulling)+
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
  
  tiff(paste0('./Plots/with_FDP/FDP_number_per_class_',temporal[t],ff,'.tiff'),units="in", width=12, height=8, res=500,compression = 'lzw')
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
          legend.text = element_text(size =18 ), # Adjust legend label text size
          legend.title = element_text(size = 18, face = "bold"))
  
  print(the_bar_plot)
  dev.off()
  
}

#Plot number of FDP per regions, overly each point with distribution over hazard classes
spect_pal<-rev(brewer.pal(5,'Spectral'))
my_pals=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])
for (t in (1:length(temporal))){
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
      region_ind<-eval(parse(text=paste0('the_reg_data[which(the_reg_data$haz_class_',temporal[t],'==cla),]')))
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
                         size = 10, show.legend = FALSE,segment.alpha = 0,segment.colour="black",force=10,force_pull =1) +
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
  
  
  tiff(paste0('./Plots/with_FDP/FDP_number_per_region_',temporal[t],ff,'.tiff'),units="in", width=8, height=17, res=500,compression = 'lzw')
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
#############################################################################################
#Plot number of FDP per regions, overly each point with distribution over risk classes (with governance)
region_col<-hcl.colors(length(all_region),palette='Dynamic')
names(region_col)<-all_region_short
temporal<-c('hist')
all_names<-c('low','moderate','high','severe','extreme') 

spect_pal<-rev(brewer.pal(5,'Spectral'))
my_pals=c('1'=spect_pal[1],'2'=spect_pal[2],'3'=spect_pal[3],'4'=spect_pal[4],'5'=spect_pal[5])
for (t in (1:length(temporal))){
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
      region_ind<-eval(parse(text=paste0('the_reg_data[which(the_reg_data$risk_class_',temporal[t],'==cla),]')))
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
                         size = 12, show.legend = FALSE,segment.alpha = 0,segment.colour="black",force=20,force_pull =1) +
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
    radius.y = 10
  )
  pie_coords <- pie_coords %>%
    mutate(xmin = center.x - radius.x,
           xmax = center.x + radius.x,
           ymin = center.y - radius.y,
           ymax = center.y + radius.y)
  
  the_bar <- data.frame(
    name=factor(all_region_short_bis,levels = all_region_short_bis ),
    value=region_count_ar/1e6)
  
  
  tiff(paste0('./Plots/with_FDP/FDP_number_gov_per_region_with_bar_',temporal[t],ff,'.tiff'),units="in", width=17, height=9, res=500,compression = 'lzw')
  the_bar_plot<-ggplot()+
    geom_bar(data=the_bar, aes(x=name, y=value),stat = "identity",fill='grey59')+
    annotation_custom_list(names(pie_list),pie_coords)+
    ylab("FDP number (million)")+
    xlab('')+
    xlim(all_region_short_bis)+
    scale_y_continuous(limits = c(0,28), expand = c(0, 0))+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text.y= element_text(size = 25),
          axis.title = element_text(size = 35),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(hjust = 0.5,size = 35),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  
  print(the_bar_plot)
  dev.off()
  
}

#######
###side visualization 
###################################
# composite index per country with governance
world_shape <- st_read(paste0(input_fold,"world_shape2.geojson"))
the_data<-readRDS(paste0(output_fold,"Intermediate/country_hazard_values.Rds"))
world_shape$hazgov_hist_class<-the_data$hazgov_hist_class
world_shape$misgov_ind<-the_data$misgov_ind
misgov_class<-cut(world_shape$misgov_ind,breaks=seq(0,1,0.2),labels=1:5)
misgov_class<-as.numeric(levels(misgov_class))[misgov_class] 
world_shape$misgov_class<-misgov_class

spect_pal<-rev(brewer.pal(5,'Spectral'))

continuous_pal<-colorRampPalette(colors = c("#2B83BA", "#ABDDA4", "#FFFFBF", "#FDAE61", "#D7191C"))(50)
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
        legend.title = element_text(size=18,face='bold',vjust=1),
        legend.text = element_text(size=15),
        legend.position = "bottom"
  )
tiff("./Plots/misgov.tiff", units="in", width=10, height=6, res=500,compression = 'lzw')
ggplot()+
  geom_sf(data=world_shape,aes(fill=misgov_ind),lwd=0.5)+
  scale_fill_gradientn(colors = continuous_pal,name='governance', limits=c(0,1),label=c('very\ngood','','','','','failed'),breaks=c(0.01,0.2,0.4,0.6,0.8,0.99),na.value="white")+
  theme_MAP

dev.off()

#highlight countries with improved hazard exposure when governance is taken into account
world_shape$haz_hist_class<-the_data$haz_hist_class
world_shape$haz_hist<-the_data$haz_hist
world_shape$hazgov_hist<-the_data$hazgov_hist

world_shape$hist_diff<-world_shape$haz_hist-world_shape$hazgov_hist

neg<-rev(carto.pal('brown.pal',n1=20))
pos<-carto.pal('blue.pal',n1=20)

my_pal<-c(neg,pos)
my_break<-seq(-0.34,0.34,0.01)

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
        legend.title = element_text(size=18,face='bold',vjust=1),
        legend.text = element_text(size=15),
        legend.key.width = unit(1, "cm"),
        legend.position = "bottom"
  )

tiff("./Plots/gov_improvement_hist.tiff", units="in", width=10, height=6, res=500,compression = 'lzw')
ggplot()+
  geom_sf(data=world_shape,aes(fill=hist_diff),lwd=0.5)+
  scale_fill_gradientn(colors = my_pal,name='shift', limits=c(-0.34,0.34),breaks=c(-0.17,0,0.17),label=c('worsened','','improved'),na.value="white")+
  theme_MAP

dev.off()

#plot the regions in a world map
world_shape$region<-the_data$region
world_shape$region[which(world_shape$region=='')]<-NA
region_col<-hcl.colors(length(all_region),palette='Dynamic')
names(region_col)<-all_region
the_lab<-names(region_col)
the_lab[which(the_lab=='East and Horn of Africa')]<-'East and Horn of Africa\nand the Great Lakes'
tiff("./Plots/regions.tiff", units="in", width=10, height=6, res=500,compression = 'lzw')
ggplot()+
  geom_sf(data=world_shape,aes(fill=region),lwd=0.5)+
  scale_fill_manual(values = region_col,name='Region',na.translate = FALSE,labels=the_lab)+
  theme_MAP
dev.off()

