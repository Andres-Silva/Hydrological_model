#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#


require(leaflet)
require(shiny)
require(leaflet.extras)
require(leaflet.extras2)
require(sp)
require(raster)
require(whitebox)
require(elevatr)
require(stars)
require(tidyverse)
require(tools)
require(vroom)
require(rgdal)
require(htmltools)
require(tidyverse)
require(patchwork)

setwd("~/Library/CloudStorage/GoogleDrive-ahsilvad@unal.edu.co/My Drive/Invest_model/model_App/Invest_model")


source(paste0(getwd(),"/functions/watershed_function2.R"))
source(paste0(getwd(),"/functions/invest_model.R"))
source(paste0(getwd(),"/functions/TWI_function.R"))

shinyApp(
 ui<-fluidPage(
  sidebarLayout(
   sidebarPanel(textInput(inputId = "value",label = "Densidad (n/ha)",),
                fileInput(inputId = "upload",
                          label =  "Upload a file",
                          accept = c(".kml")),
                downloadButton("download")),
  mainPanel(#plotOutput("myplot")
   leafletOutput("map",width = "850", height = "850")
   ))),
 
 server<-function(input,output){
  
  map<-reactive({
  req(input$upload)
  req(input$value) 
  map <- read_sf(input$upload$datapath[1])
  })
  
  
  #Reactive watershed delimitation
  invest_model<-reactive({
   
   layer<-map() 
   coords<-layer %>% st_zm() %>% 
    st_transform(crs = "epsg:4326")%>% 
    st_centroid() %>% 
    st_coordinates()
   
   dem<-get_elev_raster(locations = data.frame(x = coords[,1],
                                               y = coords[,2]),
                        z = 14,prj = 4326,
                        expand = 0.0285,clip = "bbox",src = "aws")
   
   writeRaster(dem,"dem.tif",format = "GTiff",overwrite = T)
   
   w<-whatershed_fun(lon = coords[,1],
                  lat = coords[,2],
                  "dem.tif")
  
   #Read Covariable layers
   dem<-raster("dem.tif")
   
  
   
   dir_prec<-paste0(getwd(),"/precipitacion")
   dir_pet<-paste0(getwd(),"/evapotranspiracion")
   
  
   ndvi<-raster(paste0(getwd(),"/NDVI/ndvi.tif"))
   
   pr<-raster::stack(paste0(dir_prec,"/",list.files(path =dir_prec)))
   pet<-raster::stack(paste0(dir_pet,"/",list.files(path =dir_pet)))
   landcover<-raster(paste0(getwd(),"/Landcover/landcover_2022.tif"))
   
   
   #Processing TWI variable
   TWI<-twi_fun(dem_file = "dem.tif")

   
   #Region of interest layers (ROI)
   pr_roi<-crop(pr,w)
   pet_roi<-crop(pet,w)
   
   ndvi_crop<-crop(ndvi,w)
   TWI_crop<-crop(TWI,w)
   
   ndvi_mean<-cellStats(ndvi_crop,"mean")
   
   e<-extent(w)
   lat<-(e@ymax+e@ymin)/2
   
   #Kc approximation from sentinel landcover
   
   landcover_crop<-crop(landcover,w)
   landcover_crop<-
    calc(landcover_crop,fun = function(x){ifelse(x %in% c(10,30,40,20),x,
                                                 NA)})
   
   labels<-c(10,30,40,20)
   kc_values<-c(1.25,0.8,0.7,0.8)
   
   kc_data<-data.frame(labels,kc_values)
   
   for(i in 1:nrow(kc_data)){
    values(landcover_crop)[
     values(landcover_crop)==kc_data[i,1]]<-kc_data[i,2]
   }
   
   
   kc<-resample(landcover_crop,dem,method = "ngb")
   
   pr_roi<-resample(pr_roi,dem)*0.1
   pet_roi<-resample(pet_roi,dem)*0.01

   #WWC layers
   soil_path<-paste0(getwd(),"/WWP")
   
   awc<-stack(paste0(soil_path,"/",list.files(path = soil_path,
                         pattern = "^WWP")))
   
   awc_crop<-crop(awc,w)/100
   
   awc_layer<-(1/2)*((50-0)*(awc_crop$WWP_M_sl1_250m_ll+awc_crop$WWP_M_sl2_250m_ll)+
                      (150-50)*(awc_crop$WWP_M_sl2_250m_ll+awc_crop$WWP_M_sl3_250m_ll)+
                      (300-150)*(awc_crop$WWP_M_sl3_250m_ll+awc_crop$WWP_M_sl4_250m_ll))
   
   awc_layer<-resample(awc_layer,dem)
   
   
   #Invest model 
   
   w_param<-w_fun(lat = coords[,2],ndvi_mean,TWI_crop)

   
   #Kc simulations
   #Random Kc selection
   layer2<-st_zm(layer)
   
   cells_positions<-raster::extract(kc,as_Spatial(layer2),
                                    cellnumbers=T)
   
   values(kc)[cells_positions[[1]][,1]]<-0.75
   
   sim_values<-c(0.5,0.75,1)
   kc_layers<-stack(kc)
   
   for(i in 1:length(sim_values)){
    kc_sim<-kc
    cell_target<-which(values(kc_sim) == 1.25)
    n_cells<-length(cell_target)*sim_values[[i]]
    random_cells<-sample(cell_target,size = round(n_cells))
    values(kc_sim)[random_cells] = 0.75
    kc_layers<-addLayer(kc_layers,kc_sim)
    }
   
   
   df<-NULL
   for(j in 1:nlayers(kc_layers)){
    etp_ref<-pet_fun(kc = kc_layers[[j]],ETref = pet_roi)
    
    ratio_temp<-budyko_fun(PET = etp_ref,
                            precip = pr_roi, w_param)
    
   state<-water_yield(ratio_temp,pr_roi)
   state<-state-awc_layer
   
   water_balance<-cellStats(state,"mean")
   prec_accum<-cellStats(pr_roi,"mean")
   etp_accum<-cellStats(pet_roi,"mean")
   
   df<-rbind(df,
             data.frame(water_balance = water_balance,prec_accum,
                        etp_accum,month = 1:12,
                        type =j))
   }
  
   
   df$type<-factor(df$type,levels = unique(df$type),
                   labels = c("Actual","50%","75%","100%"),ordered = T)
   
   alt_value<-raster::extract(dem,cbind(coords[,1],
                                        coords[,2]))
   
   den_value<-(1.023-0.0003*alt_value)*1000
   biomass<-0.5*round(as.numeric(input$value)*((den_value*0.25)/1000),2)
   
   biomass_df <- data.frame(
    category=c("Bosque Natural", "Aguacate"),
    count=c(96,biomass))
   
   biomass_df$ymax<-cumsum(biomass_df$count)[1:2]
   biomass_df$ymin<-c(0,96)
   biomass_df$labelPosition <- (biomass_df$ymax + biomass_df$ymin) / 2
   biomass_df$label <- paste0(biomass_df$category,":",
                              "\n", biomass_df$count)
   
   g<-ggplot(biomass_df, aes(ymax=ymax, ymin=ymin, xmax=4, 
                       xmin=3, fill=category)) +
    geom_rect() +
    geom_label( x=3.5, aes(y=labelPosition, label=label), 
                size=4) +
    ggtitle("Carbono")+
    coord_polar(theta="y")+
    scale_fill_brewer(name = "",type = "qual",palette = "Greens")+
    xlim(c(2, 4))+
    theme_void()+
    theme(
     legend.position = "none")
     
   
   df$p_ratio<-df$prec_accum/1000
   df$etp_ratio<-df$etp_accum/1000
   
   g1<-ggplot(data = df)+
    geom_bar(aes(x = month,y = water_balance),
             position = "dodge",stat = "identity",col = "black",fill = "gray80")+
    geom_segment(aes(y = prec_accum-80,yend = prec_accum-80,
                     x = -1*etp_ratio+month,xend = month,
                     group = type),
                 stat = "identity",col = "red",linewidth = 1.5)+
    geom_segment(aes(y = prec_accum-80,yend = prec_accum-80,
                     x = p_ratio+month,xend = month,
                     group = type),
                 stat = "identity",col = "blue",linewidth = 1.5)+
    scale_x_continuous(breaks = seq(1,12))+
    scale_linetype(name = "")+
   facet_wrap(~type,nrow = 2,ncol = 2)+
    ylab("Balance hídrico (mm)")+
    xlab("Mes")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "top")
   g1
   
   g1/g
  }
  )
  
  #output$myplot <- renderPlot({
   #invest_model()
   #})
  
  output$map<-renderLeaflet({
   points_df<-read.csv2("results_yield.csv")
   leaflet(data = points_df) %>% addProviderTiles(providers$Esri) %>% 
    setView(lng = -74.5,lat = 6,zoom = 7) %>% 
    addMarkers(lng = ~LONGITUD,lat = ~LATITUD,
               clusterOptions = markerClusterOptions())
    
  })
  
  output$download <- downloadHandler(
   filename = function() {
    paste0("water_yiled", ".pdf")
   },
   content = function(file) {
    ggsave(filename = file,plot = invest_model(),
           device = "pdf",dpi = 300,width = 25,height = 25,
           units = "cm")
   })
  }
 )


# Run the application 
shinyApp(ui = ui, server = server)
