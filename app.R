#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#

packages <- c('shiny', 'leaflet', 'sp',"raster","whitebox",
              "elevatr","stars","tidyverse","tools","vroom","rgdal",
              "htmltools","patchwork","bslib","waiter","scales")

install.packages(setdiff(packages,
                         rownames(installed.packages())))

require(leaflet)
require(shiny)
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
require(patchwork)
require(bslib)
require(waiter)
require(scales)

source(paste0(getwd(),"/functions/watershed_function2.R"))
source(paste0(getwd(),"/functions/invest_model.R"))
source(paste0(getwd(),"/functions/TWI_function.R"))



shinyApp(
 ui<-fluidPage(
  waiter::use_waiter(),
  theme = bs_theme(bootswatch = "flatly"), 
  titlePanel("Invest Model",windowTitle = "Invest Model"),
  tabsetPanel(
   tabPanel("Inicio",
            
    sidebarLayout(
     sidebarPanel(
      textInput(inputId = "value",label = "Densidad (n/ha)"),
      fileInput(inputId = "upload", label =  "Seleccionar archivo",
                accept = c(".kml")),
      sliderInput(inputId = "nsim",
                  label =  "Porcentaje de cambio de cobertura", 
                  min = 0, max = 100, value = 1,step = 10)),
    mainPanel(
     leafletOutput("map",width = "1000", height = "850")
   )
    )
   ),
  tabPanel("Resultados",
           mainPanel(
           plotOutput(outputId = "myplot",width =600, 
                      height = 650),
           p(
            outputId = "text",
             p("La herramienta para la cuantificación de servicios ecosistémicos 
        utilizados por el conjunto de modelos InVest aborda el balance hídrico 
        de una cuenca hidrográfica  mediante la función de Budyko,
        la cual relaciona el índice evaporativo (radio entre 
        la evapotranspiración real y precipitación) con el índice de aridez 
        (radio entre la evapotranspiración potencial y 
        evapotranspiración real). En este aplicativo Web se hace una 
        aproximación a partir de fuentes secundarias de información sobre
        la incidencia de los cultivos de aguacate en las demanda hídrica de 
        una cuenca y en el almacenamiento de carbono en la biomasa aérea viva.
        El aplicativo entrega información sobre el balance hídrico en escalas de
        tiempo mensuales para la cobertura actual y una simulación en la cual un
        porcentaje de bosque es reemplazada cultivos de aguacate.Los parámetros 
        de entrada hacen referencia a la densidad de siembra del cultivo de interés 
        por hectárea. La delimitación del área del cultivo (formato .kml) y 
        el porcentaje de bosque que será reemplazado cultivos de aguacate"), )
           ),
           downloadButton("download")
           
  )
  )
  ),
 
 server<-function(input,output){
  
  map<-reactive({
  req(input$upload)
  req(input$value) 
  req(input$nsim)
  
  map <- read_sf(input$upload$datapath[1])
  })
  
  
  #Reactive watershed delimitation
  invest_model<-reactive({
   waiter <- waiter::Waiter$new()
   waiter::Waiter$new(id = "myplot")$show()
   
   
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
    calc(landcover_crop,fun = function(x){
     ifelse(x %in% c(10,30,40,20),x,  NA)})
   
   labels<-c(10,30,40,20)
   kc_values<-c(1.25,0.8,0.7,0.8)
   
   landcover_crop<-resample(landcover_crop,dem,
                            method = 'ngb')
   
   kc_data<-data.frame(labels,kc_values)
   
   kc<-reclassify(landcover_crop,
                              rcl = kc_data)
   
   
   pr_roi<-resample(pr_roi,dem)*0.1
   pet_roi<-resample(pet_roi,dem)*0.01

   #WWC layers
   soil_path<-paste0(getwd(),"/WWP")
   
   awc<-stack(paste0(soil_path,"/",list.files(path = soil_path,
                         pattern = "^WWP")))
   
   awc_crop<-crop(awc,w)/100
   
   awc_layer<-(1/2)*((50-0)*(awc_crop$layer.1+awc_crop$layer.2)+
                      (150-50)*(awc_crop$layer.2+awc_crop$layer.3)+
                      (300-150)*(awc_crop$layer.3+awc_crop$layer.4))
   
   awc_layer<-resample(awc_layer,dem)
   
   
   #Invest model 
   
   w_param<-w_fun(lat = coords[,2],ndvi_mean,TWI_crop)

   
   #Kc simulations
   #Random Kc selection
   layer2<-st_zm(layer)
   
   
   area<-as.numeric(st_area(layer2))/10000
   
   cells_positions<-raster::extract(kc,as_Spatial(layer2),
                                    cellnumbers=T)
   
   values(kc)[cells_positions[[1]][,1]]<-0.75
   
   sim_values<-as.numeric(input$nsim)/100

   kc_sim<-kc
   cell_target<-which(values(kc_sim) == 1.25)
   n_cells<-length(cell_target)*sim_values
   random_cells<-sample(cell_target,size = round(n_cells))
   values(kc_sim)[random_cells] = 0.75
   kc_layers<-addLayer(kc,kc_sim)
   
   
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
                   labels = c("Actual",paste0(input$nsim,"%")),
                   ordered = T)
   
   alt_value<-raster::extract(dem,cbind(coords[,1],
                                        coords[,2]))
   
   den_value<-(1.023-0.0003*alt_value)*1000
   biomass_ha<-(0.5*round(as.numeric(500)*(den_value*0.25),2)/1000)
   
   current_biomass<-biomass_ha*area
  
   
   forest_area<-(sum(values(kc)==1.25,na.rm = T)*(10^2))/10000
   
   biomass_sim<-biomass_ha*(area+forest_area*sim_values)
   
   forest_biomass<-(forest_area*60)
   
   forest_sim<-(forest_area-forest_area*sim_values)*60
   
   
   biomass_df<-data.frame(cat = c("Aguacate","Bosque","Aguacate","Bosque"),
                          biomasa = c(current_biomass,forest_biomass,
                                      biomass_sim,
                                      forest_sim),
                          level = c("Actual","Actual","sim","sim"))
   
   biomass_df$level<-factor(biomass_df$level,
                            levels = unique(biomass_df$level),
                            labels = c("Actual",paste0(input$nsim,"%")))
   
   
   g<-ggplot(biomass_df, aes(x = cat,y = biomasa))+
    geom_bar(stat = "identity",col = "black",
             fill = "gray80")+
    ggtitle("Carbono almacenado")+
    xlab("Cobertura")+
    ylab("Carbono (Mg)")+
    facet_wrap(~level,nrow = 2,ncol = 2)+
    scale_y_continuous(
     labels = function(x) format(x, scientific = FALSE))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "top")
     
   g1<-ggplot(data = df)+
    geom_bar(aes(x = month,y = water_balance),
             position = "dodge",stat = "identity",col = "black",
             fill = "gray80")+
    ggtitle("Balance hídrico")+
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
  
  output$myplot <- renderPlot({
   g<-invest_model()
   g
   },res = 96)
  
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
