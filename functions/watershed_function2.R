whatershed_fun<-function(lon,lat,dem_file){
 #Fill holes in DEM 
 breached <- file.path(tempdir(), "breached.tif")
 wbt_breach_depressions_least_cost(
  dem = dem_file,
  output = breached,
  dist = 5,
  min_dist = T,
  fill = TRUE)
 
 dem_fill <- file.path(tempdir(), "dem_fill.tif")
 
 wbt_fill_depressions_wang_and_liu(
  dem = breached,
  output = dem_fill)
 
 #Flow accumulation
 
 flow_accum <- file.path(tempdir(), "flow_accum.tif")
 
 wbt_d8_flow_accumulation(input = dem_fill,
                          output = flow_accum)
 
 #Pointer
 d8_pointer <- file.path(tempdir(), "d8_pointer.tif")
 
 wbt_d8_pointer(dem = dem_fill,
                output = d8_pointer)
 
 #Pour points
 
 flow<-raster(flow_accum)
 min_value<-cellStats(flow,"min")
 
 raster_streams<-file.path(tempdir(), "raster_streams.tif")
 wbt_extract_streams(flow_accum = flow_accum,
                     output = raster_streams,
                     #Threshold hyperparameter: User setting??
                     threshold =min_value)
 
 basins<-file.path(tempdir(), "basins.tif")
 wbt_basins(d8_pntr = d8_pointer,output = basins)
 
 total_basins<-raster(basins)
 
 basin_select<-raster::extract(total_basins,cbind(lon,lat))
 
 w<-calc(total_basins,fun = function(x){ifelse(x==basin_select[[1]],1,NA)})
 return(w)}
