
whatershed_fun<-function(lon,lat,dem_file){
 
 #Fill holes in DEM 
 breached <- file.path(tempdir(), "breached.tif")
 wbt_breach_depressions_least_cost(
  dem = dem_file,
  output = breached,
  dist = 5,
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
 min_value<-cellStats(flow,"mean")
 
 raster_streams<-file.path(tempdir(), "raster_streams.tif")
 wbt_extract_streams(flow_accum = flow_accum,
                     output = raster_streams,
                     #Threshold hyperparameter: User setting??
                     threshold =min_value)
 
 outlet_distance<-file.path(tempdir(),"outlet_distance.tif")
 
 wbt_farthest_channel_head(
  d8_pntr = d8_pointer, 
  streams = raster_streams, 
  output = outlet_distance) 
 
 pour_raster<-raster(outlet_distance)
 
 coords_points<-xyFromCell(pour_raster,
                           which.max(values(pour_raster)))
 
 ppoints <- data.frame(lon = coords_points[1,1],
                       lat = coords_points[1,2])
 
 ppointsSP <- SpatialPoints(ppoints, 
                            proj4string = 
                             CRS("+proj=longlat +datum=WGS84"))
 
 pourpoints <- file.path(tempdir(), "pourpoints.shp")
 
 shapefile(ppointsSP, filename = pourpoints, 
           overwrite = TRUE)
 
 snappedpp <- file.path(tempdir(), "snappedpp.shp")
 
 #Snapped points
 wbt_jenson_snap_pour_points(pour_pts =  pourpoints,
                             streams = raster_streams,
                             output = snappedpp,
                             #Snapped distance
                             snap_dist =0.005)
 
 #Whatershed limits
 
 whatershed <- file.path(tempdir(), "whatershed.tif")
 
 wbt_watershed(d8_pntr = d8_pointer,
               pour_pts = snappedpp,
               output = whatershed)
 
 w<-raster(whatershed)
 return(w)}








