#TWI function
twi_fun<-function(dem_file){
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
 
 As_area <- file.path(tempdir(), "SCA.tif")
 
wbt_d_inf_flow_accumulation(input = dem_fill,
                            output = As_area,
                            out_type = "Specific Contributing Area")

slope <- file.path(tempdir(), "slope.tif")

wbt_slope(dem = dem_fill,
          output = slope,
          units = "degrees")

TWI_file <- file.path(tempdir(), "twi.tif")
wbt_wetness_index(sca = As_area,
                  slope = slope,
                  output = TWI_file)

TWI<-raster(TWI_file)
return(TWI)}
