#Fuh model
water_yield<-function(ratio,prec){
 return(
  (1-ratio)*prec
  )
}

#Budyko curve

budyko_fun<-function(PET,w,precip){
 ratio <- 1+(PET/precip)-(1+(PET/precip)^w)^(1/w)
 return(ratio)
 }

#Potential evaporation model

pet_fun<-function(kc,ETref){
 return(kc*ETref)
}

#Kc estimation
kc_fun<-function(tree_cover){
 return(1.000901-0.009009*tree_cover)
}

#W function
w_fun<-function(lat,ndvi,twi){
 w_param<-0.69386-0.01042*lat+2.81063*ndvi+
  0.146186*twi
 return(w_param)
 }


#Climate-Soil parameter

w_fun2<-function(AWC,prec,z){
 return(
  (z*(AWC/prec))+1.25
 )}

#Volumetric plant available content
AWC_fun<-function(root_depth,PAWC){
 return(
  PAWC*root_depth
 )}









