library(rgbif)
library(CoordinateCleaner)
library(ConR)
library(rCAT)
library(randomForest)
library(raster)


etwd("C:/Users/wieringa.3/Desktop/Rasters for TMC")

IUCN_data=read.csv("./data_deficient_bats.csv")
IUCN_test=matrix(nrow=nrow(IUCN_data),ncol=27)
colnames(IUCN_test)=c("species","IUCN_cat","rcat_area_cat","rcat_extent_cat","conr_area_cat","conr_extent_cat","conr_cat","simp_iucn_cat","binary_cat","bio1m","bio12m","bio4m","bio15m","elevm","elevmax","elevmin","latm","popm","popmin","gdpm","genus","family","trend","realm","aoo_area","eoo_area",bio2m)
IUCN_test=as.data.frame(IUCN_test)

pop=raster("./Population.tif")
gdp=raster("./GDP.tif")
bio1w = raster("./wc2.0_bio_5m_01.tif")
bio12w = raster("./wc2.0_bio_5m_12.tif")
elevation = raster("GDEM-10km-BW.tif")
bio4w = raster("./wc2.0_bio_5m_04.tif")
bio15w = raster("./wc2.0_bio_5m_15.tif")
r1=raster(ext=extent(bio1w),res=1)
bio2w=raster("./wc2.0_bio_5m_02.tif")

###For loop starting here
for(i in 1:1034){

  ####Read in gbif data
  sp=as.character(IUCN_data[i,3])
  IUCN_test[i,1]=sp
  IUCN_test[i,2]=as.character(IUCN_data[i,4])
  t1=rgbif::occ_search(scientificName = sp, hasCoordinate = TRUE,limit=200000, fields = c('decimalLatitude','decimalLongitude','basisOfRecord','name'),hasGeospatialIssue = FALSE)
  d=as.data.frame(t1$data)
  ####Remove fossils or unknown records
  d=d[!d$basisOfRecord=="UNKNOWN",]
  d=d[!d$basisOfRecord=="FOSSIL_SPECIMEN",]

  colnames(d)=c("species","basisofrecord","decimallongitude","decimallatitude")
  ####Remove problem data points
  flags <- clean_coordinates(x = d, tests = c("capitals","centroids", "equal", "gbif", "institutions","seas","zeros","outl"))
  d_clean=flags[!flags$.summary=="FALSE",]

  d_clean_coords=d_clean[,3:4] 

  ####Removing outliers on 4 STDev
  min_lat=median(d_clean_coords$decimallatitude)-(4*sd(d_clean_coords$decimallatitude))
  min_lon=median(d_clean_coords$decimallongitude)-(4*sd(d_clean_coords$decimallongitude))

  max_lat=median(d_clean_coords$decimallatitude)+(4*sd(d_clean_coords$decimallatitude))
  max_lon=median(d_clean_coords$decimallongitude)+(4*sd(d_clean_coords$decimallongitude))
  df=d_clean_coords
  df[,1][df[,1] < min_lon] <- NA
  df=na.omit(df)

  df[,1][df[,1] > max_lon] <- NA
  df=na.omit(df)


  df[,2][df[,2] < min_lat] <- NA
  df=na.omit(df)

  df[,2][df[,2] > max_lat] <- NA
  df=na.omit(df)

  ##Thinning data to make this computationally possible (5 per 1deg)
  try12a <- gridSample(df, r1, n=5)
  d_clean_coords1=try12a

  #reorder columns for analysis
  d_clean_coords2=matrix(nrow=nrow(d_clean_coords1),ncol=3)
  
  d_clean_coords2=d_clean_coords1[,c(2,1)]
  d_clean_coords2[,3]=sp
  d_clean_coords2=as.data.frame(d_clean_coords2)
  
  EOO_area=EOO.computing(d_clean_coords2)
  AOO_area=AOO.computing(d_clean_coords2)
  ####rCAT
  risk_area_rcat=AOORating(AOO_area)
  IUCN_test[i,3]=risk_area_rcat
  
  risk_extent_rcat=EOORating(EOO_area)
  IUCN_test[i,4]=risk_extent_rcat
  ####ConR
  eval_conr=IUCN.eval(d_clean_coords2)
  IUCN_test[i,5]=eval_conr$Category_AOO
  IUCN_test[i,6]=eval_conr$Category_EOO
  IUCN_test[i,7]=eval_conr$Category_code
  ####RandomForest
  #EOO/AOO
  df1=data.frame(d_clean_coords1)
  IUCN_test[i,8]=as.character(IUCN_data[i,5])
  IUCN_test[i,9]=IUCN_data[i,6]
  
  e1<-raster::extract(bio1w, df1)
  e1<-as.numeric(unlist(e1))
  bio1m<-mean(e1, na.rm=TRUE)
  IUCN_test[i,10]=bio1m
  
  e12<-extract(bio12w, df1)
  e12<-as.numeric(unlist(e12))
  bio12m<-mean(e12, na.rm=TRUE)
  IUCN_test[i,11]=bio12m
  
  e2<-extract(bio2w, df1)
  e2<-as.numeric(unlist(e2))
  bio2m<-mean(e2, na.rm=TRUE)
  IUCN_test[i,27]=bio2m
  
  e4<-extract(bio4w, df1)
  e4<-as.numeric(unlist(e4))
  bio4m<-mean(e4, na.rm=TRUE)
  IUCN_test[i,12]=bio4m
  
  e15<-extract(bio15w, df1)
  e15<-as.numeric(unlist(e15))
  bio15m<-mean(e15, na.rm=TRUE)
  IUCN_test[i,13]=bio15m
  
  e_elev<-extract(elevation, df1)
  e_elev<-as.numeric(unlist(e_elev))
  elevm<-mean(e_elev, na.rm=TRUE)
  elevmin<-min(e_elev, na.rm=TRUE) 
  elevmax<-max(e_elev, na.rm=TRUE)
  IUCN_test[i,14]=elevm
  IUCN_test[i,15]=elevmax
  IUCN_test[i,16]=elevmin
  
  latm=mean(df1$decimallatitude)
  IUCN_test[i,17]=latm
  
  e_pop<-raster::extract(pop, df1)
  e_pop<-as.numeric(unlist(e_pop))
  popm<-mean(e_pop, na.rm=TRUE)
  popmin<-min(e_pop, na.rm=TRUE) 
  IUCN_test[i,18]=popm
  IUCN_test[i,19]=popmin
  
  e_gdp<-extract(gdp, df1)
  e_gdp<-as.numeric(unlist(e_gdp))
  gdpm<-mean(e_gdp, na.rm=TRUE)
  IUCN_test[i,20]=gdpm
  
  IUCN_test[i,21]=as.character(IUCN_data[i,29])
  IUCN_test[i,22]=as.character(IUCN_data[i,28])
  IUCN_test[i,23]=as.character(IUCN_data[i,13])
  IUCN_test[i,24]=as.character(IUCN_data[i,18])
  IUCN_test[i,25]=as.numeric(AOO_area)
  IUCN_test[i,26]=as.numeric(EOO_area)
  
}



#####################################################Running Random Forest ###################
####rf1=occurence only
t_rf1=matrix(nrow=500,ncol=27)
t_rf1_cat=matrix(nrow=500,ncol=26)

colnames(t_rf1_cat)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","family","Realm","AOO_Area","bio2m","NA","NA","OOB","false_pos","false_neg","TN","FP","FN","TN")
colnames(t_rf1)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","trend","family","Realm","AOO_area","bio2m","NA","NA","OOB","false_pos","false_neg","TN","FP","FN","TN")

####rf2=occurence and genetic
t_rf2=matrix(nrow=500,ncol=28)
t_rf2_cat=matrix(nrow=500,ncol=27)

colnames(t_rf2_cat)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","eoo_area","realm","family","div_avg","area","bio2m","OOB","false_pos","false_neg","TN","FP","FN","TN")
colnames(t_rf2)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","eoo_area","trend","realm","family","div_avg","area","bio2m","OOB","false_pos","false_neg","TN","FP","FN","TN")

####rf3=all variables
t_rf3=matrix(nrow=500,ncol=31)
t_rf3_cat=matrix(nrow=500,ncol=30)

colnames(t_rf3_cat)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","eoo_area","realm","family","div_avg","forearm_mm","AET","PET","area","bio2m","OOB","false_pos","false_neg","TN","FP","FN","TN")
colnames(t_rf3)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","eoo_area","trend","realm","family","div_avg","forearm_mm","AET","PET","area","bio2m","OOB","false_pos","false_neg","TN","FP","FN","TN")

####rf4=occurence and trait
t_rf4=matrix(nrow=500,ncol=30)
t_rf4_cat=matrix(nrow=500,ncol=29)

colnames(t_rf4_cat)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","eoo_area","realm","family","forearm_mm","AET","PET","area","bio2m","OOB","false_pos","false_neg","TN","FP","FN","TN")
colnames(t_rf4)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","eoo_area","trend","realm","family","forearm_mm","AET","PET","area","bio2m","OOB","false_pos","false_neg","TN","FP","FN","TN")


####All realm/trend/family variables converted to number manually (see )
iucn_train=read.csv("C:/Users/wieringa.3/Desktop/Rasters for TMC/IUCN_testing_known_nucdiv_trait.csv")
iucn_train_res=iucn_train$binary_cat
iucn_train_vars=iucn_train

##REad in unknown bats
##Again cats had to be switched to numbers
iucn_pred=read.csv("C:/Users/wieringa.3/Desktop/Rasters for TMC/IUCN_testing_dd_nucdiv_trait2.csv")

########cat doesn't include population trend
###rf1=only occurence derived vars
pred_rf1_cat=as.data.frame(iucn_pred$species)
pred_rf1=as.data.frame(iucn_pred$species)

###rf2=occurence and genetic
pred_rf2_cat=as.data.frame(iucn_pred$species)
pred_rf2=as.data.frame(iucn_pred$species)

###rf3=all variables
pred_rf3_cat=as.data.frame(iucn_pred$species)
pred_rf3=as.data.frame(iucn_pred$species)

###rf4=occurence and trait
pred_rf4_cat=as.data.frame(iucn_pred$species)
pred_rf4=as.data.frame(iucn_pred$species)


species=as.data.frame(iucn_train$species)
colnames(species)="species"

pred_knw=as.data.frame(species)
pred_knw=cbind(pred_knw,iucn_train$iucn_group,iucn_train$X)
colnames(pred_knw)=c("species","iucn_group","X")

for(i in 1:500){
  
  
  iucn_train_vars_tr=iucn_train_vars[iucn_train_vars$binary_cat %in% 1,]
  iucn_train_vars_nt=iucn_train_vars[iucn_train_vars$binary_cat %in% 0,]
  
  ####Subsample non-threatened to have equal sample sizes
  iucn_train_vars_nt2=iucn_train_vars_nt[sample(nrow(iucn_train_vars_nt),107), ]
  iucn_train_vars2=rbind(iucn_train_vars_nt2,iucn_train_vars_tr)
  
  
  #############################################################################
  ####rf1=only occurence
  rf1_cat=randomForest(as.factor(binary_cat)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+family_num+realm_num+aoo_area+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  rf1=randomForest(as.factor(binary_cat)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+trend_num+family_num+realm_num+aoo_area+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  
  t_rf1[i,1]=i
  t_rf1[i,2:18]=rf1$importance
  t_rf1[i,21]=rf1$err.rate[1000,1]
  t_rf1[i,22]=rf1$confusion[1,3]
  t_rf1[i,23]=rf1$confusion[2,3]
  t_rf1[i,24]=rf1$confusion[1,1]
  t_rf1[i,25]=rf1$confusion[1,2]
  t_rf1[i,26]=rf1$confusion[2,1]
  t_rf1[i,27]=rf1$confusion[2,2]
  
  
  t_rf1_cat[i,1]=i
  t_rf1_cat[i,2:17]=rf1_cat$importance
  t_rf1_cat[i,20]=rf1_cat$err.rate[1000,1]
  t_rf1_cat[i,21]=rf1_cat$confusion[1,3]
  t_rf1_cat[i,22]=rf1_cat$confusion[2,3]
  t_rf1_cat[i,23]=rf1_cat$confusion[1,1]
  t_rf1_cat[i,24]=rf1_cat$confusion[1,2]
  t_rf1_cat[i,25]=rf1_cat$confusion[2,1]
  t_rf1_cat[i,26]=rf1_cat$confusion[2,2]
  
  
  ##predict for unknown (DD) Bats
  pred_try=predict(rf1,iucn_pred)
  pred_try2=predict(rf1_cat,iucn_pred)
  
  #pred_test2=as.data.frame(iucn_train_vars2$species)
  #pred_test2=cbind(pred_test2,iucn_train_vars2$X)
  pred_test=as.data.frame(rf1$predicted)
  pred_test$X <- rownames(pred_test)

  pred_rf1=cbind(pred_rf1,pred_try)
  pred_rf1_cat=cbind(pred_rf1_cat,pred_try2)
  
  pred_knw=merge(pred_knw,pred_test,by="X",all=TRUE)
  
  #############################################################################
  ####rf2=occurence and genetic
  rf2_cat=randomForest(as.factor(binary_cat)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+family_num+realm_num+aoo_area+eoo_area+div_avg+area_sqkm+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  rf2=randomForest(as.factor(binary_cat)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+aoo_area+eoo_area+trend_num+realm_num+family_num+div_avg+area_sqkm+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  
  t_rf2[i,1]=i
  t_rf2[i,2:21]=rf2$importance
  t_rf2[i,22]=rf2$err.rate[1000,1]
  t_rf2[i,23]=rf2$confusion[1,3]
  t_rf2[i,24]=rf2$confusion[2,3]
  t_rf2[i,25]=rf2$confusion[1,1]
  t_rf2[i,26]=rf2$confusion[1,2]
  t_rf2[i,27]=rf2$confusion[2,1]
  t_rf2[i,28]=rf2$confusion[2,2]
  
  
  t_rf2_cat[i,1]=i
  t_rf2_cat[i,2:20]=rf2_cat$importance
  t_rf2_cat[i,21]=rf2_cat$err.rate[1000,1]
  t_rf2_cat[i,22]=rf2_cat$confusion[1,3]
  t_rf2_cat[i,23]=rf2_cat$confusion[2,3]
  t_rf2_cat[i,24]=rf2_cat$confusion[1,1]
  t_rf2_cat[i,25]=rf2_cat$confusion[1,2]
  t_rf2_cat[i,26]=rf2_cat$confusion[2,1]
  t_rf2_cat[i,27]=rf2_cat$confusion[2,2]
  
  
  ##predict for unknown (DD) Bats
  pred_try=predict(rf2,iucn_pred)
  pred_try2=predict(rf2_cat,iucn_pred)
  
  pred_rf2=cbind(pred_rf2,pred_try)
  pred_rf2_cat=cbind(pred_rf2_cat,pred_try2)
  
  
  #############################################################################
  ####rf3=all variables
  rf3_cat=randomForest(as.factor(binary_cat)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+family_num+realm_num+aoo_area+eoo_area+div_avg+forearm_length+X30.1_AET_Mean_mm+X30.2_PET_Mean_mm+area_sqkm+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  rf3=randomForest(as.factor(binary_cat)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+aoo_area+eoo_area+trend_num+realm_num+family_num+div_avg+forearm_length+X30.1_AET_Mean_mm+X30.2_PET_Mean_mm+area_sqkm+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  
  t_rf3[i,1]=i
  t_rf3[i,2:24]=rf3$importance
  t_rf3[i,25]=rf3$err.rate[1000,1]
  t_rf3[i,26]=rf3$confusion[1,3]
  t_rf3[i,27]=rf3$confusion[2,3]
  t_rf3[i,28]=rf3$confusion[1,1]
  t_rf3[i,29]=rf3$confusion[1,2]
  t_rf3[i,30]=rf3$confusion[2,1]
  t_rf3[i,31]=rf3$confusion[2,2]
  
  
  t_rf3_cat[i,1]=i
  t_rf3_cat[i,2:23]=rf3_cat$importance
  t_rf3_cat[i,24]=rf3_cat$err.rate[1000,1]
  t_rf3_cat[i,25]=rf3_cat$confusion[1,3]
  t_rf3_cat[i,26]=rf3_cat$confusion[2,3]
  t_rf3_cat[i,27]=rf3_cat$confusion[1,1]
  t_rf3_cat[i,28]=rf3_cat$confusion[1,2]
  t_rf3_cat[i,29]=rf3_cat$confusion[2,1]
  t_rf3_cat[i,30]=rf3_cat$confusion[2,2]
  
  
  ##predict for unknown (DD) Bats
  pred_try=predict(rf3,iucn_pred)
  pred_try2=predict(rf3_cat,iucn_pred)
  
  pred_rf3=cbind(pred_rf3,pred_try)
  pred_rf3_cat=cbind(pred_rf3_cat,pred_try2)
  
  
  #############################################################################
  ####rf4=occurence and trait
  rf4_cat=randomForest(as.factor(binary_cat)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+family_num+realm_num+aoo_area+eoo_area+forearm_length+X30.1_AET_Mean_mm+X30.2_PET_Mean_mm+area_sqkm+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  rf4=randomForest(as.factor(binary_cat)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+aoo_area+eoo_area+trend_num+realm_num+family_num+forearm_length+X30.1_AET_Mean_mm+X30.2_PET_Mean_mm+area_sqkm+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  
  t_rf4[i,1]=i
  t_rf4[i,2:23]=rf4$importance
  t_rf4[i,24]=rf4$err.rate[1000,1]
  t_rf4[i,25]=rf4$confusion[1,3]
  t_rf4[i,26]=rf4$confusion[2,3]
  t_rf4[i,27]=rf4$confusion[1,1]
  t_rf4[i,28]=rf4$confusion[1,2]
  t_rf4[i,29]=rf4$confusion[2,1]
  t_rf4[i,30]=rf4$confusion[2,2]
  
  
  t_rf4_cat[i,1]=i
  t_rf4_cat[i,2:22]=rf4_cat$importance
  t_rf4_cat[i,23]=rf4_cat$err.rate[1000,1]
  t_rf4_cat[i,24]=rf4_cat$confusion[1,3]
  t_rf4_cat[i,25]=rf4_cat$confusion[2,3]
  t_rf4_cat[i,26]=rf4_cat$confusion[1,1]
  t_rf4_cat[i,27]=rf4_cat$confusion[1,2]
  t_rf4_cat[i,28]=rf4_cat$confusion[2,1]
  t_rf4_cat[i,29]=rf4_cat$confusion[2,2]
  
  
  ##predict for unknown (DD) Bats
  pred_try=predict(rf4,iucn_pred)
  pred_try2=predict(rf4_cat,iucn_pred)
  
  pred_rf4=cbind(pred_rf4,pred_try)
  pred_rf4_cat=cbind(pred_rf4_cat,pred_try2)
  
  
  
  
  print(i)
}

write.csv(pred_knw,"./rf_output_redo/known_pred_redo.csv")

write.csv(t_rf1,"./rf_output_redo/known_rf1.csv")
write.csv(t_rf1_cat,"./rf_output_redo/known_rf1_cat.csv")
write.csv(t_rf2,"./rf_output_redo/known_rf2.csv")
write.csv(t_rf2_cat,"./rf_output_redo/known_rf2_cat.csv")
write.csv(t_rf3,"./rf_output_redo/known_rf3.csv")
write.csv(t_rf3_cat,"./rf_output_redo/known_rf3_cat.csv")
write.csv(t_rf4,"./rf_output_redo/known_rf4.csv")
write.csv(t_rf4_cat,"./rf_output_redo/known_rf4_cat.csv")


