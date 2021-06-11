#####################################################Running Random Forest ###################
####rf1=occurence only
t_rf1=matrix(nrow=500,ncol=23)
t_rf1_cat=matrix(nrow=500,ncol=22)

colnames(t_rf1_cat)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","realm","family","bio2m","NA","NA","OOB","false_pos","false_neg")
colnames(t_rf1)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","trend","realm","family","bio2m","NA","NA","OOB","false_pos","false_neg")

####rf2=occurence and genetic
t_rf2=matrix(nrow=500,ncol=24)
t_rf2_cat=matrix(nrow=500,ncol=23)

colnames(t_rf2_cat)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","eoo_area","realm","family","div_avg","area","bio2m","OOB","false_pos","false_neg")
colnames(t_rf2)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","eoo_area","trend","realm","family","div_avg","area","bio2m","OOB","false_pos","false_neg")

####rf3=all variables
t_rf3=matrix(nrow=500,ncol=27)
t_rf3_cat=matrix(nrow=500,ncol=26)

colnames(t_rf3_cat)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","eoo_area","realm","family","div_avg","forearm_mm","AET","PET","area","bio2m","OOB","false_pos","false_neg")
colnames(t_rf3)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","eoo_area","trend","realm","family","div_avg","forearm_mm","AET","PET","area","bio2m","OOB","false_pos","false_neg")

####rf4=occurence and trait
t_rf4=matrix(nrow=500,ncol=26)
t_rf4_cat=matrix(nrow=500,ncol=25)

colnames(t_rf4_cat)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","eoo_area","realm","family","forearm_mm","AET","PET","area","bio2m","OOB","false_pos","false_neg")
colnames(t_rf4)=c("run","rcat","conr","bio1m","bio12m","bio4m","bio15m","elevm","elevmin","latm","popm","popmin","gdpm","aoo_area","eoo_area","trend","realm","family","forearm_mm","AET","PET","area","bio2m","OOB","false_pos","false_neg")


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


####Starting matrices for confusion matrix
t_testing=matrix(ncol = 5,nrow=5,0)
t_testing_rf1=t_testing
t_testing_rf1_cat=t_testing
t_testing_rf2=t_testing
t_testing_rf2_cat=t_testing
t_testing_rf3=t_testing
t_testing_rf3_cat=t_testing
t_testing_rf4=t_testing
t_testing_rf4_cat=t_testing

####in for loop
###rf1
for(i in 1:500){
  
  
  iucn_train_vars_lc=iucn_train_vars[iucn_train_vars$iucn_group %in% 1,]
  iucn_train_vars_ce=iucn_train_vars[iucn_train_vars$iucn_group %in% 5,]
  iucn_train_vars_nt=iucn_train_vars[iucn_train_vars$iucn_group %in% 2,]
  iucn_train_vars_vu=iucn_train_vars[iucn_train_vars$iucn_group %in% 3,]
  iucn_train_vars_en=iucn_train_vars[iucn_train_vars$iucn_group %in% 4,]
  
  ####Subsample non-threatened to have equal sample sizes
  iucn_train_vars_nt2=iucn_train_vars_nt[sample(nrow(iucn_train_vars_nt),30), ]
  iucn_train_vars_lc2=iucn_train_vars_lc[sample(nrow(iucn_train_vars_lc),30), ]
  iucn_train_vars_vu2=iucn_train_vars_vu[sample(nrow(iucn_train_vars_vu),30), ]
  iucn_train_vars2=rbind(iucn_train_vars_nt2,iucn_train_vars_lc2,iucn_train_vars_vu2,iucn_train_vars_ce,iucn_train_vars_en)
  
  #############################################################################
  ####rf1=only occurence
  rf1_cat=randomForest(as.factor(iucn_group)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+family_num+realm_num+aoo_area+realm_num+family_num+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  rf1=randomForest(as.factor(iucn_group)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+aoo_area+trend_num+realm_num+family_num+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  
  t_rf1[i,1]=i
  t_rf1[i,2:17]=rf1$importance
  t_rf1[i,20]=rf1$err.rate[1000,1]
  t_rf1[i,21]=rf1$confusion[1,3]
  t_rf1[i,22]=rf1$confusion[2,3]
  
  
  t_rf1_cat[i,1]=i
  t_rf1_cat[i,2:16]=rf1_cat$importance
  t_rf1_cat[i,19]=rf1_cat$err.rate[1000,1]
  t_rf1_cat[i,20]=rf1_cat$confusion[1,3]
  t_rf1_cat[i,21]=rf1_cat$confusion[2,3]
  
  
  ##predict for unknown (DD) Bats
  pred_try=predict(rf1,iucn_pred)
  pred_try2=predict(rf1_cat,iucn_pred)
  
  pred_rf1=cbind(pred_rf1,pred_try)
  pred_rf1_cat=cbind(pred_rf1_cat,pred_try2)
  
  
  #############################################################################
  ####rf2=occurence and genetic
  rf2_cat=randomForest(as.factor(iucn_group)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+family_num+realm_num+aoo_area+eoo_area+div_avg+area_sqkm+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  rf2=randomForest(as.factor(iucn_group)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+aoo_area+eoo_area+trend_num+realm_num+family_num+div_avg+area_sqkm+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  
  t_rf2[i,1]=i
  t_rf2[i,2:21]=rf2$importance
  t_rf2[i,22]=rf2$err.rate[1000,1]
  t_rf2[i,23]=rf2$confusion[1,3]
  t_rf2[i,24]=rf2$confusion[2,3]
  
  
  t_rf2_cat[i,1]=i
  t_rf2_cat[i,2:20]=rf2_cat$importance
  t_rf2_cat[i,21]=rf2_cat$err.rate[1000,1]
  t_rf2_cat[i,22]=rf2_cat$confusion[1,3]
  t_rf2_cat[i,23]=rf2_cat$confusion[2,3]
  
  
  ##predict for unknown (DD) Bats
  pred_try=predict(rf2,iucn_pred)
  pred_try2=predict(rf2_cat,iucn_pred)
  
  pred_rf2=cbind(pred_rf2,pred_try)
  pred_rf2_cat=cbind(pred_rf2_cat,pred_try2)
  
  
  #############################################################################
  ####rf3=all variables
  rf3_cat=randomForest(as.factor(iucn_group)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+family_num+realm_num+aoo_area+eoo_area+div_avg+forearm_length+X30.1_AET_Mean_mm+X30.2_PET_Mean_mm+area_sqkm+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  rf3=randomForest(as.factor(iucn_group)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+aoo_area+eoo_area+trend_num+realm_num+family_num+div_avg+forearm_length+X30.1_AET_Mean_mm+X30.2_PET_Mean_mm+area_sqkm+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  
  t_rf3[i,1]=i
  t_rf3[i,2:24]=rf3$importance
  t_rf3[i,25]=rf3$err.rate[1000,1]
  t_rf3[i,26]=rf3$confusion[1,3]
  t_rf3[i,27]=rf3$confusion[2,3]
  
  
  t_rf3_cat[i,1]=i
  t_rf3_cat[i,2:23]=rf3_cat$importance
  t_rf3_cat[i,24]=rf3_cat$err.rate[1000,1]
  t_rf3_cat[i,25]=rf3_cat$confusion[1,3]
  t_rf3_cat[i,26]=rf3_cat$confusion[2,3]
  
  
  ##predict for unknown (DD) Bats
  pred_try=predict(rf3,iucn_pred)
  pred_try2=predict(rf3_cat,iucn_pred)
  
  pred_rf3=cbind(pred_rf3,pred_try)
  pred_rf3_cat=cbind(pred_rf3_cat,pred_try2)
  
  
  #############################################################################
  ####rf4=occurence and trait
  rf4_cat=randomForest(as.factor(iucn_group)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+family_num+realm_num+aoo_area+eoo_area+forearm_length+X30.1_AET_Mean_mm+X30.2_PET_Mean_mm+area_sqkm+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  rf4=randomForest(as.factor(iucn_group)~binary_rcat+binary_conr+bio1m+bio12m+bio4m+bio15m+elevm+elevmin+latm+popm+popmin+gdpm+aoo_area+eoo_area+trend_num+realm_num+family_num+forearm_length+X30.1_AET_Mean_mm+X30.2_PET_Mean_mm+area_sqkm+bio2m,data = iucn_train_vars2,na.action = na.omit,ntree=1000)
  
  t_rf4[i,1]=i
  t_rf4[i,2:23]=rf4$importance
  t_rf4[i,24]=rf4$err.rate[1000,1]
  t_rf4[i,25]=rf4$confusion[1,3]
  t_rf4[i,26]=rf4$confusion[2,3]
  
  
  t_rf4_cat[i,1]=i
  t_rf4_cat[i,2:22]=rf4_cat$importance
  t_rf4_cat[i,23]=rf4_cat$err.rate[1000,1]
  t_rf4_cat[i,24]=rf4_cat$confusion[1,3]
  t_rf4_cat[i,25]=rf4_cat$confusion[2,3]
  
  
  ##predict for unknown (DD) Bats
  pred_try=predict(rf4,iucn_pred)
  pred_try2=predict(rf4_cat,iucn_pred)
  
  pred_rf4=cbind(pred_rf4,pred_try)
  pred_rf4_cat=cbind(pred_rf4_cat,pred_try2)
  
  
  
  
  ####Create giant confusion matrix
  try=rf1$confusion[,1:5]
  t_testing_rf1=t_testing_rf1+try
  try=rf1_cat$confusion[,1:5]
  t_testing_rf1_cat=t_testing_rf1_cat+try
  
  try=rf2$confusion[,1:5]
  t_testing_rf2=t_testing_rf2+try
  try=rf2_cat$confusion[,1:5]
  t_testing_rf2_cat=t_testing_rf2_cat+try
  
  try=rf3$confusion[,1:5]
  t_testing_rf3=t_testing_rf3+try
  try=rf3_cat$confusion[,1:5]
  t_testing_rf3_cat=t_testing_rf3_cat+try
  
  try=rf4$confusion[,1:5]
  t_testing_rf4=t_testing_rf4+try
  try=rf4_cat$confusion[,1:5]
  t_testing_rf4_cat=t_testing_rf4_cat+try
  
  
  
  print(i)
}


############average importance and other stuff
t_rf1_avg=colMeans(t_rf1)
t_rf1_cat_avg=colMeans(t_rf1_cat)
t_rf2_avg=colMeans(t_rf2)
t_rf2_cat_avg=colMeans(t_rf2_cat)
t_rf3_avg=colMeans(t_rf3)
t_rf3_cat_avg=colMeans(t_rf3_cat)
t_rf4_avg=colMeans(t_rf4)
t_rf4_cat_avg=colMeans(t_rf4_cat)



write.csv(t_rf1_avg,"./rf_output3/group_avg_occ_trend.csv")
write.csv(t_rf1_cat_avg,"./rf_output3/group_avg_occ_notrend.csv")

write.csv(t_rf2_avg,"./rf_output3/group_avg_occ_div_trend.csv")
write.csv(t_rf2_cat_avg,"./rf_output3/group_avg_occ_div_notrend.csv")

write.csv(t_rf3_avg,"./rf_output3/group_avg_occ_all_trend.csv")
write.csv(t_rf3_cat_avg,"./rf_output3/group_avg_occ_all_notrend.csv")

write.csv(t_rf4_avg,"./rf_output3/group_avg_occ_trait_trend.csv")
write.csv(t_rf4_cat_avg,"./rf_output3/group_avg_occ_trait_notrend.csv")


write.csv(pred_rf1,"./rf_output3/group_pred_dd.csv")