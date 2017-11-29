all = as.matrix(read.csv("P:/MAPPING/Model_data/alldat.csv", header=T, row.names=1))
rownames(all)[which(rownames(all)=="Ruminococcus_bircirculans_80_3")] = "Ruminococcus_bicirculans_80_3"

allex = vector()
for(i in 1:length(specs)){
  allex = union(allex,specs[[i]]@react_id[grep("EX_",specs[[i]]@react_id)])
}
muc = c("EX_acgam(e)","EX_cspg_b(e)","EX_cspg_c(e)","EX_cspg_a_degr(e)","EX_cspg_b_degr(e)","EX_cspg_c_degr(e)","EX_cspg_ab_rest(e)","EX_cspg_c_rest(e)","EX_acgalglcur(e)","EX_acgalidour(e)","EX_acgalidour2s(e)","EX_acgal(e)","EX_hspg(e)","EX_hspg_degr_1(e)","EX_hspg_degr_2(e)","EX_hspg_degr_3(e)","EX_hspg_degr_4(e)","EX_hspg_degr_5(e)","EX_hspg_degr_6(e)","EX_hspg_degr_7(e)","EX_hspg_degr_8(e)","EX_hspg_degr_9(e)","EX_hspg_degr_10(e)","EX_hspg_degr_11(e)","EX_hspg_degr_12(e)","EX_hspg_degr_13(e)","EX_hspg_degr_14(e)","EX_hspg_degr_15(e)","EX_hspg_rest(e)","EX_homogal(e)","EX_ha(e)","EX_ha_deg1(e)","EX_ha_pre1(e)","EX_T_antigen(e)","EX_Tn_antigen(e)","EX_core2(e)","EX_core3(e)","EX_core4(e)","EX_core5(e)","EX_core6(e)","EX_core7(e)","EX_core8(e)","EX_dsT_antigen(e)","EX_f1a(e)","EX_gncore1(e)","EX_gncore2(e)","EX_sT_antigen(e)","EX_sTn_antigen(e)","EX_acnam(e)")
gas = c("EX_aso3(e)","EX_aso4(e)","EX_ca2(e)","EX_cd2(e)","EX_cl(e)","EX_co2(e)","EX_cobalt2(e)","EX_cro4(e)","EX_cu2(e)","EX_fe2(e)","EX_h(e)","EX_h2o(e)","EX_hg2(e)","EX_k(e)","EX_mg2(e)","EX_mn2(e)","EX_na1(e)","EX_nh4(e)","EX_o2(e)","EX_pb(e)","EX_pi(e)","EX_so4(e)","EX_zn2(e)","EX_h2s(e)","EX_ppi(e)","EX_h2(e)","EX_no2(e)","EX_meoh(e)","EX_no3(e)","EX_tma(e)","EX_tmao(e)","EX_tsul(e)","EX_dms(e)","EX_dmso(e)","EX_selni(e)","EX_sel(e)","EX_mobd(e)","EX_n2(e)","EX_butso3(e)","EX_etha(e)","EX_ethso3(e)")
miss = setdiff(allex,c(muc,gas))
misstable = rbind(cbind(gas,rep("Gas",length(gas))),
                  cbind(muc,rep("Mucin",length(muc))),
                  cbind(miss,rep("",length(miss))))
#write.csv(misstable,file="P:/MAPPING/Model_data/Medium_class.csv")


diets = read.csv("P:/MAPPING/Model_data/Diet_all.csv",row.names=1,header=T)
essnut = as.character(read.csv("P:/MAPPING/Model_data/essnutall_mingrow0_1_1000lb.csv",header=F)[,1])
mclass = read.csv("P:/MAPPING/Model_data/Medium_class.csv",row.names=1,header=F)

final_diet = matrix(0,nrow=nrow(mclass),ncol=5)
rownames(final_diet) = rownames(mclass)
colnames(final_diet) = c("EU","Vege","Medit","Dach","diffconstant")

final_diet[essnut,c("EU","Vege","Medit","Dach")] = 1
final_diet[intersect(rownames(final_diet),rownames(diets)),c("EU","Vege","Medit","Dach")] = 
  as.matrix(diets[intersect(rownames(final_diet),rownames(diets)),c("EU","Vegetarian","Mediterranean","DACH")])
final_diet[rownames(mclass)[which(mclass$V2=="Gas")],"diffconstant"] = 20*10^(-6)
final_diet[rownames(mclass)[which(mclass$V2=="Organic")],"diffconstant"] = 6.7*10^(-6)

#write.csv(final_diet,"P:/MAPPING/Model_data/final_diet.csv")
#write.csv(rownames(mclass)[which(mclass$V2=="Mucin")],"P:/MAPPING/Model_data/mucin.csv")


alldiet = read.csv("P:/MAPPING/Model_data/diets.csv", header=T, row.names=1)


#######################################
############## compiling diet for unhealthy diet

all = as.matrix(read.csv("P:/MAPPING/Pediatric_crohns/alldat773_pediatric_filtered.csv", header=T, row.names=1))
load("P:/MAPPING/Models773/agora773_constrain.RData")

allex = vector()
for(i in 1:length(specs)){
  allex = union(allex,specs[[i]]@react_id[grep("EX_",specs[[i]]@react_id)])
}

diets = read.csv("P:/MAPPING/Pediatric_crohns/fluxes.csv",header=T)
rownames(diets) = as.character(diets$reaction)
essnut = as.character(read.csv("P:/MAPPING/Model_data/essnutall_mingrow0_1_1000lb.csv",header=F)[,1])
mclass = read.csv("P:/MAPPING/Model_data/Medium_class.csv",row.names=1,header=F)

diet = data.frame(reaction=allex, unclass="diet", Unhealthy=0)
rownames(diet) = allex
diet$unclass = as.character(diet$unclass)
diet[essnut,"unclass"] = "essential"
diet[intersect(essnut,as.character(diets$reaction)),"unclass"] = "essential_diet"
diet[intersect(allex,as.character(diets$reaction)),"Unhealthy"] = diets[intersect(allex,as.character(diets$reaction)),"fluxValue"]

write.csv(diet,file="P:/MAPPING/Pediatric_crohns/diet_unhealthy.csv")
