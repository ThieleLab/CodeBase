
#
# DeNoPa Sample Matching and Selection Analysis
#
#
# Important: The folder location of the DeNoPa input datasets has to be specified below.
# Please replace "/Users/enrico.glaab/denopa/" with your working folder containing the data.
#
#

folder_containing_all_data = '/Users/enrico.glaab/denopa/' # need to add the '/'-symbol at the end!

# location for DeNoPa clinical data
denopa_data_location = paste(folder_containing_all_data,'DeNoPa.v3.sav',sep="")

# location for DeNoPa genetics data
genetics_data_location = paste(folder_containing_all_data,'genetikdat.Rdata',sep="")

# location for DeNoPa body height data
denopa_bodyheight_data_location = paste(folder_containing_all_data,'bodyheight.Rdata',sep="")

# location for the filtering list for samples with lacking biospecimen volumes (provided by Brit Mollenhauer)
denopa_filtering_list_location = paste(folder_containing_all_data,'filtering_list_samples_with_lacking_biospecimen_volumes_from_brit_mollenhauer.txt',sep="")



#
# Step 1: Loading and filtering the DeNoPa clinical data
#


# load R-package 'foreign' to parse SPSS data
if(!require('foreign'))
{
  print("The R-package 'foreign' has to be installed in order to use this script (Internet connection required for automated installation)")
  install.packages('foreign')
}


# Load DeNoPa clinical data
if(file.exists(denopa_data_location)){
	options(warn=-1)
	dataset = read.spss(denopa_data_location, to.data.frame=TRUE, reencode=TRUE)
} else {
	stop("You need to specify the correct folder name of the DeNoPa.v3.sav clinical dataset in order to use this script.")
}

# remove all white-spaces
datanows = apply(as.matrix(dataset), 2, function(x) gsub(" ", "", x, fixed = TRUE, useBytes = TRUE))

# replace empty values by NA
datanows[which(datanows=="")] = rep(NA, length(which(datanows=="")))

# retrieve columns containing only missing/empty values and remove them
nacols = which(apply(datanows, 2, function(x) all(is.na(x))))
datfilt = datanows[,-nacols]


# extract data for the three main sample groups: PD subjects, healthy controls and other neurodegenerative disorders (OND)
datfilt2 = datfilt[c(which(as.matrix(datfilt[,5]) == "PDsubject"),which(as.matrix(datfilt[,5]) == "HC-healthycontrol"),which(as.matrix(datfilt[,5]) == "ONDsubject")),]

# replace gender labels by numerical values (0 to encode female, 1 to encode male)
genderorg = ifelse(datfilt2[,3]=="weiblich",0,1)
datfilt2[,3] = genderorg

# determine the parts of the data that can be converted to numeric format and the non-numeric/categorical part
datfiltnum = apply(datfilt2, 2, as.numeric)

# determine the non-numerical part of the data
datnotnum = which(apply(datfiltnum, 2, function(x) all(is.na(x))))


# determine the subjects with too high numbers of missing values and remove them
# number of missing values per sample:
nummiss = apply(datfiltnum, 1, function(x) length(which(is.na(x))))

# remove the subjects that have too high numbers of missing values (cut-off 500 used)
filtsamp = which(nummiss-length(datnotnum) > 500)
datfiltnum2 = datfiltnum[-filtsamp,]

# determine the non-numerical part of the data after filtering
datnotnum_filt = which(apply(datfiltnum2, 2, function(x) any(is.na(x))))


#
# Categorical features: retain only features with less than 10 categories and convert to integer-encoded features
#
datcat = datfilt2[-filtsamp,datnotnum_filt]

# determine categorical features without missing values
catnona = which(apply(datcat, 2, function(x) length(which(is.na(x)))==0))
datcatn = datcat[,catnona]

# remove variables with more than 10 categories
catvar = which(apply(datcatn, 2, function(x) (length(unique(x)) <= 10)))
datcatfinal = datcatn[,catvar]
# remove the "Gruppe"/c_group feature and convert the categories to integers
datcatnum = apply(datcatfinal[,-2], 2, function(x) match(x, unique(x)))


#
# Numerical features: take only those without any missing values
#
datnotna_ind = which(apply(datfiltnum2, 2, function(x) all(!is.na(x))))
datfinalnum = datfiltnum2[,datnotna_ind]

# combine converted categorical features with the filtered numerical features
datfinal = cbind(datfinalnum, datcatnum)

# label the rows according to the original subject identifiers
rownames(datfinal) = datfilt2[-filtsamp,1]

# determine the clinical outcome labels and encode them numerically
yclin = rep(1,nrow(datfinal))
yclin[which(as.matrix(datfilt2[-filtsamp,5]) == "ONDsubject")] = rep(2,length(which(as.matrix(datfilt2[-filtsamp,5]) == "ONDsubject")))
yclin[which(as.matrix(datfilt2[-filtsamp,5]) == "HC-healthycontrol")] = rep(0,length(which(as.matrix(datfilt2[-filtsamp,5]) == "HC-healthycontrol")))




#
# Step 2: Extract data on relevant confounding variables: body height, body weight, age, gender etc. 
#

# Load body height data
if(file.exists(denopa_bodyheight_data_location)){
	load(file=denopa_bodyheight_data_location)
} else {
	stop("Please specify the correct folder location for the DeNoPa body height data")
}

# Load body weight data
gewicht = datfilt[,27]
names(gewicht) = datfilt[,1]

# make sample identifier prefix identical for mappings
names(bodyheight) = gsub("P0","P",gsub("P0","P",gsub("K0","K",gsub("K0","K",names(bodyheight)))))
names(gewicht) = gsub("P0","P",gsub("P0","P",gsub("K0","K",gsub("K0","K",names(gewicht)))))

# extract age and gender data
ages = datfinal[,3]
genders = datfinal[,2]

# make sample identifier prefix identical for mappings
rownames(datfinal) = gsub("P0","P",gsub("P0","P",gsub("K0","K",gsub("K0","K",rownames(datfinal)))))

# extract relevant weight and height data via sample mappings
indheight = match(rownames(datfinal), names(bodyheight))
indgewicht = match(rownames(datfinal), names(gewicht))

selheights = bodyheight[indheight]
selgewicht = as.numeric(gewicht[indgewicht])

# impute missing weight values by median value
selgewicht[which(is.na(selgewicht))] = median(selgewicht, na.rm=T)

# compute BMIs
bmis = selgewicht / (selheights/100)^2

# filter out all dead subjects
dead = which(!is.na(datfilt2[,7]))
dead_ids = as.matrix(datfilt2[dead,1])
filter_dead = match(gsub("0","",dead_ids), rownames(datfinal))

# remove the dead subjects from each of the relevant variables
yclin = yclin[-filter_dead]
genders = genders[-filter_dead]
ages = ages[-filter_dead]
bmis = bmis[-filter_dead]

# extract information on blood pressure and other symptoms (+filter out dead subjects)
# Encode missing vlue NA numerically as 0.5
bloodpressure = ifelse(is.na(datfilt2[,match("c_Scopa_AUT_26c",colnames(datfilt2))]),0.5,ifelse(datfilt2[,match("c_Scopa_AUT_26c",colnames(datfilt2))]=="ja",1,0))
othersymp = ifelse(is.na(datfilt2[,match("c_Scopa_AUT_26d",colnames(datfilt2))]),0.5,ifelse(datfilt2[,match("c_Scopa_AUT_26d",colnames(datfilt2))]=="ja",1,0))
bloodpressurefinal = (bloodpressure[-filtsamp])[-filter_dead]
othersympfinal = (othersymp[-filtsamp])[-filter_dead]


# final clinical data without entries for dead subjects
datfinal = datfinal[-filter_dead,]

# final list of subjects
subjectsfinal = rownames(datfinal)



#
# Step 3: Filter out samples with lacking biospecimen volumes according to the following list (provided by Brit Mollenhauer)
#

if(file.exists(denopa_filtering_list_location)){
	filterdat = as.matrix(read.table(denopa_filtering_list_location, header=F))
} else {
	stop("Please specify the correct folder location for the DeNoPa sample filtering list.")
}

# make sample identifier prefix identical for mappings
filterdat2 = gsub("P0","P",gsub("P0","P",gsub("K0","K",gsub("K0","K",filterdat))))

# map the subjects to be filterd onto the subject list
mapfilter = match(filterdat2, subjectsfinal)

# final filtering based on information for lacking biospecimen volumes
mapfilternona = mapfilter[which(!is.na(mapfilter))]
allsubfilter2 = subjectsfinal[mapfilternona]
allsubfilterind2 = match(allsubfilter2, subjectsfinal)



#
# Step 4: Load DeNoPa genetics data to use it for the sample selection
#
if(file.exists(genetics_data_location)){
	load(file=genetics_data_location)
} else {
	stop("Please specify the correct location for the genetics data folger.")
}

# make sample IDs identical and compute mappings
mapgenetind = match(subjectsfinal, gsub("P0","P",gsub("P0","P",gsub("K0","K",gsub("K0","K",rownames(genetikdat))))))

# extract subjects with APOE E4 haplotypes E4/E4 and E3/E4
apoe_e4sub = subjectsfinal[c(which(genetikdat[mapgenetind,10]=="E4/E4"),which(genetikdat[mapgenetind,10]=="E3/E4"))]

# add further genetics filters
rmgenetikind2 = unique(c(which(genetikdat[mapgenetind,1]=="Mutationnachweisbar"),which(genetikdat[mapgenetind,5]=="261/263"),which(genetikdat[mapgenetind,5]=="263/263"),which(!is.na(genetikdat[mapgenetind,6]))))

# combine the previous filters with the added genetics filters to obtain the 2nd stage filterings
allfiltered2 = union(allsubfilterind2, rmgenetikind2)


# determine relevant clinical covariates after 2nd stage genetics filtering
yclinallfilt2 = yclin[-allfiltered2]
gendersallfilt2 = genders[-allfiltered2]
agesallfilt2 = ages[-allfiltered2]
bmisallfilt2 = bmis[-allfiltered2]
bloodpressureallfilt2 = bloodpressurefinal[-allfiltered2]
othersympallfilt2 = othersympfinal[-allfiltered2]
subjectsallfilt2 = subjectsfinal[-allfiltered2]


# specify filtering variables for APOE E4 status and occurrence of SNCA-related risk factor SNP (rs11931074)
apoee4filt = rep(0,length(subjectsallfilt2))
sncasnpfilt = rep(0,length(subjectsallfilt2))
rs11931074_SNP = subjectsfinal[which(genetikdat[mapgenetind,3]!="G/G")]
apoee4filt[which(!is.na(match(subjectsallfilt2,apoe_e4sub)))] = rep(1, length(which(!is.na(match(subjectsallfilt2,apoe_e4sub)))))
sncasnpfilt[which(!is.na(match(subjectsallfilt2,rs11931074_SNP)))] = rep(1, length(which(!is.na(match(subjectsallfilt2,rs11931074_SNP)))))




#
# Step 5: Run function to determine the matched sample selection
#

# Function to compute the optimized sample selection
optimize_sample_selection = function(outcome, genders, ages, bmis, othersymp, bloodpressure, apoee4, sncasnp, numsel = 10, gendweight = 1, balanceweight = 1, ageweight = 1, bmiweight = 1)
{
	
	print("Starting optimization...")
	# start with initial patient selection (patsel) and control selection (contsel) containing all patients and controls from the prior filterings
	patsel = which(outcome == 1)
	contsel = which(outcome == 0)
	
	
	# score the current selection based on how balanced the genders are (balance_score), how well matched the genders are across groups (gendiff_score),
	# how similar the age and BMI distributions are (agediff_score and bmidiff_score), how many other disease symptoms occur (othersymp_score, to be minimized),
	# how similar the groups are in terms of blood pressure problems (bloodpressure_score), and how often APOE E4 and SNCA SNP variations occur (apoee4_score, sncasnp_score, to be minimized)
	
 	balance_score = as.numeric(abs(table(genders[patsel])[1] - table(genders[patsel])[2])+abs(table(genders[contsel])[1] - table(genders[contsel])[2]))/((length(patsel)+length(contsel)))
 	gendiff_score = abs(sum(genders[patsel]) - sum(genders[contsel]))/((length(patsel)+length(contsel))/2)
 	agediff_score = as.numeric((ks.test(ages[patsel],ages[contsel])$statistic))
 	bmidiff_score = as.numeric((ks.test(bmis[patsel],bmis[contsel])$statistic))
 	othersymp_score = (sum(othersymp[patsel]) + sum(othersymp[contsel]))/((length(patsel)+length(contsel)))
 	bloodpressure_score = abs(sum(bloodpressure[patsel]) - sum(bloodpressure[contsel]))/((length(patsel)+length(contsel))/2)
 	apoee4_score = (sum(apoee4[patsel]) + sum(apoee4[contsel]))/((length(patsel)+length(contsel)))
 	sncasnp_score = (sum(sncasnp[patsel]) + sum(sncasnp[contsel]))/((length(patsel)+length(contsel)))
 	
 	# compute the sum of all individual scores
 	sumscore = balanceweight*balance_score + gendweight*gendiff_score + ageweight*agediff_score + bmiweight*bmidiff_score + othersymp_score + bloodpressure_score + apoee4_score + sncasnp_score
	
	# initialize the best score with a high value and iteratively remove patients or controls from the selection to improve the score (as long as the specified selection sizes for patients and controls have not been reached) 
	bestscore = 100000
	while((length(patsel) > numsel) || length(contsel) > numsel)
	{
	   patbest = 0
	   curscore = 100
	   bestsel = -1
	   curbest = 100
	   cursel = -1
	   patcur = 0
	   
	   if(length(patsel) > numsel)
	   {
		   # remove next patient and check whether the score for the selection improves
		   for(i in 1:length(patsel))
		   {
		     balance_score = as.numeric(abs(table(genders[patsel[-i]])[1] - table(genders[patsel[-i]])[2])+abs(table(genders[contsel])[1] - table(genders[contsel])[2]))/((length(patsel[-i])+length(contsel)))   
		     	   
			   gendiff_score = abs(sum(genders[patsel[-i]]) - sum(genders[contsel]))/((length(patsel[-i])+length(contsel))/2)
			   
			   agediff_score = as.numeric((ks.test(ages[patsel[-i]],ages[contsel])$statistic))
			   
			   bmidiff_score = as.numeric((ks.test(bmis[patsel[-i]],bmis[contsel])$statistic))
			
				 othersymp_score = (sum(othersymp[patsel[-i]]) + sum(othersymp[contsel]))/((length(patsel[-i])+length(contsel)))
		   		
		   	 bloodpressure_score = abs(sum(bloodpressure[patsel[-i]]) - sum(bloodpressure[contsel]))/((length(patsel[-i])+length(contsel))/2)
			
				 apoee4_score = (sum(apoee4[patsel[-i]]) + sum(apoee4[contsel]))/((length(patsel[-i])+length(contsel)))
		   	 
		   	 sncasnp_score = (sum(sncasnp[patsel[-i]]) + sum(sncasnp[contsel]))/((length(patsel[-i])+length(contsel)))
			
				 curscore = balanceweight*balance_score + gendweight*gendiff_score + ageweight*agediff_score + bmiweight*bmidiff_score + othersymp_score + bloodpressure_score + apoee4_score + sncasnp_score
				 
				 if(curscore <= bestscore)
				 {
				   bestscore = curscore
				   curbest = bestscore
				   bestsel = i
				   patbest = 1
				 } else if(curscore <= curbest)
				 {
				   curbest = curscore
				   cursel = i
				   patcur = 1
				 }
		   }
	   }
	   if(length(contsel) > numsel)
	   {
	   
		   # remove next control and check whether the score for the selection improves
		   for(i in 1:length(contsel))
		   {
			   balance_score = as.numeric(abs(table(genders[patsel])[1] - table(genders[patsel])[2])+abs(table(genders[contsel[-i]])[1] - table(genders[contsel[-i]])[2]))/((length(patsel)+length(contsel[-i])))
			   
			   gendiff_score = abs(sum(genders[patsel]) - sum(genders[contsel[-i]]))/((length(patsel)+length(contsel[-i]))/2)
			   
			   agediff_score = as.numeric((ks.test(ages[patsel],ages[contsel[-i]])$statistic))
			   
			   bmidiff_score = as.numeric((ks.test(bmis[patsel],bmis[contsel[-i]])$statistic))
			   
			   othersymp_score = (sum(othersymp[patsel]) + sum(othersymp[contsel[-i]]))/((length(patsel)+length(contsel[-i])))
		   		
		   	 bloodpressure_score = abs(sum(bloodpressure[patsel]) - sum(bloodpressure[contsel[-i]]))/((length(patsel)+length(contsel[-i]))/2)
		   	 
		   	 apoee4_score = (sum(apoee4[patsel]) + sum(apoee4[contsel[-i]]))/((length(patsel)+length(contsel[-i])))
		   	 
		   	 sncasnp_score = (sum(sncasnp[patsel]) + sum(sncasnp[contsel[-i]]))/((length(patsel)+length(contsel[-i])))
			
				 curscore = balanceweight*balance_score + gendweight*gendiff_score + ageweight*agediff_score + bmiweight*bmidiff_score + othersymp_score + bloodpressure_score + apoee4_score + sncasnp_score
		   
		   	 if(curscore <= bestscore)
				 {
				   bestscore = curscore
				   curbest = bestscore
				   bestsel = i
				   patbest = 0
				 } else if(curscore <= curbest)
				 {
				   curbest = curscore
				   cursel = i
				   patcur = 0
				 }
		   }
	   }
	   
	   # Check if the global score has improved
	   if(bestsel == -1)
	   {
	     cat(paste('\nNo global improvement - score:',bestscore,' - ',curscore,' - ',length(contsel),',',length(patsel)))
	     
	     # remove current selection
	     if(patcur == 0)
	     {
	       contsel = contsel[-cursel]
	     } else {
	       patsel = patsel[-cursel]
	     }	     	     
	     
	   } else {
	   
	   	 # Show current improvement and remove the selected subject
	   	 cat(paste('\nImprovement - score: ',bestscore,' - no. controls: ',length(contsel),', no. patients: ',length(patsel)))   	
	   
	     if(patbest == 0)
	     {
	       contsel = contsel[-bestsel]
	     } else {
	       patsel = patsel[-bestsel]
	     }
	          
	   }
	}
	
	
	# after initial optimization check all single-sample replacements to attain possible further improvements
	bestsel = 1
	while(bestsel != -1)
	{
		patnonsel = setdiff(which(outcome == 1),patsel)
		contnonsel = setdiff(which(outcome == 0),contsel)
		bestsel = -1
		worstsel = -1
		patbest = 0
		curscore = 100
		curbest = 100
		cursel = -1
		patcur = 0
				
		for(i in 1:length(patsel))
		{
		      
		      # check and score replacements of selected with non-selected patient
		      for(j in 1:length(patnonsel))
		      {
			     balance_score = as.numeric(abs(table(genders[c(patsel[-i],patnonsel[j])])[1] - table(genders[c(patsel,patnonsel[j])])[2])+abs(table(genders[contsel])[1] - table(genders[contsel])[2]))/((length(c(patsel[-i],patnonsel[j]))+length(contsel)))   
			     	   
				   gendiff_score = abs(sum(genders[c(patsel[-i],patnonsel[j])]) - sum(genders[contsel]))/((length(patsel)+length(contsel))/2)
				   
				   agediff_score = as.numeric((ks.test(ages[c(patsel[-i],patnonsel[j])],ages[contsel])$statistic))
				   
				   bmidiff_score = as.numeric((ks.test(bmis[c(patsel[-i],patnonsel[j])],bmis[contsel])$statistic))
				
					 othersymp_score = (sum(othersymp[c(patsel[-i],patnonsel[j])]) + sum(othersymp[contsel]))/((length(patsel)+length(contsel)))
			   		
			   	 bloodpressure_score = abs(sum(bloodpressure[c(patsel[-i],patnonsel[j])]) - sum(bloodpressure[contsel]))/((length(patsel)+length(contsel))/2)
			   	 
			   	 apoee4_score = (sum(apoee4[c(patsel[-i],patnonsel[j])]) + sum(apoee4[contsel]))/((length(patsel)+length(contsel)))
		   	 
		   	 	 sncasnp_score = (sum(sncasnp[c(patsel[-i],patnonsel[j])]) + sum(sncasnp[contsel]))/((length(patsel)+length(contsel)))
				
					 curscore =  balanceweight*balance_score + gendweight*gendiff_score + ageweight*agediff_score + bmiweight*bmidiff_score + othersymp_score + bloodpressure_score + apoee4_score + sncasnp_score
					 
					 if(curscore < bestscore)
					 {
					   bestscore = curscore
					   curbest = bestscore
					   bestsel = j
					   worstsel = i
					   patbest = 1
					 }
				 }
		}
		
		
		for(i in 1:length(contsel))
		{
		      # check and score replacements of selected with non-selected controls
		      for(j in 1:length(contnonsel))
		      {
			     balance_score = as.numeric(abs(table(genders[patsel])[1] - table(genders[c(patsel,patnonsel[j])])[2])+abs(table(genders[c(contsel[-i],contnonsel[j])])[1] - table(genders[c(contsel[-i],contnonsel[j])])[2]))/((length(patsel)+length(contsel)))   
			     	   
				   gendiff_score = abs(sum(genders[patsel]) - sum(genders[c(contsel[-i],contnonsel[j])]))/((length(patsel)+length(contsel))/2)
				   
				   agediff_score = as.numeric((ks.test(ages[patsel],ages[c(contsel[-i],contnonsel[j])])$statistic))
				   
				   bmidiff_score = as.numeric((ks.test(bmis[patsel],bmis[c(contsel[-i],contnonsel[j])])$statistic))
				
					 othersymp_score = (sum(othersymp[patsel]) + sum(othersymp[c(contsel[-i],contnonsel[j])]))/((length(patsel)+length(contsel)))
			   		
			   	 bloodpressure_score = abs(sum(bloodpressure[patsel]) - sum(bloodpressure[c(contsel[-i],contnonsel[j])]))/((length(patsel)+length(contsel))/2)
			   	 
			   	 apoee4_score = (sum(apoee4[patsel]) + sum(apoee4[c(contsel[-i],contnonsel[j])]))/((length(patsel)+length(contsel)))
		   	 
		   	 	 sncasnp_score = (sum(sncasnp[patsel]) + sum(sncasnp[c(contsel[-i],contnonsel[j])]))/((length(patsel)+length(contsel)))
				
					 curscore =  balanceweight*balance_score + gendweight*gendiff_score + ageweight*agediff_score + bmiweight*bmidiff_score + othersymp_score + bloodpressure_score + apoee4_score + sncasnp_score
					 
					 if(curscore < bestscore)
					 {
					   bestscore = curscore
					   curbest = bestscore
					   bestsel = j
					   worstsel = i
					   patbest = 0
					 }
				 }
		}
		
		# Show if improvments have been attained from the final optimization step and update the selection accordingly
		if(bestsel != -1)
		{
		  	 cat(paste('\nFinal improvement - score: ',bestscore,' - no. controls: ',length(contsel),', no. patients: ',length(patsel)))   	
		   
		     if(patbest == 0)
		     {
		       contsel = c(contsel[-worstsel],contnonsel[bestsel])
		     } else {
		       patsel = c(patsel[-worstsel],patnonsel[bestsel])
		     }
		     		          
		}
	}
	
	# return final selection of patients (patsel) and controls (contsel)
	return(list(patsel=patsel,contsel=contsel))
}



# Run the sample selection function
opt30_selection = optimize_sample_selection(yclinallfilt2, gendersallfilt2, agesallfilt2, bmisallfilt2, othersympallfilt2, bloodpressureallfilt2, apoee4filt, sncasnpfilt, numsel = 30, gendweight = 10, balanceweight = 20, ageweight = 1, bmiweight = 1)

# extract the selected patients and control identifiers
patsel = opt30_selection$patsel
contsel = opt30_selection$contsel


print("Final selection of patients:")
print(paste(subjectsallfilt2[patsel],collapse=", "))
# DKP6, DKP14, DKP16, DKP18, DKP44, DKP50, DKP51, DKP60, DKP62, DKP63, DKP81, DKP92, DKP96, DKP97, DKP100, DKP110, DKP111, DKP112, DKP114, DKP116, DKP119, DKP128, DKP131, DKP132, DKP144, DKP145, DKP151, DKP161, DKP162, DKP123


print("Final selection of controls:")
print(paste(subjectsallfilt2[contsel],collapse=", "))
# DKK12, DKK14, DKK15, DKK17, DKK27, DKK28, DKK32, DKK34, DKK35, DKK36, DKK38, DKK40, DKK42, DKK43, DKK44, DKK54, DKK55, DKK62, DKK65, DKK73, DKK76, DKK82, DKK83, DKK88, DKK89, DKK100, DKK102, DKK105, DKK112, DKK113
