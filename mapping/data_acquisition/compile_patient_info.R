require(XML)
data <- xmlParse("P:/MAPPING/Pediatric_crohns/SRA252126.sample.xml")

xml_data <- xmlToList(data)
patdata = data.frame()
for(i in 1:length(xml_data)){
  patdata[i,"SRSID"] = xml_data[[i]]$IDENTIFIERS$PRIMARY_ID
  patdata[i,"PatID"] = xml_data[[i]]$IDENTIFIERS$SUBMITTER_ID$text
  for(j in 1:length(xml_data[[i]]$SAMPLE_ATTRIBUTES)){
    patdata[i,xml_data[[i]]$SAMPLE_ATTRIBUTES[[j]]$TAG] = xml_data[[i]]$SAMPLE_ATTRIBUTES[[j]]$VALUE
  }
}

plot(as.numeric(patdata$FCP),as.numeric(patdata$Distance))
plot(as.numeric(patdata$FCP),as.numeric(patdata$BristolScore))
plot(as.numeric(patdata$FCP),as.numeric(patdata$Human.Per))
plot(as.numeric(patdata$FCP),as.numeric(patdata$PCDAI))
plot(as.numeric(patdata$FCP),as.numeric(patdata$PUCAI))
plot(as.numeric(patdata$FCP),as.numeric(patdata$log.FCP))

plot(as.numeric(patdata$PCDAI),as.numeric(patdata$PUCAI),xlim=)

plot(as.numeric(patdata$Distance),as.numeric(patdata$PUCAI))
plot(as.numeric(patdata$Distance),as.numeric(as.factor(patdata$Cluster)))

write.csv(patdata,"P:/MAPPING/Pediatric_crohns/Additional_patient_info.csv")
