#Extracting info from xml on PD project
require(XML)
library("methods")

data <- xmlParse("Y:/Microbiome/PD/DeNoPa_scripts/ERS1647277-ERS1647335.xml")

xml_data <- xmlToList(data)


infMat=matrix(ncol=2,nrow =59)
#sAMPLE NAME
for (i in 1:59){
sampnam=xml_data[[i]]$IDENTIFIERS$PRIMARY_ID
status=as.matrix(unlist(xml_data[[i]]$SAMPLE_ATTRIBUTES))
status=status[20,1]
infMat[i,1]=as.character(sampnam)
infMat[i,2]=as.character(status)
}
write.csv(file="Y:/Microbiome/PD/DeNoPa_scripts/infMat.csv",infMat)
