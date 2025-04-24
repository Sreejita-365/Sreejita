install.packages("meta")
install.packages("metafor")
install.packages("Rdpack")
library(meta)
library(metafor)

#Read the study data from the csv file
study_data=read.csv("C:/Users/SREEJITA/OneDrive/Documents/Hand Hygiene +HAI/Meta Analysis Data.csv")

#Perform meta analysis using meta package

Author_name=study_data$Author.name

meta_analysis<-metagen(data=study_data,TE=RR,seTE = SE,studlab=Author_name,common = FALSE,random=TRUE,method.tau="REML")
print(meta_analysis)

#Create a forest plot
forest(meta_analysis,main="Meta-Analysis of Relative Risks of healthcare acquired infections ",xlab="Relative Risk (RR)",slab=study_data$Author_name)

baseline_p=44.16 #prevalence of exposure during the baseline period viz before using the intervention
RR=study_data$RR;RR

#Calculate PAF for each study  using Levin's formula
calculate_paf=function(RR,baseline_p){
PAF=(baseline_p*(1-1/RR))/(baseline_p*(1-1/RR)+1)
return(PAF)
}
##A reduced RR indicates a protective effect,so the PAF will represent the proportion of cases that could be prevented by the intervention
#Calculate PAF for each study
study_data$PAF<-calculate_paf(study_data$RR,baseline_p)
print(study_data$PAF)


#Print PAF for each study
print(study_data[,c("Study ","Author_name","study_data$PAF")])


#Extract the pooled RR from the meta analysis results
pooled_RR<-meta_analysis$TE.random;pooled_RR

#Calculate combined PAF
combined_PAF<-(baseline_p*(1-1/pooled_RR))/(baseline_p*(1-1/pooled_RR)+1);combined_PAF

#Print combined PAF
print(paste("Combined Population Attributable Fraction(PAF):",round(combined_PAF,2)))

intervention_p=60.93 #prevalence of exposure during the intervention period viz after using the intervention
RR=study_data$RR;RR

#Calculate PAF for each study  using Levin's formula
calculate_paf=function(RR,intervention_p){
  PAF=(intervention_p*(1-1/RR))/(intervention_p*(1-1/RR)+1)
  return(PAF)
}
##A reduced RR indicates a protective effect,so the PAF will represent the proportion of cases that could be prevented by the intervention
#Calculate PAF for each study
study_data$PAF<-calculate_paf(study_data$RR,intervention_p)
print(study_data$PAF)


#Print PAF for each study
print(study_data[,c("Study ","Author_name","study_data$PAF")])


#Extract the pooled RR from the meta analysis results
pooled_RR<-meta_analysis$TE.random;pooled_RR

#Calculate combined PAF
combined_PAF<-(intervention_p*(1-1/pooled_RR))/(intervention_p*(1-1/pooled_RR)+1);combined_PAF

#Print combined PAF
print(paste("Combined Population Attributable Fraction(PAF):",round(combined_PAF,2)))

