###################################
# Read patient biotab and make surv
###################################

# set current directory and clinical patient biotab
setwd("~/Downloads/")
clin.patient.link <- "nationwidechildrens.org_clinical_patient_brca.txt"

# read patient biotab
clin <- read.delim(clin.patient.link, h=T, stringsAsFactors=F, sep="\t")
skip <- which(clin$bcr_patient_uuid %in% c("bcr_patient_uuid", "CDE_ID:"))
clin <- clin[-skip, ]

# make data.frame surv
death.days <- suppressWarnings(as.numeric(clin$death_days_to))
contact.days <- suppressWarnings(as.numeric(clin$last_contact_days_to))

patients <- clin$bcr_patient_barcode
surv <- data.frame(patients)	
surv$contact.days = contact.days
surv$death.days = death.days
	
if(sum(contact.days > death.days, na.rm=T) != 0) 
	stop("Last contact after death?")
	
death = rep(F, length(patients))
death[which(!is.na(death.days))] = T
days <- contact.days
days[death] = death.days[death]
surv$days = days
surv$death = death
surv$months = round(days/365.2425*12, 2)

# select patient with defined er, pr and her2 status
surv$er = clin$er_status_by_ihc
surv$pr = clin$pr_status_by_ihc
surv$her2 = clin$her2_status_by_ihc
opt = c("Negative", "Positive")
select = with(surv, which(er %in% opt & pr %in% opt & her2 %in% opt))
surv = surv[select, ]

# remove patient with neither last contact days or death days
surv <- surv[!is.na(surv$days), ]

# save surv as surv.txt	
surv$label = "Other"
surv$label[with(surv, which(er == "Negative" & pr == "Negative" & her2 == "Negative"))] = "Triple Negative"
write.table(surv, "surv.txt", col.names=T, row.names=F, sep="\t", quote=F)	

###################################
# read surv and make survival plot by label
###################################
library(survival)
library(ggplot2)
# surv <- read.table("surv.txt", h=T, sep="\t")

surv2 = surv[which(surv$days != 0), ]
surv2$label <- as.factor(surv2$label)

surv.model <- with(surv2, Surv(months,death) ~ label)
ggsurv(survfit(surv.model), xlab = "Time (month)", ylab = "Survival Portion")
survdiff(surv.model)

