path = "/Volumes/T9/MRI_base"
patients = list.files(path, full.names=T)
head(patients)


one_patient = list.files(patients[100], full.names = T) 
