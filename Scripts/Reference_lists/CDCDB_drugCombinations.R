set.seed(5081)
rm(list = ls())



# Drug combinations from CDCDB (Continuous-Drug Combination DataBase)



# Download drug combinations from CDCDB
if(!dir.exists("Databases/CDCDB/"))dir.create("Databases/CDCDB/")
# if(!file.exists("Databases/CDCDB/")){
#   download.file(url = "",
#                 destfile = "Databases/CDCDB/", method = "wget")
# }
# Download manually from: https://icc.ise.bgu.ac.il/medical_ai/CDCDB/

print(warnings())