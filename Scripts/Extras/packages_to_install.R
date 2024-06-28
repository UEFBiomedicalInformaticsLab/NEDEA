

library(tidyverse)
library(devtools)


required_packages <- system('grep -r "library(" Scripts/ | cut -d: -f2 | sort | uniq | grep -v "#" ', intern = TRUE)
required_packages <- gsub("library\\((.*)\\)", "\\1", required_packages)
required_packages <- installed.packages() %>% rownames() %>% package_info() %>% filter(package %in% required_packages)


tmp1 <- required_packages[, c("package", "ondiskversion", "source")]
colnames(tmp1) <- c("Package", "Version", "Source")
write_csv(tmp1, file = "pkg_list.csv")

# List CRAN packages 

cran_packages <- required_packages %>% filter(source == "CRAN (R 4.2.2)")
cran_packages <- paste0('    R -e \'remotes::install_version(package = "', cran_packages$package, '", version = "', cran_packages$ondiskversion, '", dependencies = TRUE, repos = c("https://cloud.r-project.org/", "http://rforge.net/"), upgrade = "never")\'')
write_lines(cran_packages, "cran_pkg.txt")


# List BioConductor packages 
bioc_packages <- required_packages %>% filter(source == "Bioconductor")

