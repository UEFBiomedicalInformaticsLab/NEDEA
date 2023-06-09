# Additional functions used in data manipulation



#Load libraries
library(tidyverse)





# Function to convert list object to dataframe 
library(tidyr)
func_list_2_df <- function(query_list){
  if(is.list(query_list)){
    df <- as.data.frame(matrix(nrow = length(names(query_list)), ncol = 2))
    for(i in 1:length(names(query_list))){
      df[i,1] <- names(query_list[i])
      df[i,2] <- paste(query_list[[i]], collapse = ";")
    }
    df <- as.data.frame(separate_rows(df, V2, sep = ";"))
    return(df)
  } else {
    print("Error: Not a list")
  }
}





# Function cbind.fill
cbind.fill <- function(...){
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
                          
                          
                          
                          
                          
# Function to extract unique elements from a string
extract_unique <- function(string) {
  elements <- unlist(strsplit(string, ","))
  unique_elements <- trimws(elements)
  unique_elements <- unique(unique_elements)
  return(paste(unique_elements, collapse = ", "))
}