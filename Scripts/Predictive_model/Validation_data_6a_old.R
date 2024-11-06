
library(tidyverse)
library(ctrdata)

#####


# Generate the database with the predetermined quury
q <- ctrGetQueryUrl(url = "https://clinicaltrials.gov/search?cond=Breast%20Cancer&term=Breast%20Carcinoma&aggFilters=results:with,status:com%20ter,studyType:int")

db <- nodbi::src_sqlite(
  dbname = "ClinicalTrials_sql",
  collection = "BreastCancer"
)

ctrLoadQueryIntoDb(
  queryterm = q,
  euctrresults = TRUE,
  con = db
)



####


result <-  dbGetFieldsIntoDf(fields = c("protocolSection.armsInterventionsModule.armGroups", 
                                        "resultsSection.participantFlowModule.groups", 
                                        "resultsSection.adverseEventsModule.eventGroups"),
                             con = db)

colnames(result) <- c("NCT_ID", "interventions", "groups", "adverse_events")


#####


# Unnest adverse events columns

result$col_df <- unlist(lapply(result$adverse_events, is.data.frame))

tmp1 <- result %>% 
  filter(col_df %in% "TRUE") %>% 
  unnest(adverse_events, names_sep = "__") %>%
  # select(!col_df)  
  select(c("NCT_ID", "interventions", "groups", 
           "adverse_events__id", "adverse_events__title", "adverse_events__description",
           "adverse_events__deathsNumAffected", "adverse_events__deathsNumAtRisk",
           "adverse_events__seriousNumAffected", "adverse_events__seriousNumAtRisk"))

tmp2 <- result %>% 
  filter(col_df %in% "FALSE") %>% 
  mutate("adverse_events" = lapply(adverse_events, as_tibble) ) %>% 
  unnest(adverse_events, names_sep = "__") %>%
  # select(!col_df)
  select(c("NCT_ID", "interventions", "groups", 
           "adverse_events__id", "adverse_events__title", "adverse_events__description",
           "adverse_events__deathsNumAffected", "adverse_events__deathsNumAtRisk",
           "adverse_events__seriousNumAffected", "adverse_events__seriousNumAtRisk"))

result <- rbind(tmp1, tmp2)
rm(list = c("tmp1", "tmp2"))



#####


# Unnest groups columns

result$col_df <- unlist(lapply(result$groups, is.data.frame))

tmp1 <- result %>% 
  filter(col_df %in% "TRUE") %>% 
  unnest(groups, names_sep = "__") %>%
  # select(!col_df) %>% 
  select(c("NCT_ID", "interventions", starts_with("groups__"), starts_with("adverse_events__")))

tmp2 <- result %>% 
  filter(col_df %in% "FALSE") %>% 
  mutate("groups" = lapply(groups, as_tibble) )%>% 
  unnest(groups, names_sep = "__") %>%
  # select(!col_df) %>% 
  select(c("NCT_ID", "interventions", starts_with("groups__"), starts_with("adverse_events__")))

result <- rbind(tmp1, tmp2)
rm(list = c("tmp1", "tmp2"))


#####


# Unnest interbventions column

result$interventions <- lapply(result$interventions, function(x){
  
  if(is.data.frame(x)){
    
    if(! all(sapply(x$interventionNames, is.character))){
      
      x$interventionNames <- sapply(x$interventionNames, as.character)
      x <- x[which(sapply(x$interventionNames, length) != 0), ] # Removing arms without interventions to remove data precessing error
      
    }
    
    x <- unnest(x, cols = colnames(x)) %>% 
      group_by(label) %>% 
      summarise(across(everything(), ~ toString(unique(.x))))
    
  } else {
    x <- as_tibble(x)
    if(nrow(x) > 1){
      x <- x %>% 
        group_by(label) %>% 
        summarise(across(everything(), ~ toString(unique(.x))))
    }
  }
  x
})

result <- result %>%
  unnest(interventions, names_sep = "__")


#####


tmp1 <- result %>% 
  select(!c("interventions__value", "groups__id", "adverse_events__id")) %>% 
  group_by(NCT_ID, interventions__label, interventions__type, groups__title, adverse_events__title) %>% 
  summarise(across(everything(), ~ toString(unique(.x))))



