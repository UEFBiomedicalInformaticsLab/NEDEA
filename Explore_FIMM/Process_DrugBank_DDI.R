



find_top_k_frequent_sequences <- function(strings_list, n = 2, k = 10) {
  require(stringi)  # For string operations
  # Generate n-grams
  generate_ngrams <- function(text, n) {
    words <- unlist(strsplit(text, split = " "))
    l <- length(words)
    ngrams <- sapply(1:(l-n+1), function(i) paste(words[i:(i+n-1)], collapse = " "))
    return(ngrams)
  }
  all_ngrams <- unlist(lapply(strings_list, function(x) generate_ngrams(x, n)))
  # Count and sort n-grams
  ngram_freq <- table(all_ngrams)
  sorted_ngram_freq <- sort(ngram_freq, decreasing = TRUE)
  # Return the top k frequent n-grams and their counts
  top_k <- min(k, length(sorted_ngram_freq))
  list(ngrams = names(sorted_ngram_freq)[1:top_k], counts = as.integer(sorted_ngram_freq[1:top_k]))
}


find_top_k_sequences_with_keyword <- function(strings_list, keyword, n = 2, k = 10) {
  require(stringi)
  # Generate n-grams
  generate_ngrams <- function(text, n) {
    words <- unlist(strsplit(text, split = " "))
    l <- length(words)
    ngrams <- sapply(1:(l-n+1), function(i) paste(words[i:(i+n-1)], collapse = " "))
    return(ngrams)
  }
  # Filter sequences with the keyword
  filter_keyword <- function(ngram_list, keyword) {
    ngram_list[sapply(ngram_list, function(ngram) grepl(keyword, ngram))]
  }
  all_ngrams <- unlist(lapply(strings_list, function(x) generate_ngrams(x, n)))
  keyword_ngrams <- filter_keyword(all_ngrams, keyword)
  # If no sequences with keyword, return empty list
  if (length(keyword_ngrams) == 0) {
    return(list(ngrams = character(0), counts = integer(0)))
  }
  # Count and sort n-grams
  ngram_freq <- table(keyword_ngrams)
  sorted_ngram_freq <- sort(ngram_freq, decreasing = TRUE)
  # Return the top k frequent n-grams and their counts
  top_k <- min(k, length(sorted_ngram_freq))
  list(ngrams = names(sorted_ngram_freq)[1:top_k], counts = as.integer(sorted_ngram_freq[1:top_k]))
}


# 
# DrugBank_drugInteractions <- read.csv("Databases/DrugBank/drug_drug_interactions.csv")
# 
# 
# 
# top_freq = find_top_k_frequent_sequences(DrugBank_drugInteractions$description, n = 4, 100)
# 
# 
# risk_sev_res = find_top_k_sequences_with_keyword(DrugBank_drugInteractions$description, keyword = "The risk or severity", n = 6, k = 50)
# 
# 
# test_eff_inc = grepl("The therapeutic efficacy of", DrugBank_drugInteractions$description) & grepl("can be increased when used in combination", DrugBank_drugInteractions$description)
# test_eff_dec = grepl("The therapeutic efficacy of", DrugBank_drugInteractions$description) & grepl("can be decreased when used in combination", DrugBank_drugInteractions$description)
# test_anticinc = grepl("may increase the anticoagulant activities", DrugBank_drugInteractions$description)
# test_anticdec = grepl("may decrease the anticoagulant activities", DrugBank_drugInteractions$description)
# test_antihinc = grepl("may increase the antihypertensive activities", DrugBank_drugInteractions$description)
# test_antihdec = grepl("may decrease the antihypertensive activities", DrugBank_drugInteractions$description)
# test_nsd = grepl("(CNS depressant)", DrugBank_drugInteractions$description)
# test_scinc = grepl("The serum concentration of", DrugBank_drugInteractions$description) & grepl("can be increased", DrugBank_drugInteractions$description)
# test_scdec = grepl("The serum concentration of", DrugBank_drugInteractions$description) & grepl("can be decreased", DrugBank_drugInteractions$description)
# test_meinc = grepl("The metabolism of", DrugBank_drugInteractions$description) & grepl("can be increased", DrugBank_drugInteractions$description)
# test_medec = grepl("The metabolism of", DrugBank_drugInteractions$description) & grepl("can be decreased", DrugBank_drugInteractions$description)
# test_absde = grepl("a decrease in the absorption ", DrugBank_drugInteractions$description)
# #
# list_tox_test = list()
# for(i in 1:length(risk_sev_res$ngrams)){
#   list_tox_test[[length(list_tox_test) + 1]] = grepl(risk_sev_res$ngrams[i], DrugBank_drugInteractions$description)
# }
# names(list_tox_test) = risk_sev_res$ngrams
# #
# sum(test_eff_inc)
# sum(test_eff_dec)
# sum(test_anticinc)
# sum(test_anticdec)
# sum(test_antihinc)
# sum(test_antihdec)
# sum(test_nsd)
# sum(test_scinc)
# sum(test_scdec)
# sum(test_meinc)
# sum(test_medec)
# sum(test_absde)
# 
# tox_test_df <- data.frame(list_tox_test, check.names = FALSE)
# 
# 
# DrugBank_drugInteractions <- cbind(DrugBank_drugInteractions, tox_test_df)