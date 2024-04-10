run_2PL_model <- function(language) {
  set.seed(42)
  
  lang_data <- readRDS(glue("data/all_forms/{language}_data.rds"))
  
  d_prod <- lang_data$all_prod
  
  bad_words_prod <- c(which(colSums(d_prod, na.rm=T) == 0), 
                      which(colSums(d_prod, na.rm=T) == colSums(!is.na(d_prod))))
  if(length(bad_words_prod) != 0) d_prod = d_prod[,-bad_words_prod]
  print(glue("{length(bad_words_prod)} words with all zero/max responses removed from {language}"))
  
  bad_ppts_prod <- c(which(rowSums(d_prod, na.rm=T) == 0))
  if(length(bad_ppts_prod) != 0) d_prod = d_prod[-bad_ppts_prod,]
  print(glue("{length(bad_ppts_prod)} ppts with all zero responses removed from {language}"))
  
  print(glue("Fitting {nrow(d_prod)} subjects and {ncol(d_prod)} words in {language}"))
  mod_string <- glue("G = 1-{ncol(d_prod)},\n",
                     "LBOUND = (1-{ncol(d_prod)}, a1, 0),\n",
                     "PRIOR = (1-{ncol(d_prod)}, d, norm, 0, 3)")
  mod <- mirt.model(mod_string)
  model <- mirt(data = d_prod, model = mod, itemtype = "2PL", 
                method = "QMCEM", verbose = TRUE, 
                technical = list(NCYCLES = 3000))
  coefs <- as_tibble(coef(model, simplify = TRUE)$items,
                     rownames = "definition")
  
  saveRDS(list(model = model,
               coefs = coefs),
          glue("data/prod_models/{language}_2PL_allforms_prod_fits.rds"))
}