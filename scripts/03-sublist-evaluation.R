COR_TYPE = "spearman" # "pearson"

make_swadesh_sublist <- function(prod_sum, list_size, k) {
  prod_sum |> 
    filter(num_langs >= k) |> 
    arrange(sd_d) |> 
    slice(1:list_size)
}

make_random_sublist <- function(prod_sum, list_size, k) {
  subk <- prod_sum |> 
    filter(num_langs >= k)
  idx <- sample(1:nrow(subk), min(list_size, nrow(subk)))
  subk |> slice(idx)
}

get_difficulty_cor <- function(sublist, test_items) {
  prod_res <- sublist |>
    left_join(test_items, by = "uni_lemma",
              multiple = "first")
  if (nrow(prod_res) == 0 || 
      complete.cases(prod_res$mean_d, prod_res$d) |> sum() == 0) {
    return(c(0, NA))
  }
  c((!is.na(prod_res$d)) |> sum(),
    cor(prod_res$mean_d, prod_res$d, 
        use = "complete.obs", method = COR_TYPE))
}

get_sumscore_cor <- function(sublist, xldf, all_prod, lang) {
  xldf_sub <- sublist |> 
    left_join(xldf |> filter(language == lang),
              by = "uni_lemma",
              multiple = "first") |> 
    filter(!is.na(uid))
  if (nrow(xldf_sub) == 0) {
    return(NA) 
  } else if (nrow(xldf_sub) == 1) {
    prod_sub <- all_prod[[lang]][,xldf_sub$uid]
  } else {
    prod_sub <- rowMeans(all_prod[[lang]][,xldf_sub$uid], na.rm = TRUE)
  }
  c(nrow(xldf_sub),
    cor(rowMeans(all_prod[[lang]], na.rm = TRUE), prod_sub,
      use = "na.or.complete", method = COR_TYPE))
}

get_fscore_cor <- function(sublist, xldf, lang, full_fscores) {
  full_model <- readRDS(glue("data/prod_models/{lang}_2PL_allforms_prod_fits.rds"))
  xldf_sub <- sublist |> 
    left_join(xldf |> filter(language == lang),
              by = "uni_lemma",
              multiple = "first") |> 
    filter(!is.na(uid))
  
  params_orig <- mod2values(full_model$model)
  
  new_vals <- params_orig |> 
    filter(item %in% xldf_sub$uid) |> 
    select(item, name, value) |> 
    pivot_wider(names_from = name, values_from = value) |> 
    mutate(a1 = 1) |> 
    select(-d) |> 
    left_join(xldf_sub |> select(uid, uni_lemma, mean_d), by = c("item" = "uid")) |> 
    rename(d = mean_d) |> 
    pivot_longer(cols = c(a1, d, g, u)) |> 
    select(-uni_lemma)
  
  params_new <- params_orig |> 
    left_join(new_vals, by = c("item", "name")) |> 
    mutate(value = coalesce(value.y, value.x), .after = parnum) |> 
    mutate(est = FALSE,
           item = `class<-`(item, value = "character")) |> 
    select(-value.x, -value.y)
  
  swad_mod <- mirt(full_model$model@Data$data, 1, pars = params_new)
  
  non_swad_ids <- setdiff(colnames(full_model$model@Data$data), xldf_sub$uid)
  masked_resps <- full_model$model@Data$data
  masked_resps[, non_swad_ids] <- NA
  swad_fscores <- fscores(swad_mod, 
                          response.pattern = masked_resps, 
                          method = "MAP")[,1] 
  
  cor(full_fscores, swad_fscores, 
      use = "complete.obs", method = COR_TYPE)
}

