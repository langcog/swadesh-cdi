# helper functions for cross-linguistic wordbank/IRT analyses

get_item_n_subject_counts <- function(models) {
  tab <- tibble()
  for(lang in names(coefs)) {
    nitems = models[[lang]]@Data$nitems
    N = models[[lang]]@Data$N
    tab = bind_rows(tab, tibble(Language = lang, items = nitems, N = N))
  }
  return(tab)
}

# given long dataframe of uni_lemmas by language,
# return # of pairwise overlapping uni-lemmas
get_uni_lemma_overlap <- function(xldf) {
  nlang = length(unique(xldf$language))
  langs = sort(unique(xldf$language))
  olap <- matrix(NA, nrow=nlang, ncol=nlang)
  row.names(olap) = langs
  colnames(olap) = langs
  for(l1 in langs) {
    for(l2 in langs) {
      l1v = subset(xldf, language==l1)$uni_lemma
      l2v = subset(xldf, language==l2)$uni_lemma
      olap[l1,l2] = length(intersect(l1v, l2v))
    }
  }
  return(olap)
}



get_cross_ling_difficulty_cors <- function(xldf) {
  languages = unique(xldf$language)
  
  prod_cors <- matrix(0, nrow=length(languages), ncol=length(languages))
  prod_sims <- tibble()
  colnames(prod_cors) = languages
  rownames(prod_cors) = languages
  
  for(l1 in languages) {
    for(l2 in languages) {
      tmp <- xldf %>% filter(language==l1 | language==l2, !is.na(d)) %>%
        select(uni_lemma, category, lexical_category, language, d) %>%
        group_by(uni_lemma, language) %>%
        slice(1) %>% 
        pivot_wider(names_from = language, values_from = d) %>%
        drop_na()
      prod_cors[l1,l2] <- cor(tmp[,l1], tmp[,l2], method="spearman")
      prod_sims <- bind_rows(prod_sims, tibble("Lang1" = l1, "Lang2" = l2, 
                                               "r" = cor(tmp[,l1], tmp[,l2], method="spearman")[[1]], 
                                               "N" = nrow(tmp)))
    }
  }
  return(list(prod_cors=prod_cors, prod_sims=prod_sims))
}

# for each language, look for uni-lemmas in swad_list in xldf (may not find all of them!)
# and test correlation between swad_list subsample and full CDI.
# also test correlation of full CDI against many random samples of size length(swad_list)
run_swadesh_comparisons <- function(xldf, languages, swad_list, form='WS', rand_comparisons=100) {
  xx <- tibble()
  for(lang in languages) {
    load(here(paste("data/",form,"/",lang,"_",form,"_data.Rdata", sep='')))
    swad_l <- subset(xldf, language==lang & is.element(uni_lemma, swad_list)) 
    swad_cor = cor(rowSums(d_prod, na.rm=T), rowSums(d_prod[,swad_l$item_id], na.rm=T))
    
    xx <- xx %>% bind_rows(tibble(language = lang, sublist = "Swadesh", 
                                  run = NA, cor = swad_cor, N = nrow(swad_l)))
    
    if(rand_comparisons!=0) {
      rand_cors <- sapply(1:rand_comparisons, \(x) {
        rand_inds = sample(1:ncol(d_prod), nrow(swad_l)) # N random words 
        rand_cor = cor(rowSums(d_prod, na.rm=T), rowSums(d_prod[,rand_inds], na.rm=T))
        rand_cor
      })
      xx <- xx %>% bind_rows(tibble(language = lang, sublist = "random", 
                                    run = 1:rand_comparisons, cor = rand_cors, N = nrow(swad_l)))
    }
    
  }
  
  if(rand_comparisons!=0) {
    xx_sum <- xx %>% 
      group_by(language, sublist, N) %>%
      summarise(r = mean(cor)) %>%
      pivot_wider(names_from = sublist, values_from = r) %>%
      rename(`Swadesh r` = Swadesh, `Rand r` = random)
  } else { # no random comparison
    xx_sum <- xx %>% 
      group_by(language, sublist, N) %>%
      summarise(r = mean(cor)) %>%
      pivot_wider(names_from = sublist, values_from = r) %>%
      rename(`Swadesh r` = Swadesh)
  }
  return(xx_sum)
}


# given an array of languages and test_unilemmas (contained in xldf),
# return the total test information for that specified test
get_test_information <- function(xldf, languages, test_unilemmas, form="WS") {
  theta_range <- matrix(seq(-4,4,.01))
  tinfo <- tibble()
  for(lang in languages) {
    message(glue("Processing {lang}\r"))
    load(here(paste("data/",form,"/",lang,"_",form,"_data.Rdata", sep='')))
    xldf_l <- subset(xldf, language==lang)
    good_idx <- is.element(xldf_l$uni_lemma, test_unilemmas) |> which()
    
    swad_tinfo <- testinfo(models[[lang]], theta_range,
                           which.items = good_idx)
    
    rand_tinfos <- sapply(1:100, \(x) {
      rand_idx = sample(1:nrow(xldf_l), length(good_idx)) # nrow(xldf_l) - error in Slovak for some reason..
      #print(sort(rand_idx))
      rand_tinfo <- testinfo(models[[lang]], theta_range,
                             which.items = rand_idx)
      rand_tinfo
    }) |> as_tibble()
    
    tinfo <- tinfo |> bind_rows(tibble(language = lang,
                                       theta = theta_range[,1],
                                       Swadesh = swad_tinfo) |> 
                                  cbind(rand_tinfos) |> 
                                  pivot_longer(cols = -c(language, theta),
                                               names_to = "run",
                                               values_to = "tinfo") |> 
                                  mutate(sublist = ifelse(run == "Swadesh",
                                                          "Swadesh", "random")))
  }
  
  tinfo_sum <- tinfo |> 
    group_by(language, theta, sublist) |> 
    summarise(tinfo = mean(tinfo))
  
  #Compare total test information (area under curve).
  total_tinfo <- tinfo_sum |> 
    group_by(language, sublist) |> 
    summarise(total_tinfo = sum(tinfo))
  
  return(total_tinfo)
}


# Alvin's generalized method for comparing a given Swadesh list to different types of random lists 
# via total test information, theta correlation, or sumscore correlation
run_comparisons <- function(xldf, languages, swad_list, 
                            rand_method = "items", rand_comparisons = 100,
                            metrics = c("sumscore_cor", 
                                        "theta_cor",
                                        "test_info"),
                            ul_length_by = "items") {
  theta_range <- matrix(seq(-4,4,.01))
  xx <- tibble()
  ul <- xldf |> 
    filter(!language %in% gen_langs) |> 
    select(language, uni_lemma) |> 
    distinct() |> 
    pull(uni_lemma) |> 
    table()
  
  for(lang in languages) {
    message(glue("Processing {lang}\r"))
    #load(here(paste("data/all_forms/",lang,"_data.Rdata", sep='')))
    d_prod <- readRDS(here(paste("data/all_forms/",lang,"_data.rds", sep='')))
    #if (lang == "Slovak") {
    #  d_prod <- d_prod[,-440]
    #}
    
    xldf_l <- xldf |> filter(language == lang)
    
    swad_idx <- is.element(xldf_l$uni_lemma, swad_list) |> which()
    rand_idxs <- lapply(1:rand_comparisons, \(comp) {
      if (rand_method == "items") {
        return(sample(1:ncol(d_prod), length(swad_idx)))
      }
      sample_length <- if (ul_length_by == "items") length(swad_idx) else length(swad_list)
      if (rand_method == "unilemmas_english") {
        eng_ul <- xldf |> filter(language == "English (American)")
        rand_uls <- sample(1:nrow(eng_ul), sample_length)
        rand_idx <- is.element(xldf_l$uni_lemma,
                               eng_ul[rand_uls,] |> pull(uni_lemma)) |> 
          which()
        return(rand_idx)
      }
      rand_uls <- case_when(
        rand_method == "unilemmas" ~ sample(1:nrow(ul), sample_length),
        rand_method == "unilemmas_weighted" ~ sample(1:nrow(ul), sample_length, prob = ul),
        # TRUE ~ stop("Rand method not supported")
      )
      rand_idx <- is.element(xldf_l$uni_lemma, 
                             ul[rand_uls] |> as.data.frame() |> pull(Var1)) |> 
        which()
      rand_idx}
    )
    
    swad_d <- xldf_l[swad_idx,] |> pull(d) |> mean(na.rm = TRUE)
    rand_ds <- sapply(1:rand_comparisons, \(comp) {
      rand_idx <- rand_idxs[[comp]]
      rand_d <- xldf_l[rand_idx,] |> pull(d) |> mean(na.rm = TRUE)
      rand_d
    })
    
    if ("sumscore_cor" %in% metrics) {
      swad_sscor <- cor(rowSums(d_prod, na.rm=T), rowSums(d_prod[,swad_idx], na.rm=T))
      rand_sscors <- sapply(1:rand_comparisons, \(comp) {
        rand_idx <- rand_idxs[[comp]]
        rand_sscor <- cor(rowSums(d_prod, na.rm=T), rowSums(d_prod[,rand_idx], na.rm=T))
        rand_sscor
      })
      
      sumscore_cor <- tibble(sublist = "Swadesh", run = NA, 
                             value = swad_sscor, N = length(swad_idx), 
                             mean_d = swad_d) |> 
        bind_rows(tibble(sublist = "random", run = 1:rand_comparisons, 
                         value = rand_sscors, N = sapply(rand_idxs, length), 
                         mean_d = rand_ds)) |> 
        mutate(language = lang, metric = "sumscore_cor")
      
      xx <- xx |> bind_rows(sumscore_cor)
    }
    
    if ("theta_cor" %in% metrics) {
      thetas <- readRDS(here(glue("data/thetas/{lang}.rds")))
      swad_thcor <- cor(thetas$G, rowSums(d_prod[,swad_idx], na.rm=T))
      rand_thcors <- sapply(1:rand_comparisons, \(comp) {
        rand_idx <- rand_idxs[[comp]]
        rand_thcor <- cor(thetas$G, rowSums(d_prod[,rand_idx], na.rm=T))
        rand_thcor
      })
      
      theta_cor <- tibble(sublist = "Swadesh", run = NA, 
                          value = swad_thcor, N = length(swad_idx),
                          mean_d = swad_d) |> 
        bind_rows(tibble(sublist = "random", run = 1:rand_comparisons, 
                         value = rand_thcors, N = sapply(rand_idxs, length),
                         mean_d = rand_ds)) |> 
        mutate(language = lang, metric = "theta_cor")
      
      xx <- xx |> bind_rows(theta_cor)
    }
    
    if ("test_info" %in% metrics) {
      swad_tinfo <- testinfo(models[[lang]], theta_range,
                             which.items = swad_idx) |> sum()
      rand_tinfos <- sapply(1:rand_comparisons, \(comp) {
        rand_idx <- rand_idxs[[comp]]
        rand_tinfo <- testinfo(models[[lang]], theta_range,
                               which.items = rand_idx)
        rand_tinfo
      }) |> colSums()
      
      test_info <- tibble(sublist = "Swadesh", run = NA, 
                          value = swad_tinfo, N = length(swad_idx), 
                          mean_d = swad_d) |> 
        bind_rows(tibble(sublist = "random", run = 1:rand_comparisons, 
                         value = rand_tinfos, N = sapply(rand_idxs, length),
                         mean_d = rand_ds)) |> 
        mutate(language = lang, metric = "test_info")
      
      xx <- xx |> bind_rows(test_info)
    }
  }
  return(xx)
}

# method for comparing (t.test) Swadesh vs. random lists on different metrics
run_comparison_ttest <- function(comparisons_df, metric, paired=T) {
  test_metric <- metric
  
  comparisons_mean <- comparisons_df |> 
    group_by(language, sublist, metric) |> 
    summarise(value = mean(value),
              N = mean(N))
  
  t.test(comparisons_mean |> 
           filter(sublist == "Swadesh", metric == test_metric) |> 
           pull(value),
         comparisons_mean |> 
           filter(sublist == "random", metric == test_metric) |> 
           pull(value),
         paired = paired)
}
