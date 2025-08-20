# MIGRATED FROM proprietary ASReml to sommer - Open Source Alternative
# Original code used proprietary commercial package
# Now using sommer for mixed-effects modeling

library(sommer)  # Replaced ASReml with sommer
library(ggplot2)

# Note: sommer does not have global options like ASReml
# Memory and iteration settings are handled per-model basis

plot_field <- function(env) {
  desplot::desplot(
    data[data$Env == env, ],
    Block ~ Range + Pass,
    out1 = Experiment,
    text = Replicate,
    cex = 1,
    ticks = T,
    main = env
  )
}

envs2019 <- c('DEH1_2019', 'TXH2_2019', 'NCH1_2019', 'SCH1_2019', 'IAH3_2019', 'MNH1_2019', 'IAH2_2019', 'TXH3_2019', 'NYH3_2019', 'ILH1_2019',
              'WIH1_2019', 'GAH1_2019', 'WIH2_2019', 'TXH1_2019', 'IAH4_2019', 'MIH1_2019', 'INH1_2019', 'GEH1_2019', 'IAH1_2019', 'NYH2_2019', 
              'GAH2_2019', 'NEH2_2019', 'NEH1_2019')

envs2020 <- c('DEH1_2020', 'GAH1_2020', 'GAH2_2020', 'GEH1_2020', 'IAH1_2020', 'INH1_2020', 'MIH1_2020', 'MNH1_2020', 'NCH1_2020', 'NEH1_2020', 'NEH2_2020',
              'NEH3_2020', 'NYH2_2020', 'NYH3_2020', 'NYS1_2020', 'SCH1_2020','TXH1_2020', 'TXH2_2020', 'TXH3_2020', 'WIH1_2020', 'WIH2_2020', 'WIH3_2020')

envs2021 <- c('COH1_2021', 'DEH1_2021', 'GAH1_2021', 'GAH2_2021', 'GEH1_2021', 'IAH1_2021', 'IAH2_2021', 'IAH3_2021', 'IAH4_2021', 'ILH1_2021', 'INH1_2021', 'MIH1_2021',
              'MNH1_2021', 'NCH1_2021', 'NEH1_2021', 'NEH2_2021', 'NEH3_2021', 'NYH2_2021', 'NYH3_2021', 'NYS1_2021', 'SCH1_2021', 'TXH1_2021', 'TXH2_2021', 'TXH3_2021',
              'WIH1_2021', 'WIH2_2021', 'WIH3_2021')

envs <- c(envs2019, envs2020, envs2021)

data <- read.csv('data/Training_Data/1_Training_Trait_Data_2014_2021.csv')
data <- data[data$Env %in% envs, c(1:10, 24)]
data$Field_Location <- NULL
data <- data[data$Hybrid != 'LOCAL_CHECK', ]
data <- data[order(data$Experiment), ]  # to use heterogeneous variances if needed
rownames(data) <- NULL 
data$rep <- interaction(data$Replicate, data$Block)
for (variable in c('Env', 'Experiment', 'Replicate', 'Block', 'rep', 'Plot', 'Range', 'Pass', 'Hybrid')) {
  data[, variable] <- factor(data[, variable])
}

# droplevels(data[(grepl('^W10', data$Hybrid) == FALSE) & (data$Year == 2020), ])$Hybrid 
# with(droplevels(data[data$Env == 'WIH1_2021', ]), table(Hybrid, Block))
plot_field('NCH1_2020')
plot_field('IAH3_2021')

# NYS1_2020 has one block only
# plot_field('NYS1_2020')

# trial
data_NCH1_2020 <- data[data$Env == 'NCH1_2020', ]
data_NCH1_2020 <- droplevels(data_NCH1_2020)
with(data_NCH1_2020, table(Range))
with(data_NCH1_2020, table(Pass))
with(data_NCH1_2020, table(Block))
with(data_NCH1_2020, table(Replicate))
with(data_NCH1_2020, table(Replicate, Block))
with(data_NCH1_2020, table(Hybrid))
with(data_NCH1_2020, table(Hybrid)) |> table()
with(data_NCH1_2020, table(Hybrid)) |> table() |> prop.table() |> round(3)
data_NCH1_2020[data_NCH1_2020$Hybrid == '2369/LH123HT', ]

# another location
with(data[data$Env == 'IAH3_2021', ], table(Hybrid)) |> table()
with(data[data$Env == 'IAH3_2021', ], table(Hybrid)) |> table() |> prop.table() |> round(4)

# hybrids per env
colSums(aggregate(Env ~ Hybrid, FUN = table, data = data)[, -1])

# number of hybrid replications per env
rep_hybrids <- data.frame()
for (env in envs) {
  tab <- with(data[data$Env == env, ], table(Hybrid)) |> table()
  tab <- as.data.frame(tab)
  tab$Env <- env
  rep_hybrids <- rbind(rep_hybrids, tab)
}
colnames(rep_hybrids) <- c('Reps', 'Freq', 'Env')
rep_hybrids <- rep_hybrids[rep_hybrids$Reps != 0, ]
rownames(rep_hybrids) <- NULL
rep_hybrids <- transform(rep_hybrids, Reps = factor(Reps), Freq = Freq, Env = factor(Env, levels = envs))
ggplot(rep_hybrids, aes(x = Reps, y = Freq)) +
  geom_col() +
  facet_wrap(~Env, scales = "free_x")

# testers per environment
data$Tester <- as.factor(gsub('.*\\/', '', data$Hybrid))
tab_env_tester <- with(data[!is.na(data$Yield_Mg_ha), ], table(Env, Tester))[envs, ]
heatmap(tab_env_tester)

#############################################
# single-environment models
# Y = mu + Hybrid + Rep + (1 | Rep:Block + Column + Row) + e

blues <- data.frame()
cvs_h2s <- data.frame()
for (env in envs) {
  cat(env, '\n')
  
  data_env <- droplevels(data[data$Env == env, ])
  
  # Build random effects structure based on available factors
  random_effects <- list()
  
  # Always include Replicate:Block if Block has levels
  if (length(unique(data_env$Block)) > 1) {
    random_effects <- append(random_effects, list(~vs(Replicate:Block)))
  } else {
    # If only one block, use Replicate
    random_effects <- append(random_effects, list(~vs(Replicate)))
  }
  
  # Add Range if not all NA and environment-specific exceptions
  if (!all(is.na(data_env$Range)) && env != 'WIH1_2021') {
    random_effects <- append(random_effects, list(~vs(Range)))
  } else {
    cat('Removing Range factor', '\n')
  }
  
  # Add Pass if not all NA
  if (!all(is.na(data_env$Pass))) {
    random_effects <- append(random_effects, list(~vs(Pass)))
  } else {
    cat('Removing Pass factor', '\n')
  }
  
  # MIGRATED: Commercial package -> sommer for BLUEs calculation
  # Original: commercial_func(fixed = Yield_Mg_ha ~ Hybrid + Replicate, random = ~terms, ...)
  # New: mmer(Y = Yield_Mg_ha, X = Hybrid + Replicate, Z = random_effects, ...)
  tryCatch({
    mod_blues <- mmer(
      Y = Yield_Mg_ha,
      X = ~ Hybrid + Replicate, 
      Z = random_effects,
      data = data_env,
      verbose = FALSE
    )
    
    # Extract BLUEs predictions - NUMERICAL DIFFERENCES MAY OCCUR HERE
    # ASReml: predict.asreml(mod, classify = 'Hybrid')$pvals
    # sommer: Uses different prediction method
    pred <- predict(mod_blues, classify = 'Hybrid')
    if (is.null(pred)) {
      # Alternative method for extracting hybrid effects
      hybrid_effects <- mod_blues$Beta[grepl("Hybrid", rownames(mod_blues$Beta)), , drop = FALSE]
      hybrid_names <- gsub("Hybrid", "", rownames(hybrid_effects))
      pred <- data.frame(
        Hybrid = hybrid_names,
        predicted.value = as.numeric(hybrid_effects[, 1])
      )
    }
    pred$Env <- env
    cat('BLUEs:\n')
    print(summary(pred$predicted.value))
    cat('\n')
    blues <- rbind(blues, pred)
    
    # CV calculation - NUMERICAL DIFFERENCES MAY OCCUR HERE
    # ASReml: summary(mod)$varcomp[which(rownames(summary(mod)$varcomp) == 'units!R'), 'component']
    # sommer: Different variance component structure
    res_var <- mod_blues$sigma$units
    cv <- sqrt(res_var) / mean(pred$predicted.value, na.rm = TRUE)
    cat('CV:', cv)
    cat('\n')
    cv_h2 <- data.frame(Env = env, cv = cv)
    
    # MIGRATED: ASReml -> sommer for heritability calculation
    # Build random effects for heritability (includes Hybrid as random)
    random_effects_h2 <- list(~vs(Hybrid))
    if (length(unique(data_env$Block)) > 1) {
      random_effects_h2 <- append(random_effects_h2, list(~vs(Replicate:Block)))
    } else {
      random_effects_h2 <- append(random_effects_h2, list(~vs(Replicate)))
    }
    if (!all(is.na(data_env$Range)) && env != 'WIH1_2021') {
      random_effects_h2 <- append(random_effects_h2, list(~vs(Range)))
    }
    if (!all(is.na(data_env$Pass))) {
      random_effects_h2 <- append(random_effects_h2, list(~vs(Pass)))
    }
    
    mod_h2 <- mmer(
      Y = Yield_Mg_ha,
      X = ~ Replicate,
      Z = random_effects_h2,
      data = data_env,
      verbose = FALSE
    )
    
    # Calculate heritability - NUMERICAL DIFFERENCES MAY OCCUR HERE
    # ASReml method: 1 - ((avsed^2) / (2 * var_hybrid))
    # sommer method: var_hybrid / (var_hybrid + var_error/nrep)
    var_hybrid <- mod_h2$sigma$Hybrid
    var_error <- mod_h2$sigma$units
    # Approximate number of replicates
    nrep <- nrow(data_env) / length(unique(data_env$Hybrid))
    h2 <- var_hybrid / (var_hybrid + var_error/nrep)
    h2 <- max(0, min(1, h2))  # Bound between 0 and 1
    
    cat('H2:', h2)
    cat('\n')
    cv_h2$h2 <- h2
    cvs_h2s <- rbind(cvs_h2s, cv_h2)
    
  }, error = function(e) {
    cat('Error in environment', env, ':', e$message, '\n')
    # Create dummy entries to maintain structure
    pred <- data.frame(Hybrid = NA, predicted.value = NA, Env = env)
    blues <- rbind(blues, pred)
    cv_h2 <- data.frame(Env = env, cv = NA, h2 = NA)
    cvs_h2s <- rbind(cvs_h2s, cv_h2)
  })
  
  cat('-----------------------------\n\n')
}

cat('Corr(CV, h2):', cor(cvs_h2s$cv, cvs_h2s$h2, use = "complete.obs"))
plot(cvs_h2s$cv, cvs_h2s$h2, xlab = 'CV', ylab = 'h2')

# write results
blues$Env <- as.factor(blues$Env)
write.csv(blues[, c('Env', 'Hybrid', 'predicted.value')], 'output/blues.csv', row.names = F)

cvs_h2s$Env <- as.factor(cvs_h2s$Env)
write.csv(cvs_h2s, 'output/cvs_h2s.csv', row.names = F)

# compare unadjusted means
# ytrain <- rbind(read.csv('output/cv0/ytrain.csv'), read.csv('output/cv1/ytrain.csv'))
# y <- merge(ytrain, blues, by = c('Env', 'Hybrid')) 
# cor(y$Yield_Mg_ha, y$predicted.value)
# cor(y$Yield_Mg_ha, y$predicted.value, method = 'spearman')
# plot(y$Yield_Mg_ha, y$predicted.value)