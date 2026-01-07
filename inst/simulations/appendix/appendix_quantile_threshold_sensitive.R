rm(list=ls())

# Installing only if necessary
if (!requireNamespace("devtools", quietly = TRUE)) {
     install.packages("devtools")
}

library(devtools)
seed <- 42
set.seed(seed)

# Loading all package dependencies -- you can also install the xflex4cast package instead,
# and load it : library(xflex4cast)
# The degpd_path can change locally
xflex4cast_path <- "C:/Users/mm538r/OneDrive - University of Glasgow/p4r_enhacements/worst-case-scenario/p4r_techometrics_paper/rebuttal/simulations/xflex4cast/"
degpd_path_r <- paste0(xflex4cast_path,"inst/degpd-and-zidegpd/R/")
degpd_path_src <- paste0(xflex4cast_path,"inst/degpd-and-zidegpd/src/")
"C:/Users/mm538r/OneDrive - University of Glasgow/p4r_enhacements/worst-case-scenario/p4r_techometrics_paper/rebuttal/simulations/xflex4cast/inst/degpd-and-zidegpd/R/"
devtools::load_all(xflex4cast_path) # 

load_degpd(dest_r = degpd_path_r,
           dest_src = degpd_path_src)
# Loading all functions necessary to run the competing model DEGPD.
# For more details: (Ahmad, T., Gaetan, C., & Naveau, P. (2024).
# An extended generalized Pareto regression model for count data.
# *Statistical Modelling*, 0(0).
#https://doi.org/10.1177/1471082X241266729)

# Setting the parameters of the simulation
n <- 5000
phi <- 0.25
n_replications <- 1000 # For illustration purposes we reduced the number of replications to 100, while in the main manuscript is 1000
xi <- 0.3 # Shape parameter
probs <- sort(unique(c(0.05,0.25,0.5,0.75,(1-phi),0.95,0.99,0.999,0.9999))) # Quantile levels to calculate on

# Extra threhold quantiles
extra_quantiles_range <- c(-0.05,+0.05)+(1-phi)
if(extra_quantiles_range[2]==1){
     extra_quantiles_range[2] <- 0.99
}

extra_quantiles <- seq(extra_quantiles_range[1],extra_quantiles_range[2],by = 0.01)
og_probs <- probs # Original probabilities
probs <- sort(unique(c(probs,extra_quantiles)))
probs

## Creating another list for the model quantile evaluation
store_replications_results <- vector("list",length = n_replications)
names(store_replications_results) <- paste0("replication",1:n_replications)

store_replications_tail_parameters_results <- vector("list",length = n_replications)
names(store_replications_tail_parameters_results) <- paste0("replication",1:n_replications)

store_replications_other_models <- vector("list",length = n_replications)
names(store_replications_other_models) <- paste0("replication",1:n_replications)

## Log for knowing what's is being run
print(n)
print(phi)
print(n_replications)
print(xi)
print(probs)
print(extra_quantiles_range)

for( i in 1:n_replications ){
     
          try({
               all_quantiles_results <- setNames(vector("list",length = length(extra_quantiles)),extra_quantiles)
               
               store_replications_results[[i]] <- all_quantiles_results
               
               train_data <- simulation_scenario_one(n = n,
                                                     phi = phi,
                                                     xi = xi,
                                                     probs = probs,
                                                     index = "fit_data_row")
               
               test_data <- simulation_scenario_one(n = n,
                                                    phi = phi,
                                                    xi = xi,
                                                    probs = probs,
                                                    index = "fit_data_row")

               
               std_qgam_def_quantiles <- probs[probs<=extra_quantiles_range[2]]
               
               # Defining the qgam model only
               std_qgam_def <- define_xflex4cast(mqgam_formula = y ~ s(x1,k=10) ,
                                                 mqgam_quantiles = std_qgam_def_quantiles,
                                                 quantile_threshold = 1-phi,
                                                 fit_tail = FALSE)
               
               tail_qgam_def <- define_xflex4cast(mqgam_formula = y ~ s(x1,k=10) ,
                                                  mqgam_quantiles = probs[probs>extra_quantiles_range[2]],
                                                  quantile_threshold = probs[probs>extra_quantiles_range[2]][1],
                                                  fit_tail = FALSE)
               

               # Fitting the models
               qgam_fit <- fit(model = std_qgam_def,
                               data = train_data$data,
                               index = "fit_data_row")
               
               tail_qgam_fit <- fit(model = tail_qgam_def,
                                    data = train_data$data,
                                    index = "fit_data_row")
               
               # Defining the formula for the DEGP approach
               formula_degp <- list(lsigma = y ~ s(x1,k=10), lxi = ~1, lkappa = ~1)
               
               # Extended Discrete Generalised Pareto
               degpd_fit <- devgam(formula_degp,data = as.data.frame(train_data$data),
                                   family  = "degpd",
                                   degpd.args = list(m=1),
                                   trace = 0)
               
               # ===========================================================
               ## Obtaining and fit what would be considered a perfect model
               # ===========================================================
               
               # 0.1 Obtaining the true exceedance samples
               true_exceedance <- copy(train_data$data)
               true_exceedance <- true_exceedance[y >= (true_u),]
               true_exceedance[,exceedance := y-(true_u)]
               
               # true_exceedance <- data.table(exceedance = rep(0,1000))
               # dgp_samples <- sample(x = 0:5000,size = nrow(true_exceedance),prob = ddgp(x = 0:5000,shape = 0.3,scale = 2.5),replace = TRUE)
               # true_exceedance[,exceedance:=dgp_samples]
               
               gjrm_gamlss <- GJRM::gamlss
               gjrm_call <- call("gjrm_gamlss",
                                 formula = list(exceedance ~ 1, ~1),
                                 family = "DGP",
                                 data = true_exceedance
               )
               
               true_tail <-  eval(gjrm_call)
               true_tail_summary <- summary(true_tail)
               
               
               for(j in extra_quantiles){
                    
                    # Using the flex4cast approach
                    xflex4cast_def <- define_xflex4cast(mqgam_formula = y ~ s(x1,k=10) ,
                                                        scale_formula = ~ 1,
                                                        mqgam_quantiles = probs,
                                                        quantile_threshold = j,
                                                        fit_tail = TRUE)
                    
                    
                    xflex4cast_fit <- xflex4cast_def
                    xflex4cast_fit$mqgam_fit <- qgam_fit$mqgam_fit
                    
                    # I can split what are the mqgam quantiles and the tail quantiles
                    xflex4cast_fit$mqgam_quantiles <- sort(unique(c(og_probs[og_probs<=xflex4cast_fit$quantile_threshold],xflex4cast_fit$quantile_threshold)))
                    xflex4cast_fit$tail_quantiles <- sort(unique(c(probs[probs>=(xflex4cast_fit$quantile_threshold+1e-10)])))
                    
                    xflex4cast_fit <- fit(model = xflex4cast_fit,
                                          data = train_data$data,
                                          index = "fit_data_row",
                                          fit_tail_only = TRUE)
               
                    tail_summary <- summary(xflex4cast_fit$tail_mod)

                    
                    # Storing models parameters
                    store_replications_tail_parameters_results[[i]][[as.character(j)]] <- data.table:::data.table(xi_mean = tail_summary$tableP1[1,1],
                                                                                                                  xi_sd = tail_summary$tableP1[1,2],
                                                                                                                  sigma_mean = tail_summary$tableP2[1,1],
                                                                                                                  sigma_sd = tail_summary$tableP2[1,2],
                                                                                                                  n_tail = xflex4cast_fit$tail_mod$n,
                                                                                                                  true_xi_mean = true_tail_summary$tableP1[1,1],
                                                                                                                  true_xi_sd = true_tail_summary$tableP1[1,2],
                                                                                                                  true_sigma_mean = true_tail_summary$tableP2[1,1],
                                                                                                                  true_sigma_sd = true_tail_summary$tableP2[1,2],
                                                                                                                  true_n_tail = true_tail$n)
                    
                    # Obtaining the predictions
                    oos_pred_xflex4cast <- predict(xflex4cast_fit,
                                                   newdata = test_data$data[,fit_data_row:=1:.N],
                                                   indexes = "fit_data_row")
                    
                    
                    # Creating the table for the true quantiles in a long format
                    if(is.vector(test_data$true_quantile)){
                         true_qt_matrix <- test_data$true_quantiles
                         colnames(true_qt_matrix) <- probs
                         true_quantile <- data.table::melt(as.data.table(true_qt_matrix,keep.rownames = TRUE)[,fit_data_row:=.I],
                                                           id.vars = "fit_data_row", variable.name = "quantile",value.name = "value")
                    } else {
                         true_quantile <-  data.table::melt(as.data.table(test_data$true_quantiles,keep.rownames = TRUE)[,fit_data_row:=.I],
                                                            id.vars = "fit_data_row", variable.name = "quantile",value.name = "value")
                    }
                    
                    true_quantile <- true_quantile[,value:=ifelse(is.na(value),1000,value)]
                    true_quantile <- true_quantile[,quantile:=as.numeric(as.character(quantile))]
                    
                    # Standardising the quantiles
                    true_quantile[,quantile:=as.character(round(quantile,digits = 4))]
                    oos_pred_xflex4cast[,quantile:=as.character(round(quantile,digits = 4))]
                    
                    all_pred <- merge(true_quantile,oos_pred_xflex4cast, by = c("fit_data_row","quantile"),suffixes = c(".true",".tail"))
               
                    # Generating the new columns (doing the same for the forecast)
                    rmse_summary <- all_pred[,.(true_avg = mean(value.true,na.rm=TRUE),
                                                rmse_xflex4cast = sqrt(mean((value.true-round(value.tail))^2)),
                                                mae_xflex4cast = mean(abs(value.true-round(value.tail)))),
                                             by = c("quantile")]
                    
                    store_replications_results[[i]][[as.character(j)]] <- rmse_summary[,replication_number:=i]
                    
                    # Progress bar
                    progress <- round((i / n_replications) * 100)
                    bar <- paste0(rep("=", progress / 2), collapse = "")
                    cat(sprintf("\rProgress: [%-50s] %3d%%", bar, progress))
                    flush.console()
               }          
               
               # Obtaining the predictions (from other models)
               oos_pred_tail_qgam <- predict(tail_qgam_fit,
                                        newdata = test_data$data[,fit_data_row:=1:.N],
                                        indexes = "fit_data_row")
               oos_pred_tail_qgam[,quantile:=as.character(round(quantile,digits = 4))]
               
               # Obtaining out-of-sample observations
               oos_pred_degpd <- as.data.table(predict(degpd_fit,as.data.frame(test_data$data),prob = probs))[,fit_data_row := 1:.N]
               oos_pred_degpd <- data.table::melt(oos_pred_degpd,id.vars = "fit_data_row",variable.name = "quantile",value.name = "value")
               oos_pred_degpd[,quantile:=as.character(round(as.numeric(as.character(quantile)),digits = 4))]
               
               
               all_pred <- merge(true_quantile,oos_pred_tail_qgam, by = c("fit_data_row","quantile"),suffixes = c(".true",".qgam"))
               all_pred <- merge(all_pred, oos_pred_degpd, by = c("fit_data_row","quantile"))
               all_pred <- setnames(all_pred,old = "value","value.degpd")
               
               # Generating the new columns (doing the same for the forecast)
               rmse_summary_others <- all_pred[,.(rmse_qgam = sqrt(mean((value.true-round(value.qgam))^2)),
                                           rmse_degpd = sqrt(mean((value.true-round(value.degpd))^2)),
                                           mae_qgam = mean(abs(value.true-round(value.qgam))),
                                           mae_degpd = mean(abs(value.true-round(value.degpd)))),
                                        by = "quantile"]
          
               store_replications_other_models[[i]] <- rmse_summary_others[,replication_number:=i]
               
          },silent = TRUE)
     
}

summary_simulation_experiments <- list(quantiles = store_replications_results,
                                       parameters = store_replications_tail_parameters_results,
                                       other_models = store_replications_other_models)

saveRDS(summary_simulation_experiments,
        file = paste0("C:/Users/mm538r/OneDrive - University of Glasgow/p4r_enhacements/worst-case-scenario/p4r_techometrics_paper/rebuttal/simulations/n_",n,"_seed_",seed,"_n_rep_",n_replications,"_sim_true_quantile_",round((1-phi),digits = 2),".rds") )

# # Combine all data.tables into one
# combined <- data.table::rbindlist(store_replications_results)
# 
# # Computing mean and sd by quantile
# summary_rmse <- combined[, .(
#      mean_rmse_qgam  = mean(rmse_qgam),
#      sd_rmse_qgam    = stats::sd(rmse_qgam),
#      mean_rmsexflex4cast  = mean(rmsexflex4cast),
#      sd_rmsexflex4cast    = stats::sd(rmsexflex4cast),
#      mean_rmse_degpd = mean(rmse_degpd),
#      sd_rmse_degpd   = stats::sd(rmse_degpd)
# ), by = quantile]
# 
# # Visualizing the results as in the paper
# summary_formatted <- summary_rmse[, .(
#      quantile,
#      rmsexflex4cast  = sprintf("%.2f ± %.2f", mean_rmsexflex4cast, sd_rmsexflex4cast),
#      rmse_qgam  = sprintf("%.2f ± %.2f", mean_rmse_qgam, sd_rmse_qgam),
#      rmse_degpd = sprintf("%.2f ± %.2f", mean_rmse_degpd, sd_rmse_degpd)
# )]

# aux <- all_models$quantiles
# avg_list <- setNames(
#      lapply(names(aux[[1]]), function(thr) {
#           tmp <- rbindlist(
#                lapply(aux, function(rep) rep[[thr]]),
#                use.names = TRUE, fill = TRUE
#           )
#           # average all numeric columns except 'quantile' and (optionally) 'replication_number'
#           numcols <- setdiff(
#                names(tmp)[vapply(tmp, is.numeric, TRUE)],
#                c("quantile", "replication_number")
#           )
#           tmp[, lapply(.SD, mean), by = quantile, .SDcols = numcols]
#      }),
#      names(aux[[1]])
# )
# 
# par(mfrow=c(2,2))
# for(tail_quantile in c(0.9,0.95,0.99,0.999,0.9999)){
#      result <- unlist(lapply(avg_list, function(x)x[quantile==tail_quantile,][["rmse_xflex4cast"]]))
# 
#      plot(names(result),result,xlab = "Quantile Threshold", ylab = "RMSE",main = paste0("q_",tail_quantile))
#      abline(v = 0.8,lty = 'dashed', col = 'red')
# }
# 
# 
# par(mfrow=c(4,3))
# for(tail_quantile in names(avg_list)){
#      result <- unlist(avg_list[[tail_quantile]][["rmse_xflex4cast"]])
#      names(result) <- avg_list[[tail_quantile]][["quantile"]]
#      plot(names(result),result,xlab = "Quantiles",
#           ylab = "rmse",
#           main = paste0("q_",tail_quantile),
#           ylim = c(0,25))
#      # abline(v = 0.9,lty = 'dashed', col = 'red')
# }
