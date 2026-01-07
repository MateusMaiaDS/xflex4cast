rm(list=ls())

library(data.table)
library(ggplot2)
# Defining the simulation parameters
n <- 5000
seed <- 42 
n_replications <- 1000
phi <- 0.05
set.seed(seed)

## Loading the results
simulation_results <- readRDS(paste0("C:/Users/mm538r/OneDrive - University of Glasgow/p4r_enhacements/worst-case-scenario/p4r_techometrics_paper/rebuttal/simulations/n_",n,"_seed_",seed,"_n_rep_",n_replications,"_sim_true_quantile_",round((1-phi),digits = 2),".rds"))

## Plotting the quantiles only
quantile_results <- simulation_results$quantiles

get_avg_list <- function(quantile_results){
     avg_list <- setNames(
          lapply(names(quantile_results[[1]]), function(thr) {
               tmp <- rbindlist(
                    lapply(quantile_results, function(rep) rep[[thr]]),
                    use.names = TRUE, fill = TRUE
               )
               # average all numeric columns except 'quantile' and 'replication_number'
               numcols <- setdiff(
                    names(tmp)[vapply(tmp, is.numeric, TRUE)],
                    c("quantile", "replication_number")
               )
               tmp[, {
                    means <- lapply(.SD, mean)
                    sds   <- lapply(.SD, sd)
                    names(means) <- paste0(numcols, "_mean")
                    names(sds)   <- paste0(numcols, "_sd")
                    c(means, sds)
               }, by = quantile, .SDcols = numcols]
               
               
          }),
          names(quantile_results[[1]])
     )
     
     return(avg_list)
}

if(phi==0.05){
     # Dropping the NULL last element
     simulation_results$quantiles <- lapply(
          simulation_results$quantiles,
          function(rep_list) {
               rep_list[setdiff(names(rep_list), "0.99")]
          }
     )
     
     # Dropping the NULL last element
     simulation_results$other_models <- lapply(
          simulation_results$other_models,
          function(rep_list) {
               rep_list[setdiff(names(rep_list), "0.99")]
          }
     )
     
}

avg_list <- get_avg_list(simulation_results$quantiles)


# avg_dt_qgam <- do.call(rbind,simulation_results$other_models)[,.(rmse_qgam_mean = mean(rmse_qgam),
#                                                                    rmse_qgam_sd = sd(rmse_qgam)),
#                                                                 by = c("quantile")]

# par(mfrow=c(2,2))
# for(tail_quantile in c(0.95,0.99,0.999,0.9999)){
#      result <- unlist(lapply(avg_list, function(x)x[quantile==tail_quantile,][["rmse_xflex4cast_mean"]]))
#      
#      plot(names(result),result,xlab = "Quantile Threshold", ylab = "RMSE",main = paste0("q_",tail_quantile))
#      if(tail_quantile %in% avg_dt_qgam[["quantile"]]){
#           points(0.9,avg_dt_qgam[quantile == tail_quantile,][["rmse_qgam_mean"]], col = "blue", pch = 4)
#      }
#      abline(v = 0.9,lty = 'dashed', col = 'red')
# }


## writing tidy version to the paper
## building the df
tail_q <- c(0.95,0.99,0.999,0.9999)

df_rmse <- rbindlist(lapply(names(avg_list), function(th)
     avg_list[[th]][quantile %in% tail_q,
                    .(quantile,
                      quantile_threshold = as.numeric(th),
                      rmse_mean = rmse_xflex4cast_mean/true_avg_mean,
                      rmse_sd   = rmse_xflex4cast_sd)]
))

# df_qgam <- avg_dt_qgam[quantile %in% tail_q,
#                        .(quantile,
#                          quantile_threshold = 1-phi,
#                          rmse_mean = rmse_qgam_mean,
#                          rmse_sd   = rmse_qgam_sd)]

df_rmse[, quantile_label :=
             paste0("alpha == ", quantile)]


p <- ggplot(df_rmse, aes(x = quantile_threshold, y = rmse_mean)) +
     geom_line() +
     geom_point() +
     # geom_errorbar(aes(ymin = rmse_mean - rmse_sd,
     #                   ymax = rmse_mean + rmse_sd),
     #               width = 0.0) +
     geom_vline(xintercept = 1-phi, linetype = "dashed", colour = "red") +
     facet_wrap(~quantile_label, scale = "free_y",labeller = label_parsed) +
     labs(x = expression(alpha^(T)),
          y = expression(RMSE),
          # title = expression("RMSE for "~hat(q)[alpha](bold(x)[d*t])~
          #                         " at upper-tail probability levels")
          ) +
     theme_bw(base_size = 11)
p
ggsave(paste0("C:/Users/mm538r/OneDrive - University of Glasgow/p4r_enhacements/worst-case-scenario/p4r_techometrics_paper/rebuttal/figures/appendix_threshold_",(1-phi),"_scaled_rmse_tail_quantiles.pdf"), p, device = cairo_pdf,
       width = 7, height = 5)

stop("Not updating the other plots")

## Looking into the parameters of the model
parameters_results <- simulation_results$parameters

library(data.table)

# Preallocate a list the size of parameters_results
get_parameter_dt <- function(parameters_results){
     
     out_list <- vector("list", length(parameters_results))

     for (i in seq_along(parameters_results)) {
          
          inner <- parameters_results[[i]]
          dt_i <- rbindlist(inner, idcol = "quantile_threshold")
          
          dt_i[, quantile_threshold := as.numeric(quantile_threshold)]
          dt_i[, replication_number := i]
          
          out_list[[i]] <- dt_i
     }
     
     out <- rbindlist(out_list)
     
     return(out)
}

dt_parameters_results <- get_parameter_dt(parameters_results)
avg_col <- setdiff(names(dt_parameters_results),c("quantile_threshold","replication_number"))
summary_dt <- dt_parameters_results[,{
     means <- lapply(.SD, mean)
     sds   <- lapply(.SD, sd)
     names(means) <- paste0(avg_col, "_mean")
     names(sds)   <- paste0(avg_col, "_sd")
     c(means, sds)
     },
     by = "quantile_threshold", .SDcols = avg_col]


xi_mean <- ggplot(dt_parameters_results)+
     geom_boxplot(mapping = aes(x = paste0("",as.character(quantile_threshold)), y = xi_mean))+
     geom_hline(aes(yintercept = dt_parameters_results[["true_xi_mean"]][1], colour = "Estimation (true u)"), linetype = "dashed", linewidth = 1.0) +
     geom_hline(aes(yintercept = 0.3, colour = "True parameter value"), linetype = "dotdash", linewidth = 1.0) +
     scale_colour_manual(name = "", values = c("Estimation (true u)" = "blue", "True parameter value" = "orange")) +
     # xlab("Selected Quantile Threshold")+
     xlab("")+
     ylab(expression(bar(hat(xi))))+
     guides(colour = guide_legend(override.aes = list(linewidth = 1.), keywidth = unit(2, "cm"))) +
     theme_bw(base_size = 20)+
     theme(legend.position = "none",
           axis.text.x = element_text(angle = 90))

sigma_mean <- ggplot(dt_parameters_results)+
     geom_boxplot(mapping = aes(x = paste0("",as.character(quantile_threshold)), y = (sigma_mean)))+
     geom_hline(aes(yintercept = (dt_parameters_results[["true_sigma_mean"]])[1], colour = "Estimation (true u)"), linetype = "dashed", linewidth = 1.0) +
     geom_hline(aes(yintercept = log(2.5), colour = "True parameter value"), linetype = "dotdash", linewidth = 1.0) +
     scale_colour_manual(name = "", values = c("Estimation (true u)" = "blue", "True parameter value" = "orange")) +
     xlab(expression(alpha^(T)))+
     ylab(expression(bar(log(hat(sigma)))))+
     guides(colour = guide_legend(override.aes = list(linewidth = 1.), keywidth = unit(2, "cm"),direction = 'horizontal')) +
     theme_bw(base_size = 20)+
     theme(legend.position = "right",
           axis.text.x = element_text(angle = 90))

# --- create versions WITH legend ---
xi_mean_l <- xi_mean + theme(legend.position = "bottom",legend.direction = "horizontal")
sigma_mean_l <- sigma_mean + theme(legend.position = "bottom",legend.direction = "horizontal")

# extract the legend
leg <- cowplot::get_legend(sigma_mean)

xi_mean_nl <- xi_mean + theme(
     legend.position = "none",
     plot.margin = margin(t = 5, r = 5, b = 0, l = 5)   # <- no bottom margin
)

sigma_mean_nl <- sigma_mean + theme(
     legend.position = "none",
     plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
)

# combine
mean_plot <- cowplot::plot_grid(
     cowplot::plot_grid(xi_mean_nl, sigma_mean_nl, ncol = 2),
     leg,
     ncol = 1,
     rel_heights = c(1, 0.10)
)

mean_plot

ggsave(paste0("C:/Users/mm538r/OneDrive - University of Glasgow/p4r_enhacements/worst-case-scenario/p4r_techometrics_paper/rebuttal/figures/appendix_threshold_",(1-phi),"_mean_parameters.pdf"), mean_plot, device = cairo_pdf,
       width = 12*1.1, height = 5*1.4)

### Getting the sd
xi_sd <- ggplot(dt_parameters_results)+
     geom_boxplot(mapping = aes(x = paste0("",as.character(quantile_threshold)), y = xi_sd))+
     geom_hline(aes(yintercept = dt_parameters_results[["true_xi_sd"]][1], colour = "Estimation (true u)"), linetype = "dashed", linewidth = 1.0) +
     scale_colour_manual(name = "", values = c("Estimation (true u)" = "blue", "True parameter value" = "orange")) +
     # xlab("Selected Quantile Threshold")+
     xlab("")+
     ylab(expression(sd(hat(xi))))+
     guides(colour = guide_legend(override.aes = list(linewidth = 1.), keywidth = unit(2, "cm"))) +
     theme_bw(base_size = 20)+
     theme(legend.position = "none",
           axis.text.x = element_text(angle = 90))

sigma_sd <- ggplot(dt_parameters_results)+
     geom_boxplot(mapping = aes(x = paste0("",as.character(quantile_threshold)), y = (sigma_sd)))+
     geom_hline(aes(yintercept = (dt_parameters_results[["true_sigma_sd"]])[1], colour = "Estimation (true u)"), linetype = "dashed", linewidth = 1.0) +
     scale_colour_manual(name = "", values = c("Estimation (true u)" = "blue", "True parameter value" = "orange")) +
     xlab(expression(alpha^(T)))+
     ylab(expression(sd(log(hat(sigma)))))+
     guides(colour = guide_legend(override.aes = list(linewidth = 1.), keywidth = unit(2, "cm"),
                                  direction = 'horizontal')) +
     theme_bw(base_size = 20)+
     theme(legend.position = "right",
           axis.text.x = element_text(angle = 90))


# --- create versions WITH legend ---
xi_sd_l <- xi_sd + theme(legend.position = "bottom")
sigma_sd_l <- sigma_sd + theme(legend.position = "bottom")

# extract the legend
leg <- cowplot::get_legend(sigma_sd)

xi_sd_nl <- xi_sd + theme(
     legend.position = "none",
     plot.margin = margin(t = 5, r = 5, b = 0, l = 5)   # <- no bottom margin
)

sigma_sd_nl <- sigma_sd + theme(
     legend.position = "none",
     plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
)

# combine
sd_plot <- cowplot::plot_grid(
     cowplot::plot_grid(xi_sd_nl, sigma_sd_nl, ncol = 2),
     leg,
     ncol = 1,
     rel_heights = c(1, 0.10)
)

sd_plot

ggsave(paste0("C:/Users/mm538r/OneDrive - University of Glasgow/p4r_enhacements/worst-case-scenario/p4r_techometrics_paper/rebuttal/figures/appendix_threshold_",(1-phi),"_sd_parameters.pdf"), sd_plot, device = cairo_pdf,
       width = 12*1.1, height = 5*1.4)

sample_size_plot <- ggplot(dt_parameters_results)+
     geom_boxplot(mapping = aes(x = paste0("",as.character(quantile_threshold)), y = n_tail),
                  linewidth = 0.3,outlier.size=0.5)+
     geom_hline(yintercept = dt_parameters_results[["true_n_tail"]][1], col = "blue", linetype = 'dashed', lwd = 1.0)+
     xlab(expression(alpha^(T)))+
     ylab(expression(n))+
     theme_bw(base_size = 20)

ggsave(paste0("C:/Users/mm538r/OneDrive - University of Glasgow/p4r_enhacements/worst-case-scenario/p4r_techometrics_paper/rebuttal/figures/appendix_threshold_",(1-phi),"_sample_size.pdf"), sample_size_plot, device = cairo_pdf,
       width = 7, height = 4)


