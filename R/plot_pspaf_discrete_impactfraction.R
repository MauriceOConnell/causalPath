plot_pspaf_discrete_impactfraction <- function(response_model= response_model, mediator_models = mediator_models, riskfactor = "phys",refval = 0,data = newd, prev = 0.0035/0.9965, ci=TRUE,boot_rep=100,ci_level=0.95,ci_type=c("norm"), PS_impactFraction = TRUE, percent = c(seq(from=0.05, to=1, by=0.05)), method = c("predict","observed"), response_name ="case"){

      PSIF_store_predict <- data.frame()
      PSIF_store_predict <- data.frame(matrix(nrow = length(percent), ncol = 3))

      PSIF_store_observed <- data.frame()
      PSIF_store_observed <- data.frame(matrix(nrow = length(percent), ncol = 3))

      row_i <- 1
      for(j in percent ){

            if( "predict" %in% c("predict","observed")){

                  store_predict <- ps_paf_impactfraction(response_model= response_model, mediator_models = mediator_models, riskfactor = "phys",refval = 0,data = newd, prev = 0.0035/0.9965, ci=TRUE,boot_rep=100,ci_level=0.95,ci_type=c("norm"), PS_impactFraction = TRUE, percent = j, method = "predict", response_name ="case" )

                  PSIF_store_predict[row_i,1] <- store_predict$bias_corrected_estimate[2]
                  PSIF_store_predict[row_i,2] <- store_predict$norm_lower[2]
                  PSIF_store_predict[row_i,3] <- store_predict$norm_upper[2]
            }

            if( "observed" %in% c("predict","observed")){

                  store_observed <- ps_paf_impactfraction(response_model= response_model, mediator_models = mediator_models, riskfactor = "phys",refval = 0,data = newd, prev = 0.0035/0.9965, ci=TRUE,boot_rep=100,ci_level=0.95,ci_type=c("norm"), PS_impactFraction = TRUE, percent = j, method = "observed", response_name ="case" )

                  PSIF_store_observed[row_i,1] <- store_observed$bias_corrected_estimate[2]
                  PSIF_store_observed[row_i,2] <- store_observed$norm_lower[2]
                  PSIF_store_observed[row_i,3] <- store_observed$norm_upper[2]
            }

            row_i <- row_i + 1

      }

      library(ggplot2)

      # a<-data.frame(percent=rep(c("5","10","15","20","25","30","35","40","45","50","55","60","65","70","75","80","85","90","95","100"),2),
      a<-data.frame(percent=rep(percent*100,2),
               bias_corrected_estimate=c(PSIF_store_predict[,1],PSIF_store_observed[,1]),
               method=c(rep(c("predict"),length(percent)),rep(c(rep(c("observed"),length(percent))))),
               norm_lower=c(PSIF_store_predict[,2],PSIF_store_observed[,2]),
               norm_upper=c(PSIF_store_predict[,3],PSIF_store_observed[,3]))

      # quartz() is only for mac apple computers
      quartz()
      k <- ggplot(a,aes(x=percent,y=bias_corrected_estimate))

      k <- k + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

      k

}
