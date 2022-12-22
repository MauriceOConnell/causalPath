plot_joint_PSPAF_impactfraction <- function(data=data, model_list = model_list, parent_list=parent_list, node_vec=node_vec, prev = 0.0035/0.9965, nsim=NULL, correct_order = 3,
                  vars = c("phys","subhtn","apob_apoa","whr"), exact = TRUE, response_model = response_model, mediator_models = mediator_models,
                  riskfactor = "phys", refval=0, calculation_method = "D", ci=TRUE,boot_rep=100, ci_type=c("norm"),ci_level=0.95, ci_level_ME=0.95,
                  PS_impactFraction = TRUE, percent = 0.3, method = "predict", response_name = "case",
                  plot = c("all_exclTotalPAF", "average_PSIFs_inclTotalPAF", "all_inclTotalPAF", "Average PAF phys", "Average PAF subhtn", "Average PAF apob_apoa", "Average PAF whr", "TotalPAF", "phys", "subhtn", "apob_apoa", "whr") ){


      dat_percent <- data.frame()
      dat_percent <- data.frame(matrix(nrow = length(percent)*( length(vars)*( length(vars) + 1) + 1) , ncol = 1))

      row_j = 1
      for(j in percent*100){

        dat_percent[(1 + (length(vars)*( length(vars) + 1) + 1)*(row_j-1) ):((length(vars)*( length(vars) + 1) + 1)*row_j),1] <- j

        row_j = row_j + 1
      }


      # length(vars)*( length(vars) + 1) + 1 is the number of rows in the output of 1 call of average_pspaf_impactfraction()
      # Then length(percent)*(length(vars)*( length(vars) + 1) + 1) is the total number of row for length(percent) calls of average_pspaf_impactfraction()
      store_predict <- data.frame()
      store_predict <- data.frame(matrix(nrow = length(percent)*(length(vars)*( length(vars) + 1) + 1) , ncol = 7))

      # store_observed <- data.frame()
      # store_observed <- data.frame(matrix(nrow = length(percent)*(length(vars)*( length(vars) + 1) + 1) , ncol = 7))

      colnames(store_predict) <- c("raw_estimate", "estimated_bias","bias_corrected_estimate","norm_lower", "norm_upper", "method","percent")
      # colnames(store_observed) <- c("raw_estimate", "estimated_bias","bias_corrected_estimate","norm_lower", "norm_upper", "method","percent")

      t_vector = c(paste0(rep(node_vec[node_vec %in% vars],times=rep(length(vars),length(vars))),'_',rep(1:length(vars),length(vars))),paste0("Average PAF ", node_vec[node_vec %in% vars]),'TotalPAF')
      row_i <- 1
      for(j in percent ){

            # if( "predict" %in% c("predict","observed")){

                  # store_predict <- ps_paf_impactfraction(response_model= response_model, mediator_models = mediator_models, riskfactor = "phys",refval = 0,data = newd, prev = 0.0035/0.9965, ci=TRUE,boot_rep=100,ci_level=0.95,ci_type=c("norm"), PS_impactFraction = TRUE, percent = j, method = "predict", response_name ="case" )

                  store_predict[  (1 + (length(vars)*( length(vars) + 1) + 1)*(row_i-1) ):((length(vars)*( length(vars) + 1) + 1)*row_i) , 1:5] <- average_pspaf_impactfraction(data=newd, model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev = 0.0035/0.9965, nsim=NULL, correct_order = 3,
                                                                                                                                                        vars = c("phys","subhtn","apob_apoa","whr"), exact = TRUE, response_model = response_model, mediator_models = mediator_models,
                                                                                                                                                        riskfactor = "phys", refval=0, calculation_method = "D", ci=TRUE,boot_rep=100, ci_type=c("norm"),ci_level=0.95, ci_level_ME=0.95,
                                                                                                                                                        PS_impactFraction = TRUE, percent =j, method = "predict", response_name = "case")

                  # t_vector should coincide with rownames
                  store_predict[(1 + (length(vars)*( length(vars) + 1) + 1)*(row_i-1) ):((length(vars)*( length(vars) + 1) + 1)*row_i),6] <- t_vector

                  store_predict[(1 + (length(vars)*( length(vars) + 1) + 1)*(row_i-1) ):((length(vars)*( length(vars) + 1) + 1)*row_i),7] <- j*100

            # }

            # if( "observed" %in% c("predict","observed")){
            #
            #       store_observed <- average_pspaf_impactfraction(data=newd, model_list=model_list, parent_list=parent_list, node_vec=node_vec, prev = 0.0035/0.9965, nsim=NULL, correct_order = 3,
            #                         vars = c("phys","subhtn","apob_apoa","whr"), exact = TRUE, response_model = response_model, mediator_models = mediator_models,
            #                         riskfactor = "phys", refval=0, calculation_method = "D", ci=TRUE,boot_rep=100, ci_type=c("norm"),ci_level=0.95, ci_level_ME=0.95,
            #                         PS_impactFraction = TRUE, percent =j, method = "observed", response_name = "case")
                  #
                  # # t_vector should coincide with rownames
                  # store_observed[(1 + (length(vars)*( length(vars) + 1) + 1)*(row_i-1) ):((length(vars)*( length(vars) + 1) + 1)*row_i),6] <- t_vector
            #
            # }

            row_i <- row_i + 1

      }

      library(ggplot2)
      ########
      #########
      # dat_percent <- data.frame()
      # dat_percent <- data.frame(matrix(nrow = length(percent)*( length(vars)*( length(vars) + 1) + 1) , ncol = 1))
      #
      # row_j = 1
      # for(j in percent*100){
      #
      #   dat_percent[(1 + (length(vars)*( length(vars) + 1) + 1)*(row_j-1) ):((length(vars)*( length(vars) + 1) + 1)*row_j),1] <- j
      #
      #   row_j = row_j + 1
      # }

      library(dplyr)

      # 1. "all_exclTotalPAF"
      # 2. "average_PSIFs_inclTotalPAF"
      # 3. "all_inclTotalPAF"
      # 4. "Average PAF phys"
      # 5. "Average PAF subhtn"
      # 6. "Average PAF apob_apoa"
      # 7. "Average PAF whr"
      # 8. "TotalPAF"
      # 9. "phys"
      # 10. "subhtn"
      # 11. "apob_apoa"
      # 12. "whr"

      ########
              ########
              ## Plot all
              ########
              ########
              # dat_all <-data.frame(percent=rep(percent*100,length(vars) + 1),
              #                          bias_corrected_estimate = c( store_predict[,3] ),
              #                          method=c( rep(c("predict_1"),length(percent)) ),
              #                          norm_lower = c(store_predict[,4] ),
              #                          norm_upper = c( store_predict[,5] ) )
              ## CHECK IF REMOVE c() does it keep row names and have row names will same name
              ## then might be able to assign method by row name
              ## if not need to create column for row names at each iteration of for loop with creating
              dat_all <- data.frame()
              dat_all <- data.frame(matrix(nrow = length(percent)*( length(vars)*( length(vars) + 1) + 1) , ncol = 5))

              dat_all <- data.frame(percent = dat_percent[,1],
                                    bias_corrected_estimate = c( store_predict[,3] ),
                                    # method=c( rep(c("predict_1"),length(percent)) ),
                                    method=c( store_predict[,6] ),
                                    norm_lower = c(store_predict[,4] ),
                                    norm_upper = c( store_predict[,5] ) )

      if( plot == "all_exclTotalPAF"){

            ###################
            ###################
            ## Just plotting average pathway-specific impact-fractions (PS-IF) for phys, whr, apob_apoa and subhtn and each of the PS-IF for each path. Not plotting the total PAF here.
            ###################
            ###################
            dat <- data.frame()
            dat <- data.frame(matrix(nrow = length(percent)*(length(vars) + 1)*length(vars), ncol = 5))
            colnames(dat) <- c("percent","bias_corrected_estimate", "method", "norm_lower", "norm_upper")

            dat_check <- data.frame(matrix(nrow = length(percent)*(length(vars) + 1)*length(vars), ncol = 6))
            colnames(dat_check) <- c("percent","bias_corrected_estimate", "method", "norm_lower", "norm_upper", "percent_check")

            for( k in 1:length(percent) ){

                  # dat_all_alternative <- filter(store_predict, method == "phys_1" | method ==  "phys_2" | method == "phys_3" | method == "phys_4" | method == "Average PAF phys" |
                  #                                              method == "subhtn_1" | method ==  "subhtn_2" | method == "subhtn_3" | method == "subhtn_4" | method == "Average PAF subhtn" |
                  #                                              method== "apob_apoa_1"| method== "apob_apoa_2" | method == "apob_apoa_3" | method == "apob_apoa_4" | method =="Average PAF apob_apoa" |
                  #                                              method ==  "whr_1" | method ==  "whr_2" | method == "whr_3" | method == "whr_4" | method =="Average PAF whr")

                  # filter(dat_all, method=="Average PAF phys")

                  dat[( (1 + (k-1)*length(percent)):(length(percent)*k) ),] <- data.frame(percent = dat_percent[ c( ( 1 + ((k - 1)*length(vars)*(length(vars) + 1) + 1*(k-1)) ):( length(vars)*(length(vars) + 1)*k + 1*(k-1)) ) ,1],
                                                                                 bias_corrected_estimate = c(store_predict[ c( ( 1 + ((k - 1)*length(vars)*(length(vars) + 1) + 1*(k-1)) ):( length(vars)*(length(vars) + 1)*k + 1*(k-1)) ),3] ),
                                                                                 method=c( store_predict[ c( ( 1 + ((k - 1)*length(vars)*(length(vars) + 1) + 1*(k-1)) ):( length(vars)*(length(vars) + 1)*k + 1*(k-1)) ),6] ),
                                                                                 norm_lower = c(store_predict[ c( ( 1 + ((k - 1)*length(vars)*(length(vars) + 1) + 1*(k-1)) ):( length(vars)*(length(vars) + 1)*k + 1*(k-1)) ),4] ),
                                                                                 norm_upper = c(store_predict[ c( ( 1 + ((k - 1)*length(vars)*(length(vars) + 1) + 1*(k-1)) ):( length(vars)*(length(vars) + 1)*k + 1*(k-1)) ),5] ) )

                  # dat_check[( (1 + (k-1)*length(percent)):(length(percent)*k) ),] <- data.frame(percent = dat_percent[ c( ( 1 + ((k - 1)*length(vars)*(length(vars) + 1) + 1*(k-1)) ):( length(vars)*(length(vars) + 1)*k + 1*(k-1)) ) ,1],
                  #                                                                bias_corrected_estimate = c(store_predict[ c( ( 1 + ((k - 1)*length(vars)*(length(vars) + 1) + 1*(k-1)) ):( length(vars)*(length(vars) + 1)*k + 1*(k-1)) ),3] ),
                  #                                                                method=c( store_predict[ c( ( 1 + ((k - 1)*length(vars)*(length(vars) + 1) + 1*(k-1)) ):( length(vars)*(length(vars) + 1)*k + 1*(k-1)) ),6] ),
                  #                                                                norm_lower = c(store_predict[ c( ( 1 + ((k - 1)*length(vars)*(length(vars) + 1) + 1*(k-1)) ):( length(vars)*(length(vars) + 1)*k + 1*(k-1)) ),4] ),
                  #                                                                norm_upper = c(store_predict[ c( ( 1 + ((k - 1)*length(vars)*(length(vars) + 1) + 1*(k-1)) ):( length(vars)*(length(vars) + 1)*k + 1*(k-1)) ),5] ),
                  #                                                                percent_check = c(store_predict[ c( ( 1 + ((k - 1)*length(vars)*(length(vars) + 1) + 1*(k-1)) ):( length(vars)*(length(vars) + 1)*k + 1*(k-1)) ),7] ) )

            }

            quartz()
                  s <- ggplot(dat,aes(x=percent,y=bias_corrected_estimate))

                  s <- s + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

                  s

           # quartz()
           #        s_alt <- ggplot(dat_all_alternative,aes(x=percent,y=bias_corrected_estimate))
           #
           #        s_alt <- s_alt + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))
           #
           #        s_alt

      } else if( plot == "average_PSIFs_inclTotalPAF"){

            ########
            ########
            ## Plot just Average PS-IF's and the Total PAF
            ########
            #########


            indices_plot_average <- data.frame()
            indices_plot_average <- data.frame(matrix(nrow = length(percent)*( length(vars) + 1) , ncol = 1))


            for(i in 1:length(percent) ){

                  indices_plot_average[ ( (1 + ( length(vars) + 1)*(i-1) ):( ( length(vars) + 1)*i ) ),1] <- c(  ( ((length(vars)^2) + 1)*i +length(vars)*(i-1) ) :( ( (length(vars)^2) + length(vars) + 1 )*i)  )

            }

            dat_average <- data.frame()
            dat_average <- data.frame(matrix(nrow = length(percent)*( length(vars) + 1), ncol = 5))
            colnames(dat_average) <- c("percent","bias_corrected_estimate", "method", "norm_lower", "norm_upper")

            ## NOTE: indices_plot_average[,1] converts the dataframe in a column vector
            dat_average <-data.frame(percent = dat_percent[indices_plot_average[,1],1],
                                     # bias_corrected_estimate = c( store_predict[c( ( (length(vars)^2) + 1):( (length(vars)^2) + length(vars) + 1 ) ),3] ),
                                     bias_corrected_estimate = c(store_predict[ indices_plot_average[,1],3] ) ,
                                     # method=c( rep(c("predict_1"),length(percent)) ),
                                     method= c( store_predict[indices_plot_average[,1],6] ),
                                     norm_lower = c(store_predict[indices_plot_average[,1],4]) ,
                                     norm_upper = c(store_predict[indices_plot_average[,1],5])  )

                  quartz()
                  t <- ggplot(dat_average,aes(x=percent,y=bias_corrected_estimate))

                  t <- t + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

                  t



                 # dat_average_check <- filter(store_predict, method == "Average PAF phys" | method == "Average PAF subhtn" | method == "Average PAF apob_apoa" | method == "Average PAF whr" | method == "TotalPAF" )
                 #
                 #  quartz()
                 #  t_check <- ggplot(dat_average_check,aes(x=percent,y=bias_corrected_estimate))
                 #
                 #  t_check <- t_check + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))
                 #
                 #  t_check

      } else if( plot == "all_inclTotalPAF" ){

              quartz()
              w <- ggplot(dat_all,aes(x=percent,y=bias_corrected_estimate))

              w <- w + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

              w

      } else if( plot == "Average PAF phys"){

            ########
            dat_Average_PAF_phys <- filter(dat_all, method=="Average PAF phys")

            quartz()
            plot_Average_PAF_phys <- ggplot(dat_Average_PAF_phys,aes(x=percent,y=bias_corrected_estimate))

            plot_Average_PAF_phys <- plot_Average_PAF_phys + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

           plot_Average_PAF_phys


            ########

      } else if( plot == "Average PAF subhtn"){

            ########

            dat_Average_PAF_subhtn <- filter(dat_all, method=="Average PAF subhtn")

            quartz()
            plot_Average_PAF_subhtn <- ggplot(dat_Average_PAF_subhtn,aes(x=percent,y=bias_corrected_estimate))

            plot_Average_PAF_subhtn <- plot_Average_PAF_subhtn + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

           plot_Average_PAF_subhtn

            ########

      } else if( plot == "Average PAF apob_apoa"){

            ########


            dat_Average_PAF_apob_apoa <- filter(dat_all, method=="Average PAF apob_apoa")

            quartz()
            plot_Average_PAF_apob_apoa <- ggplot(dat_Average_PAF_apob_apoa,aes(x=percent,y=bias_corrected_estimate))

            plot_Average_PAF_apob_apoa <- plot_Average_PAF_apob_apoa + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

           plot_Average_PAF_apob_apoa


            ########

      } else if(plot == "Average PAF whr"){

            ########

            dat_Average_PAF_whr <- filter(dat_all, method=="Average PAF whr")

            quartz()
            plot_Average_PAF_whr <- ggplot(dat_Average_PAF_whr,aes(x=percent,y=bias_corrected_estimate))

            plot_Average_PAF_whr <- plot_Average_PAF_whr + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

           plot_Average_PAF_whr


            ########

      } else if( plot == "TotalPAF" ){
            ########

            dat_TotalPAF <- filter(dat_all, method=="TotalPAF")

            quartz()
            plot_TotalPAF <- ggplot(dat_TotalPAF,aes(x=percent,y=bias_corrected_estimate))

            plot_TotalPAF <- plot_TotalPAF + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

           plot_TotalPAF

            ########

      } else if( plot == "phys"){

           ########


           dat_phys <- filter(dat_all, method=="phys_1" | method=="phys_2" | method=="phys_3" | method=="phys_4" | method=="Average PAF phys")

           # dim(dat_phys)

            quartz()
            plot_dat_phys <- ggplot(dat_phys,aes(x=percent,y=bias_corrected_estimate))

            plot_dat_phys <- plot_dat_phys + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

           plot_dat_phys

           ###########

      } else if( plot == "subhtn"){

           ########


           dat_subhtn <- filter(dat_all, method=="subhtn_1" | method=="subhtn_2" | method=="subhtn_3" | method=="subhtn_4" | method=="Average PAF subhtn")

           # dim(dat_subhtn)

            quartz()
            plot_dat_subhtn <- ggplot(dat_subhtn,aes(x=percent,y=bias_corrected_estimate))

            plot_dat_subhtn <- plot_dat_subhtn + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

           plot_dat_subhtn

           ###########

      } else if( plot == "apob_apoa"){

            ########


           dat_apob_apoa <- filter(dat_all, method=="apob_apoa_1" | method=="apob_apoa_2" | method=="apob_apoa_3" | method=="apob_apoa_4" | method=="Average PAF apob_apoa")

           # dim(dat_apob_apoa)

            quartz()
            plot_dat_apob_apoa <- ggplot(dat_apob_apoa,aes(x=percent,y=bias_corrected_estimate))

            plot_dat_apob_apoa <- plot_dat_apob_apoa + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

           plot_dat_apob_apoa

           ###############

      } else if( plot == "whr" ){

            ########

            dat_whr <- filter(dat_all, method=="whr_1" | method=="whr_2" | method=="whr_3" | method=="whr_4" | method=="Average PAF whr")

           # dim(dat_whr)

            quartz()
            plot_dat_whr <- ggplot(dat_whr,aes(x=percent,y=bias_corrected_estimate))

            plot_dat_whr <- plot_dat_whr + geom_line(aes(group=method,color=method)) + geom_errorbar(aes(ymin=norm_lower,ymax=norm_upper,color=method,width=0.2))

           plot_dat_whr

           ###########
           ########
      } else{
              stop("Must assign value to 'plot' for type of plot required.")
      }


}
