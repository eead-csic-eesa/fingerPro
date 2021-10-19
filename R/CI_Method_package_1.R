  # Updated by Ivan Lizaga 30/09/2021 to include the boolean Means = T/F
  # Modified by Ivan Lizaga to extrac the Median of the CI interval p25-p26 17/01/2020
  # Created by Borja Latorre and Ivan Lizaga 12/11/2019
  #'
  #' Loop of the CI method. It executes the CI method for all your mixtures in a loop
  #'
  #' The function quantify the predictions of each individual tracer in your dataset and calculate the CI value of each tracer.
  
  CI_Method <- function(data, points = 2000, Means = F, seed = 123456L)
  {
    set.seed(seed)
    
    sources <- nrow(inputSource(data))
    
    if (sources == 3) {
      system.time({
        
        # It does everything after saving the df in a list because it is prepared for run it in a loop of several df
        # Check when using lists instead of individual datasets
        data_l <- list()
        data_l[[1]] <- data
        
        results_list <- list()
        for (ix in 1:length(data_l)) {
          data <- data_l[[ix]]
          if (Means == T) {
          ######### Uncheck when working with mean and SD ##################################
          sources1 <- data[-nrow(data),c(3:ncol(data))] 
          id_sources <- (1:(nrow(data)-1))
          id <- paste('S',id_sources, sep='')
          sources <- cbind(id, sources1)
          mixtures <- inputSample(data[,c(1:((ncol(data)/2)+1))])
          
          results <- list()
          for (i0 in 0:(((ncol(data)+1)/2) - 3)) {
            # results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, 2000, 123456)
            results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, points, seed)
            }
          } else {
          # ###########################################
          sources <- inputSource(data)
          mixtures <- inputSample(data)
          ###########################################
                  results <- list()
          for (i0 in 0:(ncol(data)-3)) {
            # results[[i0+1]] <- triangles_random_c(sources, mixtures, i0, 2000, 123456)
            results[[i0+1]] <- triangles_virtual_c(sources, mixtures, i0, points, seed)
            }
          # All tracers, X(e.g. 1000) puntos, escogiendo al azar/virtualmente el resto de trazadores
  
          }
          # tranforma todas las columnas en numeric data
          for (i1 in 1:length(results)) {
            results[[i1]][, 3] <- as.numeric(as.character(results[[i1]][, 3]))
            results[[i1]][, 4] <- as.numeric(as.character(results[[i1]][, 4]))
            results[[i1]][, 5] <- as.numeric(as.character(results[[i1]][, 5]))
          }
          results_list[[ix]] <- results
        }
        # realiza la ecuación de distancia al centro  para cada trazador
        ################################################################
        d <-list()
        d1 <-list()
        for (ix in 1:length(data_l)) {
          for (i2 in 1:length(results)) {
            results_list1 <- results_list[[ix]] 
            d1[[i2]] <- (results_list1[[i2]]$w.S1-1/3)*(results_list1[[i2]]$w.S1-1/3)+
              (results_list1[[i2]]$w.S2-1/3)*(results_list1[[i2]]$w.S2-1/3)+
              (results_list1[[i2]]$w.S3-1/3)*(results_list1[[i2]]$w.S3-1/3)#+
          }
          d[[ix]] <-  d1
        } 
        # Añade la columna de distancia a los resultados de las coordenadas  
        ################################################################
        results_d <- list()
        results_d1 <- list()
        for (ix in 1:length(data_l)) {
          for (i3 in 1:length(results)){ 
            results_list1 <- results_list[[ix]]
            results_d1[[i3]] <- cbind(results_list1[[i3]],d[[ix]][[i3]])
            names(results_d1[[i3]])[(ncol(results_list1[[i3]])+1)] <- "d"
          }
          results_d[[ix]] <- results_d1
        }
        # Extrae los puntos del percentil 25-26 y hace la mediana
        ################################################################
        x0 <- list()
        x1 <- list()
        x2 <- list()
        x3 <- list()
        x4 <- list()
        for (ix in 1:length(data_l)) {
          for (i4 in 1:length(results)) {
            x0[[i4]]<- results_d[[ix]][[i4]][order(results_d[[ix]][[i4]]$d),]
            x1[[i4]] <- as.integer(nrow(x0[[i4]])*0.25)
            x2[[i4]] <- as.integer(nrow(x0[[i4]])*0.26)
            x3[[i4]] <- x0[[i4]][c(x1[[i4]]:x2[[i4]]),]
          }
          x4[[ix]] <- x3
        }
        # Calcula el CI usando la mediana de los puntos del percentil 25-26
        ################################################################ 
        number_interval_25_26 <- sapply(x4[[1]][1], nrow)
        CI_T <- list() 
        CI_L <- list()
        CI_p25_p26 <- list()
        for (ix in 1:length(data_l)) {
          for (i5 in 1:length(results)) {
            for (i6 in 1:number_interval_25_26) { x4[[ix]][i6]
              CI = 0
              nc = 0
              w = x4[[ix]][[i5]][i6,]$w.S1
              if (w < 0) {
                nc = -w
              } else if (w > 1) {
                nc = w - 1
              }
              CI = CI + nc ^ 2
              
              nc = 0
              w = x4[[ix]][[i5]][i6,]$w.S2
              if (w < 0) {
                nc = -w
              } else if (w > 1) {
                nc = w - 1
              }
              CI = CI + nc ^ 2
              
              nc = 0
              w = x4[[ix]][[i5]][i6,]$w.S3
              if (w < 0) {
                nc = -w
              } else if (w > 1) {
                nc = w - 1
              }
              
              CI_p25_p26[[i6]] = CI + nc ^ 2 
            } 
            CI_T[[i5]] <- CI_p25_p26
          }
          CI_L[[ix]] <- CI_T
        }
        
        CI_def  <- list()
        CI_def1 <- list()
        CI_def2 <- list()
        CI_df   <- list()
        CI_med  <- list()
        for(i7 in 1:length(data_l)){
          for(i8 in 1:length(results)){
            for(i9 in 1:number_interval_25_26){
              CI_def[[i9]] = -sqrt(as.numeric(CI_L[[i7]][[i8]][i9]))
            }
            CI_def2[[i8]] = as.numeric(CI_def)
          }
          CI_def1[[i7]] <- CI_def2
          CI_df <- do.call("cbind",CI_def1[[i7]])
          #después de calculas el CI para el intervalo de percentiles se hace la mediana
          CI_med[[i7]] <- apply(CI_df,2,median)
        }   
        
        ################################################################
        CI_R2 <- list()
        CI_R <- list()
        for (ix in 1:length(data_l)) {
          CI_R0 <- as.data.frame(do.call("cbind", CI_med))
          if (Means == T) {
            names_CI_R <- names(data[3:((ncol(data)+1)/2)])
            names_CI_R <- as.data.frame(names_CI_R)
            CI_R2 <- cbind(names_CI_R, CI_R0)
          
          } else {
            names_CI_R <- names(data[3:ncol(data)])
            names_CI_R <- as.data.frame(names_CI_R)
            CI_R2 <- cbind(names_CI_R, CI_R0)
          }
        }
        
        results$CI <- CI_R2
        
        return(results)
      })
     
    } else if (sources == 4)	{
      system.time({
        data_l <- list()
        data_l[[1]] <- data
        
        results_list <- list()
        for (ix in 1:length(data_l)) {
          data <- data_l[[ix]]
          if (Means == T) {
            ######### When working with mean and SD ##################################
            sources <- inputSource(data)
            sources <- sources[,c(1:(ncol(sources)/2))]
            mixtures <- inputSample(data)
            mixtures <- mixtures[,c(1:(ncol(mixtures)/2))]
            
            # All tracers, X(e.g. 1000) puntos, escogiendo al azar/virtualmente el resto de trazadores
            results <- list()
            for (i0 in 0:((ncol(data)/2)-2)) {
              # results[[i0+1]] <- triangles_random_c(sources, mixtures, i0, 2000, 123456)
              results[[i0 + 1]] <-
                triangles_virtual_c(sources, mixtures, i0, points, seed)
            }
            
          } else if (Means == F) {
            ######### When working with raw data ##################################
            sources <- inputSource(data)
            mixtures <- inputSample(data)
  
            results <- list()
            for (i0 in 0:(ncol(data) - 3)) {
              # results[[i0+1]] <- triangles_random_c(sources, mixtures, i0, 2000, 123456)
              results[[i0 + 1]] <-
                triangles_virtual_c(sources, mixtures, i0, points, seed)
            }
          }
          
  
          # tranforma todas las columnas en numeric data
          for (i1 in 1:length(results)) {
            results[[i1]][, 3] <- as.numeric(as.character(results[[i1]][, 3]))
            results[[i1]][, 4] <- as.numeric(as.character(results[[i1]][, 4]))
            results[[i1]][, 5] <- as.numeric(as.character(results[[i1]][, 5]))
            #Uncheck when using 4 sources
            results[[i1]][, 6] <- as.numeric(as.character(results[[i1]][, 6]))
          }
          results_list[[ix]] <- results
        }
        # realiza la ecuación de distancia al centro  para cada trazador
        ################################################################
        d <- list()
        d1 <- list()
        for (ix in 1:length(data_l)) {
          for (i2 in 1:length(results)) {
            results_list1 <- results_list[[ix]]
            d1[[i2]] <-
              (results_list1[[i2]]$w.S1 - 1 / 4) * (results_list1[[i2]]$w.S1 - 1 / 4) +
              (results_list1[[i2]]$w.S2 - 1 / 4) * (results_list1[[i2]]$w.S2 - 1 / 4) +
              (results_list1[[i2]]$w.S3 - 1 / 4) * (results_list1[[i2]]$w.S3 - 1 / 4)#+
            #Uncheck when using 4 sources
            (results_list1[[i2]]$w.S4-1/4)*(results_list1[[i2]]$w.S4-1/4)
          }
          d[[ix]] <-  d1
        }
        # Añade la columna de distancia a los resultados de las coordenadas
        ################################################################
        results_d <- list()
        results_d1 <- list()
        for (ix in 1:length(data_l)) {
          for (i3 in 1:length(results)) {
            results_list1 <- results_list[[ix]]
            results_d1[[i3]] <- cbind(results_list1[[i3]], d[[ix]][[i3]])
            names(results_d1[[i3]])[(ncol(results_list1[[i3]]) + 1)] <-
              "d"
          }
          results_d[[ix]] <- results_d1
        }
        # Extrae los puntos del percentil 25-26 y hace la mediana
        ################################################################
        x0 <- list()
        x1 <- list()
        x2 <- list()
        x3 <- list()
        x4 <- list()
        for (ix in 1:length(data_l)) {
          for (i4 in 1:length(results)) {
            x0[[i4]] <- results_d[[ix]][[i4]][order(results_d[[ix]][[i4]]$d), ]
            x1[[i4]] <- as.integer(nrow(x0[[i4]]) * 0.25)
            x2[[i4]] <- as.integer(nrow(x0[[i4]]) * 0.26)
            x3[[i4]] <- x0[[i4]][c(x1[[i4]]:x2[[i4]]), ]
            # x3[[i4]]<- results_d[[ix]][[i4]][order(results_d[[ix]][[i4]]$d),]
            #  x3[[i4]] <- as.integer(nrow(x3[[i4]])*0.5)
            #   x4[[i4]] <- x0[[i4]][x3[[i4]],]
          }
          x4[[ix]] <- x3
        }
        # Calcula el CI usando la mediana de los puntos del percentil 25-26
        ################################################################
        number_interval_25_26 <- sapply(x4[[1]][1], nrow)
        CI_T <- list()
        CI_L <- list()
        CI_p25_p26 <- list()
        for (ix in 1:length(data_l)) {
          for (i5 in 1:length(results)) {
            for (i6 in 1:number_interval_25_26) {
              x4[[ix]][i6]
              CI = 0
              nc = 0
              w = x4[[ix]][[i5]][i6, ]$w.S1
              if (w < 0) {
                nc = -w
              } else if (w > 1) {
                nc = w - 1
              }
              CI = CI + nc ^ 2
              
              nc = 0
              w = x4[[ix]][[i5]][i6, ]$w.S2
              if (w < 0) {
                nc = -w
              } else if (w > 1) {
                nc = w - 1
              }
              CI = CI + nc ^ 2
              
              nc = 0
              w = x4[[ix]][[i5]][i6, ]$w.S3
              if (w < 0) {
                nc = -w
              } else if (w > 1) {
                nc = w - 1
              }
              CI = CI + nc ^ 2
              
              #Uncheck when using 4 sources
              nc = 0
              w = x4[[ix]][[i5]][i6,]$w.S4
              if (w < 0) {
                nc = -w
              } else if (w > 1) {
                nc = w - 1
              }
              
              CI_p25_p26[[i6]] = CI + nc ^ 2
            }
            CI_T[[i5]] <- CI_p25_p26
          }
          CI_L[[ix]] <- CI_T
        }
        
        CI_def  <- list()
        CI_def1 <- list()
        CI_def2 <- list()
        CI_df   <- list()
        CI_med  <- list()
        for (i7 in 1:length(data_l)) {
          for (i8 in 1:length(results)) {
            for (i9 in 1:number_interval_25_26) {
              CI_def[[i9]] = -sqrt(as.numeric(CI_L[[i7]][[i8]][i9]))
            }
            CI_def2[[i8]] = as.numeric(CI_def)
          }
          CI_def1[[i7]] <- CI_def2
          CI_df <- do.call("cbind", CI_def1[[i7]])
          #después de calculas el CI para el intervalo de percentiles se hace la mediana
          CI_med[[i7]] <- apply(CI_df, 2, median)
        }
        
        ################################################################
        CI_R2 <- list()
        CI_R <- list()
        for (ix in 1:length(data_l)) {
          CI_R0 <- as.data.frame(do.call("cbind", CI_med))
          
          if (Means == T) {
            names_CI_R <- names(data[3:round((ncol(data)/2))])
            }
          else if (Means == F) {
          names_CI_R <- names(data[3:ncol(data)])
            }
          
             names_CI_R <- as.data.frame(names_CI_R)
          CI_R2 <- cbind(names_CI_R, CI_R0)
        }
        results$CI <- CI_R2
        
        return(results)
      })
    } else {
      print("ERROR: not implemented, please, contact the code developers")
    }
  }
  
