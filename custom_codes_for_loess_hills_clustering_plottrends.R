# trial calling home brewed r functions



# heatmap and clustering ------
cust_heatmap<- function(data_file, plot_hm = FALSE, plot_dm = FALSE,
                        calc_dist_matrix = TRUE, dist_method = "pearson",
                        linkage_method = "complete", scale_r = "row", if_col_cluster = FALSE,
                        dist_corplot = FALSE, tree_height = 50, view_genename = FALSE, 
                        cp_tree_height = 50, cp_view_genes = FALSE){
  
  library(factoextra)
  library(pheatmap)
  library(dendextend)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  
  if(calc_dist_matrix == FALSE){
    for(dist in dist_method){
      
      if(dist == "pearson"){
        dist <- "correlation"
      }
      try(pheatmap(data_file, scale = scale_r, clustering_distance_rows = dist,
                   clustering_method = linkage_method, cluster_cols = if_col_cluster,
                   treeheight_row = tree_height, show_rownames = view_genename,
                   main = paste(dist,linkage_method, sep = "_")))
    }
  }
  
  if(calc_dist_matrix == TRUE){
    
    
    
    dist_matrix <- vector("list")
    dend_object <- vector("list")
    hclust_object <- vector("list")
    
    
    for(dist in dist_method){
      
      
      if(dist == "euclidean"){
        data_scale <- t(scale(t(data_file), scale = TRUE, center = TRUE))
        dist_row <- get_dist(data_scale, dist)
        
      } else {dist_row <- get_dist(data_file, dist)}
      
      dist_matrix[[paste("dist", dist, sep = "_")]] <- dist_row
      
      cluster_heir <- hclust(dist_row, method = linkage_method)
      cluster_dendo <- as.dendrogram(cluster_heir)
      
      dend_object[[paste("h_dend", dist, linkage_method, sep = "_")]] <- cluster_dendo
      hclust_object[[paste("hclust", dist, linkage_method, sep = "_")]] <- cluster_heir
      print("clustering is done")
      
      if(plot_hm == TRUE){
        pheatmap(data_file, scale = scale_r, clustering_distance_rows = dist_row,
                 clustering_method = linkage_method, cluster_cols = if_col_cluster,
                 treeheight_row = tree_height, show_rownames = view_genename,
                 main = paste("Manual Dist", dist, linkage_method, sep = "_"))
        print("HM done")
      }
      
      
      
      
      if(plot_dm == TRUE){
        
        plot(cluster_dendo, leaflab = "perpendicular",
             main = paste("Manual Dist", dist, linkage_method, sep = "_"))
        print("Dm done")
        
      }
      
      
      if(dist_corplot == TRUE){
        d_mat <- as.matrix(dist_row)
        pheatmap(d_mat, scale = "none", cluster_rows = TRUE, cluster_cols = TRUE,
                 fontsize_row = 5, fontsize_col = 3,
                 clustering_distance_rows = dist_row, clustering_distance_cols = dist_row,
                 clustering_method = linkage_method, treeheight_row = cp_tree_height,
                 treeheight_col = cp_tree_height, show_rownames = cp_view_genes,
                 show_colnames = cp_view_genes)
        
        print("corplot_done")
        
      }
      
      
    }
    detail_list <- list("dist_dt" = dist_matrix,
                        "dend_dt" =  dend_object, "hclust_dt" = hclust_object)
    
    
    
    return(detail_list)
    
    
    
  }
  
}


# cutree and boxplot ---- 
cust_cutree <- function(object_dend, data_file, dend_view = FALSE,
                        dend_col = FALSE, cut_k = NULL, cut_h = NULL, dend_cut = FALSE,
                        text_size = 0.7, text_orient = "perpendicular", 
                        plot_boxplot = FALSE, y_axis_lim = c(0,9)){
  library(dendextend)
  library(reshape2)
  library(ggplot2)
  library(RColorBrewer)
  par(cex = text_size)
  clrs <- brewer.pal(n = cut_k, "Paired")
  
  
  if(dend_view == TRUE){
    plot(object_dend, leaflab = text_orient,
         main = paste(deparse(substitute(object_dend)), "w_cutree", sep = "_"))
  }
  
  if(dend_col == TRUE){
    col_dend <- color_branches(object_dend, k = cut_k, h = cut_h,
                               col = clrs)
    plot(col_dend, leaflab = text_orient, main = deparse(substitute(object_dend))) }
  
  if(dend_cut == TRUE){
    cut_dend <- cutree(object_dend, k = cut_k, h=cut_h, order_clusters_as_data = TRUE)
    cluster_dt <- cbind.data.frame(cut_dend, data_file)
    cluster_split <- split(cluster_dt, cluster_dt$cut_dend)
    detail_list <- list("cluster_dt_comb" = cluster_dt, "cluster_list" = cluster_split)
    
    
    
    clust_melt <- melt(cluster_dt, id = "cut_dend")
    clust_melt$cut_dend <- as.factor(clust_melt$cut_dend)
    clustwise_boxplot <- ggplot(data = clust_melt, aes(x = variable, y = value))+
      geom_boxplot(aes(fill = cut_dend), position = position_dodge(0.9), outlier.shape = NA)+
      theme(axis.text.x = element_text(angle = 90)) + coord_cartesian(ylim = y_axis_lim)+
      scale_fill_manual(values = clrs)
    
    detail_list[["cluster_boxplot"]] <- clustwise_boxplot 
    
    
    if(plot_boxplot == TRUE){
      
      print(clustwise_boxplot)
      # ggplot(data = clust_melt, aes(x = variable, y = value))+
      #   geom_boxplot(aes(fill = cut_dend), position = position_dodge(0.9), outlier.shape = NA)+
      #   theme(axis.text.x = element_text(angle = 90)) +scale_y_continuous(limits = y_axis_lim)
      
      print("hello there")
    }
    
    return(detail_list)
  } 
  
  
}

# plot list trendline and individual------

plot_list <- function(dt_list, list_num = "all", x_lim = c(0,4), y_lim = c(0,20), x_vector,
                      pallet = "Dark2", trendline = TRUE, plot_ind = TRUE,
                      plot_labels = list(tittle = "Plot", x_lab = "[cAMP]", y_lab = "FC")){
  library(RColorBrewer)
  library(pracma)
  clrs <- brewer.pal(n = length(dt_list), name = pallet)
  trend_clr <- brewer.pal(n = length(dt_list), name = "Paired")
  plot(NULL, xlim = x_lim, ylim = y_lim,
       xlab = plot_labels$x_lab, ylab = plot_labels$y_lab,
       main = plot_labels$tittle)
  
  print("NULL plot done")
  
  ifelse(list_num == "all", 
         clus_num <- 1:length(dt_list),
         clus_num <- list_num)
  print("clus_num set")
  
  if(plot_ind == TRUE){
    
    for(num in clus_num){
      clus_dt <- dt_list[[num]]
      for(row_i in 1:nrow(clus_dt)){
        lines(x_vector, clus_dt[row_i,-1], col = clrs[num], type = "o")
        y_txt <- clus_dt[row_i, (ncol(clus_dt))]
        x_txt <- x_vector[length(x_vector)]
        text(x = (x_txt+1), y = y_txt, rownames(clus_dt)[row_i], cex = 1.5, pos = 1)}
      legend("topleft", as.character(clus_num), col = clrs[clus_num], lwd = 2, cex = 1)}
    
    
  }
  if(trendline == TRUE){
    for(num in clus_num){
      clus_dt <- dt_list[[num]]
      lines(x_vector, colMeans(clus_dt[ ,-1]), col = trend_clr[num], type = "o", lwd = 5)}
    legend("topleft", as.character(clus_num), col = trend_clr[clus_num], lwd = 2, cex = 1)}
  
}


# plot loess trends to a particular folder  --- 
plot_trends <- function(data_list, filepath, folder_name, x_vector){
  k = length(data_list)
  foldername <- paste(rep(folder_name, k), c(1:k), sep = "_")
  setwd(filepath)
  sapply(1:k, function(x) dir.create(foldername[x]))
  
  file_lists <- list.files(filepath)
  file_lists <- file_lists[-(grep("png", file_lists))]
  print(file_lists)
  filename <- sapply(1:k, function(x) paste(filepath, foldername[x],"/", sep = ""))
  
  for(i in 1:k){
    #i = 1
    dt <- data_list[[i]][ ,-1]
    
    #head(dt)
    for(j in 1:nrow(dt)){
      #j  = 1
      x <- x_vector
      y <- as.numeric(dt[j, ])
      
      geneid <- rownames(dt)[j]
      plot_name <- paste(filename[i], geneid, "_", j, ".png", sep = "")
      plot_name
      png(plot_name, width = 470, height = 404, units = "px", pointsize = 15)
      
      plot(x,y,type = "l", xlab = "[cAMP]", ylab = "Fold Change", lty = 3, lwd = 2)
      
      fit_loess <- loess(y ~ x, span = 0.7)
      y_pred <- predict(fit_loess, x)
      
      lines(x, y_pred, lwd = 2)
      
      dev.off()
      
      
    }
    
  }
  
}


#hills fitting function and also save the fitted files to particular folder ---- 
hills_fitting <- function(data_file, x_vector, plot_ind = F, foldername){
  
  param_list <- NULL
  y_pred_table <- NULL
  stats_table <- NULL
  rse_fit_table <- NULL
  
  for(i in 1:nrow(data_file)){
    print(i)
    
    dt <- data_file
    x <- x_vector
    y <- as.numeric(dt[i , ])
    
    df <- as.data.frame(cbind(x,y))
    
    functn <- function(b, x, vmax, k, n)b + ((vmax-b) * x^n / (k^n + x^n))
    
    result <- NULL
    
    try({
      
      result <- nls(y ~ functn(b, x, vmax, k, n), data = df, start = c(b = mean(y[1:3]), vmax = y[8], 
                                                                       k = 0.6, n=3), algorithm = "port")
      y_pred <- predict(result, newdata = list(df$x))
      
      y_pred_gn <- c(rownames(dt)[i], as.numeric(y_pred))
      y_pred_table <- rbind(y_pred_table, y_pred_gn)
      
      stats_lm <- summary(lm(y_pred ~ df$y))
      adjusted_r <- stats_lm$adj.r.squared
      r_pval <- stats_lm$coefficients[8]
      slope <- stats_lm$coefficients[2]
      
      stats_r <- cbind.data.frame(as.character(rownames(dt)[i]), as.numeric(adjusted_r), as.numeric(r_pval))
      stats_table <- rbind(stats_table, stats_r)
      
      
      params <- c(as.character(rownames(dt)[i]), as.numeric(summary(result)$coefficients[c(1:8, 13:16)]))
      param_list <- rbind(param_list, params)
      
      rse_fit <- c(as.character(rownames(dt)[i]), as.numeric(summary(result)$sigma))
      rse_fit_table <- rbind(rse_fit_table, rse_fit)
      
      
      if(plot_ind == TRUE){
        
        
        geneid <- rownames(dt)[i]
        plot_name <- paste(foldername, geneid, "_", i, ".png", sep = "")
        
        png(plot_name, width = 470, height = 404, units = "px", pointsize = 15)
        
        plot(df$x, df$y, xlab = "[cAMP]", ylab = "FC", type = "o", lty = 2)
        lines(df$x, y_pred, lty = 1, type = "o")
        legend("bottomright", c("raw", "fitted", adjusted_r), lty = c(2,1,0), lwd = 1.5, cex = 0.75)
        dev.off()
        print("plotted")
        
      }
    })
  }
  
  colnames(param_list) <- c("genename", "b", "vmax","k", "n", "se_b", "se_vmax", "se_k", "se_n",
                            "p_b", "p_vmax", "p_k", "p_n")
  
  param_list <- as.data.frame(param_list)
  param_list[ ,2:ncol(param_list)] <- sapply(2:ncol(param_list),
                                             function(x) as.numeric(as.character(param_list[ ,x])))
  rownames(param_list) <- param_list$genename
  
  
  relative_se <- param_list[ ,c(6:9)]/param_list[ ,c(2:5)]*100
  colnames(relative_se) <- paste("r", colnames(relative_se), sep = "")
  relative_se$genename <- rownames(relative_se)
  

  
  colnames(stats_table) <- c("genename", "adj_r2", "r2_pval")
  
  colnames(y_pred_table) <- c("genename", x)
  y_pred_table <- as.data.frame(y_pred_table)
  y_pred_table[ ,2:ncol(y_pred_table)] <- apply(y_pred_table[ ,2:ncol(y_pred_table)], 2, as.numeric)
  rownames(y_pred_table) <- y_pred_table$genename
  
  colnames(rse_fit_table) <- c("genename", "rse")
  rse_fit_table <- as.data.frame(rse_fit_table)
  rownames(rse_fit_table) <- rse_fit_table$genename
  
  list_params <- list(param_list, stats_table, rse_fit_table, relative_se)
  parameter_table <- Reduce(function(d1,d2) merge(d1,d2, by = "genename"), list_params)
  parameter_table <- as.data.frame(parameter_table)
  parameter_table[ ,2:ncol(parameter_table)] <- apply(parameter_table[ ,2:ncol(parameter_table)], 2, as.numeric)
  rownames(parameter_table) <- parameter_table$genename
  
  data_file$genename <- rownames(data_file)
  
  return_list <- list("y_pred_table" = y_pred_table, "predicted_parameters" = param_list,
                      "r_statistic" = stats_table, "rse_fit" = rse_fit_table,"parameter_table" = parameter_table,
                      "input_file" = data_file)
  
  return(return_list)
  
}

#---plotting both raw and hills fitted plots-----


plot_fitted_reads <- function(raw_dt, fitted_data, x_vector, foldername, r2_table){
  for(gene_i in 1:nrow(fitted_data)){
    
    geneid <- rownames(fitted_data)[gene_i]
    
    plot_name <- paste(foldername, geneid, "_", gene_i, ".png", sep = "")
    
    png(plot_name, width = 470, height = 404, units = "px", pointsize = 15)
    
    
    y_fitted <- fitted_data[gene_i, ]
    genename <- rownames(fitted_data)[gene_i]
    
    y_raw <- raw_dt[rownames(raw_dt) == genename, ]
    genename_raw <- rownames(y_raw)
    
    r2_dt <- r2_table[r2_table$genename == genename,2]
    genename_r2 <- r2_table[r2_table$genename == genename,1]
    
    plot(x_vector, y_raw, lty = 2, type = "o", xlab = "[cAMP]", ylab = "FC",
         main = paste("check", genename, genename_raw, genename_r2, sep = "_"))
    lines(x_vector, y_fitted, lty = 1, type = "o")
    
    legend("topleft", c("raw", "fitted", round(r2_dt, digits = 3)), lty = c(2,1,0), lwd = 1.5, cex = 0.75)
    
    dev.off()
    
    
  }
  
  
}
# fit loess ---- 
fit_loess <- function(dt_file, x_vector,
                      plot_ind_graphs = FALSE,
                      filepath = NULL, l_span = 0.7){
  predicted_val_tab <- NULL
  
  for(row_i in 1:nrow(dt_file)){
    
    x <- x_vector
    y <- as.numeric(dt_file[row_i, ])
    fit_loess <- loess(y~x, span = l_span)
    predict_fit <- predict(fit_loess, new_data = data.frame(x))
    predicted_val <- c(rownames(dt_file)[row_i], predict_fit)
    predicted_val_tab <- rbind(predicted_val_tab, predicted_val)
    
    
    
    
    if(plot_ind_graphs == TRUE){
      gene_name <- rownames(dt_file)[row_i]
      plot_name <- paste(filepath, gene_name, "_", row_i, ".png", sep = "")
      
      png(plot_name, width = 470, height = 404, units = "px", pointsize = 15)
      
      plot(x,y, xlab = "[cAMP]", ylab = "FC", type = "o", lty = 2)
      lines(x, predict_fit, type = "o", lty = 1)
      legend("bottomright", c("raw", "fitted", l_span), lty = c(2,1,0), lwd = 1.5, cex = 0.75)
      dev.off()
      print(paste("plotted",row_i, sep = ""))
      
      
    }
    
    
  }
  print("loess fitting done")
  
  colnames(predicted_val_tab) <- c("genename", x_vector)
  predicted_df <- as.data.frame(predicted_val_tab)
  predicted_df[ ,2:ncol(predicted_df)] <- apply(predicted_df[ ,2:ncol(predicted_df)], 2, as.numeric)
  rownames(predicted_df) <- predicted_df[ ,1]
  print("file manipulation done")
  
  return(predicted_df)
}



# add comment-----

add_comment <- function(dt_list, genename = FALSE){
  manipulated_dt <- vector("list")
  
  for(list_i in 1:length(dt_list)){
    dt_list[[list_i]]$cluster <- c(rep(list_i), nrow(dt_list))
    if(genename == TRUE){
      dt_list[[list_i]]$genename <- rownames(dt_list[[list_i]])}
    
  }
  
  dt_list_table <- do.call(rbind, dt_list)
  
  
  dt_list[["dt_table"]] <- dt_list_table
  
  return(dt_list)
  }

#fit LM ----

LM_fit_plot <- function(data, x_reads, filename, plot = FALSE){
  lm_r2 <- NULL
  lm_pval <- NULL
  lm_rse <- NULL
  lm_slope <- NULL
  genename <- NULL
  
  for(i in 1:nrow(data)){
    x <- x_reads
    y <- as.numeric(data[i, -1])
    gi <- as.character(data[i,1])
    
    fit_lm <- summary(lm(y ~ x))
    
    r2_lm <- fit_lm$adj.r.squared
    pval_lm <- fit_lm$coefficients[8]
    rse_lm <- fit_lm$sigma
    slope_lm <-fit_lm$coefficients[5]
    
    lm_r2 <- c(lm_r2, r2_lm)
    lm_pval <- c(lm_pval, pval_lm)
    genename <- c(genename, gi)
    lm_rse <- c(lm_rse, rse_lm)
    lm_slope <- c(lm_slope, slope_lm)
    
    if(plot == TRUE){
      plot_name <- paste(filename, gi, "_", i, ".png", sep = "")
      
      png(plot_name, width = 470, height = 404, units = "px", pointsize = 15)
      
      plot(x, y, xlab = "[cAMP]", ylab = "FC", type = "o", lty = 2)
      abline(lm(y~x), col = "red")
      legend("bottomright", c("raw", "fitted", round(r2_lm,3), pval_lm),
             col = c("black", "red"), lty = c(2,1,0,0), lwd = 1.5, cex = 0.75)
      dev.off()
      print(paste("plotted",i, sep = " "))
      
    }
    
  }
  LM_fit_values <- cbind.data.frame(genename, lm_r2, lm_rse, lm_pval, lm_slope)
}    
#----LM fit using nls-----
LM_fit_nls <- function(data_file, x_vector, plot_ind = F, foldername){
  
  param_list <- NULL
  y_pred_table <- NULL
  stats_table <- NULL
  rse_fit_table <- NULL
  
  for(i in 1:nrow(data_file)){
    
    dt <- data_file
    x <- x_vector
    y <- as.numeric(dt[i , ])
    
    df <- as.data.frame(cbind(x,y))
    
    
    functn <- function(x, m, c) m *x + c
    
    result <- NULL
    
    try({
      
      result <- nls(y ~ functn(x, m , c), data = df, start = c(m = 1, c = 1), algorithm = "port")
      
      y_pred <- predict(result, newdata = list(df$x))
      
      y_pred_gn <- c(rownames(dt)[i], as.numeric(y_pred))
      y_pred_table <- rbind(y_pred_table, y_pred_gn)
      
      stats_lm <- summary(lm(y_pred ~ df$y))
      adjusted_r <- stats_lm$adj.r.squared
      r_pval <- stats_lm$coefficients[8]
      
      
      stats_r <- cbind.data.frame(as.character(rownames(dt)[i]), as.numeric(adjusted_r), as.numeric(r_pval))
      stats_table <- rbind(stats_table, stats_r)
      
      
      params <- c(as.character(rownames(dt)[i]), as.numeric(summary(result)$coefficients[c(1:4, 7:8)]))
      param_list <- rbind(param_list, params)
      
      
      rse_fit <- c(as.character(rownames(dt)[i]), as.numeric(summary(result)$sigma))
      rse_fit_table <- rbind(rse_fit_table, rse_fit)
      rse_fit_table
      
      # if(plot_ind == TRUE){
      #   
      #   
      #   geneid <- rownames(dt)[i]
      #   plot_name <- paste(foldername, geneid, "_", i, ".png", sep = "")
      #   
      #   png(plot_name, width = 470, height = 404, units = "px", pointsize = 15)
      #   
      #   plot(df$x, df$y, xlab = "[cAMP]", ylab = "FC", type = "o", lty = 2)
      #   lines(df$x, y_pred, lty = 1, type = "o")
      #   legend("bottomright", c("raw", "fitted", adjusted_r), lty = c(2,1,0), lwd = 1.5, cex = 0.75)
      #   dev.off()
      #   print("plotted")
      #   
      # }
    })
  }
  
  colnames(param_list) <- c("genename", "m", "c", "se_m", "se_c", "p_m", "p_c")
  
  param_list <- as.data.frame(param_list)
  param_list[ ,2:ncol(param_list)] <- sapply(2:ncol(param_list),
                                             function(x) as.numeric(as.character(param_list[ ,x])))
  rownames(param_list) <- param_list$genename
  
  
  # relative_se <- param_list[ ,c(6:9)]/param_list[ ,c(2:5)]*100
  # colnames(relative_se) <- paste("r", colnames(relative_se), sep = "")
  # relative_se$genename <- rownames(relative_se)
  # 
  
  
  colnames(stats_table) <- c("genename", "adj_r2", "r2_pval")
  
  colnames(y_pred_table) <- c("genename", x)
  y_pred_table <- as.data.frame(y_pred_table)
  y_pred_table[ ,2:ncol(y_pred_table)] <- apply(y_pred_table[ ,2:ncol(y_pred_table)], 2, as.numeric)
  rownames(y_pred_table) <- y_pred_table$genename
  
  
  colnames(rse_fit_table) <- c("genename", "rse")
  rse_fit_table <- as.data.frame(rse_fit_table)
  rownames(rse_fit_table) <- rse_fit_table$genename
  rse_fit_table
  
  list_params <- list(param_list, stats_table, rse_fit_table)#, relative_se)
  parameter_table <- Reduce(function(d1,d2) merge(d1,d2, by = "genename"), list_params)
  parameter_table <- as.data.frame(parameter_table)
  parameter_table[ ,2:ncol(parameter_table)] <- apply(parameter_table[ ,2:ncol(parameter_table)], 2, as.numeric)
  rownames(parameter_table) <- parameter_table$genename
  parameter_table
  
  data_file$genename <- rownames(data_file)
  
  return_list <- list("y_pred_table" = y_pred_table, "predicted_parameters" = param_list,
                      "r_statistic" = stats_table, "rse_fit" = rse_fit_table,"parameter_table" = parameter_table,
                      "input_file" = data_file)
  
  
  return(return_list)
  
}







#---- Fishers exact test with multiple categories-------
#input file should be columns : DE vs nonDE and rows: FUNCAT 

custom_fisher_test <- function(count_table, print_mat = TRUE, print_summary = TRUE, alt_hyp = "greater" ){
  
  test_pval <- NULL
  test_oddsratio <- NULL
  funcat_name <- NULL
  
  
  for(i in 1:nrow(count_table)){
    
    a <- count_table[i,1]
    b <- count_table[i,2]
    c <- sum(count_table[ ,1]) - a
    d <- sum(count_table[ ,2]) - b
    
    test_mat <- matrix(c(a,c,b,d), nrow = 2, dimnames = list(FC = c("In", "Not In"),
                                                             DE = c("In", "Not In")))
    
    
    if(print_mat == TRUE){
      print(paste("Functional Cat", rownames(count_table)[i]))
      print(test_mat)
    }
    
    
    f_test <- fisher.test(test_mat, alternative = alt_hyp)
    
    if(print_summary == TRUE){
      print(paste("FC", i))
      print(f_test)
    }
    
    test_pval <- c(test_pval, f_test$p.value)
    test_oddsratio <- c(test_oddsratio, f_test$estimate)
    funcat_name <- c(funcat_name, rownames(count_table)[i])
    
    
  }
  
  summary_table <- cbind.data.frame("MF_ Cat" = funcat_name, "pval" = test_pval, "odds_ratio" = test_oddsratio)
  
  
}

#Slope calculator 1(uses k and k+d) and 2(uses k-d,k,k+d) : ---- 

slope_calculator_nls <- function(input_dt, x_vector,d){
  functn <- function(b, x, vmax, k, n)b + ((vmax-b) * x^n / (k^n + x^n))
  
  dt <- input_dt
  slope_dt <- NULL
  for(r_i in 1:nrow(input_dt)){
    gene_name <- rownames(input_dt)[r_i]
    print(c(r_i, gene_name))
    y <- as.numeric(dt[r_i, ])
    x <- x_vector
    df <- as.data.frame(cbind(x,y))
    result <- NULL
    result <- nls(y ~ functn(b, x, vmax, k, n), data = df, 
                  start = c(b = mean(y[1:3]), vmax = y[8],
                            k = 0.6, n=3), algorithm = "port")
    
    S <- summary(result)
    k <- S$coefficients[3]
    k_plus <- k + d
    new_camp <- c(k, k_plus)
    predicted_y <- predict(result, newdata = list(x = new_camp))
    delta_y <- predicted_y[2] - predicted_y[1]
    slope <- delta_y/d
    
    gene_name <- rownames(input_dt)[r_i]
    
    gene_info <- c(gene_name, k, slope)
    slope_dt <- rbind(slope_dt, gene_info)
    
  }
  slope_dt <- as.data.frame(slope_dt)
  colnames(slope_dt) <- c("genename", "k", "slope")
  rownames(slope_dt) <- slope_dt$genename
  slope_dt[ ,2:3] <- apply(slope_dt[ ,2:3], 2, as.numeric)
  
  return(slope_dt)
}

slope_calculator_nls2 <- function(input_dt, x_vector,d){
  functn <- function(b, x, vmax, k, n)b + ((vmax-b) * x^n / (k^n + x^n))
  
  dt <- input_dt
  slope_dt <- NULL
  for(r_i in 1:nrow(input_dt)){
    gene_name <- rownames(input_dt)[r_i]
    print(c(r_i, gene_name))
    y <- as.numeric(dt[r_i, ])
    x <- x_vector
    df <- as.data.frame(cbind(x,y))
    result <- NULL
    result <- nls(y ~ functn(b, x, vmax, k, n), data = df, 
                  start = c(b = mean(y[1:3]), vmax = y[8],
                            k = 0.6, n=3), algorithm = "port")
    
    S <- summary(result)
    k <- S$coefficients[3]
    k_plus <- k + d
    k_minus <- k-d
    new_camp <- c(k_minus, k, k_plus)
    predicted_y <- predict(result, newdata = list(x = new_camp))
    delta_y <- predicted_y[3] - predicted_y[1]
    delta_x <- k_plus - k_minus
    slope <- delta_y/delta_x
    
    gene_info <- c(gene_name, k, slope)
    slope_dt <- rbind(slope_dt, gene_info)
    
  }
  slope_dt <- as.data.frame(slope_dt)
  colnames(slope_dt) <- c("genename", "k", "slope")
  rownames(slope_dt) <- slope_dt$genename
  slope_dt[ ,2:3] <- apply(slope_dt[ ,2:3], 2, as.numeric)
  
  return(slope_dt)
}

# Descriptive stats for vectors ----
desc_stat <- function(input_data){
  
  
  dist_sum <- summary(input_data)
  median_i <- as.numeric(dist_sum[3])
  mean_i <- as.numeric(dist_sum[4])
  q1_i <- as.numeric(dist_sum[2])
  q3_i <- as.numeric(dist_sum[5])
  iqr_i <- iqr(input_data)
  cv_i <- (cv(input_data)) * 100
  sd_i <- sd(input_data)
  range_low <- as.numeric(dist_sum[1])
  range_high <- as.numeric(dist_sum[6])
  QOD_i <- iqr_i/median_i
  MAD_i <- mad(input_data, constant = 1)
  MADM_i <- MAD_i/median_i
  
  desc_stats <- cbind.data.frame(median_i, iqr_i,mean_i, sd_i,
                                 q1_i, q3_i, range_high, range_low,
                                 cv_i, QOD_i, MAD_i, MADM_i)
  
  return(desc_stats)
  
  
  
}


#FUNCTION FOR CREATING SIMULATED DATA FOR SIGMOID GENES ----------
create_sigmoid_table <- function(mean_k, mean_vmax, mean_n,
                                 x_val, sd_k, sd_v, sd_n, sd_b){
  hills <- function(b, k, vmax, n,x)b + ((vmax-b) * x^n / (k^n + x^n))
  k_rand <- rnorm(100, mean = mean_k, sd = sd_k)
  vmax_rand <- rnorm(100, mean = mean_vmax, sd = sd_v)
  n_rand <- rnorm(100, mean = mean_n, sd = sd_n)
  b0_rand <- rnorm(100, mean = 1, sd = sd_b)
  
  sig_par <- cbind.data.frame(b0_rand, k_rand, vmax_rand, n_rand) #table of simulated par
  sig_dt <- NULL #making simulated data from the randomly generated parameters
  for(r_i in 1:nrow(sig_par)){
    f <- as.numeric(sig_par[r_i, ])
    y_pred <- hills(f[1], f[2], f[3], f[4], xval)
    sig_dt <- rbind(sig_dt, y_pred)
  }
  
  sig_dt <- as.data.frame(sig_dt)
  colnames(sig_dt) <- xval
  rownames(sig_dt) <- 1:nrow(sig_dt)
  sig_dt$genenum <- 1:nrow(sig_dt)
  
  stats_sig <- NULL#calculating r2 to know if the clustering is happening in a good manner
  for(r_i in 1:nrow(sig_dt)){
    sum_stats <- summary(lm(as.numeric(sig_dt[r_i, -(ncol(sig_dt))]) ~ xval))
    r2 <- sum_stats$adj.r.squared
    p <- sum_stats$coefficients[8]
    r_name <- sig_dt[r_i, ncol(sig_dt)]
    r_stats <- c(r_name,r2,p)
    stats_sig <- rbind(stats_sig, r_stats)
    
  }
  stats_sig <- as.data.frame(stats_sig)
  colnames(stats_sig) <- c("genenum","r2", "pv")
  rownames(stats_sig) <- 1:nrow(stats_sig)
  
  sig_dt_f <- merge(sig_dt, stats_sig, by = "genenum", all.x = TRUE)
  sig_dt_f <- as.data.frame(sig_dt_f)
  
  return(sig_dt_f)
  
}


# FUNCTION FOR SIMULATED LM DATA ----------
create_lin_table <- function(mean_m, mean_c, x_val){
  lin_f <-  function(m,c,x)(m*x) + c
  m_rand <- rnorm(100, mean_m, 1)
  c_rand <- rnorm(100, mean_c, 0.01)
  
  lin_par <- cbind.data.frame(m_rand, c_rand)
  lin_dt <- NULL #making simulated data from the randomly generated parameters
  for(r_i in 1:nrow(lin_par)){
    f <- as.numeric(lin_par[r_i, ])
    y_pred <- lin_f(f[1], f[2], xval)
    lin_dt <- rbind(lin_dt, y_pred)
  }
  
  lin_dt <- as.data.frame(lin_dt)
  colnames(lin_dt) <- xval
  rownames(lin_dt) <- 1:nrow(lin_dt)
  lin_dt$genenum <- 1:nrow(lin_dt)
  
  return(lin_dt)
  
}

#modified index wrapping function: -----
mod_mod <- function(array, denominator){
  new_array <- (array) - 1
  looped_array <- new_array %% denominator
  new_looped_array <- (looped_array) + 1
  return(new_looped_array)
}


# binning readsperbase function (for both terminal and Rstudio) ------
binner <- function(n, input_file, filename = NULL){ #input file format: pruned sam file -- col1 - position, col2 - reads 
  bin_seq <- seq(1, nrow(input_file),n)
  bin_tab <- NULL
  for(i in 1:length(bin_seq)){
    bin_start <- bin_seq[i]
    bin_stop <- (bin_seq[i] + (n-1))
    bin_sum <- sum(input_file[mod_mod(bin_start:bin_stop, nrow(input_file)),2])
    bin_row <- c(bin_start, bin_stop, bin_sum)
    bin_tab <- rbind(bin_tab, bin_row)
  }
  print("done")
  bin_table <- as.data.frame(bin_tab)
  colnames(bin_table) <- c("bin_start", "bin_stop", "bin_reads")
  rownames(bin_table) <- 1:nrow(bin_table)
  if(!is.null(filename)){write.table(bin_table, file = paste(filename, "binned",
                                                             sep = "."), 
                                     sep = "\t", row.names =FALSE)}
  return(bin_table)
  
  
}


#Scaling rows ---------
scaler <- function(mat, if_center = T){
  scaled_dt <- t(scale(t(mat), center = if_center))
  return(scaled_dt)
}

#Pairwise_looper ------
pairwise_looper <- function(dt){ #dt should be a matrix with only numbers and rowsnames as genenames
  cor_dt <- NULL
  for(gene_i in 1:(nrow(dt)-1)){
    for(gene_j in (gene_i+1):nrow(dt)){
      cor_test <- cor.test(x = dt[gene_i, ], y = dt[gene_j, ], method = "pearson")
      gene_pair <- paste(rownames(dt)[gene_i], rownames(dt)[gene_j], sep = "_")
      cor_result <- c(gene_pair, cor_test$estimate, cor_test$p.value)
      cor_dt <- rbind(cor_dt, cor_result)
      
      
    }
  }
  cor_dt <- as.data.frame(cor_dt)
  colnames(cor_dt) <- c("gene_pair", "cor", "pvalue")
  return(cor_dt)
}

# kmeans -----
Optim_cluster <- function(raw_dt, scaling = TRUE){ #raw_dt has first column as genename
  library(ggplot2)
  library(factoextra)
  
  if(scaling == TRUE){
    clus_input_dt <- t(scale(t(raw_dt[ ,-1])))
  } else {clus_input_dt <- raw_dt}
  
  plot(fviz_nbclust(clus_input_dt, kmeans, method = "wss", nstart = 15))#4
  plot(fviz_nbclust(clus_input_dt, kmeans, method = "silhouette", nstart = 30))#3
  plot(fviz_nbclust(clus_input_dt, kmeans, method = "gap_stat",
                    nstart = 25, verbose = FALSE))#2
  print("Cluster optimisation plots done")
  
} 

cluster_kmeans <- function(raw_dt, x_vector, scaling = TRUE, kclus = 3){#raw dt is numerical matrix with genenames as 1st row
  
  if(scaling == TRUE){
    input_dt <- t(scale(t(raw_dt[ ,-1])))
    
  } else{input_dt <- raw_dt[ ,-1]}
  
  km <- kmeans(input_dt, centers = kclus, nstart = 100, iter.max = 900)
  print(table(km$cluster))
  
  clus_plot <- fviz_cluster(km, data = input_dt,
                            palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07","#FFB5C5", "#BF87B3", "#7F5AA2", "#3F2D91", "#000080"), 
                            ellipse.type = "euclid", # Concentration ellipse
                            star.plot = TRUE, # Add segments from centroids to items
                            #repel = TRUE, # Avoid label overplotting (slow)
                            
  ) + custom_theme
  
  plot(clus_plot)
  
  pty <- c(15:18,0:2)
  clrs <- c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", "#FFB5C5", "#BF87B3", "#7F5AA2", "#3F2D91", "#000080")
  plot(NULL, xlim = c(0, 5), ylim = c(-1, 2), xlab = "[cAMP]", ylab = "Centers")
  for(i in 1:kclus){
    lines(x_vector, km$centers[i, ], col = clrs[i], pch = pty[i])
  }
  legend("bottomright", legend, col = clrs[1:kclus], c(1:kclus), lwd = 2, cex = 0.75)
  
  cluster_list <- as.data.frame(km$cluster)
  colnames(cluster_list) <- "cluster"
  cluster_list$genename <- rownames(cluster_list)
  dt_clus <- merge(raw_dt, cluster_list, by = "genename")
  rownames(dt_clus) <- dt_clus$genename
  
  if(scaling == TRUE){
    output_list <- list("dt_clus" = dt_clus, "scaled_dt" = input_dt)
  }else{output_list <- dt_clus}
  
  return(output_list)
  
}

#remove outliers ----
outlier_finder <- function(input_vector, mother_matrix = NULL){#rownames as genes
  library(stats)
  library(magrittr)
  input_vector <- as.numeric(input_vector)
  inter_qr <- IQR(input_vector)
  quant <- quantile(input_vector, probs = c(0.25,0.75))
  
  upper_val <- as.numeric(quant[2]) + (1.5 * inter_qr)
  lower_val <- as.numeric(quant[1]) - (1.5 * inter_qr)
  
  outlier_up <- which(input_vector > upper_val)
  outlier_low <- which(input_vector < lower_val)
  
  outlier_rownums <- list("out_h" = outlier_up, "out_l" =  outlier_low)
  
  if(!is.null(mother_matrix)){
    mat_outlier_up <- cbind(rownames(mother_matrix[outlier_up, ]), 
                            rep("outlier_h", nrow(mother_matrix[outlier_up, ]))) %>% 
      `colnames<-`(c("genename", "outlier_status")) %>% as.data.frame()
    mat_outlier_low <- cbind(rownames(mother_matrix[outlier_low, ]), 
                             rep("outlier_l", nrow(mother_matrix[outlier_low, ]))) %>%
      `colnames<-`(c("genename", "outlier_status")) %>% as.data.frame()
    mat_outlier <- rbind.data.frame(mat_outlier_up, mat_outlier_low)
    
    outlier_dt_list <- list("outlier_mat" = mat_outlier ,
                            "out_h_nums" = outlier_up, "out_l_nums" = outlier_low)
    
    
    return(outlier_dt_list)
    
    
  } else {
    return(outlier_rownums)
    
  }
}

#finding mode ------
my_mode <- function(x){ #a vector
  density(x)$x[which.max(density(x)$y)]
}

#CPM hills fitting ----
hills_fitting_cpm <- function(data_file, x_vector, plot_ind = F, foldername){
  
  param_list <- NULL
  y_pred_table <- NULL
  stats_table <- NULL
  rse_fit_table <- NULL
  
  for(i in 1:nrow(data_file)){
    print(i)
    
    dt <- data_file
    x <- x_vector
    y <- as.numeric(dt[i , ])
    
    df <- as.data.frame(cbind(x,y))
    
    functn <- function(b, x, vmax, k, n)b + ((vmax-b) * x^n / (k^n + x^n))
    
    result <- NULL
    
    try({
      
      result <- nls(y ~ functn(b, x, vmax, k, n), data = df, start = c(b = mean(y[1:2]), vmax = mean(y[9:10]), 
                                                                       k = 0.8, n=3), algorithm = "port")
      y_pred <- predict(result, newdata = list(df$x))
      
      y_pred_gn <- c(rownames(dt)[i], as.numeric(y_pred))
      y_pred_table <- rbind(y_pred_table, y_pred_gn)
      
      stats_lm <- summary(lm(y_pred ~ df$y))
      adjusted_r <- stats_lm$adj.r.squared
      r_pval <- stats_lm$coefficients[8]
      
      stats_r <- cbind.data.frame(as.character(rownames(dt)[i]), as.numeric(adjusted_r), as.numeric(r_pval))
      stats_table <- rbind(stats_table, stats_r)
      
      
      params <- c(as.character(rownames(dt)[i]), as.numeric(summary(result)$coefficients[c(1:8, 13:16)]))
      param_list <- rbind(param_list, params)
      
      rse_fit <- c(as.character(rownames(dt)[i]), as.numeric(summary(result)$sigma))
      rse_fit_table <- rbind(rse_fit_table, rse_fit)
      
      
      if(plot_ind == TRUE){
        
        
        geneid <- rownames(dt)[i]
        plot_name <- paste(foldername, geneid, "_", i, ".png", sep = "")
        
        png(plot_name, width = 470, height = 404, units = "px", pointsize = 15)
        
        plot(df$x, df$y, xlab = "[cAMP]", ylab = "FC", type = "o", lty = 2)
        lines(df$x, y_pred, lty = 1, type = "o")
        legend("bottomright", c("raw", "fitted", adjusted_r), lty = c(2,1,0), lwd = 1.5, cex = 0.75)
        dev.off()
        print("plotted")
        
      }
    })
  }
  
  colnames(param_list) <- c("genename", "b", "vmax","k", "n", "se_b", "se_vmax", "se_k", "se_n",
                            "p_b", "p_vmax", "p_k", "p_n")
  
  param_list <- as.data.frame(param_list)
  param_list[ ,2:ncol(param_list)] <- sapply(2:ncol(param_list),
                                             function(x) as.numeric(as.character(param_list[ ,x])))
  rownames(param_list) <- param_list$genename
  
  
  relative_se <- param_list[ ,c(6:9)]/param_list[ ,c(2:5)]*100
  colnames(relative_se) <- paste("r", colnames(relative_se), sep = "")
  relative_se$genename <- rownames(relative_se)
  
  
  
  colnames(stats_table) <- c("genename", "adj_r2", "r2_pval")
  
  colnames(y_pred_table) <- c("genename", x)
  y_pred_table <- as.data.frame(y_pred_table)
  y_pred_table[ ,2:ncol(y_pred_table)] <- apply(y_pred_table[ ,2:ncol(y_pred_table)], 2, as.numeric)
  rownames(y_pred_table) <- y_pred_table$genename
  
  colnames(rse_fit_table) <- c("genename", "rse")
  rse_fit_table <- as.data.frame(rse_fit_table)
  rownames(rse_fit_table) <- rse_fit_table$genename
  
  list_params <- list(param_list, stats_table, rse_fit_table, relative_se)
  parameter_table <- Reduce(function(d1,d2) merge(d1,d2, by = "genename"), list_params)
  parameter_table <- as.data.frame(parameter_table)
  parameter_table[ ,2:ncol(parameter_table)] <- apply(parameter_table[ ,2:ncol(parameter_table)], 2, as.numeric)
  rownames(parameter_table) <- parameter_table$genename
  
  data_file$genename <- rownames(data_file)
  
  return_list <- list("y_pred_table" = y_pred_table, "predicted_parameters" = param_list,
                      "r_statistic" = stats_table, "rse_fit" = rse_fit_table,"parameter_table" = parameter_table,
                      "input_file" = data_file)
  
  return(return_list)
  
}

#chip signal calculator -----
chip_signal_calculator <- function(sam, feat, strt_col, stop_col, smooth_type = "none", average){
  #extracting the sequence segment with signal frequency
  seg_i <- sam[(which(sam$position == feat[[strt_col]])) : (which(sam$position == feat[[stop_col]])), ] 
  if(any(is.na(seg_i$frequency))){seg_i <- seg_i[-(which(is.na(seg_i$frequency))), ]}
  if(any(seg_i$frequency == "Inf")) {seg_i <- seg_i[-(which(seg_i$frequency == "Inf")), ]}
  if (missing(smooth_type)) {
    smooth_type = "none"}
  if(missing(average)){
    average = F
  }
  if(nrow(seg_i) < 5){max_i <- "NA"} else {  
    if(smooth_type == "none"){
      if(average == F){max_i <- max(seg_i$frequency)
      }
      if(average == T){max_i = mean(seg_i$frequency)} }
    if (smooth_type == "loess"){
      l <- loess(seg_i$frequency ~ seg_i$position, span = 0.1) #smoothening 
      if(average == F){ max_i <- max(l$fitted)} else if (average == T){max_i <- mean(l$fitted)}} }
  
  return(max_i)
  
}



