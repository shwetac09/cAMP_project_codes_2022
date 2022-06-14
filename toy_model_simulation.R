hills <- function(b, x, vmax, k, n)b + ((vmax-b) * x^n / (k^n + x^n))

ind_given_camp <- function(camp, k_ind, n_ind, k_dir, n_dir){
  conc_dir <- hills(b = 0,x = camp, vmax = 1,k= k_dir,n= n_dir)
  conc_ind1 <- hills(b = 0, x = conc_dir, vmax = 1, k = k_ind, n = n_ind)
}

par_fit <- function(x, conc_vector){
  df <- cbind.data.frame(x, conc_vector)
  colnames(df) <- c("x", "y")
  ind_conc_fit <- nls(y ~ hills(b, x, vmax, k, n), data = df,
                      start = c(b = 0, vmax = 1, k = 0.6, n = n_ind),
                      algorithm = "port")
  pred_k <- coefficients(ind_conc_fit)[3]
  pred_n <- coefficients(ind_conc_fit)[4]
  return(c(pred_k, pred_n))
  
}


app_n_k_list <- NULL

for(i in 1:100){
k_sample_iter <- sample(k_dist,3)
k_dir <- k_sample_iter[1]
k_ind1 <- k_sample_iter[2]
k_ind2 <- k_sample_iter[3]

conc_dir <- hills(b = 0, x = camp, vmax = 1, k = k_dir, n = n_dir)
conc_ind1 <- hills(b = 0, x = conc_dir, vmax = 1, k = k_ind1, n = n_ind)
conc_ind2 <- hills(b = 0, x = conc_ind1, vmax = 1, k = k_ind2, n = n_ind)

try({app_n_k <- c(k_dir, par_fit(x = camp, conc_vector = conc_ind1),
             par_fit(x = camp, conc_vector = conc_ind2))

app_n_k_list <- rbind(app_n_k_list, app_n_k)})
}

