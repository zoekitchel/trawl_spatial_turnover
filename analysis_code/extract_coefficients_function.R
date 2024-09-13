lm_interaction_coefficients_se <- function(mod_name, model_avg = F, SST_fish = F, SBT_fish = F){
  coef_table_length <- (length(coef(mod_name)))
  
  if(coef_table_length == 68){
    coefs.r <- coef(mod_name)
    cov_mat <- vcov(mod_name) #get covariance matrix
  D <- matrix(c(0,1,rep(0,33),0,rep(0,32),
                   0,1,rep(0,33),1,rep(0,32),
                   0,1,rep(0,34),1,rep(0,31),
                   0,1,rep(0,35),1,rep(0,30),
                   0,1,rep(0,36),1,rep(0,29),
                   0,1,rep(0,37),1,rep(0,28),
                   0,1,rep(0,38),1,rep(0,27),
                   0,1,rep(0,39),1,rep(0,26),
                   0,1,rep(0,40),1,rep(0,25),
                   0,1,rep(0,41),1,rep(0,24),
                   0,1,rep(0,42),1,rep(0,23),
                   0,1,rep(0,43),1,rep(0,22),
                   0,1,rep(0,44),1,rep(0,21),
                   0,1,rep(0,45),1,rep(0,20),
                   0,1,rep(0,46),1,rep(0,19),
                   0,1,rep(0,47),1,rep(0,18),
                   0,1,rep(0,48),1,rep(0,17),
                   0,1,rep(0,49),1,rep(0,16),
                   0,1,rep(0,50),1,rep(0,15),
                   0,1,rep(0,51),1,rep(0,14),
                   0,1,rep(0,52),1,rep(0,13),
                   0,1,rep(0,53),1,rep(0,12),
                   0,1,rep(0,54),1,rep(0,11),
                   0,1,rep(0,55),1,rep(0,10),
                   0,1,rep(0,56),1,rep(0,9),
                   0,1,rep(0,57),1,rep(0,8),
                   0,1,rep(0,58),1,rep(0,7),
                   0,1,rep(0,59),1,rep(0,6),
                   0,1,rep(0,60),1,rep(0,5),
                   0,1,rep(0,61),1,rep(0,4),
                   0,1,rep(0,62),1,rep(0,3),
                   0,1,rep(0,63),1,rep(0,2),
                   0,1,rep(0,64),1,rep(0,1),
                   0,1,rep(0,65),1,rep(0,0)),
                 ncol = length(coef(mod_name)),
                 byrow = TRUE)
  }
  
  if(coef_table_length == 48){
    coefs.r <- coef(mod_name)
    cov_mat <- vcov(mod_name) #get covariance matrix
    D <- matrix(c(0,1,rep(0,23),0,rep(0,22),
                  0,1,rep(0,23),1,rep(0,22),
                  0,1,rep(0,24),1,rep(0,21),
                  0,1,rep(0,25),1,rep(0,20),
                  0,1,rep(0,26),1,rep(0,19),
                  0,1,rep(0,27),1,rep(0,18),
                  0,1,rep(0,28),1,rep(0,17),
                  0,1,rep(0,29),1,rep(0,16),
                  0,1,rep(0,30),1,rep(0,15),
                  0,1,rep(0,31),1,rep(0,14),
                  0,1,rep(0,32),1,rep(0,13),
                  0,1,rep(0,33),1,rep(0,12),
                  0,1,rep(0,34),1,rep(0,11),
                  0,1,rep(0,35),1,rep(0,10),
                  0,1,rep(0,36),1,rep(0,9),
                  0,1,rep(0,37),1,rep(0,8),
                  0,1,rep(0,38),1,rep(0,7),
                  0,1,rep(0,39),1,rep(0,6),
                  0,1,rep(0,40),1,rep(0,5),
                  0,1,rep(0,41),1,rep(0,4),
                  0,1,rep(0,42),1,rep(0,3),
                  0,1,rep(0,43),1,rep(0,2),
                  0,1,rep(0,44),1,rep(0,1),
                  0,1,rep(0,45),1,rep(0,0)),
                ncol = length(coef(mod_name)),
                byrow = TRUE)
  }
  
  if(model_avg == T & SST_fish == T){
    coefs.r <- coef(mod_name)[c(7,39:69)]
    
    cov_mat <- vcov(mod_name)[c(7,39:69),c(7,39:69)] #get covariance matrix
    
    D <- matrix(c(#fishing
                  1,rep(0,0),rep(0,31),
                  1,rep(0,0),1,rep(0,30),
                  1,rep(0,1),1,rep(0,29),
                  1,rep(0,2),1,rep(0,28),
                  1,rep(0,3),1,rep(0,27),
                  1,rep(0,4),1,rep(0,26),
                  1,rep(0,5),1,rep(0,25),
                  1,rep(0,6),1,rep(0,24),
                  1,rep(0,7),1,rep(0,23),
                  1,rep(0,8),1,rep(0,22),
                  1,rep(0,9),1,rep(0,21),
                  1,rep(0,10),1,rep(0,20),
                  1,rep(0,11),1,rep(0,19),
                  1,rep(0,12),1,rep(0,18),
                  1,rep(0,13),1,rep(0,17),
                  1,rep(0,14),1,rep(0,16),
                  1,rep(0,15),1,rep(0,15),
                  1,rep(0,16),1,rep(0,14),
                  1,rep(0,17),1,rep(0,13),
                  1,rep(0,18),1,rep(0,12),
                  1,rep(0,19),1,rep(0,11),
                  1,rep(0,20),1,rep(0,10),
                  1,rep(0,21),1,rep(0,9),
                  1,rep(0,22),1,rep(0,8),
                  1,rep(0,23),1,rep(0,7),
                  1,rep(0,24),1,rep(0,6),
                  1,rep(0,25),1,rep(0,5),
                  1,rep(0,26),1,rep(0,4),
                  1,rep(0,27),1,rep(0,3),
                  1,rep(0,28),1,rep(0,2),
                  1,rep(0,29),1,rep(0,1),
                  1,rep(0,30),1,rep(0,0)),
                  ncol = length(coefs.r),
                  byrow = TRUE)
    
  }
  
  if(model_avg == T & SBT_fish == T){
    coefs.r <- coef(mod_name)[c(7,39:101)]
    
    cov_mat <- vcov(mod_name)[c(7,39:101),c(7,39:101)] #get covariance matrix
    
    D <- matrix(c(#fishing
      1,0,0,rep(0,61),
      1,0,1,rep(0,61),
      1,rep(0,2),1,rep(0,60),
      1,rep(0,3),1,rep(0,59),
      1,rep(0,4),1,rep(0,58),
      1,rep(0,5),1,rep(0,57),
      1,rep(0,6),1,rep(0,56),
      1,rep(0,7),1,rep(0,55),
      1,rep(0,8),1,rep(0,54),
      1,rep(0,9),1,rep(0,53),
      1,rep(0,10),1,rep(0,52),
      1,rep(0,11),1,rep(0,51),
      1,rep(0,12),1,rep(0,50),
      1,rep(0,13),1,rep(0,49),
      1,rep(0,14),1,rep(0,48),
      1,rep(0,15),1,rep(0,47),
      1,rep(0,16),1,rep(0,46),
      1,rep(0,17),1,rep(0,45),
      1,rep(0,18),1,rep(0,44),
      1,rep(0,19),1,rep(0,43),
      1,rep(0,20),1,rep(0,42),
      1,rep(0,21),1,rep(0,41),
      1,rep(0,22),1,rep(0,40),
      1,rep(0,23),1,rep(0,39),
      1,rep(0,24),1,rep(0,38),
      1,rep(0,25),1,rep(0,37),
      1,rep(0,26),1,rep(0,36),
      1,rep(0,27),1,rep(0,35),
      1,rep(0,28),1,rep(0,34),
      1,rep(0,29),1,rep(0,33),
      1,rep(0,30),1,rep(0,32),
      1,rep(0,31),1,rep(0,31),
      #temp
      0,1,0,rep(0,61),
      0,1,rep(0,31),1,rep(0,30),
      0,1,rep(0,32),1,rep(0,29),
      0,1,rep(0,33),1,rep(0,28),
      0,1,rep(0,34),1,rep(0,27),
      0,1,rep(0,35),1,rep(0,26),
      0,1,rep(0,36),1,rep(0,25),
      0,1,rep(0,37),1,rep(0,24),
      0,1,rep(0,38),1,rep(0,23),
      0,1,rep(0,39),1,rep(0,22),
      0,1,rep(0,40),1,rep(0,21),
      0,1,rep(0,41),1,rep(0,20),
      0,1,rep(0,42),1,rep(0,19),
      0,1,rep(0,43),1,rep(0,18),
      0,1,rep(0,44),1,rep(0,17),
      0,1,rep(0,45),1,rep(0,16),
      0,1,rep(0,46),1,rep(0,15),
      0,1,rep(0,47),1,rep(0,14),
      0,1,rep(0,48),1,rep(0,13),
      0,1,rep(0,49),1,rep(0,12),
      0,1,rep(0,50),1,rep(0,11),
      0,1,rep(0,51),1,rep(0,10),
      0,1,rep(0,52),1,rep(0,9),
      0,1,rep(0,53),1,rep(0,8),
      0,1,rep(0,54),1,rep(0,7),
      0,1,rep(0,55),1,rep(0,6),
      0,1,rep(0,56),1,rep(0,5),
      0,1,rep(0,57),1,rep(0,4),
      0,1,rep(0,58),1,rep(0,3),
      0,1,rep(0,59),1,rep(0,2),
      0,1,rep(0,60),1,rep(0,1),
      0,1,rep(0,61),1,rep(0,0)),
      ncol = length(coefs.r),
      byrow = TRUE)
    
  }
  
  
  stopifnot(coef_table_length %in% c(48,68) | model_avg == T)
  
  #Get point estimates
  estimates <- D %*% coefs.r
  estimates <- estimates[, 1]
  
  #covariance
  vcov_est <- D %*% tcrossprod(cov_mat, D)
  se_estimates <- sqrt(diag(vcov_est))
  
  t_scores <- estimates / se_estimates
  
  if(model_avg == F){
  p_values <-
    2 * pt(abs(t_scores),
           ## degrees of freedom are taken from the model
           df = mod_name$df.residual,
           lower.tail = FALSE)
  }
  ## putting it all together:
  
  ## making D into "logical" to facilitate getting the names of the variables
  bool_d <- ifelse(D == 1, TRUE, FALSE)
  
  term <-
    apply(bool_d,
          1, function(id) {
            paste(names(coefs.r)[id],
                  collapse = " + ")
          })
  
  
  if(model_avg == F){
  output <- data.frame(term = term,
                       estimate = estimates,
                       se = se_estimates,
                       t_value = t_scores,
                       p_value = p_values)
  }
  
  if(model_avg == T){
    output <- data.frame(term = term,
                         estimate = estimates,
                         se = se_estimates,
                         t_value = t_scores)
  }
  
}
