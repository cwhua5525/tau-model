# ============================================================================
# This file contains AMA and ADMM functions for a single lambda
# runtime and initial(warm start control added).
# dual included in input and output
# ============================================================================
# ============================================================================
# 1.network_lasso_ama_single
# Description: Solve Network Lasso problem via AMA (single lambda setting) with Dual Gap Calculation and Time Tracking
# ============================================================================

network_lasso_ama_single <- function(
    X_raw_upper,
    rho_matrix_adjusted,
    lambda,
    max_iter = 10000,
    tol = 1e-5,
    quiet = FALSE,
    Pri.init = NULL,
    Dual.init = NULL
) {
  # -----------------------
  # 基本設定
  # -----------------------
  N <- nrow(X_raw_upper)
  d <- ncol(X_raw_upper)
  lambda_penalty <- lambda / 2
  
  degrees <- rowSums(rho_matrix_adjusted > 0)
  Lmax <- max(outer(degrees, degrees, "+")[rho_matrix_adjusted > 0])
  nu <- 1.9 / Lmax
  
  # -----------------------
  # 初始化 U
  # -----------------------
  if (is.null(Pri.init)) {
    U <- X_raw_upper
  } else {
    U <- Pri.init
  }
  
  # -----------------------
  # 初始化 edge list 與 Lambda
  # -----------------------
  edge_list <- list()
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      if (rho_matrix_adjusted[i, j] > 0) {
        edge_list[[length(edge_list) + 1]] <- c(i, j)
      }
    }
  }
  
  
  # -----------------------
  # 初始化 Lambda
  # -----------------------
  if (is.null(Dual.init)) {
    Lambda <- vector("list", length(edge_list))
    for (k in seq_along(edge_list)) {
      i <- edge_list[[k]][1]
      j <- edge_list[[k]][2]
      u_diff <- U[i, ] - U[j, ]
      norm_u_diff <- sqrt(sum(u_diff^2))
      if (norm_u_diff > 0) {
        Lambda[[k]] <- lambda * rho_matrix_adjusted[i, j] * (u_diff / norm_u_diff)
      } else {
        Lambda[[k]] <- rep(0, d)
      }
    }
  } else {
    Lambda <- Dual.init
  }
  

  
  # -----------------------
  # 建立記錄器
  # -----------------------
  fit_loss_vec <- numeric(0)
  clu_loss_vec <- numeric(0)
  total_loss_vec <- numeric(0)
  dual_obj_vec <- numeric(0)
  gap_vec <- numeric(0)
  iter_time_vec <- numeric(0)
  
  # -----------------------
  # AMA 主迴圈
  # -----------------------
  converged <- FALSE
  for (iter in 1:max_iter) {
    start_time <- Sys.time()
    
    # Step 1: 更新 U
    delta <- matrix(0, nrow = N, ncol = d)
    for (k in seq_along(edge_list)) {
      i <- edge_list[[k]][1]
      j <- edge_list[[k]][2]
      delta[i, ] <- delta[i, ] + Lambda[[k]]
      delta[j, ] <- delta[j, ] - Lambda[[k]]
    }
    U_new <- X_raw_upper + delta
    
    # Step 2: 更新 Lambda
    for (k in seq_along(edge_list)) {
      i <- edge_list[[k]][1]
      j <- edge_list[[k]][2]
      g_ij <- U_new[i, ] - U_new[j, ]
      z <- Lambda[[k]] - nu * g_ij
      norm_z <- sqrt(sum(z^2))
      radius <- lambda_penalty * rho_matrix_adjusted[i, j]
      if (norm_z > radius) {
        z <- z * (radius / norm_z)
      }
      Lambda[[k]] <- z
    }
    
    # Step 3: 計算 Primal Loss
    loss_fitting <- sum(rowSums((X_raw_upper - U_new)^2))
    loss_clustering <- 0
    for (k in seq_along(edge_list)) {
      i <- edge_list[[k]][1]
      j <- edge_list[[k]][2]
      diff_ij <- U_new[i, ] - U_new[j, ]
      loss_clustering <- loss_clustering + rho_matrix_adjusted[i, j] * sqrt(sum(diff_ij^2))
    }
    loss_total <- loss_fitting + lambda * loss_clustering
    
    # Step 4: 計算 Dual Loss
    delta_dual <- matrix(0, nrow = N, ncol = d)
    for (k in seq_along(edge_list)) {
      i <- edge_list[[k]][1]
      j <- edge_list[[k]][2]
      delta_dual[i, ] <- delta_dual[i, ] + Lambda[[k]]
      delta_dual[j, ] <- delta_dual[j, ] - Lambda[[k]]
    }
    dual_obj <- -sum(rowSums(delta_dual^2)) - 2 * sum(sapply(seq_along(edge_list), function(k) {
      i <- edge_list[[k]][1]
      j <- edge_list[[k]][2]
      sum(Lambda[[k]] * (X_raw_upper[i, ] - X_raw_upper[j, ]))
    }))
    
    gap <- loss_total - dual_obj
    
    end_time <- Sys.time()
    iter_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Step 5: 記錄
    fit_loss_vec <- c(fit_loss_vec, loss_fitting)
    clu_loss_vec <- c(clu_loss_vec, lambda * loss_clustering)
    total_loss_vec <- c(total_loss_vec, loss_total)
    dual_obj_vec <- c(dual_obj_vec, dual_obj)
    gap_vec <- c(gap_vec, gap)
    iter_time_vec <- c(iter_time_vec, iter_time)
    
    # Step 6: 印出進度
    if (!quiet) {
      cat(sprintf("Iter %4d | Obj: %.4f | Fit: %.4f | Clu: %.4f | Dual: %.4f | Gap: %.6f\n",
                  iter, loss_total, loss_fitting, lambda * loss_clustering, dual_obj, gap))
      flush.console()
    }
    
    # Step 7: 收斂判斷 (改成 Gap-based)
    if (gap < tol) {
      converged <- TRUE
      break
    }
    
    U <- U_new
  }
  
  # -----------------------
  # 組裝輸出
  # -----------------------
  result <- list(
    U = U,
    Lambda = Lambda,
    final_loss_all = loss_total,
    final_loss_fit = loss_fitting,
    final_loss_clu = lambda * loss_clustering,
    n_iter = iter,
    converged = converged,
    loss = list(
      fit = fit_loss_vec,
      clu = clu_loss_vec,
      all = total_loss_vec,
      time = iter_time_vec
    ),
    res = list(
      dual_obj = dual_obj_vec,
      gap = gap_vec
    )
  )
  
  return(result)
}



# ============================================================================
# 2.network_lasso_admm_single
# Description: Solve Network Lasso problem via ADMM (single lambda setting)
# ============================================================================

network_lasso_admm_single <- function(
    X_raw_upper,
    rho_matrix_adjusted,
    lambda,
    max_iter = 1000,
    tol_primal = 1e-4,
    tol_dual = 1e-4,
    rho_admm = 1.0,
    quiet = FALSE,
    Pri.init = NULL,
    Aux.init = NULL,
    Dual.init = NULL
) {
  # -----------------------
  # 基本設定
  # -----------------------
  N <- nrow(X_raw_upper)
  d <- ncol(X_raw_upper)
  lambda_penalty <- lambda / 2
  
  # -----------------------
  # 初始化 X
  # -----------------------
  if (is.null(Pri.init)) {
    X <- X_raw_upper
  } else {
    X <- Pri.init
  }
  
  edge_list <- list()
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      if (rho_matrix_adjusted[i, j] > 0) {
        edge_list[[length(edge_list) + 1]] <- c(i, j)
      }
    }
  }
  
  
  # -----------------------
  # 初始化 Z (Auxillary)
  # -----------------------
  if (is.null(Aux.init)) {
    Z_list <- vector("list", length(edge_list))
    for (k in seq_along(edge_list)) {
      i <- edge_list[[k]][1]
      j <- edge_list[[k]][2]
      Z_list[[k]] <- list(zij = X[i, ], zji = X[j, ])
    }
  } else {
    Z_list <- Aux.init
  }
  
  # -----------------------
  # 初始化 U (Dual)
  # -----------------------
  if (is.null(Dual.init)) {
    U_list <- vector("list", length(edge_list))
    for (k in seq_along(edge_list)) {
      i <- edge_list[[k]][1]
      j <- edge_list[[k]][2]
      U_list[[k]] <- list(uij = rep(0, d), uji = rep(0, d))
    }
  } else {
    U_list <- Dual.init
  }
  
  
  # Z_list <- vector("list", length(edge_list))
  # U_list <- vector("list", length(edge_list))
  # for (k in seq_along(edge_list)) {
  #   i <- edge_list[[k]][1]
  #   j <- edge_list[[k]][2]
  #   Z_list[[k]] <- list(zij = X[i, ], zji = X[j, ])
  #   U_list[[k]] <- list(uij = rep(0, d), uji = rep(0, d))
  # }
  
  
  
  
  
  # -----------------------
  # 建立記錄器
  # -----------------------
  fit_loss_vec <- numeric(0)
  clu_loss_vec <- numeric(0)
  total_loss_vec <- numeric(0)
  primal_res_vec <- numeric(0)
  dual_res_vec <- numeric(0)
  time_vec <- numeric(0)
  
  # -----------------------
  # 主迴圈
  # -----------------------
  for (iter in 1:max_iter) {
    iter_start_time <- Sys.time()
    
    X_old <- X
    
    # Step 1: x-update
    delta_x <- matrix(0, nrow = N, ncol = d)
    count_x <- rep(0, N)
    
    for (k in seq_along(edge_list)) {
      i <- edge_list[[k]][1]
      j <- edge_list[[k]][2]
      zij <- Z_list[[k]]$zij
      uij <- U_list[[k]]$uij
      zji <- Z_list[[k]]$zji
      uji <- U_list[[k]]$uji
      
      delta_x[i, ] <- delta_x[i, ] + rho_admm * (zij - uij)
      delta_x[j, ] <- delta_x[j, ] + rho_admm * (zji - uji)
      count_x[i] <- count_x[i] + rho_admm
      count_x[j] <- count_x[j] + rho_admm
    }
    
    for (i in 1:N) {
      if (count_x[i] > 0) {
        X[i, ] <- (X_raw_upper[i, ] + delta_x[i, ]) / (1 + count_x[i])
      } else {
        X[i, ] <- X_raw_upper[i, ]
      }
    }
    
    # Step 2: z-update
    for (k in seq_along(edge_list)) {
      i <- edge_list[[k]][1]
      j <- edge_list[[k]][2]
      rho_ij <- rho_matrix_adjusted[i, j]
      
      v1 <- X[i, ] + U_list[[k]]$uij
      v2 <- X[j, ] + U_list[[k]]$uji
      diff <- v1 - v2
      norm_diff <- sqrt(sum(diff^2))
      
      if (norm_diff == 0) {
        theta <- 0.5
      } else {
        theta <- max(1 - lambda_penalty * rho_ij / (rho_admm * norm_diff), 0.5)
      }
      
      Z_list[[k]]$zij <- theta * v1 + (1 - theta) * v2
      Z_list[[k]]$zji <- (1 - theta) * v1 + theta * v2
    }
    
    # Step 3: u-update and residuals
    primal_residual <- 0
    dual_residual <- 0
    
    for (k in seq_along(edge_list)) {
      i <- edge_list[[k]][1]
      j <- edge_list[[k]][2]
      
      uij_old <- U_list[[k]]$uij
      uji_old <- U_list[[k]]$uji
      
      U_list[[k]]$uij <- uij_old + (X[i, ] - Z_list[[k]]$zij)
      U_list[[k]]$uji <- uji_old + (X[j, ] - Z_list[[k]]$zji)
      
      primal_residual <- primal_residual + sum((X[i, ] - Z_list[[k]]$zij)^2) + sum((X[j, ] - Z_list[[k]]$zji)^2)
      dual_residual <- dual_residual + sum((rho_admm * (X[i, ] - X_old[i, ]))^2) + sum((rho_admm * (X[j, ] - X_old[j, ]))^2)
    }
    
    primal_residual <- sqrt(primal_residual)
    dual_residual <- sqrt(dual_residual)
    
    # 記錄 loss
    loss_fitting <- sum(rowSums((X_raw_upper - X)^2))
    loss_clustering <- 0
    for (k in seq_along(edge_list)) {
      i <- edge_list[[k]][1]
      j <- edge_list[[k]][2]
      diff_ij <- X[i, ] - X[j, ]
      loss_clustering <- loss_clustering + rho_matrix_adjusted[i, j] * sqrt(sum(diff_ij^2))
    }
    loss_total <- loss_fitting + lambda * loss_clustering
    
    fit_loss_vec <- c(fit_loss_vec, loss_fitting)
    clu_loss_vec <- c(clu_loss_vec, lambda * loss_clustering)
    total_loss_vec <- c(total_loss_vec, loss_total)
    primal_res_vec <- c(primal_res_vec, primal_residual)
    dual_res_vec <- c(dual_res_vec, dual_residual)
    
    iter_end_time <- Sys.time()
    time_vec <- c(time_vec, as.numeric(difftime(iter_end_time, iter_start_time, units = "secs")))
    
    # Optional: print進度
    if (!quiet) {
      cat(sprintf(
        "Iter %4d | Obj: %.4f | Fit: %.4f | Clu: %.4f | Primal: %.6f | Dual: %.6f\n",
        iter, loss_total, loss_fitting, lambda * loss_clustering, primal_residual, dual_residual
      ))
      flush.console()
    }
    
    # Check convergence
    if (primal_residual < tol_primal && dual_residual < tol_dual) {
      converged <- TRUE
      break
    } else {
      converged <- FALSE
    }
  }
  
  # -----------------------
  # 最後組裝輸出
  # -----------------------
  result <- list(
    X = X,
    Z_list = Z_list,
    U_list = U_list,
    final_loss_all = loss_total,
    final_loss_fit = loss_fitting,
    final_loss_clu = lambda * loss_clustering,
    n_iter = iter,
    converged = converged,
    loss = list(
      fit = fit_loss_vec,
      clu = clu_loss_vec,
      all = total_loss_vec,
      time = time_vec
    ),
    res = list(
      primal = primal_res_vec,
      dual = dual_res_vec
    )
  )
  
  return(result)
}
