# ==============================================
# (0) 套件
# ==============================================
library(ggplot2)
library(Rtsne)
library(stats)
library(tidyverse)

# ==============================================
# (1) 基礎工具函數（圖檢查與 SI 模型）
# ==============================================

# 1. 檢查 p x p 的 adjacency matrix 是否為連通圖 (簡易 BFS)

is_connected <- function(adj) {
  p <- nrow(adj)
  visited <- rep(FALSE, p)
  queue   <- 1
  visited[1] <- TRUE
  
  while (length(queue) > 0) {
    cur <- queue[1]
    queue <- queue[-1]
    neighbors <- which(adj[cur, ] == 1)
    for (n in neighbors) {
      if (!visited[n]) {
        visited[n] <- TRUE
        queue <- c(queue, n)
      }
    }
  }
  return(all(visited))
}

# 2. SI 模型: 每一步是「隨機選取連結(edge)」而非「隨機選取節點(node)」
#   - G: p x p adjacency (0/1)
#   - s0: 傳播初始節點(單一整數)
#   - T : 執行步數

SI_infection <- function(G, s0, T) {
  p <- nrow(G)
  infected <- rep(0, p)
  infected[s0] <- 1
  
  for (step in seq_len(T)) {
    # 找出所有 "已感染->未感染" 的邊集合
    # 注: G 為對稱或無向圖，需小心避免重複計算；以下假設無向圖
    inf_nodes    <- which(infected == 1)
    uninf_nodes  <- which(infected == 0)
    # 建立所有可能的 (i, j)；i在inf_nodes, j在uninf_nodes
    # 然後只保留 G[i,j] == 1
    possible_edges <- c()
    for (i in inf_nodes) {
      neighbors <- uninf_nodes[G[i, uninf_nodes] == 1]
      if (length(neighbors) > 0) {
        # 將邊 (i,neighbor) 記錄
        possible_edges <- rbind(possible_edges,
                                cbind(i, neighbors))
      }
    }
    if (is.null(possible_edges)) {
      # 沒有邊可傳染 => 直接結束
      break
    }
    # 在所有可能邊中隨機抽 1 條邊
    sel_idx <- sample(seq_len(nrow(possible_edges)), 1)
    edge_chosen <- possible_edges[sel_idx, ]
    
    # 感染邊上的未感染者
    new_infected <- edge_chosen[2]  # (i, j) 第2個就是未感染點
    infected[new_infected] <- 1
  }
  
  return(infected)
}

# ==============================================
# (2) 資料模擬與生成函數
# ==============================================

generate_subtype_data <- function(
    d = 3,       # latent space 維度
    p = 100,     # ROI 數量
    target_density = 0.2,  # 目標邊密度
    mu = 0.2*p,  # Poisson參數
    sigma2 = 0.1,# 誤差變異量
    n_data = 100 # 每個亞型要生成幾筆 (h(Γ), y)
) {
  # -------------------------------------------------
  # 1. 生成 Omega
  #    u_i ~ N(0, I_d)  => Omega_{ij} = <u_i, u_j>, Omega_ii = 0
  # -------------------------------------------------
  #set.seed(8054)  # 可自行固定seed
  u_mat <- matrix(rnorm(p * d, mean = 0, sd = 1), nrow = p, ncol = d)
  Omega <- tcrossprod(u_mat)  # p x p, 其中 Omega[i,j] = u_i · u_j
  diag(Omega) <- 0       # Omega_ii = 0
  
  # -------------------------------------------------
  # 2. 生成 G (使 G 連通 & density ~ 0.2)
  #    G_{ij} ~ Bernoulli( sigmoid( Omega_{ij} - rho ) )
  #    如果不連通，重新生成
  # -------------------------------------------------
  # 這裡示範簡易隨機搜尋或二分法來調整rho，直到density ~ 0.2 且 G 連通
  # 為示範，可能寫一個小型的 binary search
  logistic <- function(x) 1 / (1 + exp(-x))
  
  lower <- -5
  upper <- 5
  for (iter in 1:50) {
    mid <- (lower + upper)/2
    # 產生臨時 G
    prob_mat <- logistic(Omega - mid)
    # 下三角先產生, 再對稱
    G_tmp <- matrix(rbinom(p*p, size=1, prob=as.vector(prob_mat)), nrow=p)
    # 做成無向圖: 令對角=0, 並對稱化
    diag(G_tmp) <- 0
    G_tmp[lower.tri(G_tmp)] <- t(G_tmp)[lower.tri(G_tmp)]
    
    dens <- mean(G_tmp[upper.tri(G_tmp)] == 1)
    
    if (dens > target_density) {
      # 邊太多 => rho要調高 (讓 logistic(Omega - rho) 變小)
      lower <- mid
    } else {
      upper <- mid
    }
    
    # 檢查是否連通 & 密度是否接近
    if (abs(dens - target_density) < 0.01 && is_connected(G_tmp)) {
      # 若同時符合「近似目標密度」且「連通」，就跳出
      G <- G_tmp
      break
    } else if (iter == 50) {
      # 最後一次迭代仍無法同時滿足，就強制以目前為準
      # 或者再來一個 while 也可以
      if (!is_connected(G_tmp)) {
        # 若不連通，就簡單 repeat 整個函式, 這裡省略, 實務可再遞迴等
        message("Warning: cannot find connected G with desired density. Using final G.")
      }
      G <- G_tmp
    }
  }
  
  # -------------------------------------------------
  # 3. 選一個 s0 (可固定或隨機)
  # -------------------------------------------------
  s0 <- sample.int(p, 1)
  
  # -------------------------------------------------
  # 4. 生成 n_data 筆 (h(Γ), y)
  #    - B=1, C=1, f=0 (投影片設定)
  #    - T ~ Pois(mu), y 由 SI(G, s0, T) 生成
  #    - h(Γ_{ij}) = Omega_{ij} + |y_i - y_j| + |y_i + y_j| + ε_{ij}
  # -------------------------------------------------
  data_list <- vector("list", n_data)
  
  for (n_idx in 1:n_data) {
    # (a) 產生 y via SI
    Tval <- rpois(1, lambda = mu)
    y_vec <- SI_infection(G = G, s0 = s0, T = Tval)
    
    # (b) 產生 誤差矩陣 eps ~ N(0, sigma2)
    eps_mat <- matrix(rnorm(p*p, mean=0, sd=sqrt(sigma2)), nrow=p)
    
    # (c) 根據模型計算 h(Γ_{ij})
    #     h(Γ_{ij}) = Omega_{ij} + |y_i - y_j| + |y_i + y_j| + eps
    #   其中 y_i, y_j ∈ {0,1}
    #   注意: 令 diag=NA 或 0 視需求而定
    hGamma <- matrix(0, nrow=p, ncol=p)
    # 如果只有 0,1, 则 |y_i+y_j| 是 0(兩端都 0)或 2(兩端都 1)，或 1(一端 1,一端 0)
    # 但簡單做法直接用 abs() 即可
    for (i in 1:p) {
      for (j in 1:p) {
        tmp_val <- Omega[i,j] +
          abs(y_vec[i] - y_vec[j]) +
          abs(y_vec[i] + y_vec[j]) +
          eps_mat[i,j]
        hGamma[i,j] <- tmp_val
      }
    }
    # (d) 儲存
    data_list[[n_idx]] <- list(
      hGamma = hGamma,
      y      = y_vec,
      Tval   = sum(y_vec)   # 也可存下 T, s0 等資訊看需要
    )
  }
  
  # 回傳: 包含 (Omega, G, s0, 及對應資料列表)
  # 註: Omega, G, s0 是「整個亞型共同」；(hGamma, y) 是「每筆病人資料」
  return(list(
    Omega = Omega,
    G     = G,
    s0    = s0,
    data  = data_list
  ))
}

# ==============================================
# (3) 組內/組間距離分析
# ==============================================

compute_within_between <- function(all_subtypes) {
  J <- length(all_subtypes)
  within_d <- numeric(0)
  
  for (j in seq_len(J)) {
    Omega_j <- all_subtypes[[j]]$Omega
    data_j  <- all_subtypes[[j]]$data
    
    for (i in seq_along(data_j)) {
      hGamma_ji <- data_j[[i]]$hGamma  # 已是 logit(h(Γ))
      y_ji <- data_j[[i]]$y            # 長度為 p 的 0/1 向量
      
      if (is.null(y_ji)) next  # 若 y 缺失就跳過
      
      # 計算 |y_i - y_j| 與 |y_i + y_j|
      y_diff <- abs(outer(y_ji, y_ji, "-"))
      y_sum  <- abs(outer(y_ji, y_ji, "+"))
      
      expected_hGamma <- Omega_j + y_diff + y_sum  # 假設 B = 1, C = 1
      
      d <- norm(hGamma_ji - expected_hGamma, type = "F")
      if (!is.nan(d) && !is.na(d)) {
        within_d <- c(within_d, d)
      }
    }
  }
  
  # 組間距離: baseline 之間的距離
  between_d <- numeric(0)
  if (J > 1) {
    for (j in 1:(J - 1)) {
      Omega_j <- all_subtypes[[j]]$Omega
      for (k in (j + 1):J) {
        Omega_k <- all_subtypes[[k]]$Omega
        dist_jk <- norm(Omega_j - Omega_k, type = "F")
        between_d <- c(between_d, dist_jk)
      }
    }
  } else {
    between_d <- NA
  }
  
  within_mean  <- mean(within_d, na.rm = TRUE)
  between_mean <- mean(between_d, na.rm = TRUE)
  ratio        <- within_mean / between_mean
  
  list(
    within_mean = within_mean,
    between_mean = between_mean,
    ratio = ratio,
    within_d = within_d,
    between_d = between_d
  )
}


# ==============================================
# (4) 建圖與鄰接矩陣函數
# ==============================================

build_graph <- function(rho_matrix, type = c("full", "topk", "block", "thres"), 
                        k = 10, thres = 0.5, num_block = 2) {
  N <- nrow(rho_matrix)
  edges <- list(); weights <- c()
  
  if (type == "block") {
    # 依 num_block 均分節點（最後一組可能稍大）
    group_sizes <- rep(floor(N / num_block), num_block)
    remainder <- N %% num_block
    if (remainder > 0) {
      group_sizes[1:remainder] <- group_sizes[1:remainder] + 1
    }
    group_indices <- split(1:N, rep(1:num_block, times = group_sizes))
    
    for (group in group_indices) {
      for (i in group) {
        for (j in group) {
          if (i < j) {
            edges[[length(edges) + 1]] <- c(i, j)
            weights <- c(weights, 1)
          }
        }
      }
    }
    return(list(E = do.call(rbind, edges), w = weights))
  }
  
  if (type == "full") {
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        edges[[length(edges) + 1]] <- c(i, j)
        weights <- c(weights, rho_matrix[i, j])
      }
    }
  } else if (type == "topk") {
    for (i in 1:N) {
      neighbors <- order(rho_matrix[i, ], decreasing = TRUE)[2:(k + 1)]
      for (j in neighbors) {
        if (i < j) {
          edges[[length(edges) + 1]] <- c(i, j)
          weights <- c(weights, rho_matrix[i, j])
        }
      }
    }
  } else if (type == "thres") {
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        if (rho_matrix[i, j] > thres) {
          edges[[length(edges) + 1]] <- c(i, j)
          weights <- c(weights, rho_matrix[i, j])  # 保留原本的 rho 作為權重
        }
      }
    }
  }
  
  return(list(E = do.call(rbind, edges), w = weights))
}

# -----------------------------
# 將結果返回矩陣形式
# -----------------------------
get_adjacency_matrix <- function(graph_result, N) {
  adj_mat <- matrix(0, nrow = N, ncol = N)
  E <- graph_result$E
  w <- graph_result$w
  
  for (i in seq_len(nrow(E))) {
    u <- E[i, 1]
    v <- E[i, 2]
    adj_mat[u, v] <- w[i]
    adj_mat[v, u] <- w[i]  # 因為是無向圖
  }
  return(adj_mat)
}

# ==============================================
# (4.2) 建圖與鄰接矩陣函數
# ==============================================

build_graph_matrix <- function(rho_matrix, type = c("full", "topk", "block", "thres"), 
                               k = 10, thres = 0.5, num_block = 2) {
  N <- nrow(rho_matrix)
  W <- matrix(0, nrow = N, ncol = N)
  
  if (type == "block") {
    group_sizes <- rep(floor(N / num_block), num_block)
    remainder <- N %% num_block
    if (remainder > 0) {
      group_sizes[1:remainder] <- group_sizes[1:remainder] + 1
    }
    group_indices <- split(1:N, rep(1:num_block, times = group_sizes))
    
    for (group in group_indices) {
      for (i in group) {
        for (j in group) {
          if (i < j) {
            W[i, j] <- 1
            W[j, i] <- 1
          }
        }
      }
    }
    
  } else if (type == "full") {
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        W[i, j] <- rho_matrix[i, j]
        W[j, i] <- rho_matrix[i, j]
      }
    }
    
  } else if (type == "topk") {
    for (i in 1:N) {
      neighbors <- order(rho_matrix[i, ], decreasing = TRUE)[2:(k + 1)]
      for (j in neighbors) {
        W[i, j] <- rho_matrix[i, j]
      }
    }
    # ✅ 對稱化處理
    W <- 0.5 * (W + t(W))
    
  } else if (type == "thres") {
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        if (rho_matrix[i, j] > thres) {
          W[i, j] <- rho_matrix[i, j]
          W[j, i] <- rho_matrix[i, j]
        }
      }
    }
  }
  
  return(W)
}

# ============================================================================
# (5) 上三角前處理
# ============================================================================

# 將矩陣取出上三角（不含對角）
get_upper_triangle <- function(mat) {
  mat[lower.tri(mat, diag = TRUE)] <- NA
  return(as.vector(na.omit(as.vector(mat))))
}

# 將上三角向量還原成對稱矩陣，並強制對角為 0
restore_full_matrix <- function(upper_vec, p) {
  mat <- matrix(0, p, p)
  mat[upper.tri(mat)] <- upper_vec
  mat <- mat + t(mat)
  diag(mat) <- 0
  return(mat)
}

# 計算 true 與 est 的距離矩陣（上三角 Frobenius）
compute_omega_distance <- function(Omega_true_list, Omega_est_list) {
  num_true <- length(Omega_true_list)
  num_est  <- length(Omega_est_list)
  dist_mat <- matrix(0, nrow = num_true, ncol = num_est)
  
  for (i in 1:num_true) {
    u_true <- get_upper_triangle(Omega_true_list[[i]])
    for (j in 1:num_est) {
      u_est <- get_upper_triangle(Omega_est_list[[j]])
      dist_mat[i, j] <- sqrt(sum((u_true - u_est)^2))
    }
  }
  
  rownames(dist_mat) <- paste0("True_", seq_len(num_true))
  colnames(dist_mat) <- paste0("Est_", seq_len(num_est))
  return(dist_mat)
}
