setwd("/Users/mapinguari/Dropbox/0_1Rebirth/Papers/Submetidos/Oikos/revisao")
## ==========================================================
## Simulação de paisagens fragmentadas + distâncias entre fragmentos
## ==========================================================

## -------------------------------
## 1) Função de simulação
## -------------------------------
simular_fragmentacao <- function(X, H, F, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (H < 0 || H > X^2) stop("H deve estar entre 0 e X^2.")
  if (F <= 0) stop("F deve ser positivo.")
  if (H < F) stop("H deve ser >= F.")
  
  habitat     <- matrix(0, X, X)
  fragment_id <- matrix(0, X, X)
  
  # tamanhos alvo dos fragmentos
  tamanhos <- rep(1, F)
  restante <- H - F
  if (restante > 0) {
    extra <- as.vector(rmultinom(1, restante, prob = rep(1, F)))
    tamanhos <- tamanhos + extra
  }
  
  # helpers
  idx_to_rc <- function(idx, X) {
    list(
      r = ((idx - 1) %% X) + 1,
      c = ((idx - 1) %/% X) + 1
    )
  }
  rc_to_idx <- function(r, c, X) (c - 1) * X + r
  
  vizinhos_8 <- function(cells, X) {
    rc <- idx_to_rc(cells, X)
    r  <- rc$r
    c  <- rc$c
    
    rs <- c(r-1, r-1, r-1, r,   r,   r+1, r+1, r+1)
    cs <- c(c-1, c,   c+1, c-1, c+1, c-1, c,   c+1)
    ok <- which(rs >= 1 & rs <= X & cs >= 1 & cs <= X)
    if (length(ok) == 0) return(integer(0))
    rc_to_idx(rs[ok], cs[ok], X)
  }
  
  encosta_outro_fragmento <- function(idx, f, fragment_id, X) {
    rc <- idx_to_rc(idx, X)
    r  <- rc$r
    c  <- rc$c
    
    rs <- c(r-1, r-1, r-1, r,   r,   r+1, r+1, r+1)
    cs <- c(c-1, c,   c+1, c-1, c+1, c-1, c,   c+1)
    ok <- which(rs >= 1 & rs <= X & cs >= 1 & cs <= X)
    if (length(ok) == 0) return(FALSE)
    
    viz <- fragment_id[cbind(rs[ok], cs[ok])]
    any(viz != 0 & viz != f)
  }
  
  frag_cells <- vector("list", F)
  
  for (f in seq_len(F)) {
    alvo <- tamanhos[f]
    if (alvo <= 0) next
    
    ## --- semente ---
    cel_vaz <- which(habitat == 0)
    if (length(cel_vaz) == 0) break
    
    cand <- cel_vaz[
      !vapply(
        cel_vaz,
        function(i) encosta_outro_fragmento(i, f, fragment_id, X),
        logical(1L)
      )
    ]
    if (length(cand) == 0) next  # não conseguiu semente para esse fragmento
    
    seed <- sample(cand, 1)
    rc   <- idx_to_rc(seed, X)
    
    habitat[rc$r, rc$c]     <- 1
    fragment_id[rc$r, rc$c] <- f
    frag_cells[[f]]        <- seed
    
    ## --- crescimento ---
    while (length(frag_cells[[f]]) < alvo) {
      neigh <- vizinhos_8(frag_cells[[f]], X)
      cel_vaz <- which(habitat == 0)
      fronteira <- intersect(neigh, cel_vaz)
      if (length(fronteira) == 0) break
      
      fronteira_ok <- fronteira[
        !vapply(
          fronteira,
          function(i) encosta_outro_fragmento(i, f, fragment_id, X),
          logical(1L)
        )
      ]
      if (length(fronteira_ok) == 0) break
      
      novo <- sample(fronteira_ok, 1)
      rc2  <- idx_to_rc(novo, X)
      
      habitat[rc2$r, rc2$c]     <- 1
      fragment_id[rc2$r, rc2$c] <- f
      frag_cells[[f]]          <- c(frag_cells[[f]], novo)
    }
  }
  
  list(
    habitat = habitat,
    fragment_id = fragment_id
  )
}

## -------------------------------
## 2) Distância média ao fragmento mais próximo
## -------------------------------
dist_media_vizinho_mais_proximo <- function(fragment_id) {
  ids <- sort(setdiff(unique(as.vector(fragment_id)), 0))
  n   <- length(ids)
  if (n < 2) return(NA_real_)
  
  coords <- lapply(ids, function(id) which(fragment_id == id, arr.ind = TRUE))
  names(coords) <- ids
  
  dist_min_matrix <- matrix(Inf, n, n)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      A <- coords[[i]]
      B <- coords[[j]]
      
      d1 <- outer(A[, 1], B[, 1], "-")
      d2 <- outer(A[, 2], B[, 2], "-")
      dist <- sqrt(d1^2 + d2^2)
      
      dmin <- min(dist)
      dist_min_matrix[i, j] <- dmin
      dist_min_matrix[j, i] <- dmin
    }
  }
  
  viz_prox <- apply(dist_min_matrix, 1, min)
  mean(viz_prox)
}

## -------------------------------
## 3) Loop de simulação
##    - 100 réplicas válidas
##    - Máx. 1000 tentativas
## -------------------------------
X <- 50
H <- 250
F <- 8

N_VALIDAS      <- 100
MAX_TENTATIVAS <- 1000

dist_nn      <- numeric(0)
n_frag_real  <- integer(0)
n_hab_real   <- integer(0)

cont       <- 0   # réplicas válidas
tentativas <- 0   # tentativas totais

while (cont < N_VALIDAS && tentativas < MAX_TENTATIVAS) {
  tentativas <- tentativas + 1
  cat("Tentativa", tentativas,
      "- réplicas válidas:", cont, "de", N_VALIDAS, "\n")
  
  res  <- simular_fragmentacao(X, H, F)
  frag <- res$fragment_id
  hab  <- res$habitat
  
  ids <- sort(setdiff(unique(as.vector(frag)), 0))
  
  ## --- checar se gerou todos os fragmentos (1:F) ---
  if (length(ids) != F) {
    next  # descarta simulação
  }
  
  ## --- checar se algum fragmento ficou vazio (por garantia extra) ---
  ok <- TRUE
  for (f_id in ids) {
    if (sum(frag == f_id) == 0) {
      ok <- FALSE
      break
    }
  }
  if (!ok) next
  
  ## --- simulação é válida ---
  cont <- cont + 1
  
  dist_nn[cont]     <- dist_media_vizinho_mais_proximo(frag)
  n_frag_real[cont] <- length(ids)
  n_hab_real[cont]  <- sum(hab)
}

if (cont < N_VALIDAS) {
  warning(paste(
    "Foram obtidas apenas", cont, "réplicas válidas em",
    MAX_TENTATIVAS, "tentativas."
  ))
}

## -------------------------------
## 4) Montar data frame e salvar CSV
## -------------------------------
df <- data.frame(
  simulacao                        = 1:cont,
  dist_media_vizinho_mais_proximo  = dist_nn,
  n_fragmentos                     = n_frag_real,
  n_habitat_real                   = n_hab_real
)

nome_arquivo <- paste0(
  "resultados_fragmentacao_X", X,
  "_H", H,
  "_F", F,
  ".csv"
)

# Colunas separadas por ;  e decimais com ,
write.table(
  df,
  file      = nome_arquivo,
  sep       = ";",   # separador de colunas
  dec       = ",",   # separador decimal
  row.names = FALSE,
  col.names = TRUE,
  qmethod   = "double"
)

cat("\nFinalizado!\n")
cat("Réplicas válidas obtidas:", cont, "de", N_VALIDAS, "\n")
cat("Arquivo salvo como:", nome_arquivo, "\n")

## ----------------------------------------------------
## 5) Estatísticas finais das simulações válidas
## ----------------------------------------------------

media_dist <- mean(dist_nn, na.rm = TRUE)
desvio_dist <- sd(dist_nn, na.rm = TRUE)
n_validas_final <- cont

cat("\n================ RESULTADOS FINAIS ================\n")
cat("Média da distância mínima (100 simulações): ", media_dist, "\n", sep = "")
cat("Desvio padrão das distâncias: ", desvio_dist, "\n", sep = "")
cat("Número de simulações bem-sucedidas: ", n_validas_final, "\n", sep = "")
cat("====================================================\n")



