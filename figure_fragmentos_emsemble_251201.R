## Simulação de paisagens fragmentadas + distâncias entre fragmentos
## - Só aceita simulações que atingem exatamente H e F (e tamanhos alvo)
## IMPORTANTE: combinações de valores altos de H e F podem ser impossíveis de gerar
## Distância mínima entre fragmentos é igual a 2

## Fragmented landscape simulations + distances among fragments
## - Only accepts simulations that reach exactly H and F (and target sizes)
## IMPORTANT: combinations with high values of H and F may be impossible to generate
## The minimum distance between fragments is equal to 2


## FUNÇÕES AUXILIARES DE VERIFICAÇÃO
## AUXILIARY CHECKING FUNCTIONS

## Verifica se há contato entre fragmentos diferentes (vizinhança de Moore)
## Checks whether different fragments touch each other (Moore neighborhood)

tem_contato_moore <- function(fragment_id) {
  X <- nrow(fragment_id)
  dr <- c(-1,-1,-1, 0,0, 1,1,1)
  dc <- c(-1, 0, 1,-1,1,-1,0,1)

  ids <- sort(setdiff(unique(as.vector(fragment_id)), 0))
  for (f in ids) {
    pos <- which(fragment_id == f, arr.ind = TRUE)
    for (k in seq_len(nrow(pos))) {
      r <- pos[k, 1]; c <- pos[k, 2]
      rs <- r + dr; cs <- c + dc
      ok <- rs >= 1 & rs <= X & cs >= 1 & cs <= X
      if (!any(ok)) next

      viz <- fragment_id[cbind(rs[ok], cs[ok])]
      if (any(viz != 0 & viz != f)) return(TRUE)
    }
  }
  FALSE
}

## Verifica se algum fragmento (mesmo ID) está partido em mais de um componente
## Checks whether any fragment (same ID) is split into more than one connected component

tem_fragmento_partido <- function(fragment_id) {
  X <- nrow(fragment_id)
  dr <- c(-1,-1,-1, 0,0, 1,1,1)
  dc <- c(-1, 0, 1,-1,1,-1,0,1)

  ids <- sort(setdiff(unique(as.vector(fragment_id)), 0))
  for (f in ids) {
    cells <- which(fragment_id == f)
    if (length(cells) <= 1) next

    ## Operador booleana no índice linear
    ## Boolean operator over linear indices
    cell_set <- rep(FALSE, X * X)
    cell_set[cells] <- TRUE

    ## Busca em largura (BFS) a partir da primeira célula
    ## Breadth-first search (BFS) starting from the first cell
    queue <- cells[1]
    reached <- 0L

    while (length(queue) > 0) {
      cur <- queue[1]
      queue <- queue[-1]

      if (!cell_set[cur]) next
      cell_set[cur] <- FALSE
      reached <- reached + 1L

      ## Índice linear -> (linha, coluna) 
      ## Linear index -> (row, column) 
      r <- ((cur - 1) %% X) + 1
      c <- ((cur - 1) %/% X) + 1

      rs <- r + dr
      cs <- c + dc
      ok <- rs >= 1 & rs <= X & cs >= 1 & cs <= X
      if (!any(ok)) next

      neigh <- (cs[ok] - 1) * X + rs[ok]
      neigh <- neigh[cell_set[neigh]]
      if (length(neigh) > 0) queue <- c(queue, neigh)
    }

    ## Se não alcançou todas as células do fragmento, ele está partido
    ## If not all cells were reached, the fragment is split
    if (reached != length(cells)) return(TRUE)
  }
  FALSE
}



## 1) Função de simulação /// 1) Simulation function 

simular_fragmentacao <- function(X, H, F, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  ## Validação de entrada ///  Input validation
  if (H < 0 || H > X^2) stop("H deve estar entre 0 e X^2.") # H must be between 0 and X^2.
  if (F <= 0 || F != as.integer(F)) stop("F deve ser um inteiro positivo.") # F must be a positive integer.
  if (H < F) stop("É necessário H >= F (pelo menos 1 célula por fragmento).") # H must be >= F (at least one cell per fragment).


  ## Inicializa as matrizes da paisagem
  ## Initialize landscape matrices
  ## -----------------------------
  habitat     <- matrix(0L, nrow = X, ncol = X)  # 0 = não-habitat, 1 = habitat
  ## habitat     <- matrix(0L, nrow = X, ncol = X)  # 0 = non-habitat, 1 = habitat

  fragment_id <- matrix(0L, nrow = X, ncol = X)  # 0 = sem fragmento; 1..F = IDs dos fragmentos
  ## fragment_id <- matrix(0L, nrow = X, ncol = X)  # 0 = none; 1..F = fragment IDs

  ## -----------------------------
  ## 1) Sortear tamanhos alvo dos fragmentos
  ## -----------------------------
  ## 1) Sample target fragment sizes
  ## -----------------------------
  tamanhos <- rep(1L, F)   # pelo menos 1 célula por fragmento
  ## tamanhos <- rep(1L, F)   # at least 1 cell per fragment

  restante <- H - F        # habitat restante para distribuir
  ## restante <- H - F        # remaining habitat cells to distribute

  if (restante > 0) {
    extra <- as.vector(rmultinom(1, size = restante, prob = rep(1, F)))
    tamanhos <- tamanhos + extra
  }

  ## -----------------------------
  ## Vizinhança de Moore (8-vizinhos) e bordas rígidas (não-periódicas)
  ## -----------------------------
  ## Moore neighborhood (8-neighbor) and hard-wall (non-periodic) boundaries
  ## -----------------------------
  dr <- c(-1,-1,-1, 0,0, 1,1,1)
  dc <- c(-1, 0, 1,-1,1,-1,0,1)

  ## Verifica se (r,c) encosta em fragmento de outro ID (Moore)
  ## Checks whether (r,c) touches another fragment ID (Moore)
  encosta_outro_fragmento_rc <- function(r, c, f) {
    rs <- r + dr; cs <- c + dc
    ok <- rs >= 1 & rs <= X & cs >= 1 & cs <= X
    if (!any(ok)) return(FALSE)
    viz <- fragment_id[cbind(rs[ok], cs[ok])]
    any(viz != 0L & viz != f)
  }

  ## Retorna a fronteira de Moore (8-vizinhos) de um conjunto de células (r,c)
  ## Returns the Moore frontier (8-neighbor) of a set of (r,c) cells
  fronteira_moore <- function(cells_rc) {
    if (nrow(cells_rc) == 0) return(matrix(NA_integer_, 0, 2))
    out <- matrix(NA_integer_, nrow = 0, ncol = 2)

    for (k in seq_len(nrow(cells_rc))) {
      r <- cells_rc[k, 1]; c <- cells_rc[k, 2]
      rs <- r + dr; cs <- c + dc
      ok <- rs >= 1 & rs <= X & cs >= 1 & cs <= X
      if (any(ok)) out <- rbind(out, cbind(rs[ok], cs[ok]))
    }

    if (nrow(out) == 0) return(out)
    unique(out)
  }

  ## -----------------------------
  ## Alocar e crescer os fragmentos
  ## -----------------------------
  ## Allocate and grow fragments
  ## -----------------------------
  frag_cells <- vector("list", F)

  for (f in seq_len(F)) {
    alvo <- tamanhos[f]
    if (alvo <= 0) next

    ## ---- Posicionar semente do fragmento ----
    ## ---- Seed placement ----
    vazias_rc <- which(habitat == 0L, arr.ind = TRUE)

    if (nrow(vazias_rc) == 0) {
      warning("Acabaram as células vazias antes de alocar todo o habitat.") # Empty cells ran out before all habitat was allocated.
      break
    }

    ok_seed <- apply(vazias_rc, 1, function(rc) !encosta_outro_fragmento_rc(rc[1], rc[2], f))
    cand_seed <- vazias_rc[ok_seed, , drop = FALSE]

    if (nrow(cand_seed) == 0) {
      warning(paste("Não foi possível posicionar a semente do fragmento", f,
                    "- sem espaço disponível sem encostar em outros fragmentos.")) # Could not place the seed of fragment f, no available space without touching other fragments.
      next
    }

    pick <- cand_seed[sample.int(nrow(cand_seed), 1), ]
    r0 <- pick[1]; c0 <- pick[2]

    habitat[r0, c0]     <- 1L
    fragment_id[r0, c0] <- as.integer(f)
    frag_cells[[f]]     <- matrix(c(r0, c0), ncol = 2)

    ## ---- Crescimento do fragmento ----
    ## ---- Fragment growth ----
    while (nrow(frag_cells[[f]]) < alvo) {
      front_rc <- fronteira_moore(frag_cells[[f]])
      if (nrow(front_rc) == 0) {
        warning(paste("Fragmento", f,
                      "não conseguiu atingir o tamanho alvo.",
                      "Alocado:", nrow(frag_cells[[f]]),
                      "Alvo:", alvo)) # Fragment f could not reach the target size, Allocated: ..., Target: ...
        break
      }

      ## Manter apenas células vazias
      ## Keep only empty cells
      vazias_mask <- habitat[cbind(front_rc[,1], front_rc[,2])] == 0L
      front_rc <- front_rc[vazias_mask, , drop = FALSE]
      if (nrow(front_rc) == 0) {
        warning(paste("Fragmento", f,
                      "não conseguiu atingir o tamanho alvo.",
                      "Alocado:", nrow(frag_cells[[f]]),
                      "Alvo:", alvo)) # Fragment f could not reach the target size, Allocated: ..., Target: ...
        break
      }

      ## Filtrar células vazias da fronteira que não encostam em outros fragmentos
      ## Filter empty frontier cells that do not touch other fragments
      ok_mask <- apply(front_rc, 1, function(rc) !encosta_outro_fragmento_rc(rc[1], rc[2], f))
      front_ok <- front_rc[ok_mask, , drop = FALSE]
      if (nrow(front_ok) == 0) {
        warning(paste("Fragmento", f,
                      "não conseguiu crescer sem encostar em outros fragmentos.",
                      "Alocado:", nrow(frag_cells[[f]]),
                      "Alvo:", alvo)) # Fragment f could not grow without touching other fragments, Allocated: ..., Target: ...
        break
      }

      pick2 <- front_ok[sample.int(nrow(front_ok), 1), ]
      r1 <- pick2[1]; c1 <- pick2[2]

      habitat[r1, c1]     <- 1L
      fragment_id[r1, c1] <- as.integer(f)
      frag_cells[[f]]     <- rbind(frag_cells[[f]], cbind(r1, c1))
    }
  }

  ## -----------------------------
  ## Checagens finais (estritas)
  ## -----------------------------
  ## Final checks (strict)
  ## -----------------------------
  ids <- sort(setdiff(unique(as.vector(fragment_id)), 0L))
  tamanhos_reais <- vapply(seq_len(F), function(ff) sum(fragment_id == ff), integer(1L))

  success <- (sum(habitat) == H) &&
             (length(ids) == F) &&
             all(tamanhos_reais == tamanhos)

  list(
    habitat         = habitat,
    fragment_id     = fragment_id,
    tamanhos_alvo   = tamanhos,
    tamanhos_reais  = tamanhos_reais,
    success         = success
  )
}


## -------------------------------
## 2) Distância média ao fragmento mais próximo
## -------------------------------
## 2) Mean distance to the nearest fragment
## -------------------------------
dist_media_vizinho_mais_proximo <- function(fragment_id) {
  ids <- sort(setdiff(unique(as.vector(fragment_id)), 0))
  n   <- length(ids)
  if (n < 2) return(NA_real_)

  coords <- lapply(ids, function(id) which(fragment_id == id, arr.ind = TRUE))
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
## 3) Rodar ensemble por combinação
## -------------------------------
## 3) Run ensemble per parameter combination
## -------------------------------

## Defina suas combinações aqui (pode ser 1 valor só)
## Define your parameter combinations here (can be single values)
X_vals <- 50
H_vals <- 125
F_vals <- 2

N_VALIDAS <- 100
MAX_TENTATIVAS <- 20000

## Data frame consolidado de todas as combinações
## Consolidated data frame for all combinations
resultados_todos <- data.frame()

for (X in X_vals) {
  for (H in H_vals) {
    for (F in F_vals) {

      cat("\n====================================================\n")
      cat("Combinação: X =", X, "| H =", H, "| F =", F, "\n")
      cat("Objetivo:", N_VALIDAS, "simulações válidas | Máx:", MAX_TENTATIVAS, "tentativas\n") #number of valid simulations, number of tries
      cat("====================================================\n")
      ## Parameter combination header and run limits

      dist_nn     <- numeric(0)
      n_frag_real <- integer(0)
      n_hab_real  <- integer(0)

      cont <- 0
      tentativas <- 0

      ## Observação: com a função corrigida, os fragmentos crescem apenas pela fronteira (Moore),
      ## portanto são conectados por construção; além disso, exigimos success == TRUE (H, F e tamanhos alvo).
      ## Note: with the corrected function, fragments grow only from the Moore frontier,
      ## therefore they are connected by construction; additionally, we require success == TRUE (H, F, and target sizes).

      while (cont < N_VALIDAS && tentativas < MAX_TENTATIVAS) {
        tentativas <- tentativas + 1

        cat("Tentativa", tentativas, "| válidas:", cont, "de", N_VALIDAS, "\n")
        ## Attempt number and how many valid replicates have been obtained

        res <- simular_fragmentacao(X, H, F)

        ## Só conta se atingiu exatamente H e F (e tamanhos alvo)
        ## Only count if exactly H and F were reached (and target sizes)
        if (!isTRUE(res$success)) next

        cont <- cont + 1
        frag <- res$fragment_id
        hab  <- res$habitat

        dist_nn[cont]     <- dist_media_vizinho_mais_proximo(frag)
        n_frag_real[cont] <- length(setdiff(unique(as.vector(frag)), 0))
        n_hab_real[cont]  <- sum(hab)

        cat("Simulação válida", cont, "obtida | dist_nn =", dist_nn[cont], "\n")
        ## Valid replicate obtained and its distance metric

        ## Testes diagnósticos (não usados para filtrar as simulações válidas)
        ## Diagnostic tests (not used to filter valid simulations)
        cat("Contato entre fragmentos (Moore)? ", tem_contato_moore(frag), "\n") # Do fragments touch (Moore)?
        cat("Algum fragmento está partido (Moore)? ", tem_fragmento_partido(frag), "\n") # Is any fragment split (Moore)?
      }

      if (cont < N_VALIDAS) {
        warning(paste("Foram obtidas apenas", cont, "simulações válidas para X =", X, "H =", H, "F =", F,
                      "em", MAX_TENTATIVAS, "tentativas.")) # Only cont valid simulations were obtained within MAX_TENTATIVAS attempts.
      }

      df <- data.frame(
        X = X,
        H = H,
        F = F,
        simulacao                       = 1:cont,
        dist_media_vizinho_mais_proximo = dist_nn,
        n_fragmentos                    = n_frag_real,
        n_habitat_real                  = n_hab_real,
        tentativas_total                = tentativas
      )

      resultados_todos <- rbind(resultados_todos, df)

      ## Salvar um CSV por combinação
      ## Save one CSV per combination
      nome_arquivo <- paste0("resultados_fragmentacao_X", X, "_H", H, "_F", F, ".csv")

      write.table(
        df,
        file      = nome_arquivo,
        sep       = ";",
        dec       = ",",
        row.names = FALSE,
        col.names = TRUE,
        qmethod   = "double"
      )

      cat("\nFinalizado para X =", X, "H =", H, "F =", F, "\n") # Finished this parameter combination
      cat("Réplicas válidas obtidas:", cont, "de", N_VALIDAS, "| Tentativas:", tentativas, "de", MAX_TENTATIVAS, "\n")
      cat("Arquivo salvo como:", nome_arquivo, "\n") # Output filename

      ## Estatísticas finais da combinação
      ## Final statistics for this combination
      media_dist  <- mean(dist_nn, na.rm = TRUE)
      desvio_dist <- sd(dist_nn, na.rm = TRUE)

      cat("\n================ RESULTADOS FINAIS ================\n")
      cat("Média da distância mínima:", media_dist, "\n")
      cat("Desvio padrão das distâncias:", desvio_dist, "\n")
      cat("Número de simulações bem-sucedidas:", cont, "\n")
      cat("====================================================\n")
      ## Final summary statistics for this combination, including average and standard deviation of shortest path to nearest neighbor and the number of simulations that were successful
    }
  }
}

## Salvar CSV consolidado
## Save consolidated CSV
write.table(
  resultados_todos,
  file      = "resultados_fragmentacao_TODAS_COMBINACOES.csv",
  sep       = ";",
  dec       = ",",
  row.names = FALSE,
  col.names = TRUE,
  qmethod   = "double"
)

cat("\nTudo finalizado. Arquivo consolidado salvo como: resultados_fragmentacao_TODAS_COMBINACOES.csv\n")
## All runs finished, consolidated file saved