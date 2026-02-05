## -------------------------------
## Simulação de paisagem fragmentada
## - Conectividade 8-vizinhos
## - Fragmentos não encostam entre si
## -------------------------------

simular_fragmentacao <- function(X, H, F, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Checagens
  if (H < 0 || H > X^2) {
    stop("H deve estar entre 0 e X^2.")
  }
  if (F <= 0 || F != as.integer(F)) {
    stop("F deve ser um inteiro positivo.")
  }
  if (H < F) {
    stop("É necessário H >= F (pelo menos 1 célula por fragmento).")
  }
  
  # Matrizes
  habitat     <- matrix(0, nrow = X, ncol = X)  # 0 = não-habitat, 1 = habitat
  fragment_id <- matrix(0, nrow = X, ncol = X)  # 0 = sem fragmento
  
  ## -----------------------------
  ## 1) Sortear tamanhos de fragmentos
  ## -----------------------------
  tamanhos <- rep(1, F)
  restante <- H - F
  
  if (restante > 0) {
    extra <- as.vector(rmultinom(1, size = restante, prob = rep(1, F)))
    tamanhos <- tamanhos + extra
  }
  
  ## -----------------------------
  ## Funções auxiliares
  ## -----------------------------
  
  # índice linear -> (linha, coluna)
  idx_to_rc <- function(idx, X) {
    r <- ((idx - 1) %% X) + 1
    c <- ((idx - 1) %/% X) + 1
    list(r = r, c = c)
  }
  
  # (linha, coluna) -> índice linear
  rc_to_idx <- function(r, c, X) {
    (c - 1) * X + r
  }
  
  # Vizinhos de 8 direções
  vizinhos_8 <- function(cells, X) {
    rc <- idx_to_rc(cells, X)
    r  <- rc$r
    c  <- rc$c
    
    rs <- c(r-1, r-1, r-1, r,   r,   r+1, r+1, r+1)
    cs <- c(c-1, c,   c+1, c-1, c+1, c-1, c,   c+1)
    
    dentro <- which(rs >= 1 & rs <= X & cs >= 1 & cs <= X)
    if (length(dentro) == 0) return(integer(0))
    
    idx <- rc_to_idx(rs[dentro], cs[dentro], X)
    unique(idx)
  }
  
  # Verifica se uma célula encosta em fragmento de outro ID (8-vizinhos)
  encosta_outro_fragmento <- function(idx, f, fragment_id, X) {
    rc <- idx_to_rc(idx, X)
    r  <- rc$r
    c  <- rc$c
    
    rs <- c(r-1, r-1, r-1, r,   r,   r+1, r+1, r+1)
    cs <- c(c-1, c,   c+1, c-1, c+1, c-1, c,   c+1)
    
    dentro <- which(rs >= 1 & rs <= X & cs >= 1 & cs <= X)
    if (length(dentro) == 0) return(FALSE)
    
    viz_ids <- fragment_id[cbind(rs[dentro], cs[dentro])]
    any(viz_ids != 0 & viz_ids != f)
  }
  
  ## -----------------------------
  ## 2) Crescer fragmentos
  ## -----------------------------
  frag_cells <- vector("list", F)
  
  for (f in seq_len(F)) {
    tamanho_alvo <- tamanhos[f]
    if (tamanho_alvo == 0) next
    
    ## ---- Semente do fragmento ----
    celulas_vazias <- which(habitat == 0)
    if (length(celulas_vazias) == 0) {
      warning("Acabaram as células vazias antes de alocar todo o habitat.")
      break
    }
    
    # Filtra células que não encostam em outros fragmentos
    candidatos_seed <- celulas_vazias[
      !vapply(
        celulas_vazias,
        function(idx) encosta_outro_fragmento(idx, f, fragment_id, X),
        logical(1L)
      )
    ]
    
    if (length(candidatos_seed) == 0) {
      warning(paste("Não foi possível posicionar a semente do fragmento", f,
                    "- sem espaço sem encostar em outros fragmentos."))
      next
    }
    
    seed_idx <- sample(candidatos_seed, 1)
    rc_seed  <- idx_to_rc(seed_idx, X)
    
    habitat[rc_seed$r, rc_seed$c]     <- 1
    fragment_id[rc_seed$r, rc_seed$c] <- f
    frag_cells[[f]] <- seed_idx
    
    ## ---- Crescimento do fragmento ----
    while (length(frag_cells[[f]]) < tamanho_alvo) {
      neigh <- vizinhos_8(frag_cells[[f]], X)
      
      celulas_vazias <- which(habitat == 0)
      fronteira <- intersect(neigh, celulas_vazias)
      
      if (length(fronteira) == 0) {
        warning(paste(
          "Fragmento", f,
          "não conseguiu atingir o tamanho alvo.",
          "Alocado:", length(frag_cells[[f]]),
          "Alvo:", tamanho_alvo
        ))
        break
      }
      
      # Filtra células de fronteira que não encostam em outros fragmentos
      fronteira_ok <- fronteira[
        !vapply(
          fronteira,
          function(idx) encosta_outro_fragmento(idx, f, fragment_id, X),
          logical(1L)
        )
      ]
      
      if (length(fronteira_ok) == 0) {
        warning(paste(
          "Fragmento", f,
          "não conseguiu crescer sem encostar em outros fragmentos.",
          "Alocado:", length(frag_cells[[f]]),
          "Alvo:", tamanho_alvo
        ))
        break
      }
      
      new_idx <- sample(fronteira_ok, 1)
      rc_new  <- idx_to_rc(new_idx, X)
      
      habitat[rc_new$r, rc_new$c]     <- 1
      fragment_id[rc_new$r, rc_new$c] <- f
      frag_cells[[f]] <- c(frag_cells[[f]], new_idx)
    }
  }
  
  list(
    habitat = habitat,
    fragment_id = fragment_id,
    tamanhos_fragmentos = tamanhos
  )
}

## -------------------------------
## 3) Rodar simulação e plotar paisagem
## -------------------------------

# Parâmetros da paisagem
X <- 50      # tamanho do lattice (X x X)
H <- 250     # total de células de habitat
F <- 2      # número de fragmentos

resultado <- simular_fragmentacao(X, H, F)

# Paisagem final: 0 = matriz sem habitat, 1 = habitat
paisagem <- resultado$habitat

# Plot simples da paisagem
# (0 = branco, 1 = verde)
paleta <- c("white", "darkgreen")

image(
  t(apply(paisagem, 2, rev)),  # truque para não ficar "de cabeça para baixo"
  col  = paleta,
  axes = FALSE,
 
)
box()

## -------------------------------
## Plot colorido por fragmento
## -------------------------------

fragmentos <- resultado$fragment_id

# Paleta: branco + cores distintas
paleta_fragmentos <- c("white", rainbow(max(fragmentos)))

image(
  t(apply(fragmentos, 2, rev)),
  col = paleta_fragmentos,
  axes = FALSE,
  main = "Fragmentos de habitat (cores distintas)"
)
box()

# Contagem de células de habitat
num_celulas_habitat <- sum(resultado$habitat)
cat("Total de células preenchidas:", num_celulas_habitat, "\n")