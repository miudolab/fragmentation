# The original code comments and console messages were written in Brazilian Portuguese.
# English comments were added prior to submission to Github just below the comments in Brazilian Portuguese


# Título: Simulação de paisagem fragmentada /// Fragmented landscape simulation

# Um fragmento é definido por um conjunto de células de habitat contínuas separadas por pelo menos uma célula de outras células de habitat.
# Conectividade é definida pela vizinhança de Moore

# A fragment is defined as a set of contiguous habitat cells separated by at least one cell from other habitat cells.
# Connectivity is defined using the Moore neighborhood.

## -------------------------------

## Simulação -- simulation 

# Function arguments / Argumentos da função

# X    : Tamanho da paisagem (dimensão da grade X × X) /// Landscape size (grid dimension X × X)
# H    : Número total de células de habitat na paisagem /// Total number of habitat cells in the landscape
# F    : Número de fragmentos de habitat a serem gerados /// Number of habitat fragments to be generated
# seed : Semente aleatória para reprodutibilidade /// Random seed for reproducibility
# ---------------------------------------------------------

simular_fragmentacao <- function(X, H, F, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Checagens /// Imput validation
  if (H < 0 || H > X^2) {
    stop("H deve estar entre 0 e X^2.") #/// H must be between 0 and X^2.
  }
  if (F <= 0 || F != as.integer(F)) {
    stop("F deve ser um inteiro positivo.") #/// F must be a positive integer.
  }
  if (H < F) {
    stop("É necessário H >= F (pelo menos 1 célula por fragmento).") #/// H must be >= F (at least one cell per fragment).
  }
  
  
  # Matrizes /// Matrices
  
  habitat     <- matrix(0, nrow = X, ncol = X)  # 0 = não-habitat, 1 = habitat. /// 0 = nonhabitat, 1 = habitat
  fragment_id <- matrix(0, nrow = X, ncol = X)  # 0 = sem fragmento /// 0 = no fragment
  
  ## 1) Sortear tamanhos de fragmentos /// Sample fragment sizes
  # -----------------------------
  tamanhos <- rep(1, F).  #/// fragment sizes
  restante <- H - F.    #/// amount of habitat that can be distributed across the fragments (after one cell per fragment)
  
  if (restante > 0) {
    extra <- as.vector(rmultinom(1, size = restante, prob = rep(1, F))) 
    tamanhos <- tamanhos + extra #/// fragment sizes
  }
  

 # Funções auxiliares /// Auxiliary functions

  
  # índice linear -> (linha, coluna) /// linear index --> (row, column)
  idx_to_rc <- function(idx, X) {
    r <- ((idx - 1) %% X) + 1
    c <- ((idx - 1) %/% X) + 1
    list(r = r, c = c)
  }
  
  # (linha, coluna) -> índice linear /// (row, column) --> linear index
  rc_to_idx <- function(r, c, X) {
    (c - 1) * X + r
  }
  
  # Vizinhos de 8 direções /// neighbors of 8 directions
  vizinhos_8 <- function(cells, X) {
    rc <- idx_to_rc(cells, X)
    r  <- rc$r
    c  <- rc$c
    
    rs <- c(r-1, r-1, r-1, r,   r,   r+1, r+1, r+1)
    cs <- c(c-1, c,   c+1, c-1, c+1, c-1, c,   c+1)
    
    dentro <- which(rs >= 1 & rs <= X & cs >= 1 & cs <= X) #/// Condições de contorno rígidas (não periódicas) /// hard-wall boundary conditions
    if (length(dentro) == 0) return(integer(0))
    
    idx <- rc_to_idx(rs[dentro], cs[dentro], X)
    unique(idx)
  }
  
  # Verifica se uma célula encosta em fragmento de outro ID (vizinhança de Moore) /// # Checks whether a cell touches a fragment with a different ID (Moore neighborhood)
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
  
  # 2) Crescer fragmentos. /// Fragment growth
  # -----------------------------
  frag_cells <- vector("list", F)
  
  for (f in seq_len(F)) {
    tamanho_alvo <- tamanhos[f]
    if (tamanho_alvo == 0) next
    
    ## ---- Semente do fragmento ----
    celulas_vazias <- which(habitat == 0)
    if (length(celulas_vazias) == 0) {
      warning("Acabaram as células vazias antes de alocar todo o habitat.") #///Empty cells ran out before all habitat was allocated.
      break
    }
    
    # Filtra células que não encostam em outros fragmentos #///Filters cells that do not touch other fragments
    candidatos_seed <- celulas_vazias[
      !vapply(
        celulas_vazias,
        function(idx) encosta_outro_fragmento(idx, f, fragment_id, X),
        logical(1L)
      )
    ]
    
    if (length(candidatos_seed) == 0) {
      warning(paste("Não foi possível posicionar a semente do fragmento", f,
                    "- sem espaço sem encostar em outros fragmentos.")) #///"Could not place the seed of fragment", f,"- no available space without touching other fragments."
      next
    }
    
    seed_idx <- sample(candidatos_seed, 1)
    rc_seed  <- idx_to_rc(seed_idx, X)
    
    habitat[rc_seed$r, rc_seed$c]     <- 1
    fragment_id[rc_seed$r, rc_seed$c] <- f
    frag_cells[[f]] <- seed_idx
    
    ## ---- Crescimento do fragmento ---- /// Fragment growth
    while (length(frag_cells[[f]]) < tamanho_alvo) {
      neigh <- vizinhos_8(frag_cells[[f]], X)
      
      celulas_vazias <- which(habitat == 0)
      fronteira <- intersect(neigh, celulas_vazias)
      
      if (length(fronteira) == 0) {
        warning(paste(
          "Fragmento", f,
          "não conseguiu atingir o tamanho alvo.",
          "Alocado:", length(frag_cells[[f]]),
          "Alvo:", tamanho_alvo #///"could not reach the target size.", allocated, target"
        ))
        break
      }
      
      # Filtra células de fronteira que não encostam em outros fragmentos /// filters edge cells that did not touch in other fragments
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
        )) #/// "Fragment", f, "could not grow without touching other fragments.", "Allocated:", length(frag_cells[[f]]), "Target:",
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
## 3) Rodar simulação e plotar paisagem /// run simulation and plot the landscape
## -------------------------------

# Parâmetros da paisagem /// landscape parameters
X <- 50      # tamanho do lattice (X x X). /// lattice size
H <- 1875     # total de células de habitat ///total number of habitat cells
F <- 2      # número de fragmentos /// number of fragments

resultado <- simular_fragmentacao(X, H, F)

# Paisagem final: 0 = matriz sem habitat, 1 = habitat /// final landscape 0= nonhabitat, 1=habitat
paisagem <- resultado$habitat

# Plot simples da paisagem
# (0 = branco, 1 = verde)
paleta <- c("white", "darkgreen")

image(
  t(apply(paisagem, 2, rev)),  
  col  = paleta,
  axes = FALSE,
 
)
box()

# Plot colorido por fragmento /// using colors to identify different fragments - not used in the manuscript
-------------------------------

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