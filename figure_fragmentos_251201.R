# Os comentários e as mensagens impressas ao rodar foram escritos em Português Brasileiro.
# The original code comments and console messages were written in Brazilian Portuguese. English comments were added below.

# Por favor, veja as atualizações em https://github.com/miudolab/fragmentation/tree/main
# Please, check for updates at https://github.com/miudolab/fragmentation/tree/main


# FUNÇÕES AUXILIARES DE VERIFICAÇÃO DA VIZINHANÇA DE MOORE
# AUXILIARY CHECKING FUNCTIONS OF MOORE NEIGHBORHOOD

# Verifica se há contato entre fragmentos diferentes (vizinhança de Moore)
# Checks whether different fragments touch each other (Moore neighborhood)

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

# Verifica se algum fragmento (mesmo ID) está "partido" em mais de um componente
# Checks whether any fragment (same ID) is split into more than one connected component

tem_fragmento_partido <- function(fragment_id) {
  X <- nrow(fragment_id)
  dr <- c(-1,-1,-1, 0,0, 1,1,1)
  dc <- c(-1, 0, 1,-1,1,-1,0,1)

  ids <- sort(setdiff(unique(as.vector(fragment_id)), 0))
  for (f in ids) {
    cells <- which(fragment_id == f)
    if (length(cells) <= 1) next

    # Máscara booleana no índice linear
    # Boolean mask over linear indices
    
    
    cell_set <- rep(FALSE, X * X)
    cell_set[cells] <- TRUE

    # Busca em largura (BFS) a partir da primeira célula
    # Breadth-first search (BFS) starting from the first cell
    queue <- cells[1]
    reached <- 0L

    while (length(queue) > 0) {
      cur <- queue[1]
      queue <- queue[-1]

      if (!cell_set[cur]) next
      cell_set[cur] <- FALSE
      reached <- reached + 1L

      # Índice linear -> (linha, coluna) 
      # Linear index -> (row, column) 
      
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

    # Se não alcançou todas as células do fragmento, ele está partido em dois componentes
    # If not all cells were reached, the fragment is split in two components
    if (reached != length(cells)) return(TRUE)
  }
  FALSE
}


# FUNÇÃO PRINCIPAL: SIMULAÇÃO DE PAISAGEM FRAGMENTADA
# MAIN FUNCTION: FRAGMENTED LANDSCAPE SIMULATION

# Um fragmento é definido por um conjunto de células de habitat contínuas, separado por pelo menos uma célula de matriz de outros fragmentos.
# A fragment is defined as a set of contiguous habitat cells, separated by at least one matrix cell from other fragments.



# Conectividade é definida pela vizinhança de Moore (8-vizinhos).
# Connectivity is defined using the Moore neighborhood (8-neighbor).

simular_fragmentacao <- function(X, H, F, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # X^2 = número de células na paisagem (grade X x X)
  # X^2 = number of cells in the landscape (X x X grid)
  
  # H = número total de células de habitat
  # H = total number of habitat cells
  
  # F = número de fragmentos
  # F = number of fragments

  # -----------------------------
  # Validação de entrada
  # Input validation
  # -----------------------------
  if (H < 0 || H > X^2) stop("H deve estar entre 0 e X^2.") #/// H must be in-between 0 and X^2
  if (F <= 0 || F != as.integer(F)) stop("F deve ser um inteiro positivo.") #/// F must be a positive integer
  if (H < F) stop("É necessário H >= F (pelo menos 1 célula de hábitat por fragmento).") #/// It is mandatory that H>=F (at least one habitat cell per fragment)


  # Inicializa as paisagens
  # Initialize landscapes
  habitat     <- matrix(0L, nrow = X, ncol = X)  # 0 = não-habitat, 1 = habitat   #/// 0 = non-habitat, 1 = habitat

  fragment_id <- matrix(0L, nrow = X, ncol = X)  # 0 = sem fragmento; 1..F = IDs #/// 0 = none; 1..F = IDs

  # 1) Sortear tamanhos dos fragmentos
  # 1) Sample fragment sizes
  tamanhos <- rep(1L, F)   # pelo menos 1 célula por fragmento #/// at least 1 cell per fragment

  restante <- H - F        # habitat restante a distribuir #/// restante <- H - F        # remaining habitat cells to distribute

  if (restante > 0) {
    extra <- as.vector(rmultinom(1, size = restante, prob = rep(1, F)))
    tamanhos <- tamanhos + extra
  }

  # Verificação: vizinhança de Moore e contato com outros fragmentos ///Moore neighborhood and contact check with other fragments
  # Condições não-periódicas de contorno (bordas rígidas) /// Hard-wall (non-periodic) boundary conditions
 
  encosta_outro_fragmento_rc <- function(r, c, f) {
    rs <- c(r-1, r-1, r-1, r,   r,   r+1, r+1, r+1)
    cs <- c(c-1, c,   c+1, c-1, c+1, c-1, c,   c+1)

    dentro <- rs >= 1 & rs <= X & cs >= 1 & cs <= X
    if (!any(dentro)) return(FALSE)

    viz_ids <- fragment_id[cbind(rs[dentro], cs[dentro])]
    any(viz_ids != 0L & viz_ids != f)
  }

  # Borda de Moore para um conjunto de células (r,c)
  # Moore frontier for a set of (r,c) cells
  
  fronteira_moore <- function(cells_rc) {
    dr <- c(-1,-1,-1, 0,0, 1,1,1)
    dc <- c(-1, 0, 1,-1,1,-1,0,1)

    out <- matrix(NA_integer_, nrow = 0, ncol = 2)
    if (nrow(cells_rc) == 0) return(out)

    for (k in seq_len(nrow(cells_rc))) {
      r <- cells_rc[k, 1]; c <- cells_rc[k, 2]
      rs <- r + dr; cs <- c + dc
      ok <- rs >= 1 & rs <= X & cs >= 1 & cs <= X
      if (any(ok)) out <- rbind(out, cbind(rs[ok], cs[ok]))
    }

    if (nrow(out) == 0) return(out)
    unique(out)
  }

  # Colocar os fragmentos
  # Allocate fragments
  frag_cells <- vector("list", F)

  for (f in seq_len(F)) {
    tamanho_alvo <- tamanhos[f]
    if (tamanho_alvo <= 0) next

    # Colocar a semente do fragmento
    # Seed placement
    celulas_vazias_rc <- which(habitat == 0L, arr.ind = TRUE)

    if (nrow(celulas_vazias_rc) == 0) {
      warning("Acabaram as células vazias antes de alocar todo o habitat.") #/// There is no remaining empty cells to allocate all habitat cells
      break
    }

    # Filtrar células vazias que não encostam em outros fragmentos /// Filter empty cells that do not touch other fragments
    ok_seed <- apply(celulas_vazias_rc, 1, function(rc) {
      r <- rc[1]; c <- rc[2]
      !encosta_outro_fragmento_rc(r, c, f)
    })
    candidatos_seed_rc <- celulas_vazias_rc[ok_seed, , drop = FALSE]

    if (nrow(candidatos_seed_rc) == 0) {
      warning(paste("Não foi possível posicionar a semente do fragmento", f,
                    "- sem espaço disponível sem encostar em outros fragmentos.")) #Could not place the seed of fragment, no available space that does not touch other fragments.
      next
    }

    pick <- candidatos_seed_rc[sample.int(nrow(candidatos_seed_rc), 1), ]
    r0 <- pick[1]; c0 <- pick[2]

    habitat[r0, c0]     <- 1L
    fragment_id[r0, c0] <- as.integer(f)
    frag_cells[[f]]     <- matrix(c(r0, c0), ncol = 2)

    # Crescimento do fragmento
    # Fragment growth
    while (nrow(frag_cells[[f]]) < tamanho_alvo) {
      front_rc <- fronteira_moore(frag_cells[[f]])

      if (nrow(front_rc) == 0) {
        warning(paste("Fragmento", f, 
                      "não conseguiu atingir o tamanho alvo.", # fragment f Could not reach the target size
                      "Alocado:", nrow(frag_cells[[f]]), #/// Alocated
                      "Alvo:", tamanho_alvo)). #///target
        break
      }

      # Manter apenas as células vazias
      # Keep only empty cells
      vazias_mask <- habitat[cbind(front_rc[,1], front_rc[,2])] == 0L
      front_rc <- front_rc[vazias_mask, , drop = FALSE]

      if (nrow(front_rc) == 0) {
        warning(paste("Fragmento", f,
                      "não conseguiu atingir o tamanho alvo.", # fragment f Could not reach the target size
                      "Alocado:", nrow(frag_cells[[f]]), #/// Alocated
                      "Alvo:", tamanho_alvo)) #///target
        break
      }

      # Filtrar fronteira para não encostar em outros fragmentos
      # Filter frontier so it does not touch other fragments
      ok_mask <- apply(front_rc, 1, function(rc) {
        !encosta_outro_fragmento_rc(rc[1], rc[2], f)
      })
      front_ok <- front_rc[ok_mask, , drop = FALSE]

      if (nrow(front_ok) == 0) {
        warning(paste("Fragmento", f,
                      "não conseguiu crescer sem encostar em outros fragmentos.",  # fragment f Could not reach the target size
                      "Alocado:", nrow(frag_cells[[f]]), #/// Alocated
                      "Alvo:", tamanho_alvo)) #///target
        break
      }

      pick2 <- front_ok[sample.int(nrow(front_ok), 1), ]
      r1 <- pick2[1]; c1 <- pick2[2]

      habitat[r1, c1]     <- 1L
      fragment_id[r1, c1] <- as.integer(f)
      frag_cells[[f]]     <- rbind(frag_cells[[f]], cbind(r1, c1))
    }
  }

  list(
    habitat             = habitat,
    fragment_id         = fragment_id,
    tamanhos_fragmentos = tamanhos
  )
}


# Parâmetros da paisagem
# Landscape parameters
X <- 50
H <- 125
F <- 2
seed <- 0

# Rodar simulação
# Run simulation
res <- simular_fragmentacao(X = X, H = H, F = F, seed = seed)

# Extrair matrizes
# Extract matrices
frag <- res$fragment_id
hab  <- res$habitat

# Verificações (para garantir qualidade da simulação)
# Checks (to ensure simulation quality)
cat("Há contato entre fragmentos (Moore)? ", tem_contato_moore(frag), "\n")
cat("Algum fragmento está partido (Moore)? ", tem_fragmento_partido(frag), "\n")

# Contagem de células de habitat
# Count habitat cells
n_celulas_habitat <- sum(hab)
cat("Número total de células de habitat:", n_celulas_habitat, "\n") #total number of habitat cells

# Contagem do número total de manchas (fragmentos presentes)
# Count total number of patches (fragments present)
n_manchas <- length(setdiff(unique(as.vector(frag)), 0))
cat("Número total de manchas (fragmentos):", n_manchas, "\n") # total number of patches (fragments)

# Gráfico: manchas verdes e matriz branca
# Plot: green habitat patches and white matrix

par(mar = c(2, 2, 3, 2))
image(
  t(apply(hab, 2, rev)),
  col    = c("white", "darkgreen"),
  breaks = c(-0.5, 0.5, 1.5),
  axes   = FALSE,
  main   = sprintf("Manchas de habitat (X=%d, H=%d, F=%d)", X, H, F)
)
box()