# ============================================================
# SCRIPT 04
# Simulação recorrente com DHs, testcrosses e reciclagem de elites
# ============================================================

library(dplyr)
set.seed(123)

# ------------------------------------------------------------
# 1. Carregar entradas
# ------------------------------------------------------------

sim_inputs_real <- readRDS("sim_inputs_real.rds")

geno_base <- sim_inputs_real$geno
mapa <- sim_inputs_real$mapa
MarkEff <- sim_inputs_real$MarkEff

geno_base <- as.matrix(geno_base)
MarkEff <- as.matrix(MarkEff)

cruzamentos_B_top50 <- readRDS("cruzamentos_programa_B_top50_por_cenario.rds")
testadores_A <- readRDS("testadores_grupo_A.rds")

# ------------------------------------------------------------
# 2. Parâmetros da simulação
# ------------------------------------------------------------

n_cycles <- 6
n_crosses <- 50
n_DH_por_cruzamento <- 20

n_TC1 <- 800
n_TC2 <- 200
n_TC3 <- 100
n_elite <- 2

traitNames <- c("AP", "AE", "FF", "FM")
cenarios <- c("AP", "AE", "FF", "FM", "MT", "Culling")

# Herdabilidades reais dos caracteres
h2_traits <- c(
  AP = 0.90,
  AE = 0.91,
  FF = 0.89,
  FM = 0.89
)

# Fator de precisão por etapa do funil de testcross
h2_stage_factor <- c(
  TC1 = 0.60,
  TC2 = 0.80,
  TC3 = 0.95
)

get_h2_stage <- function(scenario, stage) {
  
  if (scenario %in% c("AP", "AE", "FF", "FM")) {
    h2_base <- h2_traits[scenario]
  } else {
    h2_base <- mean(h2_traits)
  }
  
  h2_final <- h2_base * h2_stage_factor[stage]
  
  return(as.numeric(h2_final))
}

# ------------------------------------------------------------
# 3. Conferências iniciais
# ------------------------------------------------------------

cat("\nDimensão geno_base:\n")
print(dim(geno_base))

cat("\nDimensão MarkEff:\n")
print(dim(MarkEff))

cat("\nNúmero de cruzamentos iniciais por cenário:\n")
print(table(cruzamentos_B_top50$scenario))

cat("\nTestadores do Grupo A:\n")
print(testadores_A)

cat("\nAlinhamento geno x mapa:\n")
print(all(colnames(geno_base) == mapa$marker))

cat("\nAlinhamento geno x efeitos:\n")
print(all(colnames(geno_base) == rownames(MarkEff)))

if (!"gen_pos_M" %in% colnames(mapa)) {
  stop("A coluna gen_pos_M não foi encontrada no mapa. Rode o SCRIPT 00 atualizado.")
}

# ------------------------------------------------------------
# 4. Organizar mapa para recombinação
# ------------------------------------------------------------

mapa2 <- mapa %>%
  dplyr::arrange(chr, pos)

geno_base <- geno_base[, mapa2$marker]
MarkEff <- MarkEff[mapa2$marker, ]

cat("\nAlinhamento após ordenar pelo mapa:\n")
print(all(colnames(geno_base) == mapa2$marker))
print(all(colnames(geno_base) == rownames(MarkEff)))

chr_index <- split(seq_len(nrow(mapa2)), mapa2$chr)

# ------------------------------------------------------------
# 5. Funções auxiliares
# ------------------------------------------------------------

get_score <- function(df, scenario) {
  
  if (scenario == "AP") {
    return(df$AP)
  }
  
  if (scenario == "AE") {
    return(df$AE)
  }
  
  if (scenario == "FF") {
    return(df$FF)
  }
  
  if (scenario == "FM") {
    return(df$FM)
  }
  
  if (scenario %in% c("MT", "Culling")) {
    return(rowMeans(df[, traitNames]))
  }
  
  stop("Cenário não reconhecido: ", scenario)
}

add_field_error <- function(gv, h2) {
  
  vg <- var(gv, na.rm = TRUE)
  
  if (is.na(vg) || vg == 0) {
    ve <- 1
  } else {
    ve <- vg * (1 - h2) / h2
  }
  
  pheno <- gv + rnorm(
    n = length(gv),
    mean = 0,
    sd = sqrt(ve)
  )
  
  return(pheno)
}

simulate_gamete_chr <- function(parent1_chr,
                                parent2_chr,
                                gen_pos_chr) {
  
  hap1 <- parent1_chr / 2
  hap2 <- parent2_chr / 2
  
  chr_length <- max(gen_pos_chr, na.rm = TRUE)
  
  n_rec <- rpois(1, lambda = chr_length)
  
  if (n_rec == 0) {
    
    start_hap <- sample(c(1, 2), size = 1)
    
    if (start_hap == 1) {
      return(hap1)
    } else {
      return(hap2)
    }
  }
  
  breakpoints <- sort(runif(n_rec, min = 0, max = chr_length))
  n_switch <- findInterval(gen_pos_chr, breakpoints)
  
  start_hap <- sample(c(1, 2), size = 1)
  
  if (start_hap == 1) {
    use_hap1 <- n_switch %% 2 == 0
  } else {
    use_hap1 <- n_switch %% 2 == 1
  }
  
  gamete <- ifelse(use_hap1, hap1, hap2)
  
  return(gamete)
}

simulate_DH_cross <- function(parent1_id,
                              parent2_id,
                              geno_pool,
                              mapa2,
                              chr_index,
                              n_DH = 20) {
  
  if (!(parent1_id %in% rownames(geno_pool))) {
    stop("Parent1 não encontrado: ", parent1_id)
  }
  
  if (!(parent2_id %in% rownames(geno_pool))) {
    stop("Parent2 não encontrado: ", parent2_id)
  }
  
  p1 <- geno_pool[parent1_id, ]
  p2 <- geno_pool[parent2_id, ]
  
  DH_mat <- matrix(
    NA_real_,
    nrow = n_DH,
    ncol = ncol(geno_pool)
  )
  
  colnames(DH_mat) <- colnames(geno_pool)
  
  for (i in seq_len(n_DH)) {
    
    gamete_full <- numeric(ncol(geno_pool))
    
    for (chr_i in names(chr_index)) {
      
      idx <- chr_index[[chr_i]]
      
      gamete_chr <- simulate_gamete_chr(
        parent1_chr = p1[idx],
        parent2_chr = p2[idx],
        gen_pos_chr = mapa2$gen_pos_M[idx]
      )
      
      gamete_full[idx] <- gamete_chr
    }
    
    DH_mat[i, ] <- gamete_full * 2
  }
  
  rownames(DH_mat) <- paste(
    parent1_id,
    parent2_id,
    paste0("DH", seq_len(n_DH)),
    sep = "_"
  )
  
  return(DH_mat)
}

simulate_DHs_plan <- function(plan,
                              geno_pool,
                              mapa2,
                              chr_index,
                              scenario,
                              cycle,
                              n_DH_por_cruzamento = 20) {
  
  DH_list <- list()
  info_list <- list()
  
  for (i in seq_len(nrow(plan))) {
    
    if (i %% 10 == 0) {
      cat("  Cruzamento", i, "de", nrow(plan), "\n")
    }
    
    p1 <- plan$Parent1[i]
    p2 <- plan$Parent2[i]
    
    DH_i <- simulate_DH_cross(
      parent1_id = p1,
      parent2_id = p2,
      geno_pool = geno_pool,
      mapa2 = mapa2,
      chr_index = chr_index,
      n_DH = n_DH_por_cruzamento
    )
    
    rownames(DH_i) <- paste(
      scenario,
      paste0("C", cycle),
      plan$cross_id[i],
      paste0("DH", seq_len(n_DH_por_cruzamento)),
      sep = "_"
    )
    
    DH_list[[i]] <- DH_i
    
    info_list[[i]] <- data.frame(
      scenario = scenario,
      cycle = cycle,
      cross_id = plan$cross_id[i],
      Parent1 = p1,
      Parent2 = p2,
      DH_ID = rownames(DH_i),
      stringsAsFactors = FALSE
    )
  }
  
  DH_geno <- do.call(rbind, DH_list)
  DH_info <- do.call(rbind, info_list)
  
  return(
    list(
      DH_geno = DH_geno,
      DH_info = DH_info
    )
  )
}

make_cross_candidates_base_elite <- function(base_parentals,
                                             elite_parentals) {
  
  cand <- expand.grid(
    Parent1 = base_parentals,
    Parent2 = elite_parentals,
    stringsAsFactors = FALSE
  )
  
  cand <- cand %>%
    dplyr::filter(Parent1 != Parent2) %>%
    dplyr::mutate(
      cross_id = paste(Parent1, Parent2, sep = "_")
    ) %>%
    dplyr::distinct(cross_id, .keep_all = TRUE)
  
  return(cand)
}

score_and_select_crosses <- function(candidates,
                                     geno_pool,
                                     MarkEff,
                                     scenario,
                                     n_crosses = 50) {
  
  parent_ids <- unique(c(candidates$Parent1, candidates$Parent2))
  
  geno_parents <- geno_pool[parent_ids, , drop = FALSE]
  
  gv_parents <- geno_parents %*% MarkEff
  gv_parents <- as.data.frame(gv_parents)
  colnames(gv_parents) <- traitNames
  gv_parents$Parent <- rownames(geno_parents)
  
  p1 <- match(candidates$Parent1, gv_parents$Parent)
  p2 <- match(candidates$Parent2, gv_parents$Parent)
  
  mp <- (gv_parents[p1, traitNames] + gv_parents[p2, traitNames]) / 2
  
  cross_scores <- cbind(candidates, mp)
  
  if (scenario == "Culling") {
    
    filtered <- NULL
    
    for (q in c(0.50, 0.40, 0.30, 0.20, 0.10, 0.00)) {
      
      lim_AP <- quantile(cross_scores$AP, q, na.rm = TRUE)
      lim_AE <- quantile(cross_scores$AE, q, na.rm = TRUE)
      lim_FF <- quantile(cross_scores$FF, q, na.rm = TRUE)
      lim_FM <- quantile(cross_scores$FM, q, na.rm = TRUE)
      
      temp <- cross_scores %>%
        dplyr::filter(
          AP >= lim_AP,
          AE >= lim_AE,
          FF >= lim_FF,
          FM >= lim_FM
        )
      
      if (nrow(temp) >= n_crosses) {
        filtered <- temp
        break
      }
    }
    
    if (is.null(filtered)) {
      filtered <- cross_scores
    }
    
    cross_scores <- filtered
    cross_scores$score <- rowMeans(cross_scores[, traitNames])
    
  } else {
    
    cross_scores$score <- get_score(cross_scores, scenario)
  }
  
  selected_crosses <- cross_scores %>%
    dplyr::arrange(dplyr::desc(score)) %>%
    dplyr::slice_head(n = n_crosses) %>%
    dplyr::select(Parent1, Parent2, cross_id, AP, AE, FF, FM, score)
  
  return(selected_crosses)
}

evaluate_testcross_stage <- function(DH_geno,
                                     DH_info,
                                     tester_id,
                                     geno_pool,
                                     MarkEff,
                                     scenario,
                                     h2,
                                     stage_name,
                                     n_select) {
  
  if (!(tester_id %in% rownames(geno_pool))) {
    stop("Testador não encontrado na matriz genotípica: ", tester_id)
  }
  
  tester_geno <- geno_pool[tester_id, ]
  
  hybrid_geno <- sweep(
    DH_geno,
    2,
    tester_geno,
    FUN = "+"
  ) / 2
  
  gv_hybrid <- hybrid_geno %*% MarkEff
  gv_hybrid <- as.data.frame(gv_hybrid)
  colnames(gv_hybrid) <- traitNames
  
  gv_hybrid$DH_ID <- rownames(DH_geno)
  gv_hybrid$genetic_score <- get_score(gv_hybrid, scenario)
  
  gv_hybrid$phenotype_score <- add_field_error(
    gv = gv_hybrid$genetic_score,
    h2 = h2
  )
  
  eval_stage <- DH_info %>%
    dplyr::left_join(
      gv_hybrid %>%
        dplyr::select(DH_ID, AP, AE, FF, FM, genetic_score, phenotype_score),
      by = "DH_ID"
    ) %>%
    dplyr::mutate(
      tester = tester_id,
      stage = stage_name
    )
  
  selected <- eval_stage %>%
    dplyr::arrange(dplyr::desc(phenotype_score)) %>%
    dplyr::slice_head(n = n_select)
  
  DH_geno_selected <- DH_geno[selected$DH_ID, , drop = FALSE]
  
  return(
    list(
      selected_info = selected,
      selected_geno = DH_geno_selected,
      eval_all = eval_stage
    )
  )
}

# ------------------------------------------------------------
# 6. Definir testadores
# ------------------------------------------------------------

testers <- testadores_A$Linhagem

if (length(testers) < 3) {
  stop("É necessário ter 3 testadores no arquivo testadores_grupo_A.rds")
}

tester1 <- testers[1]
tester2 <- testers[2]
tester3 <- testers[3]

cat("\nTestadores usados na simulação:\n")
cat("TC1:", tester1, "\n")
cat("TC2:", tester2, "\n")
cat("TC3:", tester3, "\n")

# ------------------------------------------------------------
# 7. Rodar simulação recorrente para os 6 cenários
# ------------------------------------------------------------

all_cycle_summary <- list()
all_elites <- list()
all_TC_summary <- list()
all_crosses_used <- list()

for (sc in cenarios) {
  
  cat("\n###################################################\n")
  cat("CENÁRIO:", sc, "\n")
  cat("###################################################\n")
  
  geno_pool <- geno_base
  
  plan_initial <- cruzamentos_B_top50 %>%
    dplyr::filter(scenario == sc) %>%
    dplyr::arrange(dplyr::desc(Y)) %>%
    dplyr::slice_head(n = n_crosses) %>%
    dplyr::mutate(
      score = Y
    )
  
  base_parentals <- unique(c(plan_initial$Parent1, plan_initial$Parent2))
  
  cat("\nNúmero de parentais base no cenário", sc, ":\n")
  print(length(base_parentals))
  
  elite_accumulated <- character(0)
  
  for (cycle in seq_len(n_cycles)) {
    
    cat("\n-----------------------------------------\n")
    cat("Ciclo", cycle, "- Cenário", sc, "\n")
    cat("-----------------------------------------\n")
    
    if (cycle == 1) {
      
      plan_cycle <- plan_initial %>%
        dplyr::select(Parent1, Parent2, cross_id, score)
      
    } else {
      
      candidates <- make_cross_candidates_base_elite(
        base_parentals = base_parentals,
        elite_parentals = elite_accumulated
      )
      
      cat("Número de cruzamentos candidatos:", nrow(candidates), "\n")
      
      plan_cycle <- score_and_select_crosses(
        candidates = candidates,
        geno_pool = geno_pool,
        MarkEff = MarkEff,
        scenario = sc,
        n_crosses = n_crosses
      )
    }
    
    cat("Número de cruzamentos usados no ciclo:", nrow(plan_cycle), "\n")
    
    all_crosses_used[[paste(sc, cycle, sep = "_C")]] <- plan_cycle %>%
      dplyr::mutate(
        scenario = sc,
        cycle = cycle
      )
    
    # --------------------------------------------------------
    # Gerar DHs
    # --------------------------------------------------------
    
    DH_obj <- simulate_DHs_plan(
      plan = plan_cycle,
      geno_pool = geno_pool,
      mapa2 = mapa2,
      chr_index = chr_index,
      scenario = sc,
      cycle = cycle,
      n_DH_por_cruzamento = n_DH_por_cruzamento
    )
    
    DH_geno <- DH_obj$DH_geno
    DH_info <- DH_obj$DH_info
    
    cat("Número de DHs geradas:", nrow(DH_geno), "\n")
    
    # --------------------------------------------------------
    # TC1: 1000 -> 800
    # --------------------------------------------------------
    
    TC1 <- evaluate_testcross_stage(
      DH_geno = DH_geno,
      DH_info = DH_info,
      tester_id = tester1,
      geno_pool = geno_pool,
      MarkEff = MarkEff,
      scenario = sc,
      h2 = get_h2_stage(sc, "TC1"),
      stage_name = "TC1",
      n_select = n_TC1
    )
    
    # --------------------------------------------------------
    # TC2: 800 -> 200
    # --------------------------------------------------------
    
    TC2 <- evaluate_testcross_stage(
      DH_geno = TC1$selected_geno,
      DH_info = TC1$selected_info %>%
        dplyr::select(scenario, cycle, cross_id, Parent1, Parent2, DH_ID),
      tester_id = tester2,
      geno_pool = geno_pool,
      MarkEff = MarkEff,
      scenario = sc,
      h2 = get_h2_stage(sc, "TC2"),
      stage_name = "TC2",
      n_select = n_TC2
    )
    
    # --------------------------------------------------------
    # TC3: 200 -> 100
    # --------------------------------------------------------
    
    TC3 <- evaluate_testcross_stage(
      DH_geno = TC2$selected_geno,
      DH_info = TC2$selected_info %>%
        dplyr::select(scenario, cycle, cross_id, Parent1, Parent2, DH_ID),
      tester_id = tester3,
      geno_pool = geno_pool,
      MarkEff = MarkEff,
      scenario = sc,
      h2 = get_h2_stage(sc, "TC3"),
      stage_name = "TC3",
      n_select = n_TC3
    )
    
    # --------------------------------------------------------
    # Selecionar top 2 elites do ciclo
    # --------------------------------------------------------
    
    elite_cycle <- TC3$selected_info %>%
      dplyr::arrange(dplyr::desc(phenotype_score)) %>%
      dplyr::slice_head(n = n_elite)
    
    elite_ids <- elite_cycle$DH_ID
    
    cat("Elites selecionadas no ciclo:\n")
    print(elite_ids)
    
    geno_elites <- TC3$selected_geno[elite_ids, , drop = FALSE]
    
    new_elites <- elite_ids[!(elite_ids %in% rownames(geno_pool))]
    
    if (length(new_elites) > 0) {
      geno_pool <- rbind(
        geno_pool,
        geno_elites[new_elites, , drop = FALSE]
      )
    }
    
    elite_accumulated <- unique(c(elite_accumulated, elite_ids))
    
    # --------------------------------------------------------
    # Resumos do ciclo
    # --------------------------------------------------------
    
    TC_summary_cycle <- data.frame(
      scenario = sc,
      cycle = cycle,
      n_parentais_base = length(base_parentals),
      n_elites_acumuladas = length(elite_accumulated),
      n_parentais_pool = nrow(geno_pool),
      n_crosses = nrow(plan_cycle),
      n_DH = nrow(DH_geno),
      n_TC1 = nrow(TC1$selected_info),
      n_TC2 = nrow(TC2$selected_info),
      n_TC3 = nrow(TC3$selected_info),
      h2_TC1 = get_h2_stage(sc, "TC1"),
      h2_TC2 = get_h2_stage(sc, "TC2"),
      h2_TC3 = get_h2_stage(sc, "TC3"),
      mean_TC1 = mean(TC1$selected_info$phenotype_score),
      mean_TC2 = mean(TC2$selected_info$phenotype_score),
      mean_TC3 = mean(TC3$selected_info$phenotype_score),
      mean_elite = mean(elite_cycle$phenotype_score),
      mean_genetic_elite = mean(elite_cycle$genetic_score),
      var_genetic_TC3 = var(TC3$selected_info$genetic_score),
      stringsAsFactors = FALSE
    )
    
    all_cycle_summary[[paste(sc, cycle, sep = "_C")]] <- TC_summary_cycle
    
    all_elites[[paste(sc, cycle, sep = "_C")]] <- elite_cycle
    
    all_TC_summary[[paste(sc, cycle, sep = "_C")]] <- bind_rows(
      TC1$selected_info %>% dplyr::mutate(selection_step = "Selected_TC1"),
      TC2$selected_info %>% dplyr::mutate(selection_step = "Selected_TC2"),
      TC3$selected_info %>% dplyr::mutate(selection_step = "Selected_TC3")
    )
    
    cat("Resumo do ciclo:\n")
    print(TC_summary_cycle)
  }
}

# ------------------------------------------------------------
# 8. Consolidar resultados
# ------------------------------------------------------------

cycle_summary_all <- dplyr::bind_rows(all_cycle_summary)
elites_all <- dplyr::bind_rows(all_elites)
TC_selected_all <- dplyr::bind_rows(all_TC_summary)
crosses_used_all <- dplyr::bind_rows(all_crosses_used)

# ------------------------------------------------------------
# 9. Salvar resultados
# ------------------------------------------------------------

saveRDS(cycle_summary_all, "sim_cycle_summary_recorrente.rds")
saveRDS(elites_all, "sim_elites_recorrente.rds")
saveRDS(TC_selected_all, "sim_TC_selected_recorrente.rds")
saveRDS(crosses_used_all, "sim_crosses_used_recorrente.rds")

