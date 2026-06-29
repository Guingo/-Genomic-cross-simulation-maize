# ============================================================
# SCRIPT 00
# Organizar entradas reais para a simulação
# ============================================================
rm(list = ls())
library(dplyr)
# ------------------------------------------------------------
# 1. Carregar arquivos
# ------------------------------------------------------------
path_geno <- "C:/DadosR/Dados Milho Tropical/Dados_Milho/Matriz_dos_Marcadores_350.rds"
path_mapa <- "C:/DadosR/Dados Milho Tropical/Dados_Milho/baybnew/map_painel360.txt"

geno <- readRDS(path_geno)
mapa <- read.table( file = path_mapa, header = TRUE, stringsAsFactors = FALSE )
efeitos_combinados <- readRDS("efeitos_combinados.rds")

plan_AP_raw <- readRDS("maxGainPlan_AP.rds")
plan_AE_raw <- readRDS("maxGainPlan_AE.rds")
plan_FF_raw <- readRDS("maxGainPlan_FF.rds")
plan_FM_raw <- readRDS("maxGainPlan_FM.rds")
plan_MT_raw <- readRDS("maxGainPlan_MT.rds")
plan_Culling_raw <- readRDS("Culling_MT.rds")

# ------------------------------------------------------------
# 2. Ajustar tipos do mapa
# ------------------------------------------------------------

mapa$chr <- as.integer(mapa$chr)
mapa$pos <- as.numeric(mapa$pos)

# ------------------------------------------------------------
# 3. Identificar SNPs comuns
# ------------------------------------------------------------

snps_comuns <- Reduce(
  intersect,
  list(
    colnames(geno),
    mapa$marker,
    efeitos_combinados$Marcador
  )
)

cat("\nNúmero de SNPs comuns entre genótipo, mapa e efeitos:\n")
print(length(snps_comuns))

# ------------------------------------------------------------
# 4. Ordenar mapa e alinhar genótipos
# ------------------------------------------------------------

mapa2 <- mapa %>%
  dplyr::filter(marker %in% snps_comuns) %>%
  dplyr::arrange(chr, pos)

geno2 <- geno[, mapa2$marker]
geno2 <- as.matrix(geno2)

cat("\nAlinhamento geno x mapa:\n")
print(all(colnames(geno2) == mapa2$marker))

# ------------------------------------------------------------
# 5. Criar posição genética aproximada
# ------------------------------------------------------------
# Como não há mapa genético em cM, o mapa físico foi reescalado.
# Cada cromossomo foi assumido com 2 Morgans.

chr_length_M <- 2.0
mapa2$gen_pos_M <- NA_real_

for (chr_i in unique(mapa2$chr)) {
  
  idx <- which(mapa2$chr == chr_i)
  pos_i <- mapa2$pos[idx]
  
  if (length(unique(pos_i)) == 1) {
    
    mapa2$gen_pos_M[idx] <- 0
    
  } else {
    
    mapa2$gen_pos_M[idx] <- (
      (pos_i - min(pos_i)) /
        (max(pos_i) - min(pos_i))
    ) * chr_length_M
  }
}

cat("\nResumo do mapa genético aproximado:\n")
print(
  mapa2 %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(
      n_snp = dplyr::n(),
      min_M = min(gen_pos_M),
      max_M = max(gen_pos_M),
      .groups = "drop"
    )
)

# ------------------------------------------------------------
# 6. Alinhar efeitos ao mapa
# ------------------------------------------------------------

efeitos2 <- efeitos_combinados %>%
  dplyr::filter(Marcador %in% mapa2$marker) %>%
  dplyr::slice(match(mapa2$marker, Marcador))

# ------------------------------------------------------------
# 7. Remover heterozigotos residuais e imputar
# ------------------------------------------------------------

genov2 <- geno2
genov2[genov2 == 1] <- NA

cat("\nGenótipos antes da imputação:\n")
print(table(genov2, useNA = "ifany"))

for (j in seq_len(ncol(genov2))) {
  
  if (any(is.na(genov2[, j]))) {
    
    tab_j <- table(genov2[, j])
    
    if (length(tab_j) == 0) {
      
      genov2[is.na(genov2[, j]), j] <- 0
      
    } else {
      
      major_allele <- as.numeric(names(tab_j)[which.max(tab_j)])
      genov2[is.na(genov2[, j]), j] <- major_allele
    }
  }
}

cat("\nGenótipos após imputação:\n")
print(table(genov2, useNA = "ifany"))

# ------------------------------------------------------------
# 8. Criar matriz de efeitos dos marcadores
# ------------------------------------------------------------

MarkEff_mat <- efeitos2 %>%
  dplyr::select(Efeito_AP, Efeito_AE, Efeito_FF, Efeito_FM) %>%
  as.matrix()

rownames(MarkEff_mat) <- efeitos2$Marcador
colnames(MarkEff_mat) <- c("AP", "AE", "FF", "FM")

# ------------------------------------------------------------
# 9. extrair planos de cruzamento
# ------------------------------------------------------------

extract_plan <- function(plan_raw, scenario_name) {
  
  if (is.list(plan_raw) && "plan" %in% names(plan_raw)) {
    plan <- plan_raw$plan
  } else {
    plan <- plan_raw
  }
  
  plan <- as.data.frame(plan)
  
  plan <- plan %>%
    dplyr::select(Parent1, Parent2, Y, K) %>%
    dplyr::mutate(
      scenario = scenario_name,
      cross_id = paste(Parent1, Parent2, sep = "_")
    ) %>%
    dplyr::select(scenario, cross_id, Parent1, Parent2, Y, K)
  
  return(plan)
}

# ------------------------------------------------------------
# 10. Organizar planos por cenário
# ------------------------------------------------------------

plan_AP <- extract_plan(plan_AP, "AP")
plan_AE <- extract_plan(plan_AE, "AE")
plan_FF <- extract_plan(plan_FF, "FF")
plan_FM <- extract_plan(plan_FM, "FM")
plan_MT <- extract_plan(plan_MT, "MT")
plan_Culling <- extract_plan(plan_Culling, "Culling")

plans_all <- dplyr::bind_rows(
  plan_AP,
  plan_AE,
  plan_FF,
  plan_FM,
  plan_MT,
  plan_Culling
)

cat("\nNúmero de cruzamentos por cenário:\n")
print(table(plans_all$scenario))

# ------------------------------------------------------------
# 11. Criar objeto final
# ------------------------------------------------------------

sim_inputs_real <- list(
  geno = genov2,
  mapa = mapa2,
  efeitos = efeitos2,
  MarkEff = MarkEff_mat,
  recombination_assumption = list(
    type = "physical_map_rescaled",
    chr_length_M = chr_length_M,
    unit = "Morgan",
    note = "Physical positions were rescaled within each chromosome because a genetic map was not available."
  ),
  plans_all = plans_all,
  plans_by_scenario = list(
    AP = plan_AP,
    AE = plan_AE,
    FF = plan_FF,
    FM = plan_FM,
    MT = plan_MT,
    Culling = plan_Culling
  )
)

saveRDS(sim_inputs_real, "sim_inputs_real.rds")
