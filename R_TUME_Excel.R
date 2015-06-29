# Code for processing TUME (www.projetotume.com)
# Authors: Eric Gorgens (gorgens at usp.br); Andre Gracioso Peres da Silva (andregracioso at gmail.com)
#
# Citation:
# Görgens, E.B.; Silva, A.G.P. (2015). R_TUME: Ferramentas de apoio para processamento dos dados de inventário 
# florestal do projeto Teste de Uso Múltiplo do Eucalyptus. Versão 1.0. Disponível em: 
# https://github.com/Gorgens/TUME/edit/master/R_TUME_Excel.R
#
###############################################################################
# version 1.0 (07/2015)

### -------------------------------------------------------------------------------
### Funcoes auxiliares
# Para cada especie, ajusta modelo hipsometrico e estima as alturas
hipsometrica <- function(tume.esp){
  
  tume.temp = tume.esp
  for (i in  seq(1, length(tume.temp$Cod), 1)) {
    if (is.na(tume.temp$Cod[i])){
      tume.temp$Cod[i] = 0
    }
  }
  
  uteis = subset(tume.temp, tume.temp$H_m != "NA" & tume.temp$Cod != c(4,5,7))
  
  logH = log(uteis$H_m)
  invD = 1/uteis$DAP_cm
  rm(uteis)
  
  modelo = lm(logH ~ invD)
  rm(logH, invD)
  tume.esp = cbind(tume.esp, invD = 1/tume.esp$DAP_cm)
  
  correcao = sum(exp(modelo$residuals))/length(modelo$residuals)
  
  estH = 0
  
  for (i in seq(1,nrow(tume.esp))){
    if (is.na(tume.esp$DAP_cm[i]) == FALSE & is.na(tume.esp$H_m[i]) == TRUE) {
      
      estH[i] = exp(predict(modelo, new = tume.esp[i,])) * correcao
      
    } else {
      
      estH[i] = tume.esp$H_m[i]
      
    }
  }
  
  return(estH)
  rm(modelo)
}

### Calcula volume das arvores individuais com base no DAP, altura e fator de forma
parcVolume <- function(tume.esp){
  tempV = na.omit(cbind(tume.esp$DAP_cm, tume.esp$estH_m))
  volume = sum((tempV[,1]^2*pi/40000) * tempV[,2] * 0.5) * 10000 / tume.esp$Parc_m2[1]
  rm(tempV)
  return(volume)
}

### Calcula e cria tabela de resumo (estatisticas) por especie para tumes com idade maior que 23 meses e SEM desbaste
resumo_pos24 <- function(tume.esp, estH_m){
  
  resumo_pos <- data.frame(N_tume = tume.esp$N_tume[1],
                           Esp = as.character(tume.esp$Esp[1]),
                           I_meses = 0,
                           Parc_m2 = 0,
                           DAPmed = 0,
                           DAPsd = 0,
                           Hmed = 0,
                           Hsd = 0,
                           Hdom = 0,
                           N_fuste = 0,
                           Sobr = 0,
                           G = 0,
                           V = 0,
                           IMA = 0,
                           B = 0)
  
  tume.esp = tume.esp[, -9] # Remover Cod2
  tume.esp = cbind(tume.esp, estH_m = estH)
  
  resumo_pos$I_meses = tume.esp$I_meses[1]
  resumo_pos$Parc_m2 = round(tume.esp$Parc_m2[1], 1)  
  resumo_pos$DAPmed = round(mean(na.omit(tume.esp$DAP_cm)), 1)
  resumo_pos$DAPsd = round(sd(na.omit(tume.esp$DAP_cm)), 1)
  resumo_pos$Hmed = round(mean(na.omit(tume.esp$H_m)), 1)
  resumo_pos$Hsd = round(sd(na.omit(tume.esp$H_m)), 1)
  resumo_pos$Hdom = round(mean(na.omit(tume.esp[tume.esp$Cod == 6, names(tume.esp) %in% c("H_m")])), 1)
  resumo_pos$N_fuste = round(length(na.omit(tume.esp$DAP_cm)) * 10000 / tume.esp$Parc_m2[1], 0)
  resumo_pos$Sobr = round(((1 - ((nrow(tume.esp) - nrow(tume.esp[tume.esp$Cod != 1 & tume.esp$Cod != 5,])) / max(tume.esp$N_arv))) * 100), 1) 
  resumo_pos$G = round(sum(na.omit(tume.esp$DAP_cm)^2 * pi /40000) * 10000 / tume.esp$Parc_m2[1], 1)
  resumo_pos$V = round(parcVolume(tume.esp), 0)  
  resumo_pos$IMA = round(resumo_pos$V[1] / (tume.esp$I_meses[1]/12), 1)
  if (e %in% ESP.DENSIDADE$Esp){
    resumo_pos$B = round(resumo_pos$V[1] * ESP.DENSIDADE[ESP.DENSIDADE$Esp == e, 2],0)
  } else {
    resumo_pos$B = ""
  } # Incluir depois
  
  return(resumo_pos)
  
}

### Calcula e cria tabela de resumo (estatisticas) por especie para tumes com idade maior que 23 meses e COM desbaste
resumo_pos24desb <- function(tume.esp, estH_m){
  
  resumo_pos <- data.frame(N_tume = tume.esp$N_tume[1],
                           Esp = as.character(tume.esp$Esp[1]),
                           I_meses = 0,
                           Parc_m2 = 0,
                           DAPmed = 0,
                           DAPsd = 0,
                           Hmed = 0,
                           Hsd = 0,
                           Hdom = 0,
                           N_fuste = 0,
                           Sobr = 0,
                           G = 0,
                           V = 0, 
                           IMA = 0,
                           B = 0)
  
  tume.esp = tume.esp[, -9]
  tume.esp = cbind(tume.esp, estH_m = estH)
  
  resumo_pos$I_meses = tume.esp$I_meses[1]
  resumo_pos$Parc_m2 = round(tume.esp$Parc_m2[1], 1)
  resumo_pos$DAPmed = round(mean(na.omit(tume.esp$DAP_cm)), 1)
  resumo_pos$DAPsd = round(sd(na.omit(tume.esp$DAP_cm)), 1)
  resumo_pos$Hmed = round(mean(na.omit(tume.esp$H_m)), 1)
  resumo_pos$Hsd = round(sd(na.omit(tume.esp$H_m)), 1)
  resumo_pos$Hdom = round(mean(na.omit(tume.esp[tume.esp$Cod == 6, names(tume.esp) %in% c("H_m")])), 1)
  resumo_pos$N_fuste = round(length(na.omit(tume.esp$DAP_cm)) * 10000 / tume.esp$Parc_m2[1], 0)
  resumo_pos$Sobr = ""
  resumo_pos$G = round(sum(na.omit(tume.esp$DAP_cm)^2 * pi /40000) * 10000 / tume.esp$Parc_m2[1], 1)
  resumo_pos$V = round(parcVolume(tume.esp), 0)
  resumo_pos$IMA = ""
  if (e %in% ESP.DENSIDADE$Esp){
    resumo_pos$B = round(resumo_pos$V[1] * ESP.DENSIDADE[ESP.DENSIDADE$Esp == e, 2],0)
  } else {
    resumo_pos$B = ""
  }
  
  return(resumo_pos)
  
}

### Calcula e cria tabela de resumo (estatisticas) por especie para tumes com idade inferior à 24 meses
resumo_pre24 <- function(tume.esp){
  
  resumo_pre <- data.frame(N_tume = tume.esp$N_tume[1],
                           Esp = as.character(tume.esp$Esp[1]),
                           I_meses = 0,
                           Parc_m2 = 0,
                           Dapmed = 0,
                           Dapsd = 0,
                           Hmed = 0,
                           Hsd = 0,
                           Hdom = 0,
                           N_fuste = 0,
                           Sobr = 0,
                           G = 0,
                           V = 0,
                           IMA = 0,
                           B = 0)
  
  resumo_pre$I_meses = tume.esp$I_meses[1]
  resumo_pre$Parc_m2 = round(tume.esp$Parc_m2[1], 1)
  resumo_pre$Dapmed = ""
  resumo_pre$Dapsd = ""
  resumo_pre$Hmed = round(mean(na.omit(tume.esp$H_m)), 1)
  resumo_pre$Hsd = round(sd(na.omit(tume.esp$H_m)), 1)
  resumo_pre$Hdom = round(mean(na.omit(tume.esp[tume.esp$Cod == 6, names(tume.esp) %in% c("H_m")])), 1) # Stick 6
  resumo_pre$N_fuste = round(length(na.omit(tume.esp$H_m)) * 10000 / tume.esp$Parc_m2[1], 0)
  resumo_pre$Sobr = round(((1 - ((nrow(tume.esp) - nrow(tume.esp[tume.esp$Cod != 1 & tume.esp$Cod != 5,])) / max(tume.esp$N_arv))) * 100), 1) 
  resumo_pre$G = ""
  resumo_pre$V = ""
  resumo_pre$IMA = ""
  resumo_pre$B = ""
  
  return(resumo_pre)
  
}

### Cria grafico de barras para a variavel estoque de volume
plotVolume <- function(tabela_resumo, l){
  
  #Ordena o Volume de forma descrescente 
  tabela_resumo = tabela_resumo[with(tabela_resumo, order(-V)), ]
  
  #Distancia horizontal para ultimo rotulo do eixo x
  end_point = 0.5 + nrow(tabela_resumo) + nrow(tabela_resumo)-1
  
  if (nrow(tabela_resumo) > 17) {
  #Armazena grafico de barras do Volume em um arquivo de extensao .jpeg
  jpeg(paste(TUME.OUT, l, ".jpg", sep=""), height = 7, width = nrow(tabela_resumo)*0.6, units = "cm", res = 300)
  } else {
  jpeg(paste(TUME.OUT, l, ".jpg", sep=""), height = 7, width = 10, units = "cm", res = 300)  
  }
    
  #Adiciona espaco na margem inferior para rotulacao do eixo x
  par(mar = c(5, 4, 2, 1) + 0.2, mgp = c(2.2,1,0))
  #Fonte padrao utilizada pelo R
  par(family="sans")
  #Configuracoes do grafico de barras
  barplot(tabela_resumo$V,
          col="grey50", 
          main=paste("TUME ", tabela_resumo$N_tume[1], " - ", tabela_resumo$I_meses[1], " meses", sep=""),
          cex.main = 0.6,
          cex.axis = 0.6,
          cex.lab = 0.6,
          ylab = "Volume (m³/ha)",
          ylim = c(0, 1.1 * max(tabela_resumo$V)),
          xlab = "",
          space = 1)
  
  # Coloca e rotaciona o rotulo do eixo X saltando os espacos entre colunas
  text(seq(1.5, end_point, by=2),
       par("usr")[3]-2, 
       srt = 60, #rotaciona 60 graus no sentido horario
       adj = 1,
       xpd = TRUE,
       labels = tabela_resumo$Esp,
       cex = 0.6)
  dev.off()
}

### Cria grafico de barras para a variavel altura media
plotHmed <- function(tabela_resumo, l){
  
  #Ordena a altura média de forma descrescente 
  tabela_resumo = tabela_resumo[with(tabela_resumo, order(-Hmed)), ]
  
  #Distancia horizontal para ultimo rotulo do eixo x
  end_point = 0.5 + nrow(tabela_resumo) + nrow(tabela_resumo)-1
  
  if (nrow(tabela_resumo) > 17) {
  #Armazena grafico de barras do Volume em um arquivo de extensao .jpeg
  jpeg(paste(TUME.OUT, l, ".jpg", sep=""), height = 7, width = nrow(tabela_resumo)*0.6, units = "cm", res = 300)
  } else {
  jpeg(paste(TUME.OUT, l, ".jpg", sep=""), height = 7, width = 10, units = "cm", res = 300) 
  }
  
  #Adiciona espaco na margem inferior para rotulacao do eixo x
  par(mar = c(5, 4, 2, 1) + 0.2, mgp = c(2.2,1,0))
  #Fonte padrao utilizada pelo R
  par(family="sans")
  #Configuracoes do grafico de barras
  barplot(tabela_resumo$Hmed,
          col="grey50", 
          main=paste("TUME ", tabela_resumo$N_tume[1], " - ", tabela_resumo$I_meses[1], " meses", sep=""),
          cex.main=0.6,
          cex.axis = 0.6,
          cex.lab = 0.6,
          ylab = "Altura media (m)",
          ylim = c(0, 1.1 * max(tabela_resumo$Hmed)),
          xlab = "",
          space = 1) #largura das barras e largura entre as barras igual
  
  # Coloca e rotaciona o rotulo do eixo X saltando os espacos entre colunas
  text(seq(1.5, end_point, by=2),
       par("usr")[3]-0.1, 
       srt = 60, #rotaciona 60 graus no sentido horario
       adj = 1,
       xpd = TRUE,
       labels = tabela_resumo$Esp,
       cex = 0.6)
  dev.off()
  
}

### -------------------------------------------------------------------------------
### Variaveis globais

# Recebe parâmetros do comando batch construído no código em VBA
args = commandArgs(trailingOnly = TRUE)

# Define pasta com arquivos de medicoes (arquivos de entrada)
TUME.PATH <- args[1]

# Define pasta para armazenamento dos arquivos de saida
TUME.OUT <- args[2]

# Define pasta com arquivos de referencia (ex: lista de densidades basicas por material genetico)
TUME.REF <- args[3]

# Cria vetor com os nomes dos arquivos
TUME.FILES <- list.files(TUME.PATH)

# Importa tabela de densidade (Densidades.csv)
ESP.DENSIDADE <- read.csv(paste(TUME.REF, "Densidades.csv", sep=""))

### -------------------------------------------------------------------------------
### Inicio da analise

# Filtra o tume para uma determinada especie
for (l in TUME.FILES){
  
  # Importa arquivo de um tume
  tume = read.csv(paste(TUME.PATH, l, sep=""), sep=",")
  
  # Salva vetor com nome das especies contidas no TUME
  TUME.ESP <- levels(tume$Esp)
  
  tabela_resumo = data.frame()
  
  if (tume$I_meses[1] > 23 & 3 %in% tume$Cod){
    for (e in TUME.ESP){
      
      tume.esp = subset(tume, tume$Esp == e)
      
      if (nrow(na.omit(tume.esp[,1:6])) > 3 & is.numeric(tume.esp$DAP_cm) == TRUE){
        
        estH = hipsometrica(tume.esp)
        
        tabela_resumo = rbind(tabela_resumo, resumo_pos24desb(tume.esp, estH))
      } else {
        
        sem_dados <- data.frame(N_tume = tume.esp$N_tume[1],
                                Esp = as.character(tume.esp$Esp[1]),
                                I_meses = 0,
                                Parc_m2 = 0,
                                Dapmed = 0,
                                Dapsd = 0,
                                Hmed = 0,
                                Hsd = 0,
                                Hdom = 0,
                                N_fuste = 0,
                                Sobr = 0,
                                G = 0,
                                V = 0,
                                IMA = 0,
                                B = 0)
        tabela_resumo <- rbind(tabela_resumo, sem_dados)
        
      }
    }
    
    plotVolume(tabela_resumo, l)

    write.csv(tabela_resumo, file = paste(TUME.OUT,"saida_", l, sep = ""),row.names=FALSE)
    
    
  } else if (tume$I_meses[1] > 23 & !(3 %in% tume$Cod)){
    for (e in TUME.ESP){
      
      tume.esp = subset(tume, tume$Esp == e)
      
      if (nrow(na.omit(tume.esp[,1:6])) > 3 & is.numeric(tume.esp$DAP_cm) == TRUE){
        
        estH = hipsometrica(tume.esp)
        
        tabela_resumo = rbind(tabela_resumo, resumo_pos24(tume.esp, estH))
      } else {
        
        sem_dados <- data.frame(N_tume = tume.esp$N_tume[1],
                                Esp = as.character(tume.esp$Esp[1]),
                                I_meses = 0,
                                Parc_m2 = 0,
                                Dapmed = 0,
                                Dapsd = 0,
                                Hmed = 0,
                                Hsd = 0,
                                Hdom = 0,
                                N_fuste = 0,
                                Sobr = 0,
                                G = 0,
                                V = 0,
                                IMA = 0,
                                B = 0)
        tabela_resumo <- rbind(tabela_resumo, sem_dados)
        
      }
    }
    
    plotVolume(tabela_resumo, l)
    
    write.csv(tabela_resumo, file = paste(TUME.OUT,"saida_", l, sep = ""),row.names=FALSE)
    
    
  } else {
    
    for (e in TUME.ESP){
      
      tume.esp = subset(tume, tume$Esp == e)
      
      if (nrow(na.omit(tume.esp[,c(1:5, 7)])) > 3){
        
        tabela_resumo = rbind(tabela_resumo, resumo_pre24(tume.esp))
        
      } else {
        
        sem_dados <- data.frame(N_tume = tume.esp$N_tume[1],
                                Esp = as.character(tume.esp$Esp[1]),
                                I_meses = 0,
                                Parc_m2 = 0,
                                Dapmed = 0,
                                Dapsd = 0,
                                Hmed = 0,
                                Hsd = 0,
                                Hdom = 0,
                                N_fuste = 0,
                                Sobr = 0,
                                G = 0,
                                V = 0,
                                IMA = 0,
                                B = 0)
        tabela_resumo <- rbind(tabela_resumo, sem_dados)
        
      }
    }
    
    plotHmed(tabela_resumo, l)
    
    write.csv(tabela_resumo, file = paste(TUME.OUT, "saida_", l, sep = ""),row.names=FALSE)
    
  }
  
}
