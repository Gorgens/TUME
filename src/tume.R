# TODO:
# 
# Code for processing TUME (ESALQ/USP)
# Author: Gorgens (gorgens at usp.br)
###############################################################################

### -------------------------------------------------------------------------------
### Funcoes auxiliares
# Estima alturas nao medidas em campo por especie de um tume
hipsometrica <- function(tume.esp){
  
  tume.temp = tume.esp
  for (i in  seq(1, length(tume.temp$Cod), 1)) {
    if (is.na(tume.temp$Cod[i])){
      tume.temp$Cod[i] = 0
    }
  }
  
  uteis = subset(tume.temp, tume.temp$H_m != "NA" & tume.temp$Cod != c(4, 7))
  
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

### Calculo do volume por especie dentro de um tume
parcVolume <- function(tume.esp){
  tempV = na.omit(cbind(tume.esp$DAP_cm, tume.esp$estH_m))
  volume = sum((tempV[,1]^2*pi/40000) * tempV[,2] * 0.5) * 10000 / tume.esp$Parc_m2[1]
  rm(tempV)
  return(volume)
}

### Calcula e cria tabela de resumo por especie dentro de um tume maior que 2 anos sem desbaste
resumo_pos24 <- function(tume.esp, estH_m){
  
  resumo_pos <- data.frame(Tume = tume.esp$N_tume[1],
                           Especie = as.character(tume.esp$Esp[1]),
                           Idade = 0,
                           Area = 0,
                           DAPmed = 0,
                           DAPsd = 0,
                           Hmed = 0,
                           Hsd = 0,
                           Hdom = 0,
                           N = 0,
                           Sobre = 0,
                           G = 0,
                           V = 0,
                           IMA = 0,
                           Biomassa = 0)
  
  tume.esp = tume.esp[, -9] # Remover Cod2
  tume.esp = cbind(tume.esp, estH_m = estH)
  
  resumo_pos$Idade = tume.esp$I_meses[1]
  resumo_pos$DAPmed = round(mean(na.omit(tume.esp$DAP_cm)), 1)
  resumo_pos$DAPsd = round(sd(na.omit(tume.esp$DAP_cm)), 1)
  resumo_pos$Hmed = round(mean(na.omit(tume.esp$H_m)), 1)
  resumo_pos$Hsd = round(sd(na.omit(tume.esp$H_m)), 1)
  resumo_pos$Hdom = round(mean(na.omit(tume.esp[tume.esp$Cod == 6, names(tume.esp) %in% c("H_m")])), 1) # Stick 6
  resumo_pos$N = round(length(na.omit(tume.esp$DAP_cm)) * 10000 / tume.esp$Parc_m2[1], 0)
  resumo_pos$Sobre = round(length(na.omit(tume.esp$DAP_cm)) / max(tume.esp$N_arv) * 100, 0) # Stick 1
  resumo_pos$G = round(sum(na.omit(tume.esp$DAP_cm)^2 * pi /40000) * 10000 / tume.esp$Parc_m2[1], 0)
  resumo_pos$V = round(parcVolume(tume.esp), 0)  
  resumo_pos$IMA = round(resumo_pos$V[1] / (tume.esp$I_meses[1]/12), 0)
  resumo_pos$Area = round(tume.esp$Parc_m2[1], 1)  
  if (e %in% ESP.DENSIDADE$Esp){
    resumo_pos$Biomassa = resumo_pos$V[1] * ESP.DENSIDADE[ESP.DENSIDADE$Esp == e, 2]
  } else {
    resumo_pos$Biomassa = 0
  } # Incluir depois
  
  return(resumo_pos)
  
}

### Calcula e cria tabela de resumo por especie dentro de um tume maior que 2 anos com desbaste
resumo_pos24desb <- function(tume.esp, estH_m){
  
  resumo_pos <- data.frame(Tume = tume.esp$N_tume[1],
                           Especie = as.character(tume.esp$Esp[1]),
                           Idade = 0,
                           Area = 0,
                           DAPmed = 0,
                           DAPsd = 0,
                           Hmed = 0,
                           Hsd = 0,
                           Hdom = 0,
                           N = 0,
                           Sobre = 0,
                           G = 0,
                           V = 0, 
                           IMA = 0,
                           Biomassa = 0)
  
  tume.esp = tume.esp[, -9]
  tume.esp = cbind(tume.esp, estH_m = estH)
  
  resumo_pos$Idade = tume.esp$I_meses[1]
  resumo_pos$DAPmed = round(mean(na.omit(tume.esp$DAP_cm)), 1)
  resumo_pos$DAPsd = round(sd(na.omit(tume.esp$DAP_cm)), 1)
  resumo_pos$Hmed = round(mean(na.omit(tume.esp$H_m)), 1)
  resumo_pos$Hsd = round(sd(na.omit(tume.esp$H_m)), 1)
  resumo_pos$Hdom = round(mean(na.omit(tume.esp[tume.esp$Cod == 6, names(tume.esp) %in% c("H_m")])), 1) # Stick 6
  resumo_pos$N = round(length(na.omit(tume.esp$DAP_cm)) * 10000 / tume.esp$Parc_m2[1], 0)
  resumo_pos$Sobre = ""
  resumo_pos$G = round(sum(na.omit(tume.esp$DAP_cm)^2 * pi /40000) * 10000 / tume.esp$Parc_m2[1], 0)
  resumo_pos$V = round(parcVolume(tume.esp), 0)
  resumo_pos$IMA = ""
  resumo_pos$Area = round(tume.esp$Parc_m2[1], 1)
  if (e %in% ESP.DENSIDADE$Esp){
    resumo_pos$Biomassa = resumo_pos$V[1] * ESP.DENSIDADE[ESP.DENSIDADE$Esp == e, 2]
  } else {
    resumo_pos$Biomassa = 0
  } # Incluir depois
  
  return(resumo_pos)
  
}

### Calcula e cria tabela de resumo por especie dentro de um tume menor que 2 anos
resumo_pre24 <- function(tume.esp){
  
  resumo_pre <- data.frame(Tume = tume.esp$N_tume[1],
                           Especie = as.character(tume.esp$Esp[1]),
                           Idade = 0,
                           Area = 0,
                           Dapmed = 0,
                           Dapsd = 0,
                           Hmed = 0,
                           Hsd = 0,
                           Hdom = 0,
                           N = 0,
                           Sobre = 0,
                           G = 0,
                           V = 0,
                           IMA = 0,
                           Biomassa = 0)
  
  resumo_pre$Idade = tume.esp$I_meses[1]
  resumo_pre$Hmed = round(mean(na.omit(tume.esp$H_m)), 1)
  resumo_pre$Hsd = round(sd(na.omit(tume.esp$H_m)), 1)
  resumo_pre$Hdom = round(mean(na.omit(tume.esp[tume.esp$Cod == 6, names(tume.esp) %in% c("H_m")])), 1) # Stick 6
  resumo_pre$N = round(length(na.omit(tume.esp$H_m)) * 10000 / tume.esp$Parc_m2[1], 0)
  resumo_pre$Sobre = round(length(na.omit(tume.esp$H_m)) / max(tume.esp$N_arv) * 100, 0) # Stick 1
  resumo_pre$Area = round(tume.esp$Parc_m2[1], 1)
  resumo_pre$Dapmed = ""
  resumo_pre$Dapsd = ""
  resumo_pre$G = ""
  resumo_pre$V = ""
  resumo_pre$IMA = ""
  resumo_pre$Biomassa = ""
  
  return(resumo_pre)
  
}

### Cria gráfico de barras para o IMA
plotVolume <- function(tabela_resumo, l){
  
  #jpeg(paste(TUME.OUT, l, ".jpg", sep=""))
  #barplot(tabela_resumo$IMA, main="IMA", names.arg = tabela_resumo$Especie, cex.names=0.3)
  #dev.off()
  
  # ordena o IMA de forma descrescente 
  tabela_resumo = tabela_resumo[with(tabela_resumo, order(-V)), ]
  
  # Distância para rotulação do eixo x
  end_point = 0.5 + nrow(tabela_resumo) + nrow(tabela_resumo)-1
  
  # cria gráfico de barras do IMA
  jpeg(paste(TUME.OUT, l, ".jpg", sep=""), height = 7, width = 15, units = "cm", res = 150)
  
  # Adiciona margem inferior para rótulo do eixo x rotacionado
  par(mar = c(7, 4, 2, 2) + 0.2)
  barplot(tabela_resumo$V,
          col="grey50", 
          main=paste("TUME ", tabela_resumo$Tume[1], " - ", tabela_resumo$Idade[1], " meses", sep=""),
          cex.main = 1,
          cex.axis = 0.5,
          cex.lab = 0.8,
          ylab = "Volume (m³/ha)",
          ylim = c(0, 1.1 * max(tabela_resumo$V)),
          xlab = "",
          space = 1)
  
  # Coloca e rotaciona o rótulo do eixo X saltando os espaço entre colunas
  text(seq(1.5, end_point, by=2),
       par("usr")[3]-0.25, 
       srt = 60,
       adj = 1,
       xpd = TRUE,
       labels = tabela_resumo$Especie,
       cex = 0.5)
  dev.off()
}

### Cria gráfico de barras para a altura média
plotHmed <- function(tabela_resumo, l){
  
  #jpeg(paste(TUME.OUT, l, ".jpg", sep=""))
  #barplot(tabela_resumo$Sobre, main="Sobrevivência", names.arg = tabela_resumo$Especie, cex.names=0.3)
  #dev.off()
  
  # ordena a sobrevivência de forma descrescente 
  tabela_resumo = tabela_resumo[with(tabela_resumo, order(-Hmed)), ]
  
  # Distância para rotulação do eixo x
  end_point = 0.5 + nrow(tabela_resumo) + nrow(tabela_resumo)-1
  
  # cria gráfico de barras do IMA
  jpeg(paste(TUME.OUT, l, ".jpg", sep=""), height = 7, width = 15, units = "cm", res = 150)
  
  # Adiciona margem inferior para rótulo do eixo x rotacionado
  par(mar = c(7, 4, 4, 2) + 0.2)
  barplot(tabela_resumo$Hmed,
          col="grey50", 
          main=paste("TUME ", tabela_resumo$Tume[1], " - ", tabela_resumo$Idade[1], " meses", sep=""),
          cex.main=1,
          cex.axis = 0.5,
          cex.lab = 0.8,
          ylab = "Altura média (m)",
          ylim = c(0, 1.1 * max(tabela_resumo$Hmed)),
          xlab = "",
          space = 1)
  
  # Coloca e rotaciona o rótulo do eixo X saltando os espaço entre colunas
  text(seq(1.5, end_point, by=2),
       par("usr")[3]-0.25, 
       srt = 60,
       adj = 1,
       xpd = TRUE,
       labels = tabela_resumo$Especie,
       cex = 0.5)
  dev.off()
  
}

plotSobre <- function(tabela_resumo, l){
  
  #jpeg(paste(TUME.OUT, l, ".jpg", sep=""))
  #barplot(tabela_resumo$Sobre, main="Sobrevivência", names.arg = tabela_resumo$Especie, cex.names=0.3)
  #dev.off()
  
  # ordena a sobrevivência de forma descrescente 
  tabela_resumo = tabela_resumo[with(tabela_resumo, order(-Sobre)), ]
  
  # Distância para rotulação do eixo x
  end_point = 0.5 + nrow(tabela_resumo) + nrow(tabela_resumo)-1
  
  # cria gráfico de barras do IMA
  jpeg(paste(TUME.OUT, l, ".jpg", sep=""), res = 150)
  
  # Adiciona margem inferior para rótulo do eixo x rotacionado
  par(mar = c(7, 4, 4, 2) + 0.2)
  barplot(tabela_resumo$Sobre,
          col="grey50", 
          main=paste("Sobrevivência - ", tabela_resumo$Idade[1], " meses", sep=""),
          cex.main=1,
          ylab = "Sobrevivência (%)",
          ylim = c(0,100),
          xlab = "",
          space = 1)
  
  # Coloca e rotaciona o rótulo do eixo X saltando os espaço entre colunas
  text(seq(1.5, end_point, by=2),
       par("usr")[3]-0.25, 
       srt = 60,
       adj = 1,
       xpd = TRUE,
       labels = tabela_resumo$Especie,
       cex = 1)
  dev.off()
  
}

### -------------------------------------------------------------------------------
### Variaveis globais

# Define pasta com arquivos de medicoes
TUME.PATH <- paste(getwd(), "/in/", sep = "")

# Define pasta com arquivos de medicoes
TUME.OUT <- paste(getwd(), "/out/", sep = "")

# Define pasta com arquivos de referência
TUME.REF <- paste(getwd(), "/referencia/", sep = "")

# Cria vetor com os nomes dos arquivos
TUME.FILES <- list.files(TUME.PATH)

# Importa tabela de densidade
ESP.DENSIDADE <- read.csv(paste(TUME.REF, "Densidades.csv", sep=""))

### -------------------------------------------------------------------------------
### Inicio da analise

# Filtra o tume para uma determinada especie
for (l in TUME.FILES){
  
  # Importa arquivo de um tume
  tume = read.csv(paste(TUME.PATH, l, sep=""), sep=",")
  
  # Salva vetor com nome das espécies contidas no TUME
  TUME.ESP <- levels(tume$Esp)
  
  tabela_resumo = data.frame()
  
  if (tume$I_meses[1] > 23 & 3 %in% tume$Cod){
    for (e in TUME.ESP){
      
      tume.esp = subset(tume, tume$Esp == e)
      
      if (nrow(na.omit(tume.esp[,1:6])) > 3 & is.numeric(tume.esp$DAP_cm) == TRUE){
        
        estH = hipsometrica(tume.esp)
        
        tabela_resumo = rbind(tabela_resumo, resumo_pos24desb(tume.esp, estH))
      } else {
        
        sem_dados <- data.frame(Tume = tume.esp$N_tume[1],
                                Especie = as.character(tume.esp$Esp[1]),
                                Idade = 0,
                                Area = 0,
                                Dapmed = 0,
                                Dapsd = 0,
                                Hmed = 0,
                                Hsd = 0,
                                Hdom = 0,
                                N = 0,
                                Sobre = 0,
                                G = 0,
                                V = 0,
                                IMA = 0,
                                Biomassa = 0)
        tabela_resumo <- rbind(tabela_resumo, sem_dados)
        
      }
    }
    
    plotVolume(tabela_resumo, l)
    
    #write.csv(tabela_resumo, file = paste(TUME.OUT, l, ".csv", sep = ""))
    #unidades = c("", "meses", "m^2", "cm", "cm", "m", "m", "m", "m^2/ha", "m^3/ha", "fustes/ha", "%", "m^3/ha/ano")
    #tabela_resumo = rbind(unidades, tabela_resumo)  
    write.csv(tabela_resumo, file = paste(TUME.OUT, "out_", l, sep = ""))
    
    
  } else if (tume$I_meses[1] > 23 & !(3 %in% tume$Cod)){
    for (e in TUME.ESP){
      
      tume.esp = subset(tume, tume$Esp == e)
      
      if (nrow(na.omit(tume.esp[,1:6])) > 3 & is.numeric(tume.esp$DAP_cm) == TRUE){
        
        estH = hipsometrica(tume.esp)
        
        tabela_resumo = rbind(tabela_resumo, resumo_pos24(tume.esp, estH))
      } else {
        
        sem_dados <- data.frame(Tume = tume.esp$N_tume[1],
                                Especie = as.character(tume.esp$Esp[1]),
                                Idade = 0,
                                Area = 0,
                                Dapmed = 0,
                                Dapsd = 0,
                                Hmed = 0,
                                Hsd = 0,
                                Hdom = 0,
                                N = 0,
                                Sobre = 0,
                                G = 0,
                                V = 0,
                                IMA = 0,
                                Biomassa = 0)
        tabela_resumo <- rbind(tabela_resumo, sem_dados)
        
      }
    }
    
    plotVolume(tabela_resumo, l)
    
    #write.csv(tabela_resumo, file = paste(TUME.OUT, l, ".csv", sep = ""))
    #unidades = c("", "meses", "m^2", "cm", "cm", "m", "m", "m", "m^2/ha", "m^3/ha", "fustes/ha", "%", "m^3/ha/ano")
    #tabela_resumo = rbind(unidades, tabela_resumo)  
    write.csv(tabela_resumo, file = paste(TUME.OUT, "out_", l, sep = ""))
    
    
  } else {
    
    for (e in TUME.ESP){
      
      tume.esp = subset(tume, tume$Esp == e)
      
      if (nrow(na.omit(tume.esp[,c(1:5, 7)])) > 3){
        
        tabela_resumo = rbind(tabela_resumo, resumo_pre24(tume.esp))
        
      } else {
        
        sem_dados <- data.frame(Tume = tume.esp$N_tume[1],
                                Especie = as.character(tume.esp$Esp[1]),
                                Idade = 0,
                                Area = 0,
                                Dapmed = 0,
                                Dapsd = 0,
                                Hmed = 0,
                                Hsd = 0,
                                Hdom = 0,
                                N = 0,
                                Sobre = 0,
                                G = 0,
                                V = 0,
                                IMA = 0,
                                Biomassa = 0)
        tabela_resumo <- rbind(tabela_resumo, sem_dados)
        
      }
    }
    
    #plotSobre(tabela_resumo, l)
    plotHmed(tabela_resumo, l)
    
    #write.csv(tabela_resumo, file = paste(TUME.OUT, l, ".csv", sep = ""))
    #unidades = c("", "meses", "m", "m", "fustes/ha", "%")
    #tabela_resumo = rbind(unidades, tabela_resumo)
    write.csv(tabela_resumo, file = paste(TUME.OUT, "out_", l, sep = ""))
    
    
  }
  
}
