# TODO:
#
# Stick 7: Remover extenção dos nomes dos arquivos 
# 
# Author: Gorgens
###############################################################################

### -------------------------------------------------------------------------------
### Funcoes auxiliares
# Estima alturas nao medidas em campo por especie de um tume
hipsometrica <- function(tume.esp){
	
	uteis = na.omit(tume.esp[, -8])
	
	logH = log(uteis$H_m)
	invD = 1/uteis$DAP_cm
	rm(uteis)
	
	modelo = lm(logH ~ invD)
	rm(logH, invD)
	tume.esp = cbind(tume.esp, invD = 1/tume.esp$DAP_cm)
	
	estH = 0
	
	for (i in seq(1,nrow(tume.esp))){
		if (is.na(tume.esp$DAP_cm[i]) == FALSE & is.na(tume.esp$H_m[i]) == TRUE){
			
			estH[i] = exp(predict(modelo, new = tume.esp[i,]))
			
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

### Calcula e cria tabela de resumo por especie dentro de um tume maior que 2 anos
resumo_pos24 <- function(tume.esp, estH_m){
	
	resumo_pos <- data.frame(Especie = as.character(tume.esp$Esp[1]),
			Idade = 0,
			DAPmed = 0,
			DAPsd = 0,
			Hmed = 0,
			Hsd = 0,
			Hdom = 0,
			G = 0,
			V = 0,
			N = 0,
			Sobre = 0,
			IMA = 0)
	
	tume.esp = cbind(tume.esp, estH_m = estH)
	
	resumo_pos$Idade = tume.esp$I_meses[1]
	resumo_pos$DAPmed = mean(na.omit(tume.esp$DAP_cm))
	resumo_pos$DAPsd = sd(na.omit(tume.esp$DAP_cm))
	resumo_pos$Hmed = mean(na.omit(tume.esp$estH_m))
	resumo_pos$Hsd = sd(na.omit(tume.esp$estH_m))
	#resumo_pos$Hdom = mean(tume.esp$H_m[tume.esp$Cod == 6]) # Stick 6
	resumo_pos$Hdom = max(na.omit(tume.esp$H_m)) # Stick 6
	resumo_pos$G = sum(na.omit(tume.esp$DAP_cm)^2/40000) * 10000 / tume.esp$Parc_m2[1]
	resumo_pos$V = parcVolume(tume.esp)  
	resumo_pos$N = max(tume.esp$N_arv) * 10000 / tume.esp$Parc_m2[1] 
	#resumo_pos$Sobre = (max(tume.esp$N_arv) - nrow(tume.esp[tume.esp$Cod == 1 & tume.esp$Cod == 5,])) / max(tume.esp$N_arv) * 100 # Stick 1
	resumo_pos$Sobre = length(na.omit(tume.esp$DAP_cm)) / max(tume.esp$N_arv) * 100 # Stick 1
	resumo_pos$IMA = resumo_pos$V[1] / (tume.esp$I_meses[1]/12)
	
	return(resumo_pos)
	
}

### Calcula e cria tabela de resumo por especie dentro de um tume menor que 2 anos
resumo_pre24 <- function(tume.esp){
	
	resumo_pre <- data.frame(Especie = as.character(tume.esp$Esp[1]),
			Idade = 0,
			Hmed = 0,
			Hsd = 0,
			N = 0,
			Sobre = 0)
	
	resumo_pre$Idade = tume.esp$I_meses[1]
	resumo_pre$Hmed = mean(na.omit(tume.esp$H_m))
	resumo_pre$Hsd = sd(na.omit(tume.esp$H_m))
	resumo_pre$N = max(tume.esp$N_arv) * 10000 / tume.esp$Parc_m2[1] 
	#resumo_pre$Sobre = (max(tume.esp$N_arv) - nrow(tume.esp[tume.esp$Cod == 1 & tume.esp$Cod == 5,])) / max(tume.esp$N_arv) * 100 # Stick 1
	resumo_pre$Sobre = length(na.omit(tume.esp$H_m)) / max(tume.esp$N_arv) * 100 # Stick 1
	
	return(resumo_pre)
	
}

### Cria gráfico de barras para o IMA
plotIMA <- function(tabela_resumo, l){
	
	#jpeg(paste(TUME.OUT, l, ".jpg", sep=""))
	#barplot(tabela_resumo$IMA, main="IMA", names.arg = tabela_resumo$Especie, cex.names=0.3)
	#dev.off()
	
	# ordena o IMA de forma descrescente 
	tabela_resumo = tabela_resumo[with(tabela_resumo, order(-IMA)), ]
	
	# Distância para rotulação do eixo x
	end_point = 0.5 + nrow(tabela_resumo) + nrow(tabela_resumo)-1
	
	# cria gráfico de barras do IMA
	jpeg(paste(TUME.OUT, l, ".jpg", sep=""))
	
	# Adiciona margem inferior para rótulo do eixo x rotacionado
	par(mar = c(7, 4, 2, 2) + 0.2)
	barplot(tabela_resumo$IMA,
			col="grey50", 
			main=paste("IMA - ", tabela_resumo$Idade[1], " meses", sep=""),
			cex.main=1,
			ylab = "Volume (m³/ha.ano)",
			ylim = c(0,5+max(tabela_resumo$IMA)),
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

### Cria gráfico de barras para o sobrevivencia
plotSobre <- function(tabela_resumo, l){
	
	#jpeg(paste(TUME.OUT, l, ".jpg", sep=""))
	#barplot(tabela_resumo$Sobre, main="Sobrevivência", names.arg = tabela_resumo$Especie, cex.names=0.3)
	#dev.off()

	# ordena a sobrevivência de forma descrescente 
	tabela_resumo = tabela_resumo[with(tabela_resumo, order(-Sobre)), ]
	
	# Distância para rotulação do eixo x
	end_point = 0.5 + nrow(tabela_resumo) + nrow(tabela_resumo)-1
	
	# cria gráfico de barras do IMA
	jpeg(paste(TUME.OUT, l, ".jpg", sep=""))
	
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
TUME.PATH <- "C:\\Users\\Gorgens\\Documents\\R\\TUME\\in\\"

# Define pasta com arquivos de medicoes
TUME.OUT <- "C:\\Users\\Gorgens\\Documents\\R\\TUME\\out\\"

# Cria vetor com os nomes dos arquivos
TUME.FILES <- list.files(TUME.PATH)

### -------------------------------------------------------------------------------
### Inicio da analise

# Filtra o tume para uma determinada especie
for (l in TUME.FILES){

	# Importa arquivo de um tume
	tume = read.csv(paste(TUME.PATH, l, sep=""), sep=",")
	
	# Salva vetor com nome das espécies contidas no TUME
	TUME.ESP <- levels(tume$Esp)
	
	tabela_resumo = data.frame()
	
	if (tume$I_meses[1] > 18 ){
		for (e in TUME.ESP){
			
			tume.esp = subset(tume, tume$Esp == e)
			
			if (nrow(na.omit(tume.esp[,1:6])) > 3 & is.numeric(tume.esp$DAP_cm) == TRUE){
				
				estH = hipsometrica(tume.esp)
				
				tabela_resumo = rbind(tabela_resumo, resumo_pos24(tume.esp, estH))
			} else {
				
				sem_dados <- data.frame(Especie = as.character(tume.esp$Esp[1]),
						Idade = 0,
						DAPmed = 0,
						DAPsd = 0,
						Hmed = 0,
						Hsd = 0,
						Hdom = 0,
						G = 0,
						V = 0,
						N = 0,
						Sobre = 0,
						IMA = 0)
				tabela_resumo <- rbind(tabela_resumo, sem_dados)
				
			}
		}
	
		#write.csv(tabela_resumo, file = paste(TUME.OUT, l, ".csv", sep = ""))
		write.csv(tabela_resumo, file = paste(TUME.OUT, "out_", l, sep = ""))
		
		plotIMA(tabela_resumo, l)
	} else {
		
		for (e in TUME.ESP){
			
			tume.esp = subset(tume, tume$Esp == e)
			
			if (nrow(na.omit(tume.esp[,c(1:5, 7)])) > 3){
				
				tabela_resumo = rbind(tabela_resumo, resumo_pre24(tume.esp))
				
			} else {
				
				sem_dados <- data.frame(Especie = as.character(tume.esp$Esp[1]),
						Idade = 0,
						Hmed = 0,
						Hsd = 0,
						N = 0,
						Sobre = 0)
				tabela_resumo <- rbind(tabela_resumo, sem_dados)
				
			}
		}
		
		#write.csv(tabela_resumo, file = paste(TUME.OUT, l, ".csv", sep = ""))
		write.csv(tabela_resumo, file = paste(TUME.OUT, "out_", l, sep = ""))
		
		plotSobre(tabela_resumo, l)
		
		
	}

}


# IMPLEMENTACOES RECENTES
# Stick 4: Criar codigo para idade < 18 - Implementado em 20140220
# Stick 2: Acrescentar na tabela resumo sd.DAP e sd.Ht - Implementado em 20140221 
# Stick 5: Acrescentar coluna para idade na tabela resumo - Implementado em 20140224
#
# Implementados mas com problemas nos códigos dos arquivos originais do TUME 
# Stick 1: Calcular sobrevivência como Cont(F & M)/max(cova) - Implementado em 20140222
# Stick 3: Calculo da altura dominante para media(H|cod==6) - Implementado em 20140222 
