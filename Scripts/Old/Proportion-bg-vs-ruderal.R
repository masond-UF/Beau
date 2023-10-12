## --------------- JUST BG AND RUDERAL -----------------------------------------

bare.ground <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 2) |>
	filter(Meter < 4) |>
	dplyr::select(Plot, Direction, Species)

for(i in 1:nrow(bare.ground)){
	if(bare.ground$Species[i] != 'Bare Ground'){
		bare.ground$Species[i] <- 'Not Bare Ground'
	}
}

colnames(bare.ground)[3] <- 'Status'
	
bare.ground.sum <- bare.ground |>
	group_by(Plot, Status) |>
	summarize(Count = n())

bare.ground.sum <- bare.ground.sum |>
	pivot_wider(names_from = Status, values_from = Count)

bare.ground.sum <- bare.ground.sum |>
	mutate(Total = `Bare Ground` + `Not Bare Ground`,
				 Prop.bg = `Bare Ground`/Total)

# Fix WS8
bare.ground.sum[25, 3] <- 0
bare.ground.sum[25, 4] <- 24
bare.ground.sum[25, 5] <- 1.00

bare.ground.sum$Prop.bg <- round(bare.ground.sum$Prop.bg, 2)
bare.ground.sum <- bare.ground.sum[,c(1,5)]

rm(list=ls()[! ls() %in% c("bare.ground.sum")])

ruderal <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 5) |>
	filter(Meter < 4) |>
	dplyr::select(Plot, Direction, Functional)

for(i in 1:nrow(ruderal)){
	if(isTRUE(ruderal$Functional[i] != 'Ruderal')){
		ruderal$Functional[i] <- 'Not Ruderal'
	}
	if(is.na(ruderal$Functional[i])){
		ruderal$Functional[i] <- 'Not Ruderal'
	}
}

ruderal.sum <- ruderal |>
	group_by(Plot, Functional) |>
	summarize(Count = n())

ruderal.sum <- ruderal.sum |>
	pivot_wider(names_from = Functional, values_from = Count)

ruderal.sum$Ruderal[is.na(ruderal.sum$Ruderal)] <- 0

ruderal.sum <- ruderal.sum |>
	mutate(Total = Ruderal + `Not Ruderal`,
				 Prop.rud = Ruderal/Total)

ruderal.sum$Prop.rud <- round(ruderal.sum$Prop.rud, 2)
ruderal.sum <- ruderal.sum[,c(1,5)]

rm(list=ls()[! ls() %in% c("bare.ground.sum", 'ruderal.sum')])

d <- merge(bare.ground.sum, ruderal.sum)

d.no.zero <- d |>
	filter(Prop.rud > 0) 

plot(d.no.zero$Prop.bg, d.no.zero$Prop.rud)
cor.test(d.no.zero$Prop.bg, d.no.zero$Prop.rud) # 0.2277083

mod <- lm(Prop.rud~Prop.bg, data = d)
anova(mod)

library("ggpubr")
ggscatter(d.no.zero, x = "Prop.bg", y = "Prop.rud", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "Proportion Ruderals After", xlab = "Proportion Bare Ground Before")

## --------------- COMBINE BOTH ------------------------------------------------

# Clear the decks
rm(list=ls())

bare.ground <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 2) |>
	filter(Meter < 4) |>
	dplyr::select(Plot, Direction, Species, Functional)

for(i in 1:nrow(bare.ground)){
	if(isTRUE(bare.ground$Functional[i] == 'Ruderal')){
		bare.ground$Species[i] <- 'Bare Ground or Ruderal'
	}
	# if(isTRUE(bare.ground$Functional[i] == 'Litter')){
	# 	bare.ground$Species[i] <- 'Bare Ground or Ruderal'
	# }
	if(isTRUE(bare.ground$Species[i] == 'Bare Ground')){
		bare.ground$Species[i] <- 'Bare Ground or Ruderal'
	}
}

bare.ground <- bare.ground |> dplyr::select(-Functional)

for(i in 1:nrow(bare.ground)){
	if(bare.ground$Species[i] != 'Bare Ground or Ruderal'){
		bare.ground$Species[i] <- 'Not Bare Ground or Ruderal'
	}
}

colnames(bare.ground)[3] <- 'Status'
	
bare.ground.sum <- bare.ground |>
	group_by(Plot, Status) |>
	summarize(Count = n())

bare.ground.sum <- bare.ground.sum |>
	pivot_wider(names_from = Status, values_from = Count)

bare.ground.sum <- bare.ground.sum |>
	mutate(Total = `Bare Ground or Ruderal` + `Not Bare Ground or Ruderal`,
				 Prop.bg = `Bare Ground or Ruderal`/Total)

# Fix WS8
bare.ground.sum[25, 3] <- 0
bare.ground.sum[25, 4] <- 24
bare.ground.sum[25, 5] <- 1.00

bare.ground.sum$Prop.bg <- round(bare.ground.sum$Prop.bg, 2)
bare.ground.sum <- bare.ground.sum[,c(1,5)]

rm(list=ls()[! ls() %in% c("bare.ground.sum")])

ruderal <- read.csv('Clean-data/Feeder-veg.csv') |>
	filter(Dataset == 'Manipulative') |>
	filter(Feeder == 'Y') |>
	filter(Time == 5) |>
	filter(Meter < 4) |>
	dplyr::select(Plot, Direction, Species, Functional)

for(i in 1:nrow(ruderal)){
	if(isTRUE(ruderal$Functional[i] == 'Ruderal')){
		ruderal$Functional[i] <- 'Bare Ground or Ruderal'
	}
	# if(isTRUE(ruderal$Species[i] == 'Bare Ground')){
	# 	ruderal$Functional[i] <- 'Bare Ground or Ruderal'
	# }
}

ruderal <- ruderal |> dplyr::select(-Species)

for(i in 1:nrow(ruderal)){
	if(isTRUE(ruderal$Functional[i] != 'Bare Ground or Ruderal')){
		ruderal$Functional[i] <- 'Not Bare Ground or Ruderal'
	}
	if(is.na(ruderal$Functional[i])){
		ruderal$Functional[i] <- 'Not Bare Ground or Ruderal'
	}
}

ruderal.sum <- ruderal |>
	group_by(Plot, Functional) |>
	summarize(Count = n())

ruderal.sum <- ruderal.sum |>
	pivot_wider(names_from = Functional, values_from = Count)

ruderal.sum$`Bare Ground or Ruderal`[is.na(ruderal.sum$`Bare Ground or Ruderal`)] <- 0

ruderal.sum <- ruderal.sum |>
	mutate(Total = `Bare Ground or Ruderal` + `Not Bare Ground or Ruderal`,
				 Prop.rud = `Bare Ground or Ruderal`/Total)

ruderal.sum$Prop.rud <- round(ruderal.sum$Prop.rud, 2)
ruderal.sum <- ruderal.sum[,c(1,5)]

rm(list=ls()[! ls() %in% c("bare.ground.sum", 'ruderal.sum')])

d <- merge(bare.ground.sum, ruderal.sum)


d.no.zero <- d |>
	filter(Prop.rud > 0) 

plot(d.no.zero$Prop.bg, d.no.zero$Prop.rud)
cor.test(d.no.zero$Prop.bg, d.no.zero$Prop.rud)

library("ggpubr")
ggscatter(d.no.zero, x = "Prop.bg", y = "Prop.rud", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "Proportion Ruderals After", xlab = "Proportion Bare Ground Before")

# Burned 
# LV6 (1x), LV3 (2x), LV1 (1x)

# BINARY 
for(i in 1:nrow(d)){
	if(d$Prop.rud[i] > 0) {
		 d$Prop.rud[i] <- 1
	}
	if(d$Prop.rud[i] == 0) {
		 d$Prop.rud[i] <- 0
	}
}
