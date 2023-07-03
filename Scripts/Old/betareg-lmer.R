# Old models
total <- func.invasive |>
	filter(Functional == 'Total.inv')
total.mod <- betareg(Proportion.trans ~ Feeder*as_factor(Time) + Site,
											 data = total, link = c("logit"), 
											 link.phi = NULL, type = c("ML"))
total.mod <- lmer(logit(Proportion) ~ Feeder*as_factor(Time) + (1|Site),
											 data = total)
Anova(total.mod)
check_model(total.mod)
emmeans(total.mod, ~Feeder*Time)

browsed <- func.invasive |>
	filter(Functional == 'Browsed perennials')
browsed.mod <- betareg(Proportion.trans ~ Feeder*as_factor(Time) + Site,
											 data = browsed, link = c("logit"), 
											 link.phi = NULL, type = c("ML"))
browsed.mod <- lmer(logit(Proportion) ~ Feeder*as_factor(Time) + (1|Site),
											 data = browsed)
Anova(browsed.mod)
check_model(browsed.mod)
emmeans(browsed.mod, ~Feeder*Time)

ruderals <- func.invasive |>
	filter(Functional == 'Ruderal')
ruderals.mod <- betareg(Proportion.trans ~ Feeder*as_factor(Time) + Site,
											 data = ruderals, link = c("logit"), 
											 link.phi = NULL, type = c("ML"))
ruderals.mod <- lmer(logit(Proportion) ~ Feeder*as_factor(Time) + (1|Site),
											 data = ruderals)
Anova(ruderals.mod)
check_model(ruderals.mod)
emmeans(ruderals.mod, ~Feeder*Time)

resistant <- func.invasive |>
	filter(Functional == 'Resistant perennials')
resistant.mod <- betareg(Proportion.trans ~ Feeder*as_factor(Time) + Site,
											 data = resistant, link = c("logit"), 
											 link.phi = NULL, type = c("ML"))
resistant.mod <- lmer(logit(Proportion) ~ Feeder*as_factor(Time) + (1|Site),
											 data = resistant)
Anova(resistant.mod)
check_model(resistant.mod)
emmeans(resistant.mod, ~Feeder*Time)
