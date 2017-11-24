setwd("/home/steffen/Documents/science/papers/submitted/Habitat properties are key drivers of Borrelia burgdorferi s.l. prevalence in Ixodes ricinus populations of deciduous forest fragments_2017_SE/")

# read in helper functions
files <- list.files(path = "./analyses/helper functions", pattern = "[.]R$", recursive=T,
                    full.names = TRUE)
for(i in 1:length(files)) source(files[i])
# load in preprocessed data (the file which should be available at the location
# given in the paper)
source("./analyses/prepare_smallFOREST-data.R")

colour1 <- "#b8d3e6"
colour2 <- "#4683ae"

################################################################################
# test differences between regions
nip_aov <- aov(n_ip_gm ~ id_region, data = all)
tHSD <- TukeyHSD(nip_aov, ordered = FALSE, conf.level = 0.95)
nip_levels <- tHSD[["id_region"]][,4]
nip_levels <- multcompLetters(nip_levels)
nip_letters <- nip_levels$Letters
nip_letters <- nip_letters[order(names(nip_letters))]

nip_tukey <- ddply(all, .(id_region), summarise,
                   id_stage = "nymphs",
                   mean = mean(n_ip_gm),
                   sd = sd(n_ip_gm),
                   number = length(n_ip_gm),
                   lower = mean - qt(0.975, number)*(sd/sqrt(number)),
                   upper = mean + qt(0.975, number)*(sd/sqrt(number)),
                   pos_y = upper + 0.03) %>%
  cbind(letters = nip_letters)

aip_aov <- aov(a_ip_gm ~ id_region, data = all)
tHSD <- TukeyHSD(aip_aov, ordered = FALSE, conf.level = 0.95)
aip_levels <- tHSD[["id_region"]][,4]
aip_levels <- multcompLetters(aip_levels)
aip_letters <- aip_levels$Letters
aip_letters <- aip_letters[order(names(aip_letters))]

aip_tukey <- ddply(all, .(id_region), summarise,
                   id_stage = "adults",
                   mean = mean(a_ip_gm),
                   sd = sd(a_ip_gm),
                   number = length(a_ip_gm),
                   lower = mean - qt(0.975, number)*(sd/sqrt(number)),
                   upper = mean + qt(0.975, number)*(sd/sqrt(number)),
                   pos_y = upper + 0.03) %>%
  cbind(letters = aip_letters)

# model building

options(contrasts=c("contr.SAS","contr.SAS")) # contr.sum = Statistica method, contr.SAS = SAS method

#---|Models|--------------------------------------------------------------------
# 1) nymphal Borrelia prevalence
variables <- colnames(all)[c(7, 10, 11, 20, 21, 30, 31, 40, 41, 48, 57:dim(all)[2])] # onle select relevant variables.
variables <- variables[-c(which(variables=="ms_l1mm_mat_bn_mean"),
                          which(variables=="id_soil_wr"),
                          which(variables=="doy"),
                          which(variables=="tree_hedera_patch"),
                          which(variables=="tree_picea_dom_patch"),
                          which(variables=="tree_ulmus_dom_patch"),
                          which(variables=="tree_acer_part_patch"),
                          which(variables=="tree_acer_dom_patch"),
                          which(variables=="tree_prunus_dom_patch"))]

n_ip_lmer <- lmer(lt_n_ip_gm ~ prop_edge5_patch +
                    lg_n_mean +
                    edge.density.forest +
                    FA_abundances_8 + I(FA_abundances_8^2) +
                    tree_prunus_part_patch +
                    lg_alpha_spl_total +
                    herb_large_disp_abund +
                    herb_regular +
                    prop_cult_100 + I(prop_cult_100^2) +
                    diam_diff_mean +
                    shrub_large_disp_alpha +
                    lg_lg_PROX.5000 +
                    tree_mixed_dom + I(tree_mixed_dom^2) +
                    tree_quercus_part_patch +
                    all_pr_padus_abund +
                    n_cold_year +
                    ba_laying_large_bn_mean +
                    ff_n.cont2_lg_mean +
                    FA_landscape_4 +
                    herb_evergreen_abund +
                    FA_taxonomic_2 +
                    lg_herb_chamaephyte +
                    rich_tree_large_mean +
                    (1 | id_region),
                  subset = all$n_pa!=0, data = all, REML = T)

# df2 <- df; df <- find.effects(n_ip_lmer, variables, degree = 2, stat = "Chisq", sbst = "p_val < 0.05 & r_sq < 0.5", order = "p_val")
# View(df)


# Anova(n_ip_lmer, type = "III", test.statistic = "Chisq"); summary(n_ip_lmer)$varcor$id_region>0
# semi.residual(n_ip_lmer, c("FA_traits_5"), levels = "id_region")
n_ip_decomp <- var.decomp(n_ip_lmer, ddf="Satterthwaite",  verbose = 1, penalised = T, SS_adjust = T); n_ip_decomp

# 2) adult Borrelia prevlence
variables <- colnames(all)[c(7, 10, 11, 20, 21, 30, 31, 40, 41, 57:dim(all)[2])]
variables <- variables[-c(which(variables=="ms_l1mm_mat_bn_mean"),
                          which(variables=="id_soilprop_edge20_patch_wr"),
                          which(variables=="doy"),
                          which(variables=="ba_bitterlich_mean"),
                          which(variables=="tree_hedera_patch"))]
a_ip_lmer <- lmer(lt_a_ip_gm ~ lg_a_mean + herb_average_abund + prop_cult_1000 + I(prop_cult_1000^2) + shrub_medium_disp_abund + trait11_mean + I(trait11_mean^2) + tree_sla + FA_soil_3 + FA_landscape_5 + I(FA_landscape_5^2) + ms_ph_mean + I(ms_ph_mean^2) + ff_bulkdensity_lg_mean + diss_tree2_spl + lg_alpha_spl_tall + all_large_disp_abund + tree_small_disp_alpha + I(tree_small_disp_alpha^2) + dens_large_lg_mean + ba_laying_total_lg_mean + FA_abundances_5 + I(FA_abundances_5^2) + FA_structure_5 + (1 | id_region), subset = a_pa!=0, data = all, REML = T)
# df2 <- df; df <- find.effects(a_ip_lmer, variables, degree = 2, stat = "Chisq", sbst = "p_val < 0.05 & r_sq < 0.5", order = "p_val"); View(df)
# Anova(a_ip_lmer, type = "III", test.statistic = "Chisq"); summary(a_ip_lmer)$varcor$id_region>0
# semi.residual(a_ip_lmer, c("FA_abundances_5"), levels = "id_region", f = 2)
# otl <- outliers(a_ip_lmer, "FA_abundances_5", f = 2)
a_ip_decomp <- var.decomp(a_ip_lmer, ddf="Satterthwaite", verbose = 1, penalised = T, SS_adjust = T); a_ip_decomp

#---|response profiles|---------------------------------------------------------
# 1) nymphs
nip_res_pro <- visreg.to.ggplot(n_ip_lmer, type = "contrast", levels = "id_region")
nip_names <- names[names(names)%in%levels(nip_res_pro$variable)]
pos_x_n <- ddply(subset(nip_res_pro, type=="fit"), .(variable), summarise, pos_x = max(x) - (max(x)-min(x))*0.5)
resp_nip <- ggplot(subset(nip_res_pro, type == "fit"), aes(x, y)) +
  geom_point(data = subset(nip_res_pro, type == "res"), aes(x, y), shape = 20, size = .3) +
  theme_bw() +
  geom_ribbon(aes(ymin = Lwr, ymax = Upr), bg = colour1) +
  geom_line(size = 1) +
  theme(aspect.ratio=1.2,
        strip.text = element_text(size = 11),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black")) +
  ylab("Logit( infection prevalence of nymphs )\n") + xlab("Response variable") +
  ylim(-4, 4) +
  facet_wrap(~variable, scales = "free_x", strip.position = "bottom", drop = F, labeller = as_labeller(nip_names), ncol = 6)

# 2) adults
aip_res_pro <- visreg.to.ggplot(a_ip_lmer, type = "contrast", levels = "id_region")
aip_names <- names[names(names)%in%levels(aip_res_pro$variable)]
pos_x_a <- ddply(subset(aip_res_pro, type=="fit"), .(variable), summarise, pos_x = max(x) - (max(x)-min(x))*0.5)
resp_aip <- ggplot(subset(aip_res_pro, type == "fit"), aes(x, y)) +
  geom_point(data = subset(aip_res_pro, type == "res"), aes(x, y), shape = 20, size = .3) +
  geom_ribbon(aes(ymin = Lwr, ymax = Upr), bg = colour2) +
  geom_line(size = 1) +
  theme_bw() +
  theme(aspect.ratio=1.2,
        strip.text = element_text(size = 11),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black")) +
  ylab("Logit( infection prevalence of adults )\n") + xlab("Response variable") +
  facet_wrap(~variable, scales = "free_x", strip.position = "bottom", drop = F, labeller = as_labeller(aip_names), ncol = 6)

#---|relative importance|-------------------------------------------------------
# 1) nymphs
n_meta <- n_ip_decomp$`Fixed Part`
n_meta <- cbind(n_meta, stage = "nymph")
n_meta$var2 <- as.character(rownames(n_meta))
n_meta <- join(n_meta, grouping[, c(1:2)], by = "var2")
n_meta$var1 <- ifelse(is.na(n_meta$var1), n_meta$var2, as.character(n_meta$var1))
n_meta <- join(n_meta, grouping[, c(1, 3)], by = "var1")
n_meta <- join(n_meta, grouping[, c(1, 4)], by = "var1")
n_meta <- join(n_meta, lut_group, by = "group")
n_meta <- join(n_meta, quantiles[, c(1, 3:9)], by = "var1")
colnames(n_meta)[which(colnames(n_meta)=="var1")] <- "var"
colnames(n_meta)[which(colnames(n_meta)=="var2")] <- "var_detail"

n_var <- as.data.frame(n_meta %>% group_by(var) %>% summarise(percent = round(sum(relImp_part*100, na.rm = T), 2)))
n_var$label <- sprintf("eta**2==%.1f", n_var$percent)
colnames(n_var)[which(colnames(n_var)=="var")] <- "variable"
n_var <- join(n_var, pos_x_n, by = "variable")
n_var$variable <- factor(n_var$variable, levels(nip_res_pro$variable)); n_var <- n_var[order(n_var$variable),]

n_drivergroup <- n_meta %>%
  group_by(group) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); n_drivergroup
n_group <- n_meta %>%
  group_by(metagroup) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); n_group
n_habitat <- n_meta %>%
  filter(type == "microhabitat" | type == "macrohabitat") %>%
  group_by(type) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); n_habitat

resp_nip <- resp_nip + geom_text(data=n_var, aes(x = pos_x, y = 3.3, label = label), parse = T)
svg(filename = "resp_nip.svg", height = 20, width = 11); resp_nip; dev.off()

# 2) adults
a_meta <- a_ip_decomp$`Fixed Part`
a_meta <- cbind(a_meta, stage = "adult")
a_meta$var2 <- as.character(rownames(a_meta))
a_meta <- join(a_meta, grouping[, c(1:2)], by = "var2")
a_meta$var1 <- ifelse(is.na(a_meta$var1), a_meta$var2, as.character(a_meta$var1))
a_meta <- join(a_meta, grouping[, c(1, 3)], by = "var1")
a_meta <- join(a_meta, grouping[, c(1, 4)], by = "var1")
a_meta <- join(a_meta, lut_group, by = "group")
a_meta <- join(a_meta, quantiles[, c(1, 3:9)], by = "var1")
colnames(a_meta)[which(colnames(a_meta)=="var1")] <- "var"
colnames(a_meta)[which(colnames(a_meta)=="var2")] <- "var_detail"

# a_meta$group <- factor(a_meta$group, levels(a_meta$group)[c(3, 2, 6, 1, 7, 5, 4, 8)])
a_var <- as.data.frame(a_meta %>% group_by(var) %>% summarise(percent = round(sum(relImp_part*100, na.rm=T), 2)))
a_var$label <- sprintf("eta**2==%.1f", a_var$percent)
colnames(a_var)[which(colnames(a_var)=="var")] <- "variable"
a_var <- join(a_var, pos_x_a, by = "variable")
a_var$variable <- factor(a_var$variable, levels(aip_res_pro$variable)); a_var <- a_var[order(a_var$variable),]

a_drivergroup <- a_meta %>%
  group_by(group) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); a_drivergroup
a_group <- a_meta %>%
  group_by(metagroup) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); a_group
a_habitat <- a_meta %>%
  filter(type == "microhabitat" | type == "macrohabitat") %>%
  group_by(type) %>%
  summarise(percent = round(sum(relImp_part*100, na.rm = T), 2),
            mean = round(mean(relImp_part*100, na.rm = T), 2),
            sd = round(sd(relImp_part*100, na.rm = T), 2)); a_habitat

resp_aip <- resp_aip + geom_text(data=a_var, aes(x = pos_x, y = 3, label = label), parse = T)
svg(filename = "resp_aip.svg", height = 20, width = 11); resp_aip; dev.off()

################################################################################
# Graphs

# Fig. 2 -
tukey <- rbind(nip_tukey, aip_tukey)
rownames(tukey) <- NULL
tukey$id_region <- as.factor(tukey$id_region)
tukey$id_region <- factor(tukey$id_region, levels(tukey$id_region)[c(4, 3, 1, 6, 5, 8, 7, 2)])
tukey$id_stage <- as.factor(tukey$id_stage)
tukey$id_stage <- factor(tukey$id_stage, levels(tukey$id_stage)[c(2, 1)])
tukey <- tukey[order(tukey$id_region),]

ggplot(tukey, aes(x = id_region, colour = id_stage)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.9), width = 0.5 , size=1.5) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.9), size = 7, stroke = 1) +
  geom_text(aes(y = pos_y, label = tukey$letter), position = position_dodge(width = 0.9), size = 5) +
  theme_bw() +
  scale_colour_manual(name = "Tick stage", values = c(colour1, colour2), labels = c("Nymphs", "Adults")) +
  scale_x_discrete(labels=c("southern\nFrance", "northern\nFrance", "Belgium", "western\nGermany", "eastern\nGermany", "southern\nSweden", "central\nSweden", "Estonia")) +
  ylab(expression(paste(italic("Borrelia burgdorferi"), " s.l. prevalence (geometric mean)"))) +
  xlab("\nRegion") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.border = element_rect(colour = "black", size = 0.5),
        legend.key = element_rect(colour = "white"))

# Fig. A3 - infection prevalence (geometric mean)
tick_subset <- all[, c(1, 4, 8, 9, 11, 12, 14, 15, 17, 38, 39, 41, 42, 44, 45, 47)]
ticks_plotting <- reshape(tick_subset,
                          idvar = "id_patch",
                          varying = list(names(tick_subset)[c(3, 10)],
                                         names(tick_subset)[c(4, 11)],
                                         names(tick_subset)[c(7, 14)],
                                         names(tick_subset)[c(8, 15)],
                                         names(tick_subset)[c(6, 13)],
                                         names(tick_subset)[c(5, 12)],
                                         names(tick_subset)[c(9, 16)]),
                          v.names = c("sum_ticks", "mean_ticks", "inf_p", "inf_ab", "n_lab", "ticks_pa", "inf_pa"),
                          times = c("nymphs", "adults"),
                          timevar = "id_stage",
                          direction = "long")
ticks_plotting$id_region <- as.factor(ticks_plotting$id_region)
ticks_plotting$id_region <- factor(ticks_plotting$id_region, levels(ticks_plotting$id_region)[c(4, 3, 1, 6, 5, 8, 7, 2)])
ticks_plotting$id_stage <- as.factor(ticks_plotting$id_stage)
ticks_plotting$id_stage <- factor(ticks_plotting$id_stage, levels(ticks_plotting$id_stage)[c(2, 1)])
ticks_plotting <- ticks_plotting[order(ticks_plotting$id_patch, ticks_plotting$id_stage),]
rownames(ticks_plotting) <- NULL

ggplot(subset(ticks_plotting, ticks_pa!=0), aes(x = id_region, y = inf_p, fill = id_stage)) +
  geom_boxplot(outlier.colour = "white") +
  geom_jitter(aes(size = n_lab, colour = id_stage), alpha = 0.8, shape = 20, position = position_jitterdodge()) +
  theme_bw() +
  scale_colour_manual(name = "Tick stage", values = c(colour1, colour2), labels = c("Nymphs", "Adults")) +
  scale_fill_manual(guide = F, values = c("white", "white")) +
  scale_x_discrete(labels=c("southern\nFrance", "northern\nFrance", "Belgium", "western\nGermany", "eastern\nGermany", "southern\nSweden", "central\nSweden", "Estonia")) +
  scale_size_continuous(name = "Ticks tested") +
  guides(colour = guide_legend(override.aes = list(size = 5)),
         size = guide_legend(ncol = 2, override.aes = list(shape = 21))) +
  ylab(expression(paste(italic("Borrelia burgdorferi"), " s.l. prevalence (geometric mean)"))) +
  xlab("\nRegion") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.border = element_rect(colour = "black", size = 0.5),
        legend.key = element_rect(colour = "white"))

# Fig. 3 - fraction of patches
ticks_plotting_region <- ddply(ticks_plotting, .(id_region, id_stage), summarise,
                               ticks_col = sum(sum_ticks, na.rm = T),
                               ticks_mean = mean(mean_ticks, na.rm = T),
                               inf_ab_mean = mean(inf_ab, na.rm = T),
                               inf_ab_sd = sd(inf_ab, na.rm = T),
                               inf_ab_ub = inf_ab_mean + qnorm(.95)*inf_ab_sd/sqrt(length(inf_ab)),
                               inf_ab_lb = inf_ab_mean - qnorm(.95)*inf_ab_sd/sqrt(length(inf_ab)),
                               ticks_lab_sum = sum(n_lab, na.rm = T),
                               ticks_lab_mean = mean(n_lab, na.rm = T),
                               inf_pr = mean(inf_p, na.rm = T),
                               frac_ticks = sum(ticks_pa)/length(ticks_pa)*100,
                               frac_g_inf = sum(inf_pa)/length(inf_pa)*100)

ggplot(ticks_plotting_region, aes(x = id_region)) +
  geom_bar(stat = "identity", aes(y = frac_g_inf, colour = id_stage, fill = id_stage), position = "dodge") +
  geom_bar(stat = "identity", aes(y = frac_ticks, colour = id_stage), fill = "transparent", position = "dodge") +
  theme_bw() +
  scale_x_discrete(labels=c("southern\nFrance", "northern\nFrance", "Belgium", "western\nGermany", "eastern\nGermany", "southern\nSweden", "central\nSweden", "Estonia")) +
  scale_fill_manual(name = "Tick stage", values = c(colour1, colour2), labels = c("Nymphs", "Adults")) +
  scale_colour_manual(name = "Tick stage", values = c(colour1, colour2), labels = c("Nymphs", "Adults")) +
  scale_alpha_manual(name = "patches\nwith", values = c(0, 1), labels = c("Ticks", "Borrelia\nin Ticks")) +
  guides(alpha = guide_legend(override.aes = list(colour = "gray30"))) +
  xlab("\nRegion") + ylab("Proportion of patches [%]\n") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        panel.border = element_rect(colour = "black", size = 0.5))

# Fig. 4 - relative importance
tick_stats <- rbind(n_meta, a_meta)

# relative importance of only habitat variables
tick_stats_habitat <- subset(tick_stats, subset = tick_stats$metagroup=="habitat")
tick_stats_habitat <- droplevels(tick_stats_habitat)

# relative importance of variable meta groups
tick_metagroup <- ddply(tick_stats, .(stage, metagroup), summarise,
                        relImp_part=sum(relImp_part*100))

tick_metagroup <- reshape(tick_metagroup, v.names = "relImp_part", timevar = "metagroup", idvar = c("stage"), direction = "wide")
tick_metagroup[is.na(tick_metagroup)] <- 0
tick_metagroup <- reshape(tick_metagroup, direction = "long"); row.names(tick_metagroup) <- NULL
tick_metagroup$metagroup <- factor(tick_metagroup$metagroup, levels(tick_metagroup$metagroup)[c(3, 2, 1, 4)]); tick_metagroup <- tick_metagroup[order(tick_metagroup$stage, tick_metagroup$metagroup),]

tick_metagroup_avrg <- ddply(tick_metagroup, .(metagroup), summarise,
                             relImp_avrg = mean(relImp_part))

g1 <- ggplot(tick_metagroup, aes(metagroup, relImp_part, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge(), size = .3) +
  theme_bw() +
  theme(legend.position="none") +
  xlab("\nDriver Group") + ylab("Relative importance [%]\n") +
  expand_limits(y=c(0,70)) +
  scale_fill_manual(name = "Prevalence", labels = c("Nymphs", "Adults"), values = c(colour1, colour2)) +
  scale_x_discrete(labels = c("Macroclimate\n", "Landscape", "Habitat", "Ontogeny")) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

g1_avrg <- ggplot(tick_metagroup_avrg, aes(metagroup, relImp_avrg)) +
  geom_bar(stat = "identity", size = .3) +
  theme_bw() +
  theme(legend.position="none") +
  xlab("\nDriver group") + ylab("Average relative importance [%]\n") +
  expand_limits(y=c(0,70)) +
  scale_x_discrete(labels = c("Macroclimate\n", "Landscape", "Habitat", "Ontogeny")) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
  )

# relative importance of scale within habitat
tick_habitat_scale <- ddply(tick_stats_habitat, .(stage, type), summarise,
                            relImp_part = sum(relImp_part*100))

tick_habitat_scale$type <- factor(tick_habitat_scale$type, levels(tick_habitat_scale$type)[c(1, 2)])
tick_habitat_scale <- tick_habitat_scale[order(tick_habitat_scale$stage, tick_habitat_scale$type),]

tick_habitat_scale_avrg <- ddply(tick_habitat_scale, .(type), summarise,
                                 relImp_avrg = mean(relImp_part))

g2 <- ggplot(tick_habitat_scale, aes(type, relImp_part, fill = stage)) +
  geom_bar(stat = "identity", position = "dodge", size = .3) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("\nscale within Habitatp") +
  expand_limits(y=c(0,70)) +
  scale_fill_manual(name = "Tick Stage", values = c(colour1, colour2), labels = c("Nymphs", "Adults")) +
  scale_x_discrete(labels = c("Macrohabitat\n", "Microhabitat\n")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

g2_avrg <- ggplot(tick_habitat_scale_avrg, aes(type, relImp_avrg)) +
  geom_bar(stat = "identity", size = .3) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("\n...within Habitatp") +
  expand_limits(y=c(0,70)) +
  scale_x_discrete(labels = c("Macrohabitat\n", "Microhabitat\n")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))


# relative importance of sub-group within habitat
tick_habitat_sub <- ddply(tick_stats_habitat, .(stage, group), summarise,
                          relImp_part = sum(relImp_part*100))

tick_habitat_sub <- reshape(tick_habitat_sub, v.names = "relImp_part", timevar = "group", idvar = c("stage"), direction = "wide")
tick_habitat_sub[is.na(tick_habitat_sub)] <- 0.01
tick_habitat_sub <- reshape(tick_habitat_sub, direction = "long"); row.names(tick_habitat_sub) <- NULL

tick_habitat_sub$group <- factor(tick_habitat_sub$group, levels(tick_habitat_sub$group)[c(1, 3, 4, 2)])
tick_habitat_sub <- tick_habitat_sub[order(tick_habitat_sub$stage, tick_habitat_sub$group),]

tick_habitat_sub_avrg <- ddply(tick_habitat_sub, .(group), summarise,
                               relImp_avrg = mean(relImp_part))


g3 <- ggplot(tick_habitat_sub, aes(group, relImp_part, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge(), size = .3) +
  theme_bw() +
  xlab("\nsub group within Habitat") +
  expand_limits(y=c(0,70)) +
  scale_fill_manual(name = "Tick Stage", labels = c("Nymphs", "Adults"), values = c(colour1, colour2)) +
  scale_x_discrete(labels = c("Functional\nproperties", "Structural\nproperties", "Diversity", "Soil")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.background = element_rect(colour = "lightgray"),
        legend.position = c(.75,.85),
        legend.key = element_rect(colour = "white"))

g3_avrg <- ggplot(tick_habitat_sub_avrg, aes(group, relImp_avrg)) +
  geom_bar(stat = "identity", size = .3) +
  theme_bw() +
  xlab("\n...within Habitat") +
  expand_limits(y=c(0,70)) +
  scale_x_discrete(labels = c("Functional\nproperties", "Structural\nproperties", "Diversity", "Soil")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.background = element_rect(colour = "lightgray"),
        legend.position = c(.85,.85),
        legend.key = element_rect(colour = "white"))

multiplot(g1, g2, g3, cols = 3, layout = matrix(c(1, 1, 1, 1, 2, 2, 3, 3, 3, 3), nrow = 1))
multiplot(g1_avrg, g2_avrg, g3_avrg, cols = 3, layout = matrix(c(1, 1, 1, 1, 2, 2, 3, 3, 3, 3), nrow = 1))
