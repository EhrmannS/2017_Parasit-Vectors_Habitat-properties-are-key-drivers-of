################################################################################
# Introduction
#
# This script is used to prepare data gathered in the smallFOREST framework for
# analyses of the 16 focal patches per landscape window. Framework: smallFOREST
# Based on field-data collected in the framework by Steffen Ehrmann and
# field-assistants (Katja Leischke, Peter # Fräßdorf, Regina Hesseman, Iris
# Gutierrez), the smallFOREST site-managers and on spatial data of the
# smallFOREST geospatial database.
#
################################################################################
# License
#
# Copyright (C) 2016 Steffen Ehrmann
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful , but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.
#
#
################################################################################
# read in packages

# plotting
require("ggplot2"); require("grid"); require("pascal"); require("gridExtra")
require("multcompView")
# data
require("plyr"); require("dplyr"); require("reshape")
# modelling
require("lmerTest"); require("lme4"); require("leaps"); require("car"); require("visreg")
# ordination
require("psych"); require("GPArotation"); require("ellipse")

################################################################################
# read in data

setwd("/home/steffen/Documents/science/papers/submitted/Habitat properties are key drivers of Borrelia burgdorferi s.l. prevalence in Ixodes ricinus populations of deciduous forest fragments_2017_SE/")

# has been created by hand after all transformations were carried out.
colnames_transformations <- read.csv("./analyses/intermediate data/colnames_transformations.csv")

# this following file already contains some transformation, which were applied
# already during the calculation of stand structural and soil properties.
main_patches <- read.csv("/home/steffen/Projects/smallFOREST/output/FINAL_trans_wp2.csv")
grouping <- read.csv("./analyses/intermediate data/grouping.csv")


# select only rows which have complete cases
tick_patches <- main_patches[complete.cases(main_patches[,c(which(colnames(main_patches)=="l_mean"),
                                                            which(colnames(main_patches)=="n_mean"),
                                                            which(colnames(main_patches)=="a_mean"),
                                                            which(colnames(main_patches)=="rH_5"):which(colnames(main_patches)=="satdef"),
                                                            which(colnames(main_patches)=="trait11_mean"):which(colnames(main_patches)=="trait47_mean"))]),]

################################################################################
# View distributions of raw data

# main_patches2 <- read.csv("./Data Output/FINAL_wp2.csv"); main_patches2 <- main_patches2[,-1]
# variables <- colnames(main_patches2)
# # remove some unwanted variables
# variables <- variables[-c(which(colnames(main_patches2)=="dom_tree_small"),
#                           which(colnames(main_patches2)=="dom_tree_large"),
#                           which(colnames(main_patches2)=="pooled"),
#                           which(colnames(main_patches2)=="id_soil_wr"),
#                           which(colnames(main_patches2)=="soiltype"))]
#
# for(i in 5:length(variables)){
#   variable_of_interest <- variables[i]
#   svg(filename = paste("id_region - ", variable_of_interest, sep=""))
#
#   melt <- melt(main_patches2[, c(1:4, which(colnames(main_patches2)==variable_of_interest))],
#               id=c("id_patch", "id_window", "window", "id_region"))
#   g1 <- ggplot(melt, aes(x = value)) +
#     geom_histogram(aes(y=..density..), colour="black", fill="white") +
#     geom_density(fill="red", alpha=.1) +
#     ggtitle(paste(variable_of_interest, sep="")) +
#     facet_wrap(~variable, scales = "free")
#   g2 <- ggplot(dat = melt, aes(x = id_region, y = value, colour = id_region)) +
#     geom_boxplot() +
#     geom_jitter(width = 0.1) +
#     coord_flip() +
#     theme(legend.position = "none")
#   multiplot(g1, g2, cols = 1)
# dev.off()
# }

################################################################################
# put together data-frames

#---|ticks|---------------------------------------------------------------------
# l = larvae; n = nymphs; a = adults; f = female, m = male; doy = day of the
# year; lg = refers to data being log-transformed; pa = presence/absence; ip =
# infection prevalence; ia = abundance of inected ticks; mir = minimum infection
# rate (set to NULL, since not analysed); gm = geometric mean of prevalence
# probabilities (see manuscript); *_min = smallest recorded value; lt = refers
# to data being logit-transformed.
#-------------------------------------------------------------------------------
# 1) assemble sub-dataframe
# --->
ticks <- tick_patches[,c(which(colnames(tick_patches)=="id_patch"):which(colnames(tick_patches)=="id_region"),
                         which(colnames(tick_patches)=="l_mean"):which(colnames(tick_patches)=="doy")
)]

# 2) take out outlayers and apply transformations based on distributions of raw
# data
# --->
ticks <- data.frame(append(ticks, list(lg_l_mean=log10(ticks$l_mean+1)), after=match("l_mean", names(ticks))))
ticks <- data.frame(append(ticks, list(l_pa=ifelse(ticks$l_mean!=0, 1, 0)), after=match("lg_l_mean", names(ticks))))
ticks <- data.frame(append(ticks, list(lg_n_mean=log10(ticks$n_mean+1)), after=match("n_mean", names(ticks))))
ticks <- data.frame(append(ticks, list(n_pa=ifelse(ticks$n_mean!=0, 1, 0)), after=match("lg_n_mean", names(ticks))))
ticks <- data.frame(append(ticks, list(lg_f_mean=log10(ticks$f_mean+1)), after=match("f_mean", names(ticks))))
ticks <- data.frame(append(ticks, list(f_pa=ifelse(ticks$f_mean!=0, 1, 0)), after=match("lg_f_mean", names(ticks))))
ticks <- data.frame(append(ticks, list(lg_m_mean=log10(ticks$m_mean+1)), after=match("m_mean", names(ticks))))
ticks <- data.frame(append(ticks, list(m_pa=ifelse(ticks$m_mean!=0, 1, 0)), after=match("lg_m_mean", names(ticks))))
ticks <- data.frame(append(ticks, list(lg_a_mean=log10(ticks$a_mean+1)), after=match("a_mean", names(ticks))))
ticks <- data.frame(append(ticks, list(a_pa=ifelse(ticks$a_mean!=0, 1, 0)), after=match("lg_a_mean", names(ticks))))

ticks$n_ip_mir <- NULL
ticks$n_ia_mir <- NULL
n_min <- min(ticks$n_ip_gm[ticks$n_ip_gm!=0])*0.9
ticks$n_ip_gm <- round(ifelse(ticks$n_pa!=0, ifelse(ticks$n_ip_gm!=0, ticks$n_ip_gm, n_min/ticks$n_lab_sum), 0), 2)
ticks$lt_n_ip_gm <- log(ticks$n_ip_gm+0.009/(1.009-ticks$n_ip_gm))
ticks$n_ia_gm <- round(ticks$n_mean*ticks$n_ip_gm)
ticks$f_ip_mir <- NULL
ticks$f_ia_mir <- NULL
f_min <- min(ticks$f_ip_gm[ticks$f_ip_gm!=0])*0.9
ticks$f_ip_gm <- round(ifelse(ticks$f_pa!=0, ifelse(ticks$f_ip_gm!=0, ticks$f_ip_gm, f_min/ticks$f_lab_sum), 0), 2)
ticks$lt_f_ip_gm <- log(ticks$f_ip_gm+0.009/(1.009-ticks$f_ip_gm))
ticks$f_ia_gm <- round(ticks$f_mean*ticks$f_ip_gm)
ticks$m_ip_mir <- NULL
ticks$m_ia_mir <- NULL
m_min <- min(ticks$m_ip_gm[ticks$m_ip_gm!=0])-0.01
ticks$m_ip_gm <- round(ifelse(ticks$m_pa!=0, ifelse(ticks$m_ip_gm!=0, ticks$m_ip_gm, m_min/ticks$m_lab_sum), 0), 2)
ticks$lt_m_ip_gm <- log(ticks$m_ip_gm+0.009/(1.009-ticks$m_ip_gm))
ticks$m_ia_gm <- round(ticks$m_mean*ticks$m_ip_gm)
ticks$a_ip_mir <- NULL
ticks$a_ia_mir <- NULL
a_min <- min(ticks$a_ip_gm[ticks$a_ip_gm!=0])-0.01
ticks$a_ip_gm <- round(ifelse(ticks$a_pa!=0, ifelse(ticks$a_ip_gm!=0, ticks$a_ip_gm, a_min/ticks$a_lab_sum), 0), 2)
ticks$lt_a_ip_gm <- log(ticks$a_ip_gm+0.009/(1.009-ticks$a_ip_gm))
ticks$a_ia_gm <- round(ticks$a_mean*ticks$a_ip_gm)

ticks <- data.frame(append(ticks, list(lg_n_ia_gm=log10(ticks$n_ia_gm+1)), after=match("n_ia_gm", names(ticks))))
ticks <- data.frame(append(ticks, list(lg_f_ia_gm=log10(ticks$f_ia_gm+1)), after=match("f_ia_gm", names(ticks))))
ticks <- data.frame(append(ticks, list(lg_m_ia_gm=log10(ticks$m_ia_gm+1)), after=match("m_ia_gm", names(ticks))))
ticks <- data.frame(append(ticks, list(lg_a_ia_gm=log10(ticks$a_ia_gm+1)), after=match("a_ia_gm", names(ticks))))

# 3) analytical graphs
# --->
# pairs.panels(ticks[, c(5:dim(ticks)[2])], bg=c("yellow", "blue")[ticks$window], pch=21, main="ticks", hist.col="gray")
# cor.ticks <- cor(ticks[, c(5:dim(ticks)[2])], use="complete.obs"); plotcorr(cor.ticks, col=rgb(colorfun((cor.ticks+1)/2), maxColorValue=255), mar = c(0.1, 0.1, 0.1, 0.1))
# ticks_melt <- melt(ticks_ab, id=c("id_patch", "id_window", "window", "id_region"))
# ticks.gg <- ggplot(ticks_melt, aes(x=value)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(fill="red", alpha=.1) + ggtitle("distributions tick-data") + facet_wrap(~variable, scales = "free")
# ticks.gg
# scatter <- ggplot(dat=all, aes(x=ms_ph_mean, y=log10(n_col+1))) + geom_point(aes(colour=id_region)) + geom_point(shape=1) + geom_smooth() + geom_smooth(method=lm, se=F, colour="red") + scale_colour_manual(values=colours) + theme(legend.position=c(1,1), legend.justification=c(1,1)); g <- scatter; g
# scatter2 <- ggplot(dat=main_patches, aes(x=, y=, colour=age3)) + geom_point() + geom_smooth() + geom_smooth(method=lm, se=F, colour="red") + scale_color_gradient(low="darkkhaki", high="darkgreen"); g <- scatter2; g

# 4) factorial analysis
# --->
# not happening here, since this is the response

#---|stand structure|-----------------------------------------------------------
#
#-------------------------------------------------------------------------------
# 1) assemble sub-dataframe
# --->
structure <- tick_patches[,c(which(colnames(tick_patches)=="id_patch"):which(colnames(tick_patches)=="id_region"),
                       which(colnames(tick_patches)=="height_stand_mean"):which(colnames(tick_patches)=="ba_laying_large_bn_mean"),
                       which(colnames(tick_patches)=="species_mingling_mean"):which(colnames(tick_patches)=="biomass_total_lg_mean"),
                       which(colnames(tick_patches)=="diss_tree1_tree2"):which(colnames(tick_patches)=="diss_tree2_spl"),
                       which(colnames(tick_patches)=="herb_canopy_height"):which(colnames(tick_patches)=="herb_canopy_height_cv"),
                       which(colnames(tick_patches)=="tree_canopy_height"):which(colnames(tick_patches)=="tree_canopy_height_cv")
)]
structure <- structure[, -c(which(colnames(structure)=="dw_num_total_lg_mean"),
                            which(colnames(structure)=="dw_num_large_bn_mean"),
                            which(colnames(structure)=="gap_fraction_mean"),
                            which(colnames(structure)=="dw_diam_max_lg_mean"),
                            which(colnames(structure)=="dw_piles_bn_sum")
)]

# 2) take out outlayers and apply transformations based on distributions of raw
# data
# --->
# since I checked that all data are calculated correctly, I will replace really
# huge/small values with values similar to respective values (maintaining the
# order of values).
structure$diss_tree1_spl[is.na(structure$diss_tree1_spl)] <- min(structure$diss_tree1_spl, na.rm=T)-0.0001
structure$diss_tree2_spl[structure$diss_tree2_spl==min(structure$diss_tree2_spl, na.rm=T)] <- NA
structure$diss_tree2_spl[is.na(structure$diss_tree2_spl)] <- min(structure$diss_tree2_spl, na.rm=T)-0.0001
structure$dens_total_lg_mean[structure$dens_total_lg_mean==max(structure$dens_total_lg_mean, na.rm=T)] <- -12345
structure$dens_total_lg_mean[structure$dens_total_lg_mean==max(structure$dens_total_lg_mean, na.rm=T)] <- -23456
structure$dens_total_lg_mean[structure$dens_total_lg_mean==-23456] <- max(structure$dens_total_lg_mean, na.rm=T)+.1 # add 1 to keep the ranking
structure$dens_total_lg_mean[structure$dens_total_lg_mean==-12345] <- max(structure$dens_total_lg_mean, na.rm=T)+.1 # add 1 to keep the ranking
structure$dens_large_lg_mean[structure$dens_large_lg_mean==min(structure$dens_large_lg_mean, na.rm=T)] <- 12345
structure$dens_large_lg_mean[structure$dens_large_lg_mean==min(structure$dens_large_lg_mean, na.rm=T)] <- 23456
structure$dens_large_lg_mean[structure$dens_large_lg_mean==23456] <- min(structure$dens_large_lg_mean, na.rm=T)-.1 # add 1 to keep the ranking
structure$dens_large_lg_mean[structure$dens_large_lg_mean==12345] <- min(structure$dens_large_lg_mean, na.rm=T)-.1 # add 1 to keep the ranking
structure$slenderness_small_mean[structure$slenderness_small_mean==min(structure$slenderness_small_mean, na.rm=T)] <- NA
structure$slenderness_small_mean[is.na(structure$slenderness_small_mean)] <- min(structure$slenderness_small_mean, na.rm=T)-.01 # add 1 to keep the ranking
structure$slenderness_large_mean[structure$slenderness_large_mean==0] <- NA
structure$slenderness_large_mean[is.na(structure$slenderness_large_mean)] <- min(structure$slenderness_large_mean, na.rm=T)-.001 # since they are all 0 there is no ranking to abide to.
structure$ba_total_lg_mean[structure$ba_total_lg_mean==max(structure$ba_total_lg_mean, na.rm=T)] <- NA
structure$ba_total_lg_mean[is.na(structure$ba_total_lg_mean)] <- max(structure$ba_total_lg_mean, na.rm=T)+.1 # add 1 to keep the ranking
structure$ba_dead_bn_mean <- ifelse(structure$ba_dead_bn_mean==0, 0, 1); colnames(structure)[which(colnames(structure)=="ba_dead_bn_mean")] <- "bn_ba_dead_bn_mean"
structure$diam_iqr <- log10(structure$diam_iqr+1); colnames(structure)[which(colnames(structure)=="diam_iqr")] <- "lg_diam_iqr"
structure$diam_biggest <- log10(structure$diam_biggest+1); colnames(structure)[which(colnames(structure)=="diam_biggest")] <- "lg_diam_biggest"
structure$diam_cv_mean <- log10(structure$diam_cv_mean+1); colnames(structure)[which(colnames(structure)=="diam_cv_mean")] <- "lg_diam_cv_mean"
structure$ba_evergr_bn_mean <- ifelse(structure$ba_evergr_bn_mean==0, 0, 1); colnames(structure)[which(colnames(structure)=="ba_evergr_bn_mean")] <- "bn_ba_evergr_bn_mean"
structure$tree_canopy_height_cv[structure$tree_canopy_height_cv==max(structure$tree_canopy_height_cv, na.rm=T)] <- -12345
structure$tree_canopy_height_cv[structure$tree_canopy_height_cv==max(structure$tree_canopy_height_cv, na.rm=T)] <- -23456
structure$tree_canopy_height_cv[structure$tree_canopy_height_cv==-23456] <- max(structure$tree_canopy_height_cv, na.rm=T)+.1 # add 1 to keep the ranking
structure$tree_canopy_height_cv[structure$tree_canopy_height_cv==-12345] <- max(structure$tree_canopy_height_cv, na.rm=T)+.1 # add 1 to keep the ranking
structure$tree_canopy_height_cv <- log10(structure$tree_canopy_height_cv+1); colnames(structure)[which(colnames(structure)=="tree_canopy_height_cv")] <- "lg_tree_canopy_height_cv"

# 3) analytical graphs
# --->
# structure_melt1 <- melt(structure[,c(1:50)], id=c("id_patch", "id_window", "window", "id_region"))
# structure.gg1 <- ggplot(structure_melt1, aes(x=value)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(fill="red", alpha=.1) + ggtitle("distributions structure-data 1") + facet_wrap(~variable, scales = "free")
# structure.gg1
# structure_melt2 <- melt(structure[,c(1:5, 51:dim(structure)[2])], id=c("id_patch", "id_window", "window", "id_region"))
# structure.gg2 <- ggplot(structure_melt2, aes(x=value)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(fill="red", alpha=.1) + ggtitle("distributions structure-data 2") + facet_wrap(~variable, scales = "free")
# structure.gg2

# 4) factorial analysis
# --->
# fa.parallel(structure[c(5:dim(structure)[2])])
structure_fa <- fa(structure[c(5:dim(structure)[2])], nfactors=8, rotate="varimax", fm="ml", scores="regression")
# fa.graph(structure_fa, cut=.5, main = "factors structure")
factors_structure <- data.frame(structure_fa$scores); colnames(factors_structure) <- c(paste("FA_structure_", 1:dim(factors_structure)[2], sep=""))

#---|abundance within vegetation layers|----------------------------------------
#
#-------------------------------------------------------------------------------
# 1) assemble sub-dataframe
# --->
abundances <- tick_patches[,c(which(colnames(tick_patches)=="id_patch"):which(colnames(tick_patches)=="id_region"),                       which(colnames(tick_patches)=="herb_total_abund"):which(colnames(tick_patches)=="herb_average_abund"),
                              which(colnames(tick_patches)=="herb_small_disp_abund"):which(colnames(tick_patches)=="herb_large_disp_abund"),
                              which(colnames(tick_patches)=="herb_w_nuts_abund"):which(colnames(tick_patches)=="herb_pr_padus_abund"),
                              which(colnames(tick_patches)=="shrub_total_abund"):which(colnames(tick_patches)=="shrub_average_abund"),
                              which(colnames(tick_patches)=="shrub_small_disp_abund"):which(colnames(tick_patches)=="shrub_large_disp_abund"),
                              which(colnames(tick_patches)=="shrub_w_nuts_abund"):which(colnames(tick_patches)=="shrub_pr_padus_abund"),
                              which(colnames(tick_patches)=="tree_total_abund"):which(colnames(tick_patches)=="tree_average_abund"),
                              which(colnames(tick_patches)=="tree_small_disp_abund"):which(colnames(tick_patches)=="tree_large_disp_abund"),
                              which(colnames(tick_patches)=="tree_w_nuts_abund"):which(colnames(tick_patches)=="tree_pr_padus_abund"),
                              which(colnames(tick_patches)=="all_total_abund"):which(colnames(tick_patches)=="ratio_shrub_tree_abund"),
                              which(colnames(tick_patches)=="all_small_disp_abund"):which(colnames(tick_patches)=="all_large_disp_abund"),
                              which(colnames(tick_patches)=="all_w_nuts_abund"):which(colnames(tick_patches)=="all_pr_padus_abund")
)]

abundances <- abundances[, -c(which(colnames(abundances)=="all_w_legume_abund"),
                              which(colnames(abundances)=="shrub_w_legume_abund"),
                              which(colnames(abundances)=="tree_w_legume_abund")
)]
# 2) take out outlayers and apply transformations based on distributions of raw
# data
# --->
abundances$ratio_herb_shrub_abund <- log10(abundances$ratio_herb_shrub_abund+1); colnames(abundances)[which(colnames(abundances)=="ratio_herb_shrub_abund")] <- "lg_ratio_herb_shrub_abund"
abundances$ratio_shrub_tree_abund <- log10(abundances$ratio_shrub_tree_abund+1); colnames(abundances)[which(colnames(abundances)=="ratio_shrub_tree_abund")] <- "lg_ratio_shrub_tree_abund"

# 3) analytical graphs
# --->
# abundances_melt <- melt(abundances[,c(1:dim(abundances)[2])], id=c("id_patch", "id_window", "window", "id_region"))
# abundances.gg1 <- ggplot(abundances_melt, aes(x=value)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(fill="red", alpha=.1) + ggtitle("distributions abundance-data 1") + facet_wrap(~variable, scales = "free")
# abundances.gg1
# cor.abundances <- cor(abundances[,c(5:dim(abundances)[2])])
# plotcorr(cor.abundances, col=rgb(colorfun((cor.abundances+1)/2), maxColorValue=255), mar = c(0.1, 0.1, 0.1, 0.1))

# 4) factorial analysis
# --->
# fa.parallel(abundances[c(5:dim(abundances)[2])])
abundances_fa <- fa(abundances[c(5:dim(abundances)[2])], nfactors=8, rotate="varimax", fm="ml", scores="regression")
# fa.graph(abundances_fa, cut=.5, main = "factors abundances")
factors_abundances <- data.frame(abundances_fa$scores); colnames(factors_abundances) <- c(paste("FA_abundances_", 1:dim(factors_abundances)[2], sep=""))

#---|plant traits|--------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# 1) assemble sub-dataframe
# --->
traits <- tick_patches[,c(which(colnames(tick_patches)=="id_patch"):which(colnames(tick_patches)=="id_region"),
                          which(colnames(tick_patches)=="herb_weight_disp"),
                          which(colnames(tick_patches)=="herb_small_disp_alpha"):which(colnames(tick_patches)=="herb_large_disp_alpha"),
                          which(colnames(tick_patches)=="herb_sla"):which(colnames(tick_patches)=="herb_w_branching"),
                          which(colnames(tick_patches)=="shrub_weight_disp"),
                          which(colnames(tick_patches)=="shrub_small_disp_alpha"):which(colnames(tick_patches)=="shrub_large_disp_alpha"),
                          which(colnames(tick_patches)=="shrub_sla"),
                          which(colnames(tick_patches)=="tree_weight_disp"),
                          which(colnames(tick_patches)=="tree_small_disp_alpha"):which(colnames(tick_patches)=="tree_large_disp_alpha"),
                          which(colnames(tick_patches)=="tree_sla"),
                          which(colnames(tick_patches)=="all_weight_disp"),
                          which(colnames(tick_patches)=="all_small_disp_alpha"):which(colnames(tick_patches)=="all_large_disp_alpha"),
                          which(colnames(tick_patches)=="all_sla")
)]

# 2) take out outlayers and apply transformations based on distributions of raw
# data
# --->
traits$herb_weight_disp <- log10(traits$herb_weight_disp+1); colnames(traits)[which(colnames(traits)=="herb_weight_disp")] <- "lg_herb_weight_disp"
traits$herb_chamaephyte <- log10(traits$herb_chamaephyte+1); colnames(traits)[which(colnames(traits)=="herb_chamaephyte")] <- "lg_herb_chamaephyte"
traits$herb_therophyte <- log10(traits$herb_therophyte+1); colnames(traits)[which(colnames(traits)=="herb_therophyte")] <- "lg_herb_therophyte"
traits$herb_geophyte <- log10(traits$herb_geophyte+1); colnames(traits)[which(colnames(traits)=="herb_geophyte")] <- "lg_herb_geophyte"
traits$herb_grass <- log10(traits$herb_grass+1); colnames(traits)[which(colnames(traits)=="herb_grass")] <- "lg_herb_grass"
traits$herb_fern <- ifelse(traits$herb_fern==0, 0, 1); colnames(traits)[which(colnames(traits)=="herb_fern")] <- "bn_herb_fern"
traits$herb_prostrate <- ifelse(traits$herb_prostrate==0, 0, 1); colnames(traits)[which(colnames(traits)=="herb_prostrate")] <- "bn_herb_prostrate"
traits$herb_erect[traits$herb_erect==min(traits$herb_erect, na.rm=T)] <- NA
traits$herb_erect[is.na(traits$herb_erect)] <- min(traits$herb_erect, na.rm=T)-0.1
traits$shrub_weight_disp <- log10(traits$shrub_weight_disp+1); colnames(traits)[which(colnames(traits)=="shrub_weight_disp")] <- "lg_shrub_weight_disp"
traits$tree_weight_disp <- log10(traits$tree_weight_disp+1); colnames(traits)[which(colnames(traits)=="tree_weight_disp")] <- "lg_tree_weight_disp"
traits$all_weight_disp <- log10(traits$all_weight_disp+1); colnames(traits)[which(colnames(traits)=="all_weight_disp")] <- "lg_all_weight_disp"

# 3) analytical graphs
# --->
# traits_melt <- melt(traits[,c(1:dim(traits)[2])], id=c("id_patch", "id_window", "window", "id_region"))
# traits.gg <- ggplot(traits_melt, aes(x=value)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(fill="red", alpha=.1) + ggtitle("distributions plant traits") + facet_wrap(~variable, scales = "free")
# traits.gg

# 4) factorial analysis
# --->
# fa.parallel(traits[c(5:dim(traits)[2])])
traits_fa <- fa(traits[c(5:dim(traits)[2])], nfactors=7, rotate="varimax", fm="ml", scores="regression")
# fa.graph(traits_fa, cut=.5, main = "factors forest structure")
factors_traits <- data.frame(traits_fa$scores); colnames(factors_traits) <- c(paste("FA_traits_", 1:dim(factors_traits)[2], sep=""))

#---|taxonomic|--------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# 1) assemble sub-dataframe
# --->
taxonomic <- tick_patches[,c(which(colnames(tick_patches)=="id_patch"):which(colnames(tick_patches)=="id_region"),
                             which(colnames(tick_patches)=="rich_tree_total_mean"):which(colnames(tick_patches)=="rich_tree_sgof_patch"),
                             which(colnames(tick_patches)=="alpha_spl_tall"):which(colnames(tick_patches)=="beta_total"),
                             which(colnames(tick_patches)=="tree_mixed_dom"):which(colnames(tick_patches)=="tree_evergreen_patch"))]

# 2) take out outlayers and apply transformations based on distributions of raw
# data
# --->
taxonomic$alpha_spl_tall <- log10(taxonomic$alpha_spl_tall+1); colnames(taxonomic)[which(colnames(taxonomic)=="alpha_spl_tall")] <- "lg_alpha_spl_tall"
taxonomic$alpha_spl_short <- log10(taxonomic$alpha_spl_short+1); colnames(taxonomic)[which(colnames(taxonomic)=="alpha_spl_short")] <- "lg_alpha_spl_short"
taxonomic$alpha_spl_total <- log10(taxonomic$alpha_spl_total+1); colnames(taxonomic)[which(colnames(taxonomic)=="alpha_spl_total")] <- "lg_alpha_spl_total"

# 3) analytical graphs
# --->
# taxonomic_melt <- melt(taxonomic[,c(1:which(colnames(taxonomic)=="beta_total"))], id=c("id_patch", "id_window", "window", "id_region"))
# taxonomic.gg <- ggplot(taxonomic_melt, aes(x=value)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(fill="red", alpha=.1) + ggtitle("distributions taxonomic") + facet_wrap(~variable, scales = "free")
# taxonomic.gg

# 4) factorial analysis
# --->
# fa.parallel(taxonomic[c(5:which(colnames(taxonomic)=="beta_total"))])
taxonomic_fa <- fa(taxonomic[c(5:which(colnames(taxonomic)=="beta_total"))], nfactors=4, rotate="varimax", fm="ml", scores="regression")
# fa.graph(taxonomic_fa, cut=.5, main = "factors forest structure")
factors_taxonomic <- data.frame(taxonomic_fa$scores); colnames(factors_taxonomic) <- c(paste("FA_taxonomic_", 1:dim(factors_taxonomic)[2], sep=""))


#---|soil|----------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# 1) assemble sub-dataframe
# --->
soil <- tick_patches[,c(which(colnames(tick_patches)=="id_patch"):which(colnames(tick_patches)=="id_region"),
                        which(colnames(tick_patches)=="ff_weight_lg_mean"):which(colnames(tick_patches)=="ms_n.p_org_lg_mean"),
                        which(colnames(tick_patches)=="trait11_mean"):which(colnames(tick_patches)=="trait47_mean"),
                        which(colnames(tick_patches)=="id_soil_wr")
)]
soil <- soil[, -c(which(colnames(soil)=="ff_weight_lg_mean"),
                  which(colnames(soil)=="ms_weight_mean"):which(colnames(soil)=="ms_rest_weight_mean"),
                  which(colnames(soil)=="vol_core_mean"),
                  which(colnames(soil)=="vol_root_mean"),
                  which(colnames(soil)=="ms_l1mm_vol_mean"),
                  which(colnames(soil)=="ff_c.cont1_lg_mean"),
                  which(colnames(soil)=="ms_c.cont1_lg_mean"),
                  which(colnames(soil)=="tot_c.cont1_lg_mean"),
                  which(colnames(soil)=="ff_n.cont1_lg_mean"),
                  which(colnames(soil)=="ms_n.cont1_lg_mean"),
                  which(colnames(soil)=="tot_n.cont1_lg_mean"),
                  which(colnames(soil)=="ff_c.conc_mean"):which(colnames(soil)=="ms_c.conc_lg_mean"),
                  which(colnames(soil)=="ff_n.conc_lg_mean"):which(colnames(soil)=="ms_n.conc_lg_mean"),
                  which(colnames(soil)=="ms_p.conc_tot_lg_mean"):which(colnames(soil)=="ms_p.conc_org_lg_mean"),
                  which(colnames(soil)=="ms_p.cont_perc.inorg_lg_mean"):which(colnames(soil)=="ms_p.cont_perc.org_lg_mean"),
                  which(colnames(soil)=="tot_c.cont1_lg_mean"):which(colnames(soil)=="tot_c.cont2_lg_mean"),
                  which(colnames(soil)=="tot_n.cont1_lg_mean"):which(colnames(soil)=="tot_n.cont2_lg_mean")
)]

# 2) take out outlayers and apply transformations based on distributions of raw
# data
# --->
soil$ms_p.cont_org_lg_mean[soil$ms_p.cont_org_lg_mean==min(soil$ms_p.cont_org_lg_mean, na.rm=T)] <- NA
soil$ms_p.cont_org_lg_mean[is.na(soil$ms_p.cont_org_lg_mean)] <- min(soil$ms_p.cont_org_lg_mean, na.rm=T)-0.01
soil$trait13_mean <- soil$trait13_mean*100

# 3) analytical graphs
# --->
# pairs.panels(soil[, c(5:dim(soil)[2])], pch='.', main="soil 1", hist.col="gray")
# pairs.panels(soil[, c(14:22)], pch='.', main="soil 2", hist.col="gray")
# cor.soil <- cor(soil[5:dim(soil)[2]], use="complete.obs"); plotcorr(cor.soil, col=rgb(colorfun((cor.soil+1)/2), maxColorValue=255), mar = c(0.1, 0.1, 0.1, 0.1), cex.lab=.7)
# soil_melt <- melt(soil[,c(1:22)], id=c("id_patch", "id_window", "window", "id_region"))
# soil.gg <- ggplot(soil_melt, aes(x=value)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(fill="red", alpha=.1) + ggtitle("distributions soil-data") + facet_wrap(~variable, scales = "free")
# soil.gg # looks kind of ok.
# scatter <- ggplot(dat=soil, aes(x=ff_bulkdensity_mean, y=ff_depth_mean)) + geom_point(aes(colour=id_region)) + geom_point(shape=1) + geom_smooth() + geom_smooth(method=lm, se=F, colour="red") + scale_colour_manual(values=colours) + theme(legend.position=c(1,1), legend.justification=c(1,1)); g <- scatter; g
# scatter2 <- ggplot(dat=main_patches, aes(x=height_stand_max, y=ba_dead_mean, colour=age3)) + geom_point() + geom_smooth() + geom_smooth(method=lm, se=F, colour="red") + scale_color_gradient(low="darkkhaki", high="darkgreen"); g <- scatter2; g

# 4) factorial analysis
# --->
# fa.parallel(soil[,5:(dim(soil)[2]-1)])
soil_fa <- fa(soil[,5:(dim(soil)[2]-1)], nfactors=6, rotate="varimax", fm="ml", scores="regression")
# fa.graph(soil_fa, cut=.5, main = "factors soil")
factors_soil <- data.frame(soil_fa$scores); colnames(factors_soil) <- c(paste("FA_soil_", 1:dim(factors_soil)[2], sep=""))


#---|landscape|-----------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# 1) assemble sub-dataframe
# --->
landscape <- tick_patches[,c(which(colnames(tick_patches)=="id_patch"):which(colnames(tick_patches)=="id_region"),
                             which(colnames(tick_patches)=="cv_altitude"),
                             which(colnames(tick_patches)=="surface_edge5_patch"):which(colnames(tick_patches)=="prop_edge20_patch"),
                             which(colnames(tick_patches)=="NND"):which(colnames(tick_patches)=="prop_hedg_500"),
                             which(colnames(tick_patches)=="edge_density_250"):which(colnames(tick_patches)=="prop_edge20_500"),
                             which(colnames(tick_patches)=="crown.cover"):which(colnames(tick_patches)=="patch.density")
)]

# 2) take out outlayers and apply transformations based on distributions of raw
# data
# --->
landscape[is.na(landscape)] <- 0

landscape$cv_altitude <- log10(landscape$cv_altitude+1); colnames(landscape)[which(colnames(landscape)=="cv_altitude")] <- "lg_cv_altitude"
landscape$surface_edge5_patch <- log10(landscape$surface_edge5_patch+1); colnames(landscape)[which(colnames(landscape)=="surface_edge5_patch")] <- "lg_surface_edge5_patch"
landscape$surface_edge10_patch <- log10(landscape$surface_edge10_patch+1); colnames(landscape)[which(colnames(landscape)=="surface_edge10_patch")] <- "lg_surface_edge10_patch"
landscape$surface_edge20_patch <- log10(landscape$surface_edge20_patch+1); colnames(landscape)[which(colnames(landscape)=="surface_edge20_patch")] <- "lg_surface_edge20_patch"
landscape$NND <- log10(landscape$NND+5); colnames(landscape)[which(colnames(landscape)=="NND")] <- "lg_NND"
landscape$PROX.5000 <- log10(landscape$PROX.5000+1); colnames(landscape)[which(colnames(landscape)=="PROX.5000")] <- "lg_PROX.5000"
landscape$lg_PROX.5000 <- log10(landscape$lg_PROX.5000+1); colnames(landscape)[which(colnames(landscape)=="lg_PROX.5000")] <- "lg_lg_PROX.5000"
landscape$prop_forest_50 <- ifelse(landscape$prop_forest_50==0, 0, 1); colnames(landscape)[which(colnames(landscape)=="prop_forest_50")] <- "bn_prop_forest_50"
landscape$prop_forest_100 <- log10(landscape$prop_forest_100+1); colnames(landscape)[which(colnames(landscape)=="prop_forest_100")] <- "lg_prop_forest_100"
landscape$prop_forest_250 <- log10(landscape$prop_forest_250+1); colnames(landscape)[which(colnames(landscape)=="prop_forest_250")] <- "lg_prop_forest_250"
landscape$prop_forest_500 <- log10(landscape$prop_forest_500+1); colnames(landscape)[which(colnames(landscape)=="prop_forest_500")] <- "lg_prop_forest_500"
landscape$prop_forest_1000 <- log10(landscape$prop_forest_1000+1); colnames(landscape)[which(colnames(landscape)=="prop_forest_1000")] <- "lg_prop_forest_1000"
landscape$prop_hedg_50 <- sqrt(landscape$prop_hedg_50); colnames(landscape)[which(colnames(landscape)=="prop_hedg_50")] <- "sqrt_prop_hedg_50"
landscape$prop_hedg_100 <- sqrt(landscape$prop_hedg_100); colnames(landscape)[which(colnames(landscape)=="prop_hedg_100")] <- "sqrt_prop_hedg_100"
landscape$prop_hedg_250 <- sqrt(landscape$prop_hedg_250); colnames(landscape)[which(colnames(landscape)=="prop_hedg_250")] <- "sqrt_prop_hedg_250"
landscape$prop_hedg_500 <- sqrt(landscape$prop_hedg_500); colnames(landscape)[which(colnames(landscape)=="prop_hedg_500")] <- "sqrt_prop_hedg_500"
landscape$edge_density_250 <- log10(landscape$edge_density_250+1); colnames(landscape)[which(colnames(landscape)=="edge_density_250")] <- "lg_edge_density_250"
landscape$prop_edge5_250 <- asin(sqrt(landscape$prop_edge5_250)); colnames(landscape)[which(colnames(landscape)=="prop_edge5_250")] <- "asin_prop_edge5_250"
landscape$prop_edge10_250 <- asin(sqrt(landscape$prop_edge10_250)); colnames(landscape)[which(colnames(landscape)=="prop_edge10_250")] <- "asin_prop_edge10_250"
landscape$prop_edge20_250 <- asin(sqrt(landscape$prop_edge20_250)); colnames(landscape)[which(colnames(landscape)=="prop_edge20_250")] <- "asin_prop_edge20_250"
landscape$edge_density_500 <- log10(landscape$edge_density_500+1); colnames(landscape)[which(colnames(landscape)=="edge_density_500")] <- "lg_edge_density_500"
landscape$prop_edge5_500 <- asin(sqrt(landscape$prop_edge5_500)); colnames(landscape)[which(colnames(landscape)=="prop_edge5_500")] <- "asin_prop_edge5_500"
landscape$prop_edge10_500 <- asin(sqrt(landscape$prop_edge10_500)); colnames(landscape)[which(colnames(landscape)=="prop_edge10_500")] <- "asin_prop_edge10_500"
landscape$prop_edge20_500 <- asin(sqrt(landscape$prop_edge20_500)); colnames(landscape)[which(colnames(landscape)=="prop_edge20_500")] <- "asin_prop_edge20_500"
landscape$forest.cover <- log10(landscape$forest.cover+1); colnames(landscape)[which(colnames(landscape)=="forest.cover")] <- "lg_forest.cover"
landscape$prop_past_1000 <- landscape$prop_past_1000*100
landscape$prop_cult_50 <- landscape$prop_cult_50*100
landscape$prop_cult_1000 <- landscape$prop_cult_1000*100

# 3) analytical graphs
# --->
# pairs.panels(landscape[, c(5:14)], pch='.', main="landscape", hist.col="gray")
# pairs.panels(landscape[, c(15:24)], pch='.', main="landscape", hist.col="gray")
# pairs.panels(landscape[, c(25:34)], pch='.', main="landscape", hist.col="gray")
# pairs.panels(landscape[, c(35:44)], pch='.', main="landscape", hist.col="gray")
# pairs.panels(landscape[, c(45:48)], pch='.', main="landscape", hist.col="gray")
# cor.landscape <- cor(landscape[5:dim(landscape)[2]]); plotcorr(cor.landscape, col=rgb(colorfun((cor.landscape+1)/2), maxColorValue=255), mar = c(0.1, 0.1, 0.1, 0.1), cex.lab=.7)
# landscape_melt <- melt(landscape[,1:dim(landscape)[2]], id=c("id_patch", "id_window", "window", "id_region"))
# landscape.gg <- ggplot(landscape_melt, aes(x=value)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(fill="red", alpha=.1) + ggtitle("distributions landscape-data") + facet_wrap(~variable, scales = "free")
# landscape.gg # transformations necessary -> done

# 4) factorial analysis
# --->
# fa.parallel(landscape[,5:dim(landscape)[2]])
landscape_fa <- fa(landscape[,5:dim(landscape)[2]], nfactors=5, rotate="varimax", fm="ml", scores="regression")
# fa.graph(landscape_fa, cut=.6, main = "factors landscape")
factors_landscape <- data.frame(landscape_fa$scores); colnames(factors_landscape) <- c(paste("FA_landscape_", 1:dim(factors_landscape)[2], sep=""))


#---|microclimate|--------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# 1) assemble sub-dataframe
# --->
microclimate <- tick_patches[,c(which(colnames(tick_patches)=="id_patch"):which(colnames(tick_patches)=="id_region"),
                                which(colnames(tick_patches)=="doy"):which(colnames(tick_patches)=="satdef")
)]

# 2) take out outlayers and apply transformations based on distributions of raw
# data
# --->
microclimate <- microclimate[, -which(colnames(microclimate)=="doy")]
colnames(microclimate)[which(colnames(microclimate)=="start_time")] <- "sample_time"
microclimate$temp_soil[microclimate$temp_soil==min(microclimate$temp_soil, na.rm=T)] <- NA
microclimate$temp_soil[is.na(microclimate$temp_soil)] <- 15.5
microclimate$temp_soil[microclimate$temp_soil==min(microclimate$temp_soil, na.rm=T)] <- NA
microclimate$temp_soil[is.na(microclimate$temp_soil)] <- 15.5
microclimate$satdef <- log10(microclimate$satdef+1); colnames(microclimate)[which(colnames(microclimate)=="satdef")] <- "lg_satdef"

# 3) analytical graphs
# --->
# pairs.panels(microclimate[, c(5:(dim(microclimate)[2]))], pch='.', main="microclimate 1", hist.col="gray")
# cor.microclimate <- cor(microclimate[, c(5:(dim(microclimate)[2]-2))], use="complete.obs"); plotcorr(cor.microclimate, col=rgb(colorfun((cor.microclimate+1)/2), maxColorValue=255), mar = c(0.1, 0.1, 0.1, 0.1))
# microclimate_melt <- melt(microclimate, id=c("id_patch", "id_window", "window", "id_region"))
# microclimate.gg <- ggplot(microclimate_melt, aes(x=value)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(fill="red", alpha=.1) + ggtitle("distributions microclimate-data") + facet_wrap(~variable, scales = "free")
# microclimate.gg

# 4) factorial analysis
# --->
# fa.parallel(microclimate[, c(5:dim(microclimate)[2])])
microclimate_fa <- fa(microclimate[, c(5:dim(microclimate)[2])], nfactors=2, rotate="varimax", fm="ml", scores="regression")
# fa.graph(microclimate_fa, cut=.5, main = "factors microclimate")
factors_microclimate <- data.frame(microclimate_fa$scores); colnames(factors_microclimate) <- c(paste("FA_microclimate_", 1:dim(factors_microclimate)[2], sep=""))


#---|macroclimate|--------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# 1) assemble sub-dataframe
# --->
macroclimate <- tick_patches[,c(which(colnames(tick_patches)=="id_patch"):which(colnames(tick_patches)=="id_region"),
                                # which(colnames(tick_patches)=="bio1"):which(colnames(tick_patches)=="bio7"),
                                # which(colnames(tick_patches)=="bio11"):which(colnames(tick_patches)=="bio12"),
                                # which(colnames(tick_patches)=="bio15"),
                                which(colnames(tick_patches)=="gdd"):which(colnames(tick_patches)=="temp_max_year"),
                                which(colnames(tick_patches)=="n_snow_year"):which(colnames(tick_patches)=="n_warm_year")
)]

# 2) take out outlayers and apply transformations based on distributions of raw
# data
# --->
# not happening here as these data have been prepared already and are used "as
# is".

# 3) analytical graphs
# --->
# pairs.panels(macroclimate[, c(5:dim(macroclimate)[2])], pch='.', main="macroclim", hist.col="gray")
# cor.macroclimate <- cor(macroclimate[5:dim(macroclimate)[2]], use="complete.obs"); plotcorr(cor.macroclimate, col=rgb(colorfun((cor.macroclimate+1)/2), maxColorValue=255), mar = c(0.1, 0.1, 0.1, 0.1), cex.lab=.7)
# climate_melt <- melt(macroclimate[,1:dim(macroclimate)[2]], id=c("id_patch", "id_window", "window", "id_region"))
# macroclimate.gg <- ggplot(climate_melt, aes(x=value)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(fill="red", alpha=.1) + ggtitle("distributions macroclimate-data") + facet_wrap(~variable, scales = "free")
# macroclimate.gg # well...

# 4) factorial analysis
# --->
# fa.parallel(macroclimate[,5:dim(macroclimate)[2]])
macroclimate_fa <- fa(macroclimate[,5:dim(macroclimate)[2]], nfactors = 3, rotate="varimax", fm="ml", scores="regression")
# fa.graph(macroclimate_fa, cut=.5, main = "factors macroclimate")
factors_macroclim <- data.frame(macroclimate_fa$scores); colnames(factors_macroclim) <- c(paste("FA_macroclimate_", 1:dim(factors_macroclim)[2], sep=""))


#---|meta-data|-----------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# 1) assemble sub-dataframe
# --->
meta <- tick_patches[, c(which(colnames(tick_patches)=="id_patch"): which(colnames(tick_patches)=="id_region"),
                         which(colnames(tick_patches)=="X_long"): which(colnames(tick_patches)=="Y_lat2"),
                         which(colnames(tick_patches)=="surface"): which(colnames(tick_patches)=="age4"))]

# 2) take out outlayers and apply transformations based on distributions of raw
# data
# --->
meta$surface <- log10(meta$surface+1); colnames(meta)[which(colnames(meta)=="surface")] <- "lg.surface"
meta$perimeter <- log10(meta$perimeter+1); colnames(meta)[which(colnames(meta)=="perimeter")] <- "lg.perimeter"
meta$IndexAge_orig <- log10(meta$IndexAge_orig+1); colnames(meta)[which(colnames(meta)=="IndexAge_orig")] <- "lg.IndexAge_orig"
meta$IndexAge_buff <- log10(meta$IndexAge_buff+1); colnames(meta)[which(colnames(meta)=="IndexAge_buff")] <- "lg.IndexAge_buff"

# 3) analytical graphs
# --->
# pairs.panels(meta[, c(5:dim(meta)[2])], pch=".", main="meta", hist.col="gray")

# 4) factorial analysis
# --->
# -


################################################################################
# combine all (transformed) variables for model-building

# combined factors
factors <- data.frame(factors_microclimate, factors_macroclim, factors_structure, factors_abundances, factors_traits, factors_taxonomic, factors_soil, factors_landscape)

# overall table without ticks being scaled
all_ticks <- data.frame(meta, structure[, c(5:dim(structure)[2])], abundances[, c(5:dim(abundances)[2])], traits[, c(5:dim(traits)[2])], taxonomic[, c(5:dim(taxonomic)[2])], soil[, c(5:dim(soil)[2])], landscape[, c(5:dim(landscape)[2])], microclimate[, c(5:(dim(microclimate)[2]))], macroclimate[, c(5:dim(macroclimate)[2])])
all_ticks_scale <- scale(all_ticks[, c(6:which(colnames(all_ticks)=="trait47_mean"),
                                       which(colnames(all_ticks)=="lg_cv_altitude"):dim(all_ticks)[2])],
                         center = F, scale = F)
all <- cbind(ticks, n_points = tick_patches$n_points_wp2, all_ticks_scale, id_soil_wr = all_ticks[, c(which(colnames(all_ticks)=="id_soil_wr"))])

# restore some variables after scaling (age, surface)
all$surface2 <- as.factor(ifelse(all$surface2==0, 0, 1))
all$surface3 <- as.factor(ifelse(all$surface3==0, 0, 1))
all$age2 <- as.factor(ifelse(all$age2==0, 0, 1))
all$age3 <- as.factor(ifelse(all$age3==0, 0, 1))
all$age4 <- as.factor(ifelse(all$age4==0, 0, 1))
all$bn_ba_evergr_bn_mean <- as.factor(ifelse(all$bn_ba_evergr_bn_mean==0, 0, 1))
all$bn_ba_dead_bn_mean <- as.factor(ifelse(all$bn_ba_dead_bn_mean==0, 0, 1))
all$ms_l1mm_mat_bn_mean <- as.factor(ifelse(all$ms_l1mm_mat_bn_mean==0, 0, 1))

all <- cbind(all, factors)


################################################################################
# Meta-Data

all_num <- all[sapply(all, is.numeric)]

probs <- c(.025, .5, .975, 0, 1)

quantiles <- round(data.frame(matrix(unlist(lapply(all_num,  quantile, probs = probs, na.rm = T, name = FALSE)), ncol = 5, byrow = T)), 2)
means <- round(data.frame(matrix(unlist(lapply(all_num, mean, na.rm = T)), ncol = 1, byrow = T)), 2)
sds <- round(data.frame(matrix(unlist(lapply(all_num, sd, na.rm = T)), ncol = 1, byrow = T)), 2)
quantiles <- cbind(orig = colnames(all_num), quantiles, means, sds)
quantiles <- cbind(quantiles, var2=paste("I(", quantiles$orig, "^2)", sep=""))
quantiles <- quantiles[, c(1, 9, 2:8)]
colnames(quantiles) <- c("var1", "var2", "2.5%", "50%", "97.5%", "min", "max", "mean", "sd")
# does of course not make sense to use for all variables, but this is the
# fastest and most consistent way.


################################################################################
# Look-up table for driver groups

lut_group <- data.frame(group = c("ticks", "microclimate", "macroclimate", "landscape", "structural properties", "functional properties", "taxonomic properties", "soil", "method"),
                        metagroup = c("ticks", "habitat", "macroclimate", "landscape", "habitat", "habitat", "habitat", "habitat", "method"))

all_num <- sapply(all, is.numeric)
all_num <- all[, all_num]


################################################################################
# Names for plotting (only relevant, i.e. significant in the models, variables).
# See Tab. A4.1 in the supporting information for the meaning of the different
# name-components

names <- c(id_patch = "id_patch",
           id_window = "id_window",
           window = "window",
           id_region = "id_region",
           l_mean = "l_mean",
           lg_l_mean = "abundance \nlarvae [log]",
           l_pa = "l_pa",
           n_sum = "n_sum",
           n_mean = "n_mean",
           lg_n_mean = "abundance \nnymphs [log]",
           n_pa = "n_pa",
           n_lab_sum = "n_lab_sum",
           n_lab_mean = "n_lab_mean",
           n_ip_gm = "n_ip_gm",
           n_ia_gm = "n_ia_gm",
           lg_n_ia_gm = "lg_n_ia_gm",
           n_inf_pa = "n_inf_pa",
           f_sum = "f_sum",
           f_mean = "f_mean",
           lg_f_mean = "lg_f_mean",
           f_pa = "f_pa",
           f_lab_sum = "f_lab_sum",
           f_lab_mean = "f_lab_mean",
           f_ip_gm = "f_ip_gm",
           f_ia_gm = "f_ia_gm",
           lg_f_ia_gm = "lg_f_ia_gm",
           f_inf_pa = "f_inf_pa",
           m_sum = "m_sum",
           m_mean = "m_mean",
           lg_m_mean = "lg_m_mean",
           m_pa = "m_pa",
           m_lab_sum = "m_lab_sum",
           m_lab_mean = "m_lab_mean",
           m_ip_gm = "m_ip_gm",
           m_ia_gm = "m_ia_gm",
           lg_m_ia_gm = "lg_m_ia_gm",
           m_inf_pa = "m_inf_pa",
           a_sum = "a_sum",
           a_mean = "a_mean",
           lg_a_mean = "abundance \nadults [log]",
           a_pa = "a_pa",
           a_lab_sum = "a_lab_sum",
           a_lab_mean = "a_lab_mean",
           a_ip_gm = "a_ip_gm",
           a_ia_gm = "a_ia_gm",
           lg_a_ia_gm = "lg_a_ia_gm",
           a_inf_pa = "a_inf_pa",
           doy = "doy",
           lt_n_ip_gm = "lt_n_ip_gm",
           lt_f_ip_gm = "lt_f_ip_gm",
           lt_m_ip_gm = "lt_m_ip_gm",
           lt_a_ip_gm = "lt_a_ip_gm",
           n_points = "n_points",
           X_long2 = "X_long2",
           Y_lat = "Y_lat",
           Y_lat2 = "Y_lat2",
           lg.surface = "patch size [log]",
           surface2 = "surface2",
           surface3 = "surface3",
           lg.perimeter = "lg.perimeter",
           perim_surf = "perim_surf",
           lg.IndexAge_orig = "patch age \nindex [log]",
           lg.IndexAge_buff = "lg.IndexAge_buff",
           age2 = "age2",
           age3 = "age3",
           age4 = "age4",
           height_stand_mean = "height_stand_mean",
           height_stand_max = "height_stand_max",
           treelayers_min = "treelayers_min",
           treelayers_max = "treelayers_max",
           md_total_mean = "average distance\n~ all trees",
           md_large_mean = "average distance\n~ large trees",
           dens_total_lg_mean = "dens_total_lg_mean",
           dens_large_lg_mean = "dens. trees w/\n d130 > 30 cm",
           slenderness_total_mean = "slenderness_total_mean",
           slenderness_small_mean = "slenderness_small_mean",
           slenderness_large_mean = "slenderness_large_mean",
           ba_bitterlich_mean = "basal area \n(bitterlich)",
           ba_total_lg_mean = "ba_total_lg_mean",
           ba_total_lg_min = "ba_total_lg_min",
           ba_total_lg_max = "ba_total_lg_max",
           ba_large_lg_mean = "ba_large_lg_mean",
           ba_small_lg_mean = "ba_small_lg_mean",
           ba_decid_lg_mean = "ba_decid_lg_mean",
           bn_ba_evergr_bn_mean = "bn_ba_evergr_bn_mean",
           bn_ba_dead_bn_mean = "bn_ba_dead_bn_mean",
           ba_laying_total_lg_mean = "deadwood\n~ ground",
           ba_laying_total_lg_max = "max(deadwood)\n~ ground",
           ba_laying_large_bn_mean = "deadwood > 15cm\n~ ground",
           species_mingling_mean = "species_mingling_mean",
           diam_diff_mean = "diameter\n differentiation",
           lg_diam_cv_mean = "lg_diam_cv_mean",
           lg_diam_iqr = "lg_diam_iqr",
           diam_med = "diam_med",
           lg_diam_biggest = "lg_diam_biggest",
           dw_diam_mean_mean = "dw_diam_mean_mean",
           vol_total_lg_mean = "vol_total_lg_mean",
           vol_total_lg_min = "vol_total_lg_min",
           vol_total_lg_max = "vol_total_lg_max",
           biomass_total_lg_mean = "biomass_total_lg_mean",
           diss_tree1_tree2 = "species dissimilarity\n~ tree1 vs. tree2 layer",
           diss_tree1_spl = "diss_tree1_spl",
           diss_tree2_spl = "species dissimilarity\n~ tree2 vs. saplings",
           herb_canopy_height = "herb_canopy_height",
           herb_canopy_height_cv = "CV canopy height\n~ herb layer",
           tree_canopy_height = "tree_canopy_height",
           lg_tree_canopy_height_cv = "lg_tree_canopy_height_cv",
           herb_total_abund = "total herb abund.",
           herb_average_abund = "mean abund.\n~ herb layer",
           herb_small_disp_abund = "herb_small_disp_abund",
           herb_medium_disp_abund = "herb_medium_disp_abund",
           herb_medium_large_disp_abund = "herb_medium_large_disp_abund",
           herb_large_disp_abund = "abund. disp. > 2 g\n~ herb layer",
           herb_w_nuts_abund = "herb_w_nuts_abund",
           herb_w_nuts_e_abund = "herb_w_nuts_e_abund",
           herb_w_berries_abund = "herb_w_berries_abund",
           herb_w_legume_abund = "herb_w_legume_abund",
           herb_evergreen_abund = "abund. evergreens\n~ herb layer",
           herb_pr_padus_abund = "herb_pr_padus_abund",
           shrub_total_abund = "shrub_total_abund",
           shrub_average_abund = "shrub_average_abund",
           shrub_small_disp_abund = "shrub_small_disp_abund",
           shrub_medium_disp_abund = "abund. disp. > 0.1 g\n~ shrub layer",
           shrub_medium_large_disp_abund = "shrub_medium_large_disp_abund",
           shrub_large_disp_abund = "shrub_large_disp_abund",
           shrub_w_nuts_abund = "abund. w/ nuts\n~ shrub layer",
           shrub_w_nuts_e_abund = "abund. w/ nuts(+e)\n~ shrub layer",
           shrub_w_berries_abund = "shrub_w_berries_abund",
           shrub_evergreen_abund = "shrub_evergreen_abund",
           shrub_pr_padus_abund = "shrub_pr_padus_abund",
           tree_total_abund = "total abund.\n~ tree layer",
           tree_average_abund = "tree_average_abund",
           tree_small_disp_abund = "tree_small_disp_abund",
           tree_medium_disp_abund = "tree_medium_disp_abund",
           tree_medium_large_disp_abund = "tree_medium_large_disp_abund",
           tree_large_disp_abund = "abund. disp. > 2 g\n~ tree layer",
           tree_w_nuts_abund = "tree_w_nuts_abund",
           tree_w_nuts_e_abund = "abund. w/ nuts(+e)\n~ tree layer",
           tree_w_berries_abund = "abund. w/ berries\n~ tree layer",
           tree_evergreen_abund = "tree_evergreen_abund",
           tree_pr_padus_abund = "tree_pr_padus_abund",
           all_total_abund = "all_total_abund",
           shrub_tree_total_abund = "total abund.\n~ shrub+tree layer",
           ratio_herb_tree_abund = "ratio_herb_tree_abund",
           lg_ratio_herb_shrub_abund = "lg_ratio_herb_shrub_abund",
           lg_ratio_shrub_tree_abund = "lg_ratio_shrub_tree_abund",
           all_small_disp_abund = "all_small_disp_abund",
           all_medium_disp_abund = "all_medium_disp_abund",
           all_medium_large_disp_abund = "all_medium_large_disp_abund",
           all_large_disp_abund = "abund. disp. > 2 g\n~ all layers",
           all_w_nuts_abund = "abund. w/ nuts\n~ all layers",
           all_w_nuts_e_abund = "all_w_nuts_e_abund",
           all_w_berries_abund = "all_w_berries_abund",
           all_evergreen_abund = "all_evergreen_abund",
           all_pr_padus_abund = "abund. Pr. padus\n~ all layers",
           lg_herb_weight_disp = "lg_herb_weight_disp",
           herb_small_disp_alpha = "herb_small_disp_alpha",
           herb_medium_disp_alpha = "herb_medium_disp_alpha",
           herb_medium_large_disp_alpha = "herb_medium_large_disp_alpha",
           herb_large_disp_alpha = "herb_large_disp_alpha",
           herb_sla = "herb_sla",
           lg_herb_chamaephyte = "% chamaephyte CWM\n~ herb layer",
           lg_herb_geophyte = "lg_herb_geophyte",
           herb_hemicryptophyte = "herb_hemicryptophyte",
           herb_phanerophyte = "herb_phanerophyte",
           lg_herb_therophyte = "lg_herb_therophyte",
           herb_perennial = "herb_perennial",
           bn_herb_fern = "bn_herb_fern",
           herb_forb = "herb_forb",
           lg_herb_grass = "% grasses [log]\n~ herb layer (CWM)",
           herb_asc_to_pros = "% w/ asc./pro. hab. \n~ herb layer (CWM)",
           herb_erect = "herb_erect",
           bn_herb_prostrate = "bn_herb_prostrate",
           herb_regular = "% w/ reg. leaf-dist.\n~ herb layer (CWM)",
           herb_rosette = "herb_rosette",
           herb_semi_rosette = "herb_semi_rosette",
           herb_wo_branching = "herb_wo_branching",
           herb_w_branching = "herb_w_branching",
           lg_shrub_weight_disp = "weight disp.\n~ shrub layer",
           shrub_small_disp_alpha = "shrub_small_disp_alpha",
           shrub_medium_disp_alpha = "shrub_medium_disp_alpha",
           shrub_medium_large_disp_alpha = "shrub_medium_large_disp_alpha",
           shrub_large_disp_alpha = "richness disp. > 2 g\n~ shrub layer",
           shrub_sla = "specific leaf area\n~ shrub layer",
           lg_tree_weight_disp = "lg_tree_weight_disp",
           tree_small_disp_alpha = "richness disp. < 0.1 g\n~ tree layer",
           tree_medium_disp_alpha = "tree_medium_disp_alpha",
           tree_medium_large_disp_alpha = "tree_medium_large_disp_alpha",
           tree_large_disp_alpha = "tree_large_disp_alpha",
           tree_sla = "specific leaf area\n~ tree layer",
           lg_all_weight_disp = "lg_all_weight_disp",
           all_small_disp_alpha = "all_small_disp_alpha",
           all_medium_disp_alpha = "all_medium_disp_alpha",
           all_medium_large_disp_alpha = "all_medium_large_disp_alpha",
           all_large_disp_alpha = "richness disp. > 2 g\n~ all layers",
           all_sla = "all_sla",
           rich_tree_total_mean = "rich_tree_total_mean",
           rich_tree_small_mean = "rich_tree_small_mean",
           rich_tree_large_mean = "richness\ntrees d130 > 30 cm",
           rich_tree_short_mean = "rich_tree_short_mean",
           rich_tree_tall_mean = "rich_tree_tall_mean",
           rich_tree_sgof_mean = "rich_tree_sgof_mean",
           rich_tree_total_patch = "rich_tree_total_patch",
           rich_tree_small_patch = "rich_tree_small_patch",
           rich_tree_large_patch = "rich_tree_large_patch",
           rich_tree_short_patch = "rich_tree_short_patch",
           rich_tree_tall_patch = "rich_tree_tall_patch",
           lg_alpha_spl_tall = "richness [log]\nsaplings h > 1.3 m",
           lg_alpha_spl_short = "lg_alpha_spl_short",
           lg_alpha_spl_total = "richness [log]\nall saplings",
           gamma_herb = "gamma_herb",
           gamma_shrub = "gamma_shrub",
           gamma_tree = "gamma_tree",
           gamma_total = "gamma_total",
           mean_alpha_herb = "mean_alpha_herb",
           mean_alpha_shrub = "mean_alpha_shrub",
           mean_alpha_tree = "richness\n~ tree layer",
           mean_alpha_total = "mean_alpha_total",
           beta_herb = "beta_herb",
           beta_shrub = "beta_shrub",
           beta_tree = "beta_tree",
           beta_total = "beta_total",
           tree_mixed_dom = "% patch codominated\n~ tree layer",
           tree_fagus_dom = "tree_fagus_dom",
           tree_fagus_part_patch = "tree_fagus_part_patch",
           tree_quercus_dom_patch = "tree_quercus_dom_patch",
           tree_quercus_part_patch = "% patch with Quercus",
           tree_populus_dom_patch = "tree_populus_dom_patch",
           tree_populus_part_patch = "tree_populus_part_patch",
           tree_betula_dom_patch = "tree_betula_dom_patch",
           tree_betula_part_patch = "tree_betula_part_patch",
           tree_alnus_dom_patch = "tree_alnus_dom_patch",
           tree_alnus_part_patch = "tree_alnus_part_patch",
           tree_prunus_dom_patch = "tree_prunus_dom_patch",
           tree_prunus_part_patch = "% patch with Prunus",
           tree_fraxinus_dom_patch = "tree_fraxinus_dom_patch",
           tree_fraxinus_part_patch = "tree_fraxinus_part_patch",
           tree_acer_dom_patch = "tree_acer_dom_patch",
           tree_acer_part_patch = "tree_acer_part_patch",
           tree_tilia_dom_patch = "tree_tilia_dom_patch",
           tree_tilia_part_patch = "tree_tilia_part_patch",
           tree_ulmus_dom_patch = "tree_ulmus_dom_patch",
           tree_ulmus_part_patch = "tree_ulmus_part_patch",
           tree_carpinus_dom_patch = "tree_carpinus_dom_patch",
           tree_carpinus_part_patch = "tree_carpinus_part_patch",
           tree_larix_dom_patch = "tree_larix_dom_patch",
           tree_larix_part_patch = "tree_larix_part_patch",
           tree_picea_dom_patch = "tree_picea_dom_patch",
           tree_picea_part_patch = "tree_picea_part_patch",
           tree_pinus_dom_patch = "tree_pinus_dom_patch",
           tree_pinus_part_patch = "tree_pinus_part_patch",
           tree_ilex_part_patch = "tree_ilex_part_patch",
           tree_hedera_patch = "tree_hedera_patch",
           tree_evergreen_patch = "tree_evergreen_patch",
           ff_depth_lg_mean = "ff_depth_lg_mean",
           ff_bulkdensity_lg_mean = "bulkdensity [log]\n~ forest floor",
           ms_l1mm_mat_bn_mean = "ms_l1mm_mat_bn_mean",
           ms_bulkdensity_mean = "ms_bulkdensity_mean",
           ff_c.cont2_lg_mean = "ff_c.cont2_lg_mean",
           ms_c.cont2_lg_mean = "ms_c.cont2_lg_mean",
           ff_n.cont2_lg_mean = "Nitrogen cont. [log]\n~ forest floor",
           ms_n.cont2_lg_mean = "ms_n.cont2_lg_mean",
           ff_c.n_lg_mean = "ff_c.n_lg_mean",
           ms_c.n_lg_mean = "C/N [log]\n~ mineral soil",
           ff_ph_mean = "ph-value\n~ forest floor",
           ms_ph_mean = "ph-value\n~ mineral soil",
           ms_p.cont_tot_lg_mean = "ms_p.cont_tot_lg_mean",
           ms_p.cont_inorg_lg_mean = "ms_p.cont_inorg_lg_mean",
           ms_p.cont_org_lg_mean = "ms_p.cont_org_lg_mean",
           ms_n.p_tot_lg_mean = "N/P(tot) [log]\n~ mineral soil",
           ms_n.p_inorg_lg_mean = "ms_n.p_inorg_lg_mean",
           ms_n.p_org_lg_mean = "ms_n.p_org_lg_mean",
           trait11_mean = "specific leaf area\n~ ground",
           trait13_mean = "% leaf-litter carbon\n~ ground",
           trait14_mean = "‰ leaf nitrogen\n~ ground",
           trait15_mean = "‰ leaf phosphorous\n~ ground",
           trait47_mean = "leaf dry matter content\n~ ground",
           lg_cv_altitude = "lg_cv_altitude",
           lg_surface_edge5_patch = "lg_surface_edge5_patch",
           lg_surface_edge10_patch = "lg_surface_edge10_patch",
           lg_surface_edge20_patch = "lg_surface_edge20_patch",
           prop_edge5_patch = "% edge habitat\n~ 5 m into patch",
           prop_edge10_patch = "prop_edge10_patch",
           prop_edge20_patch = "prop_edge20_patch",
           lg_NND = "lg_NND",
           lg_lg_PROX.5000 = "proximity [log]\n~ 5000 m buffer",
           bn_prop_forest_50 = "bn_prop_forest_50",
           lg_prop_forest_100 = "lg_prop_forest_100",
           lg_prop_forest_250 = "% forest [log]\n~ 250 m buffer",
           lg_prop_forest_500 = "lg_prop_forest_500",
           lg_prop_forest_1000 = "lg_prop_forest_1000",
           prop_forest_5000 = "prop_forest_5000",
           prop_cult_50 = "% agriculture\n~ 50 m buffer",
           prop_cult_100 = "% agriculture\n ~ 100 m buffer",
           prop_cult_250 = "prop_cult_250",
           prop_cult_500 = "prop_cult_500",
           prop_cult_1000 = "% agriculture\n~ 1000 m buffer",
           prop_cult_5000 = "prop_cult_5000",
           prop_past_50 = "prop_past_50",
           prop_past_100 = "prop_past_100",
           prop_past_250 = "prop_past_250",
           prop_past_500 = "prop_past_500",
           prop_past_1000 = "% pasture\n~ 1000 m buffer",
           prop_past_5000 = "prop_past_5000",
           sqrt_prop_hedg_50 = "sqrt_prop_hedg_50",
           sqrt_prop_hedg_100 = "sqrt_prop_hedg_100",
           sqrt_prop_hedg_250 = "sqrt_prop_hedg_250",
           sqrt_prop_hedg_500 = "sqrt_prop_hedg_500",
           lg_edge_density_250 = "lg_edge_density_250",
           asin_prop_edge5_250 = "asin_prop_edge5_250",
           asin_prop_edge10_250 = "asin_prop_edge10_250",
           asin_prop_edge20_250 = "asin_prop_edge20_250",
           lg_edge_density_500 = "lg_edge_density_500",
           asin_prop_edge5_500 = "asin_prop_edge5_500",
           asin_prop_edge10_500 = "asin_prop_edge10_500",
           asin_prop_edge20_500 = "asin_prop_edge20_500",
           crown.cover = "crown.cover",
           lg_forest.cover = "% forest [log]\n~ 500 m buffer",
           edge.density.forest = "edge density\n~ 500 m buffer",
           n.patches = "nr forest patches\n~ 500 m buffer",
           patch.density = "patch.density",
           sample_time = "time of\n sampling",
           flag_condition = "flag wetness",
           temp_soil = "soil temperature",
           rH_5 = "relative humidity\n~ 5 cm height",
           rH_130 = "relative humidity\n~ 130 cm height",
           rH = "relative humidity",
           temp_air = "temp_air",
           lg_satdef = "lg_satdef",
           gdd = "gdd",
           cdd = "chilling degree days\n~ since 01.01.2013",
           precip_pre = "precip_pre",
           precip_year = "precip_year",
           temp_mean_pre = "temp_mean_pre",
           temp_mean_year = "temp_mean_year",
           temp_min_pre = "temp_min_pre",
           temp_min_year = "temp_min_year",
           temp_max_pre = "temp_max_pre",
           temp_max_year = "temp_max_year",
           n_snow_year = "n_snow_year",
           n_dew_pre = "n_dew_pre",
           n_dew_year = "n_dew_year",
           n_cold_pre = "n_cold_pre",
           n_cold_year = "days < 8 °C\n~ since 01.01.2013",
           n_med_pre = "days w/ temp > 8 °C\n~ previous 30 days",
           n_med_year = "days w/ temp > 8 °C\n~ since 01.01.2013",
           n_warm_pre = "n_warm_pre",
           n_warm_year = "n_warm_year",
           id_soil_wr = "id_soil_wr",
           FA_microclimate_1 = "FA_microclimate_1",
           FA_microclimate_2 = "FA_microclimate_2",
           FA_macroclimate_1 = "FA_macroclimate_1",
           FA_macroclimate_2 = "FA_macroclimate_2",
           FA_macroclimate_3 = "FA_macroclimate_3",
           FA_structure_1 = "FA_structure_1",
           FA_structure_2 = "FA_structure_2",
           FA_structure_3 = "FA_structure_3",
           FA_structure_4 = "deadwood (FA)",
           FA_structure_5 = "tree slenderness (FA)",
           FA_structure_6 = "FA_structure_6",
           FA_structure_7 = "FA_structure_7",
           FA_structure_8 = "FA_structure_8",
           FA_abundances_1 = "abund. disp. > 0.1 g\n~ herb layer (FA)",
           FA_abundances_2 = "FA_abundances_2",
           FA_abundances_3 = "FA_abundances_3",
           FA_abundances_4 = "FA_abundances_4",
           FA_abundances_5 = "abund. Pr. padus\n~all layers (FA)",
           FA_abundances_6 = "FA_abundances_6",
           FA_abundances_7 = "FA_abundances_7",
           FA_abundances_8 = "abund. w/ nuts\n~ shrub layer (FA)",
           FA_traits_1 = "FA_traits_1",
           FA_traits_2 = "% w/ branching\n~ herb layer (FA)",
           FA_traits_3 = "FA_traits_3",
           FA_traits_4 = "FA_traits_4",
           FA_traits_5 = "FA_traits_5",
           FA_traits_6 = "FA_traits_6",
           FA_traits_7 = "FA_traits_7",
           FA_taxonomic_1 = "FA_taxonomic_1",
           FA_taxonomic_2 = "beta-diversity\n~ all layers",
           FA_taxonomic_3 = "FA_taxonomic_3",
           FA_taxonomic_4 = "FA_taxonomic_4",
           FA_soil_1 = "FA_soil_1",
           FA_soil_2 = "FA_soil_2",
           FA_soil_3 = "inorganic Phosphor\n~ mineral soil (FA)",
           FA_soil_4 = "FA_soil_4",
           FA_soil_5 = "FA_soil_5",
           FA_soil_6 = "FA_soil_6",
           FA_landscape_1 = "FA_landscape_1",
           FA_landscape_2 = "FA_landscape_2",
           FA_landscape_3 = "FA_landscape_3",
           FA_landscape_4 = "patch density\n~ 500 m buffer (FA)",
           FA_landscape_5 = "% edge habitat\n~ patch (FA)")
