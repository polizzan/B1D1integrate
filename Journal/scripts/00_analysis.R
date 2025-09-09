#####################
# CLEAR ENVIRONMENT #
#####################

rm(list = ls())

##############
# USER INPUT #
##############

## start / end ages of reproductive period
fert.ages <- c(12, 55)

## start / end ages for plotting
plot.ages <- c(15, 45)

## main country of analysis
select.country <- "Sweden"

## other countries to be considered
countries <- c("Denmark", "England and Wales", "Finland", "USA")

#############################
# INSTALL AND LOAD PACKAGES #
#############################

## packages to be installed from cran
from.cran <- c("cowplot", "ggh4x", "here", 
               "HMDHFDplus", "patchwork", 
               "plyr", "tidyverse")

## check if installed, else install
for(i in from.cran){
  
  if(system.file(package = i) == ""){install.packages(i)}
  
}

## load packages
library(tidyverse)

############
# SET PATH #
############

here::i_am("scripts/00_analysis.R")

##################
# CREATE FOLDERS #
##################

## output folder
if(!dir.exists(here::here("output"))){dir.create(here::here("output"))}

#############
# FUNCTIONS #
#############

## life table from all-cause mortality rates
life.table <-
  function(df, radix = 1){
    
    df %>%
      mutate(
        
        Age = Age,
        mx = mx,
        px = ifelse(Age > 0, 
                    ## px = 1 - qx
                    ## qx = mx / (1 + (1 - (ax)) * mx)
                    ## for age groups > 0: ax values are assumed to be 0.5
                    1 - (mx / (1 + (1 - (0.5)) * mx)),
                    ## for age group 0: ax value is based on Table 3-2 in Andreev & Kingkade (2015)
                    ifelse(mx >= 0.06891,
                           1 - (mx / (1 + (1 - (0.31411)) * mx)),
                           ifelse(mx >= 0.01724,
                                  1 - (mx / (1 + (1 - (0.04667 + 3.88089 * mx)) * mx)),
                                  1 - (mx / (1 + (1 - (0.14903 - 2.05527 * mx)) * mx))))), 
        lx = head(cumprod(c(radix, px)), -1),
        dx = lx * (1 - px),
        Lx = ifelse(mx == 0, lx, dx / mx),
        .keep = "none"
        
      )
    
  }

#######################
# AUXILIARY DATA SETS #
#######################

## vector with all selected countries
select.countries <- c(select.country, countries) 

## data set with HMD country codes
HMDcode <-
  rbind(
    c("Australia", "AUS"), 
    c("Austria", "AUT"),
    c("Belarus", "BLR"),
    c("Belgium", "BEL"),
    c("Bulgaria", "BGR"),
    c("Canada", "CAN"),
    c("Chile", "CHL"),
    c("Croatia", "HRV"),
    c("Czech Republic", "CZE"),
    c("Denmark", "DNK"),
    c("England and Wales", "GBRTENW"),
    c("Estonia", "EST"),
    c("Finland", "FIN"),
    c("France", "FRATNP"),
    c("Germany", "DEUTNP"),
    c("Greece", "GRC"),
    c("Hong Kong", "HKG"),
    c("Hungary", "HUN"),
    c("Iceland", "ISL"),
    c("Ireland", "IRL"),
    c("Israel", "ISR"),
    c("Italy", "ITA"),
    c("Japan", "JPN"),
    c("Latvia", "LVA"),
    c("Lithuania", "LTU"),
    c("Luxembourg", "LUX"),
    c("Netherlands", "NLD"),
    c("New Zealand", "NZL_NP"),
    c("Northern Ireland", "GBR_NIR"),
    c("Norway", "NOR"),
    c("Poland", "POL"),
    c("Portugal", "PRT"),
    c("Russia", "RUS"),
    c("Scotland", "GBR_SCO"),
    c("Slovakia", "SVK"),
    c("Slovenia", "SVN"),
    c("South Korea", "KOR"),
    c("Spain", "ESP"),
    c("Sweden", "SWE"),
    c("Switzerland", "CHE"),
    c("Taiwan", "TWN"),
    c("Ukraine", "UKR"),
    c("United Kingdom", "GBR_NP"),
    c("USA", "USA")
  ) %>% 
  data.frame()

names(HMDcode) <- c("Name", "Code")

########################
# CALCULATE INDICATORS #
########################

frameworks.combined <- 
  map_df(select.countries, function(x){
    
    ####################
    # GET COUNTRY CODE #
    ####################

    code <- HMDcode %>% filter(Name == x) %>% pull(Code)
    
    #####################
    # LOAD COUNTRY DATA #
    #####################
    
    ## cohort age-specific fertility rates
    CNTRY.fx <-
      HMDHFDplus::readHFD(unz(description = here::here("data", "asfr.zip"), 
                              filename = "asfrVH.txt")) %>% 
      filter(Code == code) %>% 
      select(Cohort, Age, fx = ASFR)
    
    ## cohort age-specific mortality rates
    CNTRY.mx <-
      HMDHFDplus::readHMD(unz(description = here::here("data", "c_death_rates.zip"), 
                              filename = paste0("cMx_1x1/", code, ".cMx_1x1.txt"))) %>% 
      select(Cohort = Year, Age, mx = Female)
    
    #############################################
    # FIND COHORTS WITH COMPLETE FERTILITY DATA # 
    #############################################
    
    CNTRY.fx.complete <- 
      CNTRY.fx %>% 
      filter(!is.na(fx)) %>% 
      mutate(min = min(Age),
             max = max(Age), .by = Cohort) %>% 
      filter(min <= min(fert.ages),
             max >= max(fert.ages),
             Age %in% min(fert.ages):max(fert.ages)) %>% 
      mutate(n = n(), .by = Cohort) %>% 
      filter(n == length(min(fert.ages):max(fert.ages))) %>% 
      mutate(alpha = cumsum(fx),
             beta = rev(cumsum(rev(fx))), .by = Cohort) %>% 
      filter(alpha != 0, beta != 0) %>% 
      select(Cohort, Age, fx)
    
    #############################################
    # FIND COHORTS WITH COMPLETE MORTALITY DATA #
    #############################################
    
    CNTRY.mx.complete <-
      CNTRY.mx %>% 
      filter(!is.na(mx)) %>% 
      mutate(min = min(Age),
             max = max(Age), .by = Cohort) %>% 
      filter(min == 0,
             max >= max(fert.ages),
             Age %in% 0:max(fert.ages)) %>% 
      mutate(n = n(), .by = Cohort) %>% 
      filter(n == length(0:max(fert.ages))) %>% 
      select(Cohort, Age, mx)
    
    ###########################################################
    # FIND COHORTS WITH COMPLETE FERTILITY AND MORTALITY DATA #
    ###########################################################
    
    cohorts <-
      intersect(unique(CNTRY.fx.complete$Cohort),
                unique(CNTRY.mx.complete$Cohort))
    
    ####################
    # GET COHORT ASFRS #
    ####################
    
    CNTRY.ma <- 
      CNTRY.fx.complete %>% 
      filter(Cohort %in% cohorts) %>% 
      select(Cohort, Age, ma = fx)
    
    ########################
    # GET COHORT EXPOSURES #
    ########################
    
    CNTRY.La <-
      CNTRY.mx.complete %>% 
      filter(Cohort %in% cohorts) %>%
      group_by(Cohort) %>% 
      group_modify(~ {
        
        life.table(df = .x)
        
      }) %>% 
      ungroup() %>% 
      select(Cohort, Age, La = Lx)
    
    #################################################
    # CALCULATE INDICATORS FOR DIFFERENT FRAMEWORKS #
    ################################################# 
    
    ## born once, die once (translated)
    CNTRY.NRR <-
      CNTRY.ma %>%
      left_join(CNTRY.La) %>% 
      mutate(NRR = sum(La * ma),
             NRRa = rev(cumsum(rev(La * ma))),
             Tempo = sum((Age + 0.5) * La * ma) / NRR, .by = Cohort) %>% 
      mutate(f = (La * ma) / NRR,
             S = NRRa / NRR,
             h = (La * ma) / NRRa,
             Country = x,
             Framework = "NRR") %>% 
      select(Country, Cohort, Age, f, S, h, Tempo, Quantum = NRRa, Framework)
    
    ## conventional (extended)
    CNTRY.GRR <-
      CNTRY.ma %>%
      mutate(GRR = sum(ma),
             GRRa = rev(cumsum(rev(ma))),
             Tempo = sum((Age + 0.5) * ma) / GRR, .by = Cohort) %>% 
      mutate(f = ma / GRR,
             S = GRRa / GRR,
             h = ma / GRRa,
             Country = x,
             Framework = "GRR") %>% 
      select(Country, Cohort, Age, f, S, h, Tempo, Quantum = GRRa, Framework)
    
    ## combine frameworks
    df <- 
      CNTRY.NRR %>% 
      add_row(CNTRY.GRR) %>% 
      pivot_wider(names_from = "Framework",
                  names_glue = "{.value}.{Framework}",
                  values_from = c("f", "S", "h", "Tempo", "Quantum"))
    
    df
    
  })

#################################
# CREATE CONSTANTS FOR PLOTTING #
#################################

## cohort selection
cohort.min <- plyr::round_any(min(frameworks.combined$Cohort), f = ceiling, accuracy = 20)
cohort.max <- plyr::round_any(max(frameworks.combined$Cohort), f = floor, accuracy = 20)
cohorts <- seq(cohort.min, cohort.max, 10)

#############################
# PREPARE DATA FOR PLOTTING #
#############################

## ratios of density, survival, and hazard by country and cohort
fSh.combined <-
  frameworks.combined  %>% 
  select(Country, Cohort, Age, f.NRR, f.GRR, S.NRR, S.GRR, h.NRR, h.GRR) %>% 
  mutate(f.OMEGA = f.NRR / f.GRR,
         S.OMEGA = S.NRR / S.GRR,
         h.OMEGA = h.NRR / h.GRR) %>% 
  pivot_longer(cols = starts_with(c("f.", "S.", "h.")),
               names_to = "Indicator",
               values_to = "Value") %>% 
  separate_wider_delim(Indicator, delim = ".", names = c("Indicator", "Framework")) %>% 
  mutate(Indicator = factor(Indicator, 
                            levels = c("f", "S", "h"), 
                            labels = c("Density", "Survival", "Hazard")),
         Framework = factor(Framework, 
                            levels = c("NRR", "GRR", "OMEGA"),
                            labels = c("B1D1", "Conventional", "Ratio")))

## formulas for ratios of density, survival, and hazard
label.df <-
  data.frame(
    Indicator =  rep(c("Density", "Survival", "Hazard"), 3),
    Framework = rep(c("B1D1", "Conventional", "Ratio"), each = 3),
    Label = c("frac('\U2113(a)m(a)', 'NRR')", "frac('NRR(a)', 'NRR')", "frac('\U2113(a)m(a)', 'NRR(a)')",
              "frac('m(a)', 'GRR')", "frac('GRR(a)', 'GRR')", "frac('m(a)', 'GRR(a)')",
              "frac('\U2113(a)', '\U03A9')", "frac('\U03A9(a)', '\U03A9')", "frac('\U2113(a)', '\U03A9(a)')"),
    x = rep(Inf, 9),
    y = rep(Inf, 9)
  ) %>% mutate(hjust = 1.1,
               vjust = 1.1,
               Indicator = factor(Indicator, levels = c("Density", "Survival", "Hazard")),
               Framework = factor(Framework, levels = c("B1D1", "Conventional", "Ratio")))

## OMEGA by country and cohort
OMEGA.combined <-
  frameworks.combined %>% 
  mutate(min = min(Age), .by = c(Country, Cohort)) %>%
  filter(Age == min,
         Cohort %in% cohort.min:cohort.max) %>% 
  mutate(OMEGA = Quantum.NRR / Quantum.GRR,
         CNTRY = Country == select.country,
         Indicator = "\U03A9")

## tempo by country and cohort
tempo.combined <-
  frameworks.combined %>% 
  mutate(min = min(Age), .by = c(Country, Cohort)) %>%
  filter(Age == min,
         Cohort %in% cohort.min:cohort.max) %>% 
  mutate(tempo = Tempo.NRR / Tempo.GRR,
         CNTRY = Country == select.country,
         Indicator = "Coale's MAB / MAB")

## l(a) by cohort for main country
la.cohort <-  
  HMDHFDplus::readHMD(unz(description = here::here("data", "c_death_rates.zip"), 
                          filename = paste0("cMx_1x1/", 
                                            HMDcode$Code[which(HMDcode$Name == select.country)], 
                                            ".cMx_1x1.txt"))) %>% 
  select(Cohort = Year, Age, mx = Female) %>% 
  filter(Age <= max(plot.ages),
         Cohort %in% 
           unique(frameworks.combined$Cohort[frameworks.combined$Country == select.country]),
         Cohort %in% cohorts) %>% 
  group_by(Cohort) %>% 
  group_modify(~ {
    
    life.table(df = .x)
    
  }) %>% 
  ungroup() %>% 
  filter(Age >= min(plot.ages)) %>% 
  mutate(Indicator = "\U2113(a)") %>% 
  select(Indicator, Cohort, Age, la = lx)

## OMEGA(a) by cohort for main country
OMEGAa.cohort <-
  frameworks.combined %>% 
  filter(Country == select.country,
         Age %in% min(plot.ages):max(plot.ages),
         Cohort %in% cohorts) %>% 
  mutate(OMEGAa = Quantum.NRR / Quantum.GRR,
         Indicator = paste0("\U03A9", "(a)")) 

############
# FIGURE 1 #
############

## plot all countries in loop 
for(x in select.countries){
  
  ## determine axis ranges
  density.range <- 
    fSh.combined %>% 
    filter(Country == x, 
           Age %in% min(plot.ages):max(plot.ages),
           Cohort %in% cohorts, 
           Indicator == "Density",
           Framework != "Ratio") %>%
    arrange(-Value) %>% 
    slice_head(n = 1) %>%
    pull(Value)
  
  density.range <- plyr::round_any(density.range, f = ceiling, accuracy = 0.1)
  
  ratio.range <- 
    fSh.combined %>% 
    filter(Country == x, 
           Age %in% min(plot.ages):max(plot.ages),
           Cohort %in% cohorts, 
           Framework == "Ratio") %>%
    group_by(Value > 1) %>%
    mutate(Value = ifelse(Value > 1, Value, 1 / Value)) %>% 
    arrange(-Value) %>% 
    slice_head(n = 1) %>%
    ungroup() %>% 
    pull(Value)
  
  ratio.range[1] <- plyr::round_any(1 / ratio.range[1], f = floor, accuracy = 0.1)
  ratio.range[2] <- plyr::round_any(ratio.range[2], f = ceiling, accuracy = 0.1)
  
  ## plot
  fig.1 <-
    ggplot(fSh.combined %>% 
             filter(Country == x,
                    Age %in% min(plot.ages):max(plot.ages), 
                    Cohort %in% cohorts),
           aes(x = Age, y = Value, group = Cohort, color = Cohort)) +
    geom_line(alpha = 0.5) +
    geom_label(data = label.df, 
               aes(x = x, y = y, label = bquote(.(Label)),
                   hjust = hjust, vjust = vjust), 
               inherit.aes = FALSE, parse = TRUE, size = 3.5, label.size = NA, fill = NA, family = "serif") +
    scale_x_continuous(breaks = seq(min(plot.ages), max(plot.ages), 5)) +
    scale_color_viridis_c(limits = c(cohort.min, cohort.max),
                          breaks = c(cohort.min, (cohort.min + cohort.max) / 2, cohort.max),
                          end = 0.85, option = "inferno") +
    theme_bw() +
    labs(y = "") +
    theme(aspect.ratio = 1,
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "top",
          legend.title.position = "top",
          legend.title = element_text(hjust = 0.5),
          legend.ticks = element_blank(),
          text = element_text(family = "serif")) +
    ggh4x::facet_grid2(Framework ~ Indicator, scales = "free_y", independent = "y") +
    ggh4x::facetted_pos_scales(y = list(scale_y_continuous(limits = c(0, density.range),
                                                           breaks = seq(0, density.range, length.out = 5)),
                                        scale_y_continuous(limits = c(0, 1),
                                                           breaks = seq(0, 1, length.out = 5)),
                                        scale_y_continuous(limits = c(0, 1),
                                                           breaks = seq(0, 1, length.out = 5)),
                                        scale_y_continuous(limits = c(0, density.range),
                                                           breaks = seq(0, density.range, length.out = 5)),
                                        scale_y_continuous(limits = c(0, 1),
                                                           breaks = seq(0, 1, length.out = 5)),
                                        scale_y_continuous(limits = c(0, 1),
                                                           breaks = seq(0, 1, length.out = 5)),
                                        scale_y_continuous(limits = c(ratio.range[1], ratio.range[2]),
                                                           breaks = seq(ratio.range[1], ratio.range[2], length.out = 5)),
                                        scale_y_continuous(limits = c(ratio.range[1], ratio.range[2]),
                                                           breaks = seq(ratio.range[1], ratio.range[2], length.out = 5)),
                                        scale_y_continuous(limits = c(ratio.range[1], ratio.range[2]),
                                                           breaks = seq(ratio.range[1], ratio.range[2], length.out = 5))))
  
  ## save plot
  ## pdf
  cairo_pdf(filename = paste0(here::here("output"), "/Figure-1-", x, ".pdf"),
            height = 8, width = 8)
  
  print(fig.1)
  
  dev.off()
  
  ## svg
  ggsave(fig.1,
         filename = paste0(here::here("output"), "/Figure-1-", x, ".svg"),
         height = 8, width = 8)
  
  ## eps
  cairo_ps(filename = paste0(here::here("output"), "/Figure-1-", x, ".eps"),
           width = 8, 
           height = 8)
  
  print(fig.1)
  
  dev.off()
  
}

#########################################
# BUILD SUB-FIGURES FOR FIGURES 2 and 3 #
#########################################

## determine axis ranges
OMEGA.min <- plyr::round_any(min(OMEGA.combined$OMEGA), f = floor, accuracy = 0.1)

tempo.min <- plyr::round_any(min(tempo.combined$tempo), f = floor, accuracy = 0.01)
tempo.max <- plyr::round_any(max(tempo.combined$tempo), f = ceiling, accuracy = 0.01)

la.min <- plyr::round_any(min(la.cohort$la), f = floor, accuracy = 0.1)
OMEGAa.min <- plyr::round_any(min(OMEGAa.cohort$OMEGAa), f = floor, accuracy = 0.1)
axis.min <- min(la.min, OMEGAa.min)

## OMEGA plot
OMEGA.plot <-
  ggplot(OMEGA.combined , 
         aes(x = Cohort, y = OMEGA, group = Country, color = CNTRY)) +
  scale_x_continuous(breaks = seq(cohort.min, cohort.max, 10)) +
  scale_y_continuous(breaks = seq(OMEGA.min, 1, 0.1),
                     minor_breaks = seq(OMEGA.min, 1, 0.05)) +
  guides(y = guide_axis(minor.ticks = TRUE),
         color = guide_legend(position = "inside")) +
  scale_color_manual("", breaks = c("TRUE", "FALSE"), labels = c(paste(select.country), "Other"), values = c("black", "gray80")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title.y = element_blank()) +
  coord_cartesian(xlim = c(cohort.min, cohort.max),
                  ylim = c(OMEGA.min, 1))

## tempo plot
tempo.plot <-   
  ggplot(tempo.combined, 
         aes(x = Cohort, y = tempo, group = Country, color = CNTRY)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(cohort.min, cohort.max, 10)) +
  scale_y_continuous(breaks = seq(tempo.min, tempo.max, 0.01),
                     minor_breaks = seq(tempo.min, tempo.max, 0.005)) +
  guides(y = guide_axis(minor.ticks = TRUE),
         color = guide_legend(position = "inside")) +
  scale_color_manual("", breaks = c("TRUE", "FALSE"), labels = c(paste(select.country), "Other"), values = c("black", "gray80")) +
  facet_wrap(~ Indicator) +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title.y = element_blank(),
        strip.background = element_rect(color = "transparent", fill = "white"),
        strip.text = element_text(size = 14),
        legend.position.inside = c(0.2, 0.1),
        text = element_text(family = "serif")) +
  coord_cartesian(xlim = c(cohort.min, cohort.max),
                  ylim = c(tempo.min, tempo.max)) +
  ggh4x::force_panelsizes(total_width = unit(5, "in"),
                          total_height = unit(5, "in"))

## l(a) plot
la.plot <-
  ggplot(la.cohort,
         aes(x = Age, y = la, group = Cohort, color = Cohort)) +
  geom_line(linewidth = 0.5) +
  scale_x_continuous(breaks = seq(min(plot.ages), max(plot.ages), 5)) +
  scale_y_continuous(breaks = seq(axis.min, 1, 0.1),
                     minor_breaks = seq(axis.min, 1, 0.05)) +
  guides(y = guide_axis(minor.ticks = TRUE)) +
  scale_color_viridis_c(limits = c(cohort.min, cohort.max),
                        breaks = c(cohort.min, (cohort.min + cohort.max) / 2, cohort.max),
                        labels = NULL,
                        end = 0.85, option = "inferno") +
  facet_wrap(~ Indicator) +
  theme_bw() +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        legend.position = "bottom",
        legend.title.position = "bottom",
        legend.title = element_text(size = 14, hjust = 0.5),
        legend.ticks = element_blank(),
        legend.key.width = unit(0.909, "in"),
        legend.margin = margin(c(0, 0, 0, 0)),
        text = element_text(family = "serif")) +
  coord_cartesian(xlim = c(min(plot.ages), max(plot.ages)),
                  ylim = c(axis.min, 1)) 

## OMEGA(a) plot
OMEGAa.plot <-
  ggplot(OMEGAa.cohort, 
         aes(x = Age, y = OMEGAa, group = Cohort, color = Cohort)) +
  geom_line(linewidth = 0.5) +
  scale_x_continuous(breaks = seq(min(plot.ages), max(plot.ages), 5)) +
  scale_y_continuous(breaks = seq(axis.min, 1, 0.1),
                     minor_breaks = seq(axis.min, 1, 0.05)) +  
  guides(y = guide_axis(minor.ticks = TRUE)) +
  scale_color_viridis_c(limits = c(cohort.min, cohort.max),
                        breaks = c(cohort.min, (cohort.min + cohort.max) / 2, cohort.max),
                        end = 0.85, option = "inferno") +
  facet_wrap(~ Indicator) +
  theme_bw() +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        text = element_text(family = "serif")) +
  coord_cartesian(xlim = c(min(plot.ages), max(plot.ages)),
                  ylim = c(axis.min, 1)) 

##################
# BUILD FIGURE 2 #
##################

## get legend
g <- cowplot::get_plot_component(la.plot, "guide-box", return_all = TRUE)

## build plot
fig.2 <-
  cowplot::ggdraw(OMEGA.plot +
                    geom_line(linewidth = 1) + 
                    facet_wrap(~ Indicator) +
                    theme(strip.background = element_rect(color = "transparent", fill = "white"),
                          strip.text = element_text(size = 14),
                          axis.title.x = element_blank(),
                          axis.text = element_text(size = 12),
                          axis.ticks.x = element_blank(),
                          panel.grid = element_blank(),
                          plot.margin = margin(c(0.25, 0.25, 2.5, 0.25), unit = "cm"),
                          legend.position.inside = c(0.3, 0.1),
                          text = element_text(family = "serif"))  +
                    ggh4x::force_panelsizes(total_width = unit(5, "in"),
                                            total_height = unit(5, "in"))) +
  cowplot::draw_plot(OMEGAa.plot + 
                       theme(legend.position = "none"), 
                     x = 0.15, y = 0.57,
                     height = 0.35, width = 0.35) +
  cowplot::draw_plot(la.plot + 
                       theme(legend.position = "none"), 
                     x = 0.55, y = 0.2,
                     height = 0.35, width = 0.35) +
  cowplot::draw_plot(g[[3]]$grobs[[1]], 0.02125, -0.84)

## save plot
## pdf
cairo_pdf(filename = here::here("output", paste0("Figure-2.pdf")),
          width = 7, 
          height = 7)

print(fig.2)

dev.off()

## svg
ggsave(fig.2,
       filename = here::here("output", paste0("Figure-2.svg")),
       width = 7, 
       height = 7)

## eps
cairo_ps(filename = here::here("output", paste0("Figure-2.eps")),
         width = 7, 
         height = 7)

print(fig.2)

dev.off()

##################
# BUILD FIGURE 3 #
##################

## build plot
fig.3 <-
  cowplot::ggdraw(tempo.plot +
                    theme(axis.title.x = element_text(size = 14, vjust = -1),
                          axis.text = element_text(size = 12),
                          axis.ticks.x = element_blank(),
                          panel.grid = element_blank(),
                          plot.margin = margin(c(0.25, 0.25, 0.25, 0.25), unit = "cm"))) +
  cowplot::draw_plot(OMEGA.plot + 
                       geom_line(linewidth = 0.5, show.legend = FALSE) +
                       scale_x_continuous(breaks = seq(cohort.min, cohort.max, 20),
                                          minor_breaks = seq(cohort.min, cohort.max, 10)) +
                       guides(x = guide_axis(minor.ticks = TRUE)) +
                       facet_wrap(~ Indicator, labeller = label_bquote(Omega)) +
                       theme(strip.background = element_blank(),
                             strip.text = element_text(size = 10),
                             axis.title.x = element_text(size = 10),
                             axis.text = element_text(size = 8),
                             axis.ticks = element_blank(),
                             panel.border = element_blank(),
                             plot.background = element_blank(),
                             text = element_text(family = "serif")), 
                     x = 0.44, y = 0.175,
                     height = 0.45, width = 0.45)

## save plot
## pdf
cairo_pdf(filename = here::here("output", paste0("Figure-3.pdf")),
          width = 7, 
          height = 7)

print(fig.3)

dev.off()

## svg
ggsave(fig.3,
       filename = here::here("output", paste0("Figure-3.svg")),
       width = 7, 
       height = 7)

## eps
cairo_ps(filename = here::here("output", paste0("Figure-3.eps")),
         width = 7, 
         height = 7)

print(fig.3)

dev.off()
