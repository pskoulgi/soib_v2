source('../SoIB_v2 functions.R')
library(tidyverse)
library(ggdist)
library(ggridges)
library(ggpubr)
library(wesanderson)

#### Read trends, standardize to pre-2000 ####
rawTrends = read.csv("../assorted_trends_1.csv")
trendsStd = stdtrends(rawTrends) %>% 
  mutate(nmfreqbyspec = nmfreqbyspec - 100)

#### Expand each time group to years it covers ####
tg <- c("before 2000", "2000-2006", "2007-2010", "2011-2012", "2013", "2014", "2015", "2016", "2017", "2018")
bef2000 <- tibble(
  year = c(1994.5, 1995, 1996, 1997, 1998, 1999, 1999.5),
  timeBand = rep(c(tg[1]))
)
from2000To2006 <- tibble(
  year = c(1999.5, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2006.5),
  timeBand = rep(c(tg[2]))
)
from2007To2010 <- tibble(
  year = c(2006.5, 2007, 2008, 2009, 2010, 2010.5),
  timeBand = rep(c(tg[3]))
)
from2011To2012 <- tibble(
  year = c(2010.5, 2011, 2012, 2012.5),
  timeBand = rep(c(tg[4]))
)
for2013 <- tibble(
  year = c(2012.5, 2013, 2013.5),
  timeBand = rep(c(tg[5]))
)
for2014 <- tibble(
  year = c(2013.5, 2014, 2014.5),
  timeBand = rep(c(tg[6]))
)
for2015 <- tibble(
  year = c(2014.5, 2015, 2015.5),
  timeBand = rep(c(tg[7]))
)
for2016 <- tibble(
  year = c(2015.5, 2016, 2016.5),
  timeBand = rep(c(tg[8]))
)
for2017 <- tibble(
  year = c(2016.5, 2017, 2017.5),
  timeBand = rep(c(tg[9]))
)
for2018 <- tibble(
  year = c(2017.5, 2018, 2018.5),
  timeBand = rep(c(tg[10]))
)
hsBef2000 <- left_join(bef2000, trendsStd, by = c("timeBand" = "timegroupsf"))
hs2000To2006 <- left_join(from2000To2006, trendsStd, by = c("timeBand" = "timegroupsf"))
hs2007To2010 <- left_join(from2007To2010, trendsStd, by = c("timeBand" = "timegroupsf"))
hs2011To2012 <- left_join(from2011To2012, trendsStd, by = c("timeBand" = "timegroupsf"))
hs2013 <- left_join(for2013, trendsStd, by = c("timeBand" = "timegroupsf"))
hs2014 <- left_join(for2014, trendsStd, by = c("timeBand" = "timegroupsf"))
hs2015 <- left_join(for2015, trendsStd, by = c("timeBand" = "timegroupsf"))
hs2016 <- left_join(for2016, trendsStd, by = c("timeBand" = "timegroupsf"))
hs2017 <- left_join(for2017, trendsStd, by = c("timeBand" = "timegroupsf"))
hs2018 <- left_join(for2018, trendsStd, by = c("timeBand" = "timegroupsf"))
hsWithYear <- rbind(
  hsBef2000, hs2000To2006, hs2007To2010, hs2011To2012,
  hs2013, hs2014, hs2015, hs2016, hs2017, hs2018) %>% 
  mutate(
    lower = nmfreqbyspec - (nmsebyspec/2),
    upper = nmfreqbyspec + (nmsebyspec/2),
    yLabel = paste0(year)
  )
hsWithYear$timeBand <- factor(
  hsWithYear$timeBand, 
  levels = tg
)

#### Plot time series as line ribbons ####

# Draxing ticks before and after labels https://stackoverflow.com/a/44258401
x_tick <- seq(1994, 2018) + 0.5
len <- length(x_tick)

# lineRibbonPlot <- ggplot(hsWithYear %>% filter(species == "Indian Courser"), aes(x = year, y = nmfreqbyspec, ymin = lower, ymax = upper, color = species, fill = species)) + 
lineRibbonPlot <- ggplot(hsWithYear, aes(x = year, y = nmfreqbyspec, ymin = lower, ymax = upper, color = species, fill = species)) +
  geom_lineribbon(alpha = 0.5) +
  scale_x_continuous(breaks = c(seq(1994, 2018), x_tick),
                     labels = c("...", "...", paste0(seq(1996, 2018)), rep(c(""), len))) +
  scale_y_continuous(
    breaks = c(-100, -50, 0, 50, 100, 150, 200),
    labels = c("-100", "-50", "Pre-2000 baseline", "+50", "+100", "+150", "+200")) +
  coord_cartesian(ylim = c(-100, 200)) +
  theme_minimal() +
  theme(
    axis.ticks.x = element_line(color = c(rep(NA, len), rep("black", len-1))),
    panel.grid.major.y = element_line(linetype = "dotted", size = 0.6, color = "darkgray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()) +
  xlab("years") +
  ylab("change in eBird abundance index \n") +
  ggtitle("Line at mean, and width of ribbon as se, of estimated frequencies")
lineRibbonPlot
# ggsave("outputs/trendsAsLineribbons.png", lineRibbonPlot, width = 11, height = 7, bg = "white")

#### Plot time series as normal distributions ####

# Convert mean, se to normal distribution https://stackoverflow.com/a/56446124
nn = 5000
hsWYear_distr <- hsWithYear %>%
  mutate(low  = nmfreqbyspec - 3 * nmsebyspec, high = nmfreqbyspec + 3 * nmsebyspec) %>%
  uncount(nn, .id = "row") %>%
  mutate(x    = (1 - row/nn) * low + row/nn * high, 
         norm = 10000*dnorm(x, nmfreqbyspec, nmsebyspec))

rawTrends_distr <- trendsStd %>% 
  mutate(low  = nmfreqbyspec - 3 * nmsebyspec, high = nmfreqbyspec + 3 * nmsebyspec) %>%
  uncount(nn, .id = "row") %>%
  mutate(x    = (1 - row/nn) * low + row/nn * high, 
         norm = dnorm(x, nmfreqbyspec, nmsebyspec))
rawTrends_distr$timegroupsf <- factor(
  rawTrends_distr$timegroupsf, 
  levels = tg
)

# Drop pre-2000, and plot
rawTrends_distr1 <- rawTrends_distr %>% filter((timegroupsf %in% tg[! tg %in% c("before 2000")])) # %>% 
# filter(species == "Indian Courser")
densRidgesPlot <- ggplot(rawTrends_distr1, aes(x, timegroupsf, height = norm, fill = after_scale(alpha(colour, 0.4)), color = species)) +
  geom_density_ridges(stat = "identity", scale = 0.9) +
  geom_point(aes(x = nmfreqbyspec, y = timegroupsf, size = 2.5)) +
  scale_x_continuous(
    breaks = c(-200, -100, 0, 100, 200),
    labels = c("-200%", "-100%", "Pre-2000 baseline", "+100%", "+200%"),
    expand = c(0, 0),
    limits = c(-200, 200)) +
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.3))) +
  theme_ridges(center_axis_labels = TRUE) +
  theme(
    panel.grid.major.x = element_line(linetype = "dotted", size = 0.6, color = "darkgray"),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.ticks.x = element_line("black")) +
  geom_vline(xintercept = 0) +
  guides(size = "none") +
  ylab("time steps \n") +
  xlab("\n change in eBird abundance index") +
  ggtitle("Using mean and se of estimated frequencies as \n mean and standard deviation, respectively, of Normal distributions")
densRidgesPlot
# ggsave("outputs/trendsAsDensityRidges.png", densRidgesPlot, width = 11, height = 7, bg = "white")

x_tick_pre2000Bas <- seq(1999, 2018) + 0.5
lineribbonsDiscPlot <- rawTrends_distr1 %>% filter(species == "Indian Courser") %>% 
  group_by(timegroups, species) %>% 
  median_qi(x, .width = c(.50, .80, .95)) %>% 
  ggplot(aes(x = timegroups, y = x, ymin = .lower, ymax = .upper, color = species)) +
  geom_lineribbon() +
  geom_point(aes(size = 0.75)) +
  scale_fill_brewer() +
  scale_color_brewer(palette = "Greens") +
  geom_hline(yintercept = 0) +
  scale_x_continuous(
    breaks = c(seq(1999, 2018), x_tick_pre2000Bas),
    labels = c("", "2000", "2001", rep(c(""), 2006-2000-2), 
               "2006", "2007", rep(c(""), 2010-2006-2), 
               "2010", "2011", rep(c(""), 2012-2010-2), 
               paste0(seq(2012, 2018)), rep(c(""), length(x_tick_pre2000Bas))),
    limits = c(1999.5, 2018)) +
  scale_y_continuous(
    breaks = c(-100, 0, 100, 200),
    labels = c("-100%", "Pre-2000 baseline", "+100%", "+200%"),
    expand = c(0, 0),
    limits = c(-150, 300)) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(color = "darkgray", size = 0.6),
    panel.grid.major.y = element_line(linetype = "dotted", size = 0.6, color = "darkgray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    # axis.ticks.x = element_line(color = c(rep(NA, length(x_tick_pre2000Bas)), rep("black", length(x_tick_pre2000Bas)-1)))
    axis.ticks.x = element_line(
      color = c(rep(NA, length(x_tick_pre2000Bas)),
                NA, "black", rep(NA, 2006-2000-1),
                "black", rep(NA, 2010-2006-1),
                "black", rep(NA, 2012-2010-1),
                rep("black", 2018-2012+1)))) +
  guides(size = "none", color = "none") +
  xlab("\n time steps") +
  ylab("change in eBird abundance index \n") +
  labs(title = "Indian Courser (declining)")
lineribbonsDiscPlot
# ggsave("outputs/trendsAsLineribbonsDisc.png", lineribbonsDiscPlot, width = 11, height = 7, bg = "white")

lineribbonsDiscPlot_blockXaxis <- rawTrends_distr1 %>% filter(species == "Indian Courser") %>% 
  group_by(timegroups, species) %>% 
  median_qi(x, .width = c(.50, .80, .95)) %>% 
  ggplot(aes(x = timegroups, y = x, ymin = .lower, ymax = .upper, color = species)) +
  geom_lineribbon() +
  geom_point(aes(size = 0.75)) +
  scale_fill_brewer() +
  scale_color_brewer(palette = "Greens") +
  geom_hline(yintercept = 0) +
  geom_bracket(
    inherit.aes = FALSE, 
    xmin = c(2000, 2006, 2010, 2012, seq(2013, 2017)) + 0.5, 
    xmax = c(2006, 2010, 2012, seq(2013, 2018)) + 0.5,
    y.position = -130,
    bracket.shorten = 0.15,
    tip.length = 0.05,
    vjust = 3,
    label = tg[-1],
    label.size = 3) +
  scale_x_continuous(limits = c(1999.5, 2018.5)) +
  scale_y_continuous(
    breaks = c(-100, 0, 100, 200),
    labels = c("-100%", "Pre-2000 baseline", "+100%", "+200%"),
    expand = c(0, 0),
    limits = c(-150, 300)) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(linetype = "dotted", size = 0.6, color = "darkgray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  guides(size = "none", color = "none") +
  xlab("\n time steps") +
  ylab("change in eBird abundance index \n") +
  labs(title = "Indian Courser (declining)")
lineribbonsDiscPlot_blockXaxis
# ggsave("outputs/trendsAsLineribbonsDisc_blockXaxis.png", lineribbonsDiscPlot_blockXaxis, width = 11, height = 7, bg = "white")

lineribbonsDiscPlot_blockXaxis_noblockLabel <- rawTrends_distr1 %>% filter(species == "Indian Courser") %>% 
  group_by(timegroups, species) %>% 
  median_qi(x, .width = c(.50, .80, .95)) %>% 
  ggplot(aes(x = timegroups, y = x, ymin = .lower, ymax = .upper, color = species)) +
  geom_lineribbon() +
  geom_point(aes(size = 0.75)) +
  scale_fill_brewer() +
  scale_color_brewer(palette = "Greens") +
  geom_hline(yintercept = 0) +
  geom_bracket(
    inherit.aes = FALSE, 
    xmin = c(2000, 2006, 2010) + 0.5, 
    xmax = c(2006, 2010, 2012) + 0.5,
    y.position = -130,
    bracket.shorten = 0.15,
    # tip.length = 0.05,
    vjust = 3,
    # label = tg[-1],
    label = c("", "", "")) +
  scale_x_continuous(
    breaks = c(seq(1999, 2018), x_tick_pre2000Bas),
    labels = c("", "2000", "2001", rep(c(""), 2006-2000-2), 
               "2006", "2007", rep(c(""), 2010-2006-2), 
               "2010", "2011", rep(c(""), 2012-2010-2), 
               paste0(seq(2012, 2018)), rep(c(""), length(x_tick_pre2000Bas))),
    limits = c(1999.5, 2018)) +
  scale_y_continuous(
    breaks = c(-100, 0, 100, 200),
    labels = c("-100%", "Pre-2000 baseline", "+100%", "+200%"),
    expand = c(0, 0),
    limits = c(-150, 300)) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(color = "darkgray", size = 0.6),
    panel.grid.major.y = element_line(linetype = "dotted", size = 0.6, color = "darkgray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    # axis.ticks.x = element_line(color = c(rep(NA, length(x_tick_pre2000Bas)), rep("black", length(x_tick_pre2000Bas)-1)))
    axis.ticks.x = element_line(
      color = c(rep(NA, length(x_tick_pre2000Bas)),
                NA, "black", rep(NA, 2006-2000-1),
                "black", rep(NA, 2010-2006-1),
                "black", rep(NA, 2012-2010-1),
                rep("black", 2018-2012+1)))) +
  guides(size = "none", color = "none") +
  xlab("\n time steps") +
  ylab("change in eBird abundance index \n") +
  labs(title = "Indian Courser (declining)")
lineribbonsDiscPlot_blockXaxis_noblockLabel
# ggsave("outputs/trendsAsLineribbonsDisc_blockXaxis_noblockLabel.png", lineribbonsDiscPlot_blockXaxis_noblockLabel, width = 11, height = 7, bg = "white")

# rawTrends_distr1 %>% filter(species == "Indian Courser") %>% 
#   ggplot(aes(x = timegroups, y = x, fill = stat(.width))) +
#   stat_lineribbon(.width = ppoints(50)) +
#   scale_fill_distiller() +
#   # geom_point() +
#   # geom_point(aes(x = timegroups, y = mean(x))) +
#   # geom_point(data = medians, aes(x = timegroups, y = x)) +
#   labs(title = "Indian Courser (declining)")

medians <- rawTrends_distr1 %>% filter(species == "Indian Courser") %>% 
  group_by(timegroups, species) %>% 
  median_qi(x)
lineribbonsContPlot <- ggplot() +
  stat_lineribbon(data = filter(rawTrends_distr1, species == "Indian Courser"), aes(x = timegroups, y = x, fill = stat(.width), color = species), .width = ppoints(50), alpha = 0.2) +
  scale_fill_distiller() +
  scale_color_brewer(palette = "Greens") +
  geom_point(data = medians, aes(x = timegroups, y = x, color = species, size = 0.75)) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(
    breaks = c(-100, 0, 100, 200),
    labels = c("-100%", "Pre-2000 baseline", "+100%", "+200%"),
    expand = c(0, 0),
    limits = c(-150, 300)) +
  scale_x_continuous(
    breaks = c(seq(1999, 2018), x_tick_pre2000Bas),
    labels = c("", "2000", "2001", rep(c(""), 2006-2000-2), 
               "2006", "2007", rep(c(""), 2010-2006-2), 
               "2010", "2011", rep(c(""), 2012-2010-2), 
               paste0(seq(2012, 2018)), rep(c(""), length(x_tick_pre2000Bas))),
    limits = c(1999.5, 2018)) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(color = "darkgray"),
    panel.grid.major.y = element_line(linetype = "dotted", size = 0.6, color = "darkgray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    # legend.position = "bottom",
    legend.title = element_blank(),
    # axis.ticks.x = element_line("darkgray"),
    # axis.ticks.x = element_line(color = c(rep(NA, length(x_tick_pre2000Bas)), rep("black", length(x_tick_pre2000Bas)-1)))
    axis.ticks.x = element_line(
      color = c(rep(NA, length(x_tick_pre2000Bas)),
                NA, "black", rep(NA, 2006-2000-1),
                "black", rep(NA, 2010-2006-1),
                "black", rep(NA, 2012-2010-1),
                rep("black", 2018-2012+1)))) +
  guides(size = "none", color = "none") +
  xlab("\n time steps") +
  ylab("change in eBird abundance index \n") +
  labs(title = "Indian Courser (declining)")
lineribbonsContPlot
# ggsave("outputs/trendsAsLineribbonsCont.png", lineribbonsContPlot, width = 11, height = 7, bg = "white")

# For multiple species in lineribbon plot, refer to
# https://github.com/mjskay/tidybayes/issues/103
# https://github.com/EcoClimLab/growth_phenology/pull/39; 
# https://stackoverflow.com/a/68127360
# https://gradientdescending.com/how-to-use-multiple-color-scales-in-ggplot-with-ggnewscale/

medians_allSpecies <- rawTrends_distr1 %>%
  group_by(timegroups, species) %>% 
  median_qi(x)
lineribbonsDiscPlot_blockXaxis_multSpecies <- ggplot() +
  stat_lineribbon(
    data = rawTrends_distr1, 
    aes(
      x = timegroups, 
      y = x, 
      fill = species, 
      color = species,
      group = paste(group, after_stat(.width))),
    .width = c(.5, .8, .95), 
    alpha = 0.25) +
  scale_color_manual(values = wes_palette("BottleRocket2"), guide = FALSE) +
  scale_fill_manual(values = wes_palette("BottleRocket2")) +
  geom_point(data = medians_allSpecies, aes(x = timegroups, y = x, color = species, size = 0.75)) +
  geom_hline(yintercept = 0) +
  geom_bracket(
    inherit.aes = FALSE, 
    xmin = c(2000, 2006, 2010, 2012, seq(2013, 2017)) + 0.5, 
    xmax = c(2006, 2010, 2012, seq(2013, 2018)) + 0.5,
    y.position = -130,
    bracket.shorten = 0.15,
    vjust = 3,
    label = tg[-1],
    label.size = 3) +
  scale_x_continuous(limits = c(1999.5, 2018.5)) +
  scale_y_continuous(
    breaks = c(-100, 0, 100, 200),
    labels = c("-100%", "Pre-2000 baseline", "+100%", "+200%"),
    expand = c(0, 0),
    limits = c(-150, 300)) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(linetype = "dotted", size = 0.6, color = "darkgray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  coord_cartesian(expand = FALSE, clip = "off") +
  guides(size = "none", color = "none", fill = guide_legend(override.aes= list(alpha = 0.8))) +
  xlab("\n time steps") +
  ylab("change in eBird abundance index \n")
lineribbonsDiscPlot_blockXaxis_multSpecies
# ggsave("outputs/trendsAsLineribbonsDisc_multSpecies.png", lineribbonsDiscPlot_blockXaxis_multSpecies, width = 11, height = 7, bg = "white")

lineribbonsContPlot_blockXaxis_multSpecies <- ggplot() +
  stat_lineribbon(
    data = rawTrends_distr1, 
    aes(
      x = timegroups, 
      y = x, 
      fill = species, 
      color = species,
      group = paste(group, after_stat(.width))),
    .width = ppoints(50),
    alpha = 1/50) +
  scale_color_manual(values = wes_palette("BottleRocket2")) +
  scale_fill_manual(values = wes_palette("BottleRocket2")) +
  geom_point(data = medians_allSpecies, aes(x = timegroups, y = x, color = species, size = 0.75)) +
  geom_hline(yintercept = 0) +
  geom_bracket(
    inherit.aes = FALSE, 
    xmin = c(2000, 2006, 2010, 2012, seq(2013, 2017)) + 0.5, 
    xmax = c(2006, 2010, 2012, seq(2013, 2018)) + 0.5,
    y.position = -130,
    bracket.shorten = 0.15,
    vjust = 3,
    label = tg[-1],
    label.size = 3) +
  scale_x_continuous(limits = c(1999.5, 2018.5)) +
  scale_y_continuous(
    breaks = c(-100, 0, 100, 200),
    labels = c("-100%", "Pre-2000 baseline", "+100%", "+200%"),
    expand = c(0, 0),
    limits = c(-150, 300)) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(linetype = "dotted", size = 0.6, color = "darkgray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  coord_cartesian(expand = FALSE, clip = "off") +
  guides(size = "none", color = "none", fill = guide_legend(override.aes= list(alpha = 0.8))) +
  xlab("\n time steps") +
  ylab("change in eBird abundance index \n")
lineribbonsContPlot_blockXaxis_multSpecies
# ggsave("outputs/trendsAsLineribbonsCont_multSpecies.png", lineribbonsContPlot_blockXaxis_multSpecies, width = 11, height = 7, bg = "white")

# geom_path try
#### Plot as density gradient ridges ####

# Drop pre-2000, and plot
rawTrends_distr5 <- rawTrends_distr %>% filter((timegroupsf %in% tg[! tg %in% c("before 2000")])) %>% 
  filter(species == "Ashy Prinia")
densRidgesPlot5 <- ggplot(rawTrends_distr5, aes(x, timegroupsf, height = norm, fill = stat(x))) +
  # change scale to control how much the ridges overlap vertically
  geom_density_ridges_gradient(stat = "identity", scale = 0.8, rel_min_height = 0.01, gradient_lwd = 1, color = "white", alpha = 0.8) +
  # geom_path(aes(x = nmfreqbyspec, y = timegroupsf)) +
  # geom_path() +
  geom_point(aes(x = nmfreqbyspec, y = timegroupsf, fill = nmfreqbyspec), size = 3, shape = 21, color = "lightslategrey") +
  scale_x_continuous(
    breaks = c(-200, -100, 0, 100, 200),
    labels = c("-200%", "-100%", "Pre-2000 baseline", "+100%", "+200%"),
    expand = c(0, 0),
    limits = c(-200, 200)) +
  scale_fill_gradient2(
    low = "#f1a340", mid = "#f7f7f7", high = "#998ec3", # brewer diverging PuOr https://colorbrewer2.org/#type=diverging&scheme=PuOr&n=3
    breaks = c(-200, -100, 0, 100, 200),
    labels = c("200% decrease", "100% decrease", "Pre-2000 baseline", "100% increase", "200% increase"),
    limits = c(-200, 200)) +
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.3))) +
  theme_ridges(center_axis_labels = TRUE) +
  theme(
    panel.grid.major.x = element_line(linetype = "dotted", size = 0.6, color = "darkgray"),
    legend.title = element_blank(),
    axis.ticks.x = element_line("black")) +
  geom_vline(xintercept = 0) +
  guides(size = "none") +
  ylab("time steps \n") +
  xlab("\n change in eBird abundance index") +
  ggtitle("Ashy Prinia (increasing)")
densRidgesPlot5
# ggsave("outputs/trendsAsDensGradRidges_AshyPrinia.png", densRidgesPlot5, width = 11, height = 7, bg = "white")

#### TRYING multiple species density ridges in the same plot ####

# # Drop pre-2000, and plot
# rawTrends_distr7 <- rawTrends_distr %>% filter((timegroupsf %in% tg[! tg %in% c("before 2000")])) # %>%
#   # filter(species == "Ashy Prinia")
# densRidgesPlot7 <- ggplot(rawTrends_distr7, aes(x, timegroupsf, height = norm, fill = stat(x))) +
#   geom_density_ridges_gradient(stat = "identity", scale = 0.8, rel_min_height = 0.01, gradient_lwd = 1, color = "white", alpha = 0.8) +
#   # geom_point(aes(x = nmfreqbyspec, y = timegroupsf, fill = nmfreqbyspec), size = 3, shape = 21, color = "lightgrey") +
#   geom_point(aes(x = nmfreqbyspec, y = timegroupsf, fill = nmfreqbyspec, shape = species), size = 3, color = "lightgrey") +
#   scale_shape_manual(values = c(21, 22, 23, 24)) +
#   scale_x_continuous(
#     breaks = c(-200, -100, 0, 100, 200),
#     labels = c("-200%", "-100%", "Pre-2000 baseline", "+100%", "+200%"),
#     expand = c(0, 0),
#     limits = c(-200, 200)) +
#   scale_fill_gradient2(
#     low = "#f1a340", mid = "#f7f7f7", high = "#998ec3", # brewer diverging PuOr https://colorbrewer2.org/#type=diverging&scheme=PuOr&n=3
#     breaks = c(-200, -100, 0, 100, 200),
#     labels = c("200% decrease", "100% decrease", "Pre-2000 baseline", "100% increase", "200% increase"),
#     limits = c(-200, 200)) +
#   scale_y_discrete(expand = expansion(mult = c(0.01, 0.3))) +
#   theme_ridges(center_axis_labels = TRUE) +
#   theme(
#     panel.grid.major.x = element_line(linetype = "dotted", size = 0.6, color = "darkgray"),
#     legend.title = element_blank(),
#     axis.ticks.x = element_line("black")) +
#   geom_vline(xintercept = 0) +
#   guides(size = "none") +
#   ylab("time steps \n") +
#   xlab("\n change in eBird abundance index") +
#   ggtitle("Multiple species combined (not working)")
# densRidgesPlot7
# # ggsave("outputs/trendsAsDensGradRidges_MultipleSpecies_notWorking.png", densRidgesPlot7, width = 11, height = 7, bg = "white")

#### TRYING: place legend at the bottom and reorient labels ####
# # Drop pre-2000, and plot
# rawTrends_distr6 <- rawTrends_distr %>% filter((timegroupsf %in% tg[! tg %in% c("before 2000")])) %>% 
#   filter(species == "Red-necked Falcon")
# densRidgesPlot6 <- ggplot(rawTrends_distr6, aes(x, timegroupsf, height = norm, fill = stat(x))) +
#   geom_density_ridges_gradient(stat = "identity", scale = 0.9, rel_min_height = 0.01, gradient_lwd = 1, color = "lightgrey", alpha = 0.8) +
#   geom_point(aes(x = nmfreqbyspec, y = timegroupsf, fill = nmfreqbyspec), size = 3, shape = 21, color = "darkgrey") +
#   scale_x_continuous(
#     breaks = c(-200, -100, 0, 100, 200),
#     labels = c("-200%", "-100%", "Pre-2000 baseline", "+100%", "+200%"),
#     expand = c(0, 0),
#     limits = c(-200, 200)) +
#   scale_fill_gradient2(
#     low = 'red', mid = 'white', high = 'green',
#     breaks = c(-200, -100, 0, 100, 200),
#     labels = c("-200%", "-100%", "Pre-2000 baseline", "+100%", "+200%"),
#     limits = c(-200, 200),
#     guide = guide_legend(
#       label.position = "bottom",
#       direction = "horizontal",
#       label.theme = element_text(angle = 90))) +
#   scale_y_discrete(expand = expansion(mult = c(0.01, 0.3))) +
#   theme_ridges(center_axis_labels = TRUE) +
#   theme(
#     panel.grid.major.x = element_line(linetype = "dotted", size = 0.6, color = "darkgray"),
#     legend.title = element_blank(),
#     legend.position = "bottom",
#     axis.ticks.x = element_line("black")) +
#   geom_vline(xintercept = 0) +
#   guides(size = "none") +
#   ylab("time steps \n") +
#   xlab("\n change in eBird abundance index") +
#   ggtitle("Using mean and se of estimated frequencies as \n mean and standard deviation, respectively, of Normal distributions")
# densRidgesPlot6

#### TRYING ridgeline grey instead of white ####
# # Drop pre-2000, and plot
# rawTrends_distr4 <- rawTrends_distr %>% filter((timegroupsf %in% tg[! tg %in% c("before 2000")])) %>% 
#   filter(species == "Ashy Prinia")
# densRidgesPlot4 <- ggplot(rawTrends_distr4, aes(x, timegroupsf, height = norm, fill = stat(x))) +
#   geom_density_ridges_gradient(stat = "identity", scale = 0.9, rel_min_height = 0.01, gradient_lwd = 1, color = "darkgrey", alpha = 0.8) +
#   geom_point(aes(x = nmfreqbyspec, y = timegroupsf, fill = nmfreqbyspec, size = 1.5), shape = 21, color = "lightgrey") +
#   scale_x_continuous(
#     breaks = c(-200, -100, 0, 100, 200),
#     labels = c("-200%", "-100%", "Pre-2000 baseline", "+100%", "+200%"),
#     expand = c(0, 0),
#     limits = c(-200, 200)) +
#   scale_fill_gradient2(
#     low = 'red', mid = 'white', high = 'green',
#     breaks = c(-200, -100, 0, 100, 200),
#     labels = c("-200%", "-100%", "Pre-2000 baseline", "+100%", "+200%"),
#     limits = c(-200, 200)) +
#   scale_y_discrete(expand = expansion(mult = c(0.01, 0.3))) +
#   theme_ridges(center_axis_labels = TRUE) +
#   theme(
#     panel.grid.major.x = element_line(linetype = "dotted", size = 0.6, color = "darkgray"),
#     # legend.position = "bottom",
#     legend.title = element_blank(),
#     axis.ticks.x = element_line("black")) +
#   geom_vline(xintercept = 0) +
#   guides(size = "none") +
#   ylab("time steps \n") +
#   xlab("\n change in eBird abundance index") +
#   ggtitle("Using mean and se of estimated frequencies as \n mean and standard deviation, respectively, of Normal distributions")
# densRidgesPlot4

