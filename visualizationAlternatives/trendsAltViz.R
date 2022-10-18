source('../SoIB_v2 functions.R')
library(tidyverse)
library(ggdist)
library(ggridges)

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
rawTrends_distr <- rawTrends_distr %>% filter((timegroupsf %in% tg[! tg %in% c("before 2000")])) # %>% 
  # filter(species == "Indian Courser")
densRidgesPlot <- ggplot(rawTrends_distr, aes(x, timegroupsf, height = norm, fill = after_scale(alpha(colour, 0.4)), color = species)) +
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
