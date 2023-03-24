## install.packages("data.table")
library(data.table)
library(tidyverse)

## set the path
path <- 'E:/APEX-SWAT-GW/APEX-SWAT-GW with R'
setwd(path)
getwd()

## read ****.MSA file in APEX folder
sit1 <-
  read.table(
    'SITE1.MSA',
    skip = 10,
    comment.char = "",
    header = F,
    fill = T
  )
sit1r <-
  rename(sit1, c(
    "GSub" = "V2",
    "Year" = "V3",
    "class" = "V5",
    "value" = "V18"
  ))
sit1r <- sit1r[, c("GSub", "Year", "class", "value")]

## filter required parameters
sit1_wid <-
  pivot_wider(sit1r, names_from = class, values_from = value)
sit1_s <- select(sit1_wid, GSub, Year, PRCP, IRGA, ET, WYLD)
fwrite(sit1_s, "sit1.csv")

## scientific notation is not allowed
op <- options(scipen = 99999)
sit_start <- fread("sit1.csv")

##calculate PRK
sit1_PRK <- mutate(sit_start, PRK = PRCP + IRGA - ET - WYLD)
sit1_PRK <- select(sit1_PRK, GSub, Year, PRK)
head(sit1_PRK)

## filter required subbasins
#sit1_PRK <- filter(sit1_PRK, GSub == 3 | GSub == 4| GSub == 5| GSub == 6| GSub == 8| GSub == 9| GSub == 10| GSub == 13| GSub == 15| GSub == 17| GSub == 18| GSub == 21)
#sit1_PRK<- filter(sit1_PRK, Year != 1990 & Year != 1991)

## read **gw.CSV from SWAT-GW folder
sit1_gw <- read.table('zygw.csv', header = T, sep = ',')
zy_merge <- merge(sit1_gw, sit1_PRK,
                  by = intersect(names(sit1_gw), names(sit1_PRK)))
head(zy_merge)

## calculate the variation of groundwater level(H)
zy_S <-
  mutate(zy_merge, S = (PRK + LARCHRG - DA_RCHG - SA_IRR - Qgw - REVAP))
zy_H <- mutate(zy_S, H = (S / u) / 1000)
head(zy_H)

## output calculation results
fwrite(zy_H, "zy_H.csv")
