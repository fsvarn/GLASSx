###################################################
# Quick check to examine whether microglia are correlated with other signatures
# Updated: 2020.04.10
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(reshape)


##################################################
rm(list=ls())

# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in data
q <- 
"
SELECT * 
FROM analysis.davoli_immune_score im
UNION
SELECT aliquot_barcode, CONCAT('Muller ', signature_name)  AS signature_name, enrichment_score
FROM analysis.muller_tam_score tm 
"

dat <- dbGetQuery(con,q)

mat <- data.frame(pivot_wider(dat, names_from = aliquot_barcode, values_from = enrichment_score))
rownames(mat) <- mat[,1]
mat <- mat[,-1]
mat <- t(mat)

cor(mat,method="s")

#                        NK.cells Macrophages CD4.mature Muller Microglia
# NK.cells              1.0000000   0.3376630  0.5410595        0.2306348
# Macrophages           0.3376630   1.0000000  0.4785151        0.4977778
# CD4.mature            0.5410595   0.4785151  1.0000000        0.3575701
# Muller Microglia      0.2306348   0.4977778  0.3575701        1.0000000
# Muller Macrophages    0.5105471   0.8242213  0.6418709        0.4158560
# Macrophages.M2        0.3392596   0.9685218  0.4952228        0.4814991
# CD8.effector          0.7592786   0.4246606  0.4703657        0.2146966
# Dendritic             0.4695568   0.7255236  0.5395383        0.4571239
# B.cells               0.5361608   0.4031107  0.5385371        0.3866344
# T.reg                 0.4717043   0.5732111  0.4697823        0.1025572
# CD8.effector.NK.cells 0.7838420   0.4166863  0.5050089        0.2275690
# Macrophages.M1        0.4728293   0.9016750  0.5875026        0.4263778
#                       Muller Macrophages Macrophages.M2 CD8.effector Dendritic
# NK.cells                       0.5105471      0.3392596    0.7592786 0.4695568
# Macrophages                    0.8242213      0.9685218    0.4246606 0.7255236
# CD4.mature                     0.6418709      0.4952228    0.4703657 0.5395383
# Muller Microglia               0.4158560      0.4814991    0.2146966 0.4571239
# Muller Macrophages             1.0000000      0.8579351    0.6134387 0.8721788
# Macrophages.M2                 0.8579351      1.0000000    0.4328153 0.7480854
# CD8.effector                   0.6134387      0.4328153    1.0000000 0.5848674
# Dendritic                      0.8721788      0.7480854    0.5848674 1.0000000
# B.cells                        0.5342402      0.4148362    0.5473699 0.5534723
# T.reg                          0.6828898      0.6093651    0.5223386 0.6556318
# CD8.effector.NK.cells          0.6590326      0.4432382    0.8519542 0.5939385
# Macrophages.M1                 0.8621861      0.9003250    0.5214634 0.7242990
#                         B.cells     T.reg CD8.effector.NK.cells Macrophages.M1
# NK.cells              0.5361608 0.4717043             0.7838420      0.4728293
# Macrophages           0.4031107 0.5732111             0.4166863      0.9016750
# CD4.mature            0.5385371 0.4697823             0.5050089      0.5875026
# Muller Microglia      0.3866344 0.1025572             0.2275690      0.4263778
# Muller Macrophages    0.5342402 0.6828898             0.6590326      0.8621861
# Macrophages.M2        0.4148362 0.6093651             0.4432382      0.9003250
# CD8.effector          0.5473699 0.5223386             0.8519542      0.5214634
# Dendritic             0.5534723 0.6556318             0.5939385      0.7242990
# B.cells               1.0000000 0.4095615             0.5519157      0.4805586
# T.reg                 0.4095615 1.0000000             0.4863984      0.6373450
# CD8.effector.NK.cells 0.5519157 0.4863984             1.0000000      0.5620241
# Macrophages.M1        0.4805586 0.6373450             0.5620241      1.0000000