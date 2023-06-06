#Packages needed and data set-up

library(microbiome)
library(phyloseq)
library(tidyverse)
library(MetBrewer)
library(ade4)
library(vegan)

# Data

load("SM3-ps_Chiarello_et_al2020.rData") #Unmodified france data
ASV_florida_completo<-read.csv("~/Cosas de la maestría/PRIMER SEMESTRE/BIOGEOQUÍMICA Y ECOLOGÍA MICROBIANA MARINA/TRABAJO FINAL/ByEMM/florida/ASV_florida.csv") #Florida data 300+ samples