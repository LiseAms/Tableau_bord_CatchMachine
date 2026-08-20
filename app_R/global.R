# =============================================================================
# global.R
#
# Exécuté une seule fois au démarrage de l'app (avant ui/server), et partagé
# entre toutes les sessions utilisateurs. On y charge :
#   - les packages
#   - la fonction clean_data()
#   - les référentiels fixes (taille/poids, Info_pecheur), qui ne sont
#     PAS uploadés par l'utilisateur mais packagés avec l'app.
#
# À ADAPTER : vérifie les chemins ci-dessous selon l'emplacement réel de tes
# fichiers de référence dans le projet.
# =============================================================================
cat(">>> global.R démarré\n")

# chargement des library 
library(shiny)
library(shinycssloaders)
library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(scales)
library(sf)
library(rnaturalearth)
library(knitr)
library(kableExtra)
library(bslib)
library(DT)
library(viridis)  
library(leaflet)

source("appR/clean_data.R")

# --- Référentiels fixes, chargés une seule fois ---
cat(">>> avant read_excel\n")                                                   # sert a verifier que ce code tourne bien quanq l'app est RUN

sp_taille_poids <- read_excel("tableau_taille_max_especes.xlsx")                 # excel essentiel a avoir dans le dossier de l'app. On le charge ici. COmprend les tailles et poids maximaux des especes. Permet d'enlever les grosses aberrations
cat(">>> tableau_taille_max_especes chargé\n")                                  # sert a verifier que ce code tourne bien quanq l'app est RUN et que le document taille_poids est bien chargé

tableau_a_b_especes <- read_csv("tableau_coefficient_a_b_especes.csv")          # excel essentiel a avoir dans le dossier de l'app. On le charge ici, il comprend les coeeficient a et b des 3 sources : Fishbase, SIH et Obsbio de l'Ifremer

banyuls_shp <- st_read("carto/shp_banyuls.shp")%>%                              # shapefile de banyuls
  mutate(Nom_SousZoneDePeche = trimws(sub(" -.*", "", Id_Zone_su))) %>%
  st_transform(crs = 4326)

  
print(names(banyuls_shp))
cat(">>> banyuls_shp chargé\n")


calanque_shp <- st_read("carto/shp_PNC.shp") %>%                                # shapefile du PNC
  mutate(Nom_SousZoneDePeche = zone) %>%
  st_set_crs(2154) %>%
  st_transform(crs = 4326)   # on converti tout de suite en WGS84 pour leaflet

# Rajouter ici de nouveaux shapefile de sous zone des AMP afin de pouvoir les intégrer dans le Tableau de bord. Il faut également rajouter une ligne de code dans le script app.R (partie cartographie)
