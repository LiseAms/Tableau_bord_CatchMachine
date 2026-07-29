# =============================================================================
# R/clean_data.R
#
# Reprend toute la logique du document clean_data.qmd, encapsulée dans une
# fonction réutilisable, pour pouvoir nettoyer un fichier brut CatchMachine
# uploadé dans l'app Shiny (au lieu de passer par un script Quarto exécuté à la
# main).
#
# =============================================================================

library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(readr)
library(shinycssloaders)

clean_data <- function(data_raw,
                        sp_taille_poids,
                        tableau_a_b_especes,
                       
                        date_min = "2025-01-01") {

  # ---------------------------------------------------------------------
  # 1) Renommage des colonnes et filtre de date
  # ---------------------------------------------------------------------
  pro_data <- data_raw %>%
    rename(Taille_MM = `Taille (mm)`, Poids_g = `Poids (g)`) %>%
    filter(`Heure debut session` >= as.POSIXct(paste(date_min, "00:00:00"))) # date min a changer si besoin 

  # ---------------------------------------------------------------------
  # 2) Formatage des colonnes et transformation des valeurs
  # ---------------------------------------------------------------------
  colnames(pro_data) <- gsub(" ", "_", colnames(pro_data))

  pro_data <- pro_data %>%
    filter(!is.na(Heure_debut_session)) %>%
    mutate(
      jour_semaine = wday(Heure_debut_session, label = TRUE, abbr = FALSE, week_start = 1), # code pour obtenir le jour de la semaine de la déclaration
      mois = month(Heure_debut_session, label = TRUE, abbr = FALSE), # code pour obtenir le mois de la déclaration
      annee = year(Heure_debut_session) # code pour obtenir l'année de la déclaration
    )

  pro_data <- pro_data %>%
    mutate(Durée_de_la_session = replace(Durée_de_la_session, Durée_de_la_session == 0, NA)) %>% # on corrige les durées de session qui sont à 0; ce n'est pas possible de sortir 0 minutes donc on met NA
    mutate(Nokill = ifelse(is.na(Nokill), 0, Nokill)) # la colonne Nokill contient des NA, on les remplace en 0

  pro_data <- pro_data %>%
    mutate(Id_Prise = as.numeric(str_remove_all(Id_Prise, "[^0-9]"))) # Pour passer les Id_Prise en numerique, il faut enlever les symboles présents avant les numéros, ce qui apparait des fois devant les entrées

  pro_data <- pro_data %>%
    mutate(Poids_g = abs(Poids_g)) # des poids sont parfois renseignés en nombre négatif; alors on met tous les poids en absolu pour s'affranchir de ces cas là

  # ---------------------------------------------------------------------
  # 3) Doublons
  # ---------------------------------------------------------------------
  pro_data <- pro_data %>%
    distinct(across(-c(Id_Prise, Horodate, Latitude, Longitude, Lien_photo)), .keep_all = TRUE) # doublon considéré quand tout est égal sauf les colonnes entre parenthèses

  # ---------------------------------------------------------------------
  # 4) Zones de pêche aberrantes (-1 -> 8)
  # ---------------------------------------------------------------------
  pro_data <- pro_data %>%
    mutate(modification_main = NA_character_) # on créé une colonne modification_main pour ajouter des corrections à la main 

  pro_data <- pro_data %>%
    mutate(
      
      # certaines zones de pêche sont incohérentes et ne correspondent à rien, donc on les enlève et on met dans Méditerranée (hors zone réglementée par CatchMachine). 
      modification_main = ifelse(ID_ZoneDePeche == -1,
                                  "Zone de pêche = -1 donc on reassigne en 8",
                                  NA_character_),
      ID_ZoneDePeche  = ifelse(ID_ZoneDePeche == -1, 8, ID_ZoneDePeche),
      Nom_ZoneDePeche = ifelse(is.na(Nom_ZoneDePeche),
                               "Méditerranée (hors zone réglementée par CatchMachine)",
                               Nom_ZoneDePeche)
    )

  # ---------------------------------------------------------------------
  # 5) Corrections manuelles connues
  #    (ATTENTION : ce bloc contient une correction ponctuelle codée en dur,
  #    identique au script Quarto original — Id_Prise == 41161. Si de
  #    nouvelles corrections manuelles ponctuelles sont identifiées, il faut
  #    les ajouter ici de la même façon, avec leur justification.)
  # ---------------------------------------------------------------------
  pro_data <- pro_data %>%
    mutate(
      Nb_de_prises = if_else(Id_Prise == 41161, 1, Nb_de_prises),
      modification_main = if_else(
        Id_Prise == 41161,
        if_else(
          is.na(modification_main),
          "modification nb_prise car photo montre 1 indiv",
          paste(modification_main, "modification nb_prise car photo montre 1 indiv", sep = " ; ")
        ),
        modification_main
      )
    )

  # ---------------------------------------------------------------------
  # 6) Taille / poids max par espèce
  # ---------------------------------------------------------------------
  
  # on utilise le dataframe chargé dans global.R pour les tailles et poids max des especes. On les integre ici au jeu de données que l'on utilise
  pro_data <- pro_data %>%
    left_join(sp_taille_poids %>%
                select(NomScientifique, Taille_max_mm, Poids_max_g) %>%
                distinct(),
              by = "NomScientifique") %>%
    mutate(Taille_max_mm = as.numeric(Taille_max_mm),
           Poids_max_g = as.numeric(Poids_max_g))


  # ---------------------------------------------------------------------
  # 7) Mode de pêche (bateau / bord / sous-marin)
  # ---------------------------------------------------------------------
  # il n'y a pas de de colonne mode de pêche dans le tableau de données, on le reconstruit dans la nouvelle colonne Mode_peche 
  pro_data <- pro_data %>%
    mutate(Mode_peche = NA_character_)

  # 1.  **Technique** : certains mots-clés dans la technique déclarée (sous-marin, indienne, agachon, coulée, trou) indiquent sans ambiguïté la pêche sous-marine.
  # 2.  **Zone/sous-zone de pêche** : quand la zone est renseignée, on peut déduire le mode de pêche directement (certaines sous-zones ne sont accessibles qu'en bateau ou ne sont que du bord. On suppose ici que les zones sont correctement déclarées par les utilisateurs.
  # 3.  **Technique seule** (repli) : si on n'a toujours pas de mode de pêche, on se base sur la technique déclarée seule. Certaines techniques restent ambiguës (ex. "Lancer" peut se faire du bord ou d'un bateau) -\> classées "Indeterminable" à ce stade.

    pro_data <- pro_data %>%
    mutate(
      Mode_peche = if_else(
        str_detect(
          Technique,
          regex("sous[- ]?marin|indienne|agachon|coulée|trou", ignore_case = TRUE)
        ),
        "Sous-marin",
        Mode_peche
      )) %>%
    mutate(
      Mode_peche = if_else(
        is.na(Mode_peche),
        case_when(
          ID_ZoneDePeche == 7 & ID_SousZoneDePeche %in% c(13, 14, 15) ~ "Bord",
          ID_ZoneDePeche == 7 & ID_SousZoneDePeche %in% c(16:21) ~ "Bateau",
          ID_ZoneDePeche == 4 & ID_SousZoneDePeche %in% c(35,33,34,39,52,53,54,55,56,83,58,59,60,61,62,63,64,65,66,67,68,81,70,71,72,73,74,75,76,77,78,79,82,84,24,25) ~ "Bateau", # sous zone de pêche, a alimenter si changements 
          TRUE ~ NA_character_
        ),
        Mode_peche
      )) %>%
    mutate(
      Mode_peche = if_else(
        is.na(Mode_peche),
        case_when(
          # changer/adapter ici les techniques de pêche si besoin
          Technique == "Dérive" ~ "Bateau",
          Technique == "Broumé" ~ "Indeterminable",
          Technique == "Calée verticale (à l'ancre)" ~ "Bateau",
          Technique == "Lancer" ~ "Indeterminable",
          Technique == "Pêche du bord (rockfishing)" ~ "Bord",
          Technique == "Surfcasting" ~ "Bord",
          Technique == "Traine" ~ "Bateau",
          Technique == "À la palangrotte" ~ "Indeterminable",
          TRUE ~ NA_character_
        ),
        Mode_peche
      ))

  
  # ---------------------------------------------------------------------
  # 8) Adhésion fédération
  # ---------------------------------------------------------------------
  pro_data <- pro_data %>%
    mutate(adhesion_fede = if_else(is.na(Nom_Fédération), "non", "oui"))

  # ---------------------------------------------------------------------
  # 9) Correction des valeurs aberrantes (nb de prises / nb de pêcheurs)
  #    NB : reprend à l'identique les règles du script Quarto d'origine,
  #    y compris le point à vérifier sur les blocs Cap Corse / Bonifacio
  #    (2e ligne testant ID_ZoneDePeche == 7 au lieu de 6 / 5).
  # ---------------------------------------------------------------------
  
  
  # --- Fonction utilitaire : ajoute un message dans Modification_main sans écraser l'existant ---
  ajouter_modif <- function(existant, nouveau, condition) {
    if_else(
      condition,
      if_else(is.na(existant) | existant == "",
              nouveau,
              paste(existant, nouveau, sep = " ; ")),
      existant
    )
  }

  # les especes concernées par le seuil haut du nombre de prise par déclaration. especes a changer ou adapter dans cette liste si besoin 
  especes_exclues <- c("Maquereau commun", "Sardine commune", "Maquereau espagnol",
                       "Maquereaux", "Girelle", "Girelle paon") 
  
  pro_data <- pro_data %>%
    mutate(flag_prises = FALSE, flag_pecheurs = FALSE) %>%
    
    # les corrections sont par zone, comme ça c'est adaptable si les besoins des gestionnaires changent, mais c'est similaire pour toutes les zones initialement
    # --- PNMGL (zone 1) ---
    # on met en place les conditions et modifications. 
    mutate(
      cond = ID_ZoneDePeche == 1 & Nb_de_prises > 150, #Ligne à flagger si elle appartient à la zone 1 et dépasse 150 prises
      Nb_de_prises = if_else(cond, 150, Nb_de_prises), #  Si la condition est vraie, on remplace Nb_de_prises par 150, sinon on garde la valeur initiale
      flag_prises = flag_prises | cond,  # On ajoute cette condition au drapeau existant (pour garder en mémoire les changements)
      cond = ID_ZoneDePeche == 1 & !NomVernaculaire %in% especes_exclues & Nb_de_prises > 100, # zone 1, espèce non exclue, et plus de 100 prises
      Nb_de_prises = if_else(cond, 1, Nb_de_prises), # alors on remplace par 1; sinon on garde la valeur initiale
      flag_prises = flag_prises | cond # Mise à jour du drapeau 
    ) %>%
    
    # --- Banyuls (zone 7) ---
    mutate(
      cond = ID_ZoneDePeche == 7 & Nb_de_prises > 150,
      Nb_de_prises = if_else(cond, 150, Nb_de_prises),
      flag_prises = flag_prises | cond,
      cond = ID_ZoneDePeche == 7 & !NomVernaculaire %in% especes_exclues & Nb_de_prises > 100,
      Nb_de_prises = if_else(cond, 1, Nb_de_prises),
      flag_prises = flag_prises | cond
    ) %>%
    
    # --- Calanques (zone 4) ---
    mutate(
      cond = ID_ZoneDePeche == 4 & Nb_de_prises > 150,
      Nb_de_prises = if_else(cond, 150, Nb_de_prises),
      flag_prises = flag_prises | cond,
      cond = ID_ZoneDePeche == 4 & !NomVernaculaire %in% especes_exclues & Nb_de_prises > 100,
      Nb_de_prises = if_else(cond, 1, Nb_de_prises),
      flag_prises = flag_prises | cond
    ) %>%
    
    # --- Cap Corse (zone 6) ---
    mutate(
      cond = ID_ZoneDePeche == 6 & Nb_de_prises > 150,
      Nb_de_prises = if_else(cond, 150, Nb_de_prises),
      flag_prises = flag_prises | cond,
      cond = ID_ZoneDePeche == 6 & !NomVernaculaire %in% especes_exclues & Nb_de_prises > 100,
      Nb_de_prises = if_else(cond, 1, Nb_de_prises),
      flag_prises = flag_prises | cond
    ) %>%
    
    # --- Bonifacio (zone 5) ---
    mutate(
      cond = ID_ZoneDePeche == 5 & Nb_de_prises > 150,
      Nb_de_prises = if_else(cond, 150, Nb_de_prises),
      flag_prises = flag_prises | cond,
      cond = ID_ZoneDePeche == 5 & !NomVernaculaire %in% especes_exclues & Nb_de_prises > 100,
      Nb_de_prises = if_else(cond, 1, Nb_de_prises),
      flag_prises = flag_prises | cond
    ) %>%
    
    # --- Hors AMP (zone 8) ---
    mutate(
      cond = ID_ZoneDePeche == 8 & Nb_de_prises > 150,
      Nb_de_prises = if_else(cond, 150, Nb_de_prises),
      flag_prises = flag_prises | cond,
      cond = ID_ZoneDePeche == 8 & !NomVernaculaire %in% especes_exclues & Nb_de_prises > 100,
      Nb_de_prises = if_else(cond, 1, Nb_de_prises),
      flag_prises = flag_prises | cond
    ) %>%
    
    # --- Correction Nb_conservées / Nb_relâchées ---
    # on a les endroits ou il faut changer (cond et flag), donc on peut opérer aux modifications
    mutate(
      cond = Nb_relâchées == 0 & Nb_conservées > 150, # Si aucune remise à l'eau n'est renseignée et que le nombre conservé dépasse 150, on plafonne Nb_conservées à 150
      Nb_conservées = if_else(cond, 150, Nb_conservées),
      flag_prises = flag_prises | cond,
      cond = Nb_conservées == 0 & Nb_relâchées > 150, # Si aucune conservation n'est renseignée et que le nombre relâché dépasse 150, on plafonne Nb_relâchées à 150
      Nb_relâchées = if_else(cond, 150, Nb_relâchées),
      flag_prises = flag_prises | cond
    ) %>%
    mutate(
      cond = Nb_relâchées == 0 & !NomVernaculaire %in% especes_exclues & Nb_conservées > 100, # Pour les espèces non exclues, si aucune remise à l'eau n'est renseignée , et que le nombre conservé dépasse 100, on remplace par 1
      Nb_conservées = if_else(cond, 1, Nb_conservées),
      flag_prises = flag_prises | cond,
      cond = Nb_conservées == 0 & !NomVernaculaire %in% especes_exclues & Nb_relâchées > 100, # Pour les espèces non exclues, si aucune conservation n'est renseignée, et que le nombre relâché dépasse 100, on remplace par 1
      Nb_relâchées = if_else(cond, 1, Nb_relâchées),
      flag_prises = flag_prises | cond
    ) %>%
    
    # --- Nombre de pêcheurs ---
    mutate(
      cond = ID_ZoneDePeche == 4 & Mode_peche == "Bateau" & Nb_de_pêcheurs > 12, # Zone 4, pêche en bateau : si le nombre de pêcheurs dépasse 12, on remplace par 1
      Nb_de_pêcheurs = if_else(cond, 1, Nb_de_pêcheurs),
      flag_pecheurs = flag_pecheurs | cond,
      cond = ID_ZoneDePeche == 7 & Mode_peche == "Bateau" & Nb_de_pêcheurs > 8, # Zone 7, pêche en bateau : seuil de 8
      Nb_de_pêcheurs = if_else(cond, 1, Nb_de_pêcheurs),
      flag_pecheurs = flag_pecheurs | cond,
      cond = ID_ZoneDePeche == 1 & Mode_peche == "Bateau" & Nb_de_pêcheurs > 12,  # Zone 1, pêche en bateau : seuil de 12
      Nb_de_pêcheurs = if_else(cond, 1, Nb_de_pêcheurs),
      flag_pecheurs = flag_pecheurs | cond
    ) %>%
    mutate(
      cond = Nb_de_pêcheurs > 30 | Nb_de_pêcheurs <= 0, # Si le nombre de pêcheurs est supérieur à 30 ou inférieur ou égal à 0, on le remplace par 1
      Nb_de_pêcheurs = if_else(cond, 1, Nb_de_pêcheurs),
      flag_pecheurs = flag_pecheurs | cond
    ) %>%
    
    # --- Écriture finale dans Modification_main ---
    mutate(
      modification_main = ajouter_modif(modification_main, "Quantité de prise renseignée aberrante", flag_prises),
      modification_main = ajouter_modif(modification_main, "nombre de pêcheur renseigné aberrant", flag_pecheurs)
    ) %>%
    select(-cond, -flag_prises, -flag_pecheurs)

  # ---------------------------------------------------------------------
  # 10) Identifiant de sortie reconstruit
  # ---------------------------------------------------------------------
  pro_data <- pro_data %>%
    # CatchMachine n'indique pas correctement les sorties : on considère donc qu'il s'agit de la même sortie quand on retrouve le même pêcheur, la même date/heure de début de session, le même mode de pêche et la même zone de pêche.
    mutate(
      ID_sortie_corrige = paste(Id_Abonné, Heure_debut_session, Mode_peche, Nom_ZoneDePeche, sep = "_"),
      ID_sortie_corrige = sprintf("S%06d", as.integer(factor(ID_sortie_corrige)))
    )

  pro_data
}
