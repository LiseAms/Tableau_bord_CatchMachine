#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
source("global.R", local = FALSE)
options(shiny.maxRequestSize = 50 * 1024^2)
library(bslib)
library(DT)

# global.R doit charger : les packages, clean_data(), sp_taille_poids,
# ET pour la correction des poids : tableau_a_b_especes


# =============================================================================
# Fonction de correction des poids (nécessité d'avoir le tableau des coefficients a et b : tableau_a_b_especes)
# =============================================================================
appliquer_correction_poids <- function(data_raw, tableau_a_b_especes) {
  
  data_raw %>%
    left_join(tableau_a_b_especes, by = "NomScientifique") %>%
    mutate(
      Poids_theorique_SIH    = a_SIH * Taille_MM ^ b_SIH, # calcul des poids théoriques a partir des coefficients a et b du SIH
      Poids_theorique_FB     = a_FB * (Taille_MM / 10) ^ b_FB, # calcul des poids théoriques a partir des coefficients a et b de FishBase
      Poids_theorique_ObsBio = a_obs * (Taille_MM / 10) ^ b_obs, # calcul des poids théoriques a partir des coefficients a et b de Obsbio (IFremer)
      
      nb_theo_non_na = rowSums(!is.na(across(c(Poids_theorique_SIH, Poids_theorique_FB, Poids_theorique_ObsBio)))), # compte combien de poids théoriques sont disponibles (= combien de références)
      nb_med         = rowSums(across(c(origine_SIH, origine_FB, origine_obs)) == "med", na.rm = TRUE), # Combien de ces coefficients sont calculés à partir de données de Méditerranée
      cas_equivalent_unique = nb_theo_non_na == 1 | (nb_theo_non_na > 1 & nb_med == 1), # Note si il y n'y a qu'une seule référence théorique ou une seule de Méditerranée
      
      # condition de validité ObsBio (robuste ET taille dans la gamme), appliquée partout de façon cohérente, y compris cas_equivalent_unique
      obs_valide = !is.na(Poids_theorique_ObsBio) &
        Taille_MM >= Taille_min_ratio & Taille_MM <= Taille_max_ratio & ratio_robuste == TRUE,  # code pour valider si ObsBio est valide = True -> utilisable
      
      # on selectionne un seul poids théorique par ligne . 
      poids_theo_ref = case_when(
        cas_equivalent_unique & origine_SIH == "med"                            ~ Poids_theorique_SIH, # poids théorique = SIH car seul dispo et de Med
        cas_equivalent_unique & origine_FB  == "med"                            ~ Poids_theorique_FB,  # poids théorique = Fishbase car seul dispo et de Med
        cas_equivalent_unique & origine_obs == "med" & obs_valide               ~ Poids_theorique_ObsBio, # poids théorique = Obsbio car seul dispo et de Med
        cas_equivalent_unique & obs_valide                                      ~ Poids_theorique_ObsBio, # poids théorique = Obsbio car seul dispo 
        cas_equivalent_unique & !is.na(Poids_theorique_FB)                      ~ Poids_theorique_FB, # poids théorique = Fishbase car seul dispo
        cas_equivalent_unique & !is.na(Poids_theorique_SIH)                     ~ Poids_theorique_SIH, # poids théorique = SIH car seul dispo
        !cas_equivalent_unique & origine_obs == "med" & obs_valide              ~ Poids_theorique_ObsBio, # poids théorique = Obsbio car vient de Med 
        !cas_equivalent_unique & origine_FB == "med"                            ~ Poids_theorique_FB, # poids théorique = FB car vient de Med et Obsbio non/pas dispo
        !cas_equivalent_unique & origine_SIH == "med"                           ~ Poids_theorique_SIH, # poids théorique = SIH car car vient de Med et Obsbio et FB non/pas dispo
        !cas_equivalent_unique & !is.na(Poids_theorique_ObsBio) & obs_valide    ~ Poids_theorique_ObsBio, # poids théorique = Obsbio car ordre de priorité 
        !cas_equivalent_unique & !is.na(Poids_theorique_FB)                     ~ Poids_theorique_FB, # poids théorique = FB car ordre de priorité 
        !cas_equivalent_unique & !is.na(Poids_theorique_SIH)                    ~ Poids_theorique_SIH # poids théorique = SIH car ordre de priorité 
      ),
      
      ratio_si_total    = ifelse(Taille_MM < Taille_max_mm * 1.1, (Poids_g / Nb_de_prises) / poids_theo_ref, NA_real_), # on regarde si c'est un poids total qui est renseigné
      ratio_si_moyenne  = ifelse(Taille_MM < Taille_max_mm * 1.1, Poids_g / poids_theo_ref, NA_real_), # si c'est une déclaration d'un individu moyen du groupe 
      ratio_si_total_kg = ifelse(Taille_MM < Taille_max_mm * 1.1, (Poids_g * 1000 / Nb_de_prises) / poids_theo_ref, NA_real_), # même chose mais on regarde s'il y a une erreur d'unité
      ratio_si_moy_kg   = ifelse(Taille_MM < Taille_max_mm * 1.1, (Poids_g * 1000) / poids_theo_ref, NA_real_), # même chose mais on regarde s'il y a une erreur d'unité
      
       # diagnostique de poids
      diag_poids = case_when(
        is.na(poids_theo_ref)                                     ~ "aucun poids théorique, poids obs gardé",
        Taille_MM > Taille_max_mm * 1.2 & !is.na(Taille_max_mm)   ~ "taille aberrante car > à taille max +20%, poids obs gardé",
        Nb_de_prises == 0                                          ~ "pas de capture",
        Taille_MM == 0                                             ~ "pas de taille, poids obs gardé",
        Poids_g == 0                                               ~ "poids absent",
        Nb_de_prises == 1 & between(ratio_si_moyenne, 0.5, 1.5)   ~ "poids obs ok",
        Nb_de_prises == 1 & between(ratio_si_moy_kg, 0.5, 1.5)    ~ "erreur unité : poids saisi en kg (individu)",
        Nb_de_prises == 1                                          ~ "poids obs aberrant",
        Nb_de_prises > 1 & between(ratio_si_total, 0.5, 1.5)      ~ "TOTAL GARDÉ",
        Nb_de_prises > 1 & between(ratio_si_moyenne, 0.5, 1.5)    ~ "poids donné moyen d'un individu",
        Nb_de_prises > 1 & between(ratio_si_total_kg, 0.5, 1.5)   ~ "TOTAL GARDÉ - erreur unité : poids saisi en kg",
        Nb_de_prises > 1 & between(ratio_si_moy_kg, 0.5, 1.5)     ~ "erreur unité : poids saisi en kg (individu)",
        TRUE                                                        ~ "poids aberrant"
      ),
      
      # calcul du poids en fonction du diagnostique
      poids_final = case_when(
        is.na(poids_theo_ref)                                     ~ Poids_g,
        Taille_MM > Taille_max_mm * 1.1 & !is.na(Taille_max_mm)   ~ Poids_g,
        Nb_de_prises == 0                                          ~ NA_real_,
        Taille_MM == 0                                             ~ Poids_g,
        Poids_g == 0                                               ~ poids_theo_ref,
        Nb_de_prises == 1 & between(ratio_si_moyenne, 0.5, 1.5)   ~ Poids_g,
        Nb_de_prises == 1 & between(ratio_si_moy_kg, 0.5, 1.5)    ~ Poids_g * 1000,
        Nb_de_prises == 1                                          ~ poids_theo_ref,
        Nb_de_prises > 1 & between(ratio_si_total, 0.5, 1.5)      ~ Poids_g,
        Nb_de_prises > 1 & between(ratio_si_moyenne, 0.5, 1.5)    ~ Poids_g,
        Nb_de_prises > 1 & between(ratio_si_total_kg, 0.5, 1.5)   ~ Poids_g * 1000,
        Nb_de_prises > 1 & between(ratio_si_moy_kg, 0.5, 1.5)     ~ Poids_g * 1000,
        TRUE                                                        ~ poids_theo_ref
      )
    ) %>%
    # verification finale
    mutate(
      diag_poids  = if_else(poids_final > 1000000, "taille et poids aberrant, on met NA", diag_poids),
      poids_final = if_else(poids_final > 1000000, NA_real_, poids_final)
    ) %>%
    mutate(
      poids_final_total = case_when(
        grepl("TOTAL GARDÉ", diag_poids)                            ~ poids_final,
        is.na(poids_final)                                           ~ NA_real_,
        grepl("aucun poids théorique, poids obs gardé", diag_poids) ~ poids_final,
        grepl("pas de taille, poids obs gardé", diag_poids)         ~ poids_final,
        TRUE                                                          ~ poids_final * Nb_de_prises
      )
    )
}
# =============================================================================
# UI
# =============================================================================
# le UI c'est l'interface, ce que l'utilisateur.ice du Tableau de bord va manipuler. C'est donc les affichages, les menu déroulants pour selection de certaines données, thème, design... 

ui <- page_navbar(
  title = "Tableau de bord - CatchMachine",
  theme = bs_theme(version = 5),
  
  header = tags$style(HTML("
    .dashboard-card { border: 1px solid #444; border-radius: 4px; margin-bottom: 15px; }
    .dashboard-card .card-header { font-weight: 600; background: #fff; border-bottom: 1px solid #444; }
    .card-header-flex { display:flex; justify-content:space-between; align-items:center; }
    .year-select { max-width: 160px; }
  ")),
  
  nav_panel("Informations générales",
            layout_columns(
              col_widths = c(3, 5, 4), # largeur des colonnes pour la page 1, le total doit faire 12 
              fill = FALSE,
              
              # ================= COLONNE GAUCHE =================
              
              #chaque case est codée de la même façon : 
              
              # div(                            
              #   class = "d-flex flex-column",   # format de la case, ici on veut des bordures 
              #   
              #   card(class = "dashboard-card",  # format case avec bordures
              #        card_header("TITRE"),      # ici on met le titre de la case , c'est ce qui s'affiche en haut de la case concernée 
              #
              #                     #la suite du code concerne le contenu de la case et est spécifique à chaque type d'affichage, ici les explications : 
              #                                   # fileInput = permet de loader un fichier ; checkboxInput = permet de cocher une case ; selectInput = permet de selectionner un élément dans un menu déroulant qui va determiner ce qui est traité dans la case
              #                                   # withSpinner() = affiche le logo de chargement le temps que les codes tournent . 
              #                                   # dans withSpinner() il y a plusieurs arguments différents, selon le type d'élément affiché (uiOutput = texte ; DTOutput = tableau déroulant ; plotOutput = graphique ; tableOutput = tableau simple ;)
              #   ),
                
              div(
                class = "d-flex flex-column",
                
                card(class = "dashboard-card",
                     card_header("Choisir le CSV brut CatchMachine"),
                     fileInput("file", NULL, accept = ".csv", buttonLabel = "Browse...", placeholder = "No file selected") # permet de créer la case pour l'upload du csv
                ),
                
                card(class = "dashboard-card",
                     card_header("Options"),
                     checkboxInput("utiliser_correction_poids",
                                   "Utiliser la correction des poids (CPUE / Biomasse)", # case pour utiliser le script de correction des poids ou non
                                   value = FALSE),
                     p(class = "text-muted small mb-0",
                       "Coché : CPUE et Biomasse utilisent les poids corrigés (ObsBio / FishBase / SIH). ", # texte qui s'affiche pour expliquer la case
                       "Décoché : elles utilisent le poids déclaré brut tel quel.")
                ),
                
                card(class = "dashboard-card",
                     card_header(
                       div(class = "card-header-flex",
                           span("Informations générale"), # case avec Des informations générales, le titre de la case peut être changé ici
                           selectInput("annee_info", NULL, choices = "Toutes les années", width = "160px") # permet de choisir les années, soit toutes soit une par une, commence à 2025, les nouvelles années s'accumulent automatiquement
                       )
                     ),
                     withSpinner(uiOutput("resume_texte"), type = 4, color = "steelblue")
                ),
                
                card(class = "dashboard-card",
                     card_header("Espèces"),
                     p("Liste des espèces capturées", class = "text-muted small mb-1"), # Case pour la liste des espèces pêchées
                     withSpinner(DTOutput("table_especes"), type = 4, color = "steelblue")
                )
              ),
              
              # ================= COLONNE CENTRALE =================
              div(
                class = "d-flex flex-column",
                
                card(class = "dashboard-card",
                     card_header(
                       div(class = "card-header-flex",
                           span("Evolution temporelle annuelle des captures"), # Case pour l'évolution des captures
                           selectInput("annee_evolution", NULL, choices = "Toutes les années", width = "160px")
                       )
                     ),
                     withSpinner(plotOutput("plot_evolution", height = "500px"), type = 4, color = "steelblue") # ggplot de l'évolution, la taille peut être changée dans le height =
                ),
                
                card(class = "dashboard-card", # case pour les CPUE : contient 3 onglets 
                     card_header("CPUE"),
                     tabsetPanel(
                       id = "onglets_cpue",
                       
                       # ---- Onglet 1 : CPUE générales (mode + zone) ----
                       tabPanel("CPUE générales",
                                br(),
                                withSpinner(plotOutput("plot_cpue_nb_mode", height = "320px"), type = 4, color = "steelblue"), # si le nom du plot change dans le SERVER, changer ici, taille du plot modifiable dans le height
                                withSpinner(plotOutput("plot_cpue_poids_mode", height = "320px"), type = 4, color = "steelblue") # si le nom du plot change dans le SERVER, changer ici, taille du plot modifiable dans le height
                       ),
                       
                       
                       # ---- Onglet 2 : CPUE des 10 espèces les plus pêchées, par mode ----
                       tabPanel("CPUE espèces",
                                br(),
                                selectInput("mode_peche_cpue", "Mode de pêche", choices = NULL, width = "220px"), # menu déroulant pour les 3 modes de pêche
                                withSpinner(plotOutput("plot_cpue_espece_nb", height = "380px"), type = 4, color = "steelblue"), # si le nom du plot change dans le SERVER, changer ici, taille du plot modifiable dans le height
                                withSpinner(plotOutput("plot_cpue_espece_poids", height = "380px"), type = 4, color = "steelblue") # si le nom du plot change dans le SERVER, changer ici, taille du plot modifiable dans le height
                       ),
                       
                       # ---- Onglet 3 : Tableaux résumés CPUE espèces ----
                       tabPanel("Tableaux résumés",
                                br(),
                                selectInput("mode_peche_cpue_tableau", "Mode de pêche", choices = NULL, width = "220px"), # menu déroulant pour les 3 modes de pêche
                                h5("CPUE (nb/pêcheur/sortie) par espèce"),
                                withSpinner(tableOutput("table_cpue_nb"), type = 4, color = "steelblue"),
                                hr(),
                                h5("CPUE (g/pêcheur/sortie) par espèce"),
                                withSpinner(tableOutput("table_cpue_poids"), type = 4, color = "steelblue")
                       )
                     )
                )
              ),
              
              
              # ================= COLONNE DROITE =================
              div(
                class = "d-flex flex-column",
                
                card(class = "dashboard-card",
                     card_header(
                       div(class = "card-header-flex",
                           span("Top 5 des espèces les plus capturées"), 
                           selectInput("annee_top10", NULL, choices = "Toutes les années", width = "160px") # menu déroulant pour choisir l'année (ou toutes)
                       )
                     ),
                     withSpinner(plotOutput("plot_top10", height = "280px"), type = 4, color = "steelblue") # ggplot des 5 especes les plus capturées en nombre
                ),
                
                card(class = "dashboard-card",
                     card_header(
                       div(class = "card-header-flex",
                           span("Sorties infructueuses"),
                           selectInput("annee_infructueuses", NULL, choices = "Toutes les années", width = "160px") # menu déroulant pour choisir l'année (ou toutes)
                       )
                     ),
                     withSpinner(plotOutput("plot_infructueuses", height = "280px"), type = 4, color = "steelblue")
                ),
                
                card(class = "dashboard-card",
                     card_header(
                       div(class = "card-header-flex",
                           span("Biomasse"),
                           selectInput("annee_biomasse", NULL, choices = "Toutes les années", width = "160px")
                       )
                     ),
                     withSpinner(uiOutput("biomasse_texte"), type = 4, color = "steelblue") 
                )
              )
            )
  ),
  
  nav_panel("Informations complémentaires",  # onglet/page 2
            layout_columns(
              col_widths = c(6,6),
              fill = FALSE,
              
              # ================= COLONNE GAUCHE =================
              div(
                class = "d-flex flex-column",
                
                card(class = "dashboard-card",
                     card_header("Distribution des tailles par espèce"),
                     selectInput("espece_taille", "Choisir une espèce", # menu déroulant pour choisir l'espèce, possibilité de commencer à écrire pour choisir l'espèce
                                 choices = NULL,
                                 selected = NULL,
                                 width = "300px"),
                     withSpinner(plotOutput("plot_taille_espece", height = "450px"), type = 4, color = "steelblue")
                )
              ), 
              # ================= COLONNE CENTRALE =================
              div(
                class = "d-flex flex-column", 
              
                card(class = "dashboard-card",
                     card_header("Relation taille-poids (RTP) par espèce"), 
                     selectInput("espece_rtp", "Choisir une espèce",  # menu déroulant pour choisir l'espèce, possibilité de commencer à écrire pour choisir l'espèce
                                 choices = NULL, 
                                 width = "300px"), 
                     withSpinner(plotOutput("plot_rtp_espece", height = "450px"), type = 4, color = "steelblue")
                )
              ), 
                
              
              
           )
          
), 
  nav_panel("Cartographie", # onglet 3/page 3 : disponible pour les AMP avec des sous zones de pêche 
            layout_columns(
              fill = FALSE,
              
             
              div(
                class = "d-flex flex-column",
                
                card(class = "dashboard-card",
                     card_header("Cartographie des prises par espèce dans la réserve"),
                     selectInput("espece_cartographie", "Choisir une espèce", choices = NULL, width = "300px"),
                     uiOutput("cartographie_ui")
                ),
                
                card(class = "dashboard-card",
                     card_header("Cartographie des prises toutes espèces confondues"),
                     uiOutput("cartographie_globale_ui")
                )
              )
              )

))


# =============================================================================
# SERVER
# =============================================================================
# contenu du Tableau de Bord : partie qui code concrêtement ce qu'il y a dans le TB et ce qui est appelé par le UI pour être ensuite affiché

server <- function(input, output, session) {
  
  # -------------------------------------------------------------------
  # Lecture + nettoyage :
  # -------------------------------------------------------------------
  raw_data <- reactive({                                                        # reactive veut dire qu'on peut intéragir avec ce qui est uploadé. Ici le fichier mis dans la case correspondante
    req(input$file)
    read_delim(input$file$datapath,
               delim = ";", escape_double = FALSE, trim_ws = TRUE,              # format de lecture du csv
               skip = 2)                                                        # Si le format de csv CatchMachine change et qu'il n'y a plus 2 lignes de texte avant le réel tableau de données, alors il faut changer dans le skip=2 et mettre la nouvelle valeur
  })
  
  df <- reactive({
    req(raw_data())
    
    # verifier que les données sont dans le bon format et fonctionne bien. Si ça ne fonctionne pas, un message d'erreur s'affiche si le jeu de données n'est pas conforme avec le code
    validate(
      need(
        all(c("Heure debut session", "Taille (mm)", "Poids (g)") %in% names(raw_data())),
        "Le fichier importé ne ressemble pas à un export brut CatchMachine (colonnes attendues manquantes). Vérifiez que le CSV brut est importé, et non un fichier déjà nettoyé."
      )
    )
    
    withProgress(message = "Nettoyage des données en cours...", value = 0.3, {  # affiche un message si c'est long à charger pour que les utilisateurs sachent que ça fonctionne et que ça ne beug pas. 
      result <- tryCatch({
        clean_data(                                                             # appelle le script clean_data qui est fourni dans le dossier 
          data_raw          = raw_data(),                                       # récupère les données catchmachine chargées précedemment. 
          sp_taille_poids   = sp_taille_poids,                                  # récupère le tableau de données des tailles et poids max des espèces. Ce fichier doit être dans le dossier avec l'application pour que cela fonctionne
          date_min          = "2025-01-01"                                      # à changer ici si la date de début de l'analyse des données est à changer
        )
      }, error = function(e) {
        validate(paste("Erreur lors du nettoyage des données :", conditionMessage(e))) 
      })
      incProgress(0.7, detail = "Terminé")
      result
    })
  })
  
  # -------------------------------------------------------------------
  # Poids utilisé (brut ou corrigé) -> colonne commune poids_final_total
  # pour que CPUE et Biomasse s'écrivent une seule fois, quel que soit l'état de la case à cocher.
  # -------------------------------------------------------------------
  data_avec_poids <- reactive({              
    req(df())
    
    if (isTRUE(input$utiliser_correction_poids)) {                              # si la case de correction des poids est cochée dans l'appli, alors cette partie est utilisée
      validate(
        need(exists("tableau_a_b_especes"),                                     # vérifie que le tableau avec les coefficients a et b est bien présent dans le dossier
             "Le référentiel tableau_a_b_especes n'est pas chargé dans global.R : la correction des poids est indisponible.")
      )
      appliquer_correction_poids(df(), tableau_a_b_especes)                     # applique la fonction appliquer_correction_poids définie au début de ce script
    } else {
      # Pas de correction : le poids déclaré brut est utilisé tel quel comme "total" de la ligne
      df() %>% mutate(poids_final_total = Poids_g)
    }
  })
  
  # -------------------------------------------------------------------
  # Menus déroulants : Rajouter en dessous les nouveaux menus si besoin

  
    # Années disponibles -> alimente automatiquement TOUS les menus déroulants d'année dès qu'un nouveau CSV est chargé
  observeEvent(df(), {
    annees <- sort(unique(df()$annee))                                          # fait la liste de toutes les années disponibles dans les données 
    choix  <- c("Toutes les années", as.character(annees))
    
    updateSelectInput(session, "annee_info", choices = choix, selected = "Toutes les années")
    updateSelectInput(session, "annee_evolution", choices = choix, selected = "Toutes les années")
    updateSelectInput(session, "annee_top10", choices = choix, selected = "Toutes les années")
    updateSelectInput(session, "annee_infructueuses", choices = choix, selected = "Toutes les années")
    updateSelectInput(session, "annee_biomasse", choices = choix, selected = "Toutes les années")
    
    
    # Modes de pêche disponibles (hors "Indeterminable") -> alimente le
    # menu de l'onglet CPUE espèces ET celui de l'onglet Tableaux résumés
    modes <- sort(unique(df()$Mode_peche)) # sort les modes de pêche des données
    modes <- modes[modes != "Indeterminable"] # on enleve le mode indterminable
    updateSelectInput(session, "mode_peche_cpue", choices = modes, selected = modes[1])
    updateSelectInput(session, "mode_peche_cpue_tableau", choices = modes, selected = modes[1])
    
    # espèces disponibles pour les histogrammes de classe de taille
    especes_dispo <- sort(unique(df()$NomVernaculaire))
    updateSelectInput(session, "espece_taille", choices = especes_dispo, selected = especes_dispo[1])
    
    # RTP : especes disponibles pour faire le RTP
    especes_rtp <- sort(unique(df()$NomVernaculaire))
    updateSelectInput(session, "espece_rtp", choices = especes_rtp, selected = especes_rtp[1])
    
   
  })
  
  # Petite fonction utilitaire de filtrage par année, réutilisée partout
  filtrer_annee <- function(data, annee_choisie) {
    if (is.null(annee_choisie) || annee_choisie == "Toutes les années") {
      return(data)
    }
    dplyr::filter(data, as.character(annee) == annee_choisie)
  }
  
  # -------------------------------------------------------------------
  # Calcule la limite haute d'un axe Y à partir des VRAIES moustaches du
  # boxplot (boxplot.stats), groupe par groupe, plutôt qu'à partir d'un
  # simple quantile des données brutes : ça évite de couper une boîte ou
  # sa moustache quand les données sont étalées ou peu nombreuses.
  # -------------------------------------------------------------------
  calc_ymax_boxplot <- function(valeurs, groupes, marge = 1.05) {
    d <- data.frame(valeur = valeurs, groupe = groupes) %>%
      filter(!is.na(valeur))
    
    if (nrow(d) == 0) return(NA_real_)
    
    moustaches_sup <- d %>%
      group_by(groupe) %>%
      summarise(sup = boxplot.stats(valeur)$stats[5], .groups = "drop") %>%
      pull(sup)
    
    max(moustaches_sup, na.rm = TRUE) * marge
  }
  
  # =====================================================================
  # - Résumé (nb pêcheurs, sorties, captures) - filtré par année
  # =====================================================================
  df_info <- reactive({
    req(df())
    filtrer_annee(df(), input$annee_info)
  })
  
  output$resume_texte <- renderUI({
    req(df_info())
    
    nb_pecheurs     <- n_distinct(df_info()$Id_Abonné)
    nb_sorties      <- n_distinct(df_info()$ID_sortie_corrige)
    capture_totale  <- sum(df_info()$Nb_de_prises, na.rm = TRUE)
    date_min        <- min(df_info()$Heure_debut_session, na.rm = TRUE)
    date_max        <- max(df_info()$Heure_debut_session, na.rm = TRUE)
    
    tagList(
      p(strong(format(nb_pecheurs, big.mark = " ")), " pêcheurs déclarants sur la période sélectionnée."),
      p(strong(format(nb_sorties, big.mark = " ")), " sorties de pêche sur la période sélectionnée."),
      p(strong(format(capture_totale, big.mark = " ")), " captures déclarées sur la période sélectionnée.")
    )
  })
  
  # =====================================================================
  # d - Espèces les plus pêchées -> liste déroulante (DT), sans sélection
  # multiple : c'est un affichage, pas un outil de sélection de lignes.
  # =====================================================================
  liste_especes <- reactive({
    req(df())
    df() %>%
      group_by(NomVernaculaire) %>%
      summarise(nb_captures = sum(Nb_de_prises, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(nb_captures)) %>%
      mutate(pourcentage = round(100 * nb_captures / sum(nb_captures), 2)) %>%
      rename("Espèce" = NomVernaculaire, "Nb captures" = nb_captures, "% du total" = pourcentage)
  })
  
  output$table_especes <- renderDT({
    datatable(
      liste_especes(),
      rownames = FALSE,
      selection = "none",   # <- retire la possibilité de sélectionner des lignes (mono ou multiple)
      options = list(
        dom = "t",
        scrollY = "553px",
        scrollCollapse = TRUE,
        paging = FALSE,
        ordering = FALSE
      )
    )
  })
  
  # =====================================================================
  # e - top 5 des espèces pêchées, filtrable par année
  # =====================================================================
  df_top10 <- reactive({
    req(df())
    filtrer_annee(df(), input$annee_top10)
  })
  
  espece_top10 <- reactive({
    req(df_top10())
    
    validate(
      need(all(c("NomScientifique", "Nb_de_prises", "CodeFAO", "NomVernaculaire") %in% names(df_top10())),
           "Colonnes attendues manquantes après nettoyage (NomScientifique, Nb_de_prises, CodeFAO, NomVernaculaire)")
    )
    
    top5 <- df_top10() %>%
      group_by(NomScientifique) %>%
      summarise(total = sum(Nb_de_prises, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total)) %>%
      slice_head(n = 5)
    
    top5 %>%
      left_join(
        df_top10() %>% select(CodeFAO, NomVernaculaire, NomScientifique) %>% distinct(),
        by = "NomScientifique"
      ) %>%
      arrange(desc(total))
  })
  
  output$plot_top10 <- renderPlot({
    req(espece_top10())
    ggplot(espece_top10(), aes(x = total, y = reorder(NomVernaculaire, total))) +
      geom_col(fill = "skyblue1") +
      geom_text(aes(label = total), hjust = -0.2, size = 6) +
      labs(title = "Top 5 des espèces les plus capturées", x = "Nombre de prises", y = "Espèce") +
      theme_minimal() +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10),
                         expand = expansion(mult = c(0, 0.15)))+
      theme(
        text = element_text(size = 16), # taille de base (axes, legendes, etc.)
        plot.title = element_text(face = "bold", size = 16), # titre principal
        axis.title = element_text(size = 14), # titres des axes
        axis.text = element_text(size = 13), # textes des graduations
        axis.text.y = element_text(size = 14), # labels y (noms d'espèces)
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13)
      )
  })
  
  # =====================================================================
  # f - Évolution temporelle des captures - filtrée par année si choisie
  # =====================================================================
  df_evolution <- reactive({
    req(df())
    filtrer_annee(df(), input$annee_evolution)
  })
  
  prise_mois <- reactive({
    req(df_evolution())
    
    ordre_mois <- c("janvier", "février", "mars", "avril", "mai", "juin",
                    "juillet", "août", "septembre", "octobre", "novembre", "décembre")
    
    df_evolution() %>%
      mutate(mois = factor(mois, levels = ordre_mois)) %>%
      group_by(annee, mois) %>%
      summarise(nb_prises = sum(Nb_de_prises, na.rm = TRUE), .groups = "drop") %>%
      arrange(annee, mois)
  })
  
  output$plot_evolution <- renderPlot({
    req(prise_mois())
    
    ordre_mois <- c("janvier", "février", "mars", "avril", "mai", "juin",
                    "juillet", "août", "septembre", "octobre", "novembre", "décembre")
    
    annee_courante <- year(Sys.Date())
    mois_courant   <- month(Sys.Date())
    
    ggplot(prise_mois(), aes(x = mois, y = nb_prises,
                             color = factor(annee), group = annee)) +
      
      {if (annee_courante %in% prise_mois()$annee)
        annotate("rect", xmin = mois_courant, xmax = Inf, ymin = 0, ymax = Inf,
                 fill = "grey90", alpha = 0.4)} +
      
      geom_line(linewidth = 1.4) +
      geom_point(size = 3.5, shape = 21, fill = "white", stroke = 1.5) +
      geom_label(aes(label = format(nb_prises, big.mark = "\u202f")),
                 vjust = -0.6, size = 4, fontface = "bold",
                 label.padding = unit(0.15, "lines"),
                 label.size = 0, fill = "white", alpha = 0.75,
                 show.legend = FALSE) +
      
      scale_color_brewer(palette = "Set2") +
      scale_y_continuous(labels = label_number(big.mark = "\u202f"),
                         expand = expansion(mult = c(0, 0.15))) +
      scale_x_discrete(limits = ordre_mois) +
      
      labs(title = "Nombre de prises de pêche par mois",
           subtitle = paste0("Total : ", format(sum(df_evolution()$Nb_de_prises, na.rm = TRUE), big.mark = "\u202f"), " prises"),
           x = NULL, y = "Nombre de prises", color = "Année",
           caption = "La zone grisée (si présente) indique les mois de l'année en cours.") +
      
      theme_minimal(base_size = 12) +
      theme(
        plot.title    = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40", size = 14),
        plot.caption  = element_text(color = "grey40", size = 12, hjust = 0),
        axis.text.x   = element_text(angle = 30, hjust = 1, size = 15),
        axis.text.y   = element_text(size = 15),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 12),
        legend.position  = "right"
      )
  })
  
  # =====================================================================
  # Sorties infructueuses -> camembert, filtrable par année
  # =====================================================================
  df_infructueuses <- reactive({
    req(df())
    filtrer_annee(df(), input$annee_infructueuses)
  })
  
  sorties_bredouilles <- reactive({
    req(df_infructueuses())
    
    df_infructueuses() %>%
      group_by(Nom_ZoneDePeche, ID_sortie_corrige) %>%
      summarise(nb_captures = sum(Nb_de_prises, na.rm = TRUE), .groups = "drop") %>%
      group_by(Nom_ZoneDePeche) %>%
      summarise(
        nb_sorties     = n_distinct(ID_sortie_corrige),
        nb_bredouilles = sum(nb_captures == 0),
        .groups = "drop"
      ) %>%
      mutate(nb_avec_captures = nb_sorties - nb_bredouilles)
  })
  
  camembert_df <- reactive({
    req(sorties_bredouilles())
    
    sorties_bredouilles() %>%
      summarise(
        nb_bredouilles    = sum(nb_bredouilles),
        nb_avec_captures  = sum(nb_avec_captures)
      ) %>%
      select(nb_bredouilles, nb_avec_captures) %>%
      tidyr::pivot_longer(
        cols = everything(),
        names_to = "statut",
        values_to = "n"
      ) %>%
      mutate(pct = n / sum(n) * 100)
  })
  
  output$plot_infructueuses <- renderPlot({
    req(camembert_df())
    
    ggplot(camembert_df(), aes(x = "", y = n, fill = statut)) +
      geom_col(width = 1, color = "white", linewidth = 0.8) +
      geom_text(
        aes(label = n),
        position = position_stack(vjust = 0.5),
        color = "black",
        fontface = "bold",
        size = 5
      ) +
      coord_polar(theta = "y") +
      scale_fill_manual(
        name = "Légende",
        values = c(
          "nb_avec_captures" = "palegreen",
          "nb_bredouilles" = "#8f99fb"
        ),
        labels = c(
          "nb_avec_captures" = "Nombre de sorties avec capture",
          "nb_bredouilles" = "Nombre de sorties infructueuse"
        )
      ) +
      labs(title = "Sorties infructueuses") +
      theme_void(base_size = 16) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)
      )
  })
  
  # =====================================================================
  # ONGLET 1 : CPUE générales (mode + zone), sur poids brut ou corrigé
  # selon la case à cocher, via data_avec_poids(). Échelle Y dynamique :
  # aucune limite fixe n'est posée, elle s'adapte aux données.
  # =====================================================================
  cpue_data <- reactive({
    req(data_avec_poids())
    
    data_avec_poids() %>%
      group_by(ID_sortie_corrige, Mode_peche) %>%
      summarise(
        nb_captures    = sum(Nb_de_prises, na.rm = TRUE),
        poids_captures = sum(poids_final_total, na.rm = TRUE),
        nb_pecheurs    = first(Nb_de_pêcheurs),
        zone           = first(Nom_ZoneDePeche),
        date           = first(Heure_debut_session),
        mois           = first(mois),
        .groups = "drop"
      ) %>%
      mutate(
        cpue_sortie_nb    = nb_captures / nb_pecheurs,     # ind/pêcheur/sortie
        cpue_sortie_poids = poids_captures / nb_pecheurs   # g/pêcheur/sortie
      ) %>%
      filter(Mode_peche != "Indeterminable")
  })
  
  output$plot_cpue_nb_mode <- renderPlot({
    req(cpue_data())
    d <- cpue_data()
    ymax <- calc_ymax_boxplot(d$cpue_sortie_nb, d$Mode_peche)
    ggplot(cpue_data(), aes(x = Mode_peche, y = cpue_sortie_nb, fill = Mode_peche)) +
      geom_boxplot() +
      scale_y_continuous(
        breaks = scales::pretty_breaks(n = 7),
        labels = scales::label_number(),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(ylim = c(0, ymax)) + 
      labs(title = "CPUE (nb/pêcheur/sortie) par plateforme de pêche",
           x = "Plateforme de pêche", y = "CPUE (ind/pêcheur/sortie)") +
      theme_bw() +
      theme(
        legend.position = "none",
        plot.title    = element_text(face = "bold", size = 16),
        axis.text.x   = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.y   = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 12)
        )
  })
  
  output$plot_cpue_poids_mode <- renderPlot({
    req(cpue_data())
    d <- cpue_data()
    ymax <- calc_ymax_boxplot(d$cpue_sortie_poids, d$Mode_peche)
    ggplot(cpue_data(), aes(x = Mode_peche, y = cpue_sortie_poids, fill = Mode_peche)) +
      geom_boxplot() +
      scale_y_continuous(
        breaks = scales::pretty_breaks(n = 7),
        labels = scales::label_number(),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(ylim = c(0, ymax)) + 
      labs(title = "CPUE (g/pêcheur/sortie) par plateforme de pêche",
           x = "Plateforme de pêche", y = "CPUE (g/pêcheur/sortie)") +
      theme_bw() +
      theme(
        legend.position = "none",
        plot.title    = element_text(face = "bold", size = 16),
        axis.text.x   = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.y   = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 12)
      )
  })
  
  
  
  
  # =====================================================================
  # ONGLET 2 : CPUE des 10 espèces les plus pêchées, pour le mode de
  # pêche choisi dans le menu déroulant. Le top 10 est calculé une seule
  # fois sur l'ensemble des données (comme dans le script d'origine),
  # puis on affiche la CPUE de ces 10 espèces pour le mode sélectionné.
  # Échelle Y dynamique : pas de limite fixe, elle s'adapte aux données.
  # =====================================================================
  cpue_espece_top10 <- reactive({
    req(data_avec_poids())
    
    data_avec_poids() %>%
      group_by(NomScientifique, NomVernaculaire) %>%
      summarise(total = sum(Nb_de_prises, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total)) %>%
      slice_head(n = 10)
  })
  
  ordre_especes_cpue <- reactive({
    req(cpue_espece_top10())
    cpue_espece_top10() %>%
      arrange(desc(total)) %>%
      pull(NomVernaculaire)
  })
  
  cpue_espece_mode <- reactive({
    req(data_avec_poids(), cpue_espece_top10(), input$mode_peche_cpue)
    
    data_avec_poids() %>%
      filter(Mode_peche == input$mode_peche_cpue) %>%
      semi_join(cpue_espece_top10(), by = "NomScientifique") %>%
      group_by(ID_sortie_corrige, NomScientifique, NomVernaculaire) %>%
      summarise(
        nb_captures    = sum(Nb_de_prises, na.rm = TRUE),
        poids_captures = sum(poids_final_total, na.rm = TRUE),
        nb_pecheurs    = first(Nb_de_pêcheurs),
        zone           = first(Nom_ZoneDePeche),
        date           = first(Heure_debut_session),
        mois           = first(mois),
        .groups = "drop"
      ) %>%
      mutate(
        cpue_sortie_nb    = nb_captures / nb_pecheurs,
        cpue_sortie_poids = poids_captures / nb_pecheurs,
        NomVernaculaire   = factor(NomVernaculaire, levels = ordre_especes_cpue())
      )
  })
  
  output$plot_cpue_espece_nb <- renderPlot({
    req(cpue_espece_mode())
    d <- cpue_espece_mode()
    ymax <- calc_ymax_boxplot(d$cpue_sortie_nb, d$NomVernaculaire)
    
    ggplot(cpue_espece_mode(), aes(x = NomVernaculaire, y = cpue_sortie_nb, fill = NomVernaculaire)) +
      geom_boxplot(outlier.alpha = 0.4) +
      scale_y_continuous(
        breaks = scales::pretty_breaks(n = 7),
        labels = scales::label_number(),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(ylim = c(0, ymax)) +
      
      labs(
        title = paste0("CPUE des 10 espèces les plus capturées - ", input$mode_peche_cpue),
        x = "Espèce", y = "CPUE (ind/pêcheur/sortie)"
      ) +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        plot.title    = element_text(face = "bold", size = 16),
        axis.title.x = element_text(size = 15),
        axis.text.y   = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 12)
        )
      
  })
  
  output$plot_cpue_espece_poids <- renderPlot({
    req(cpue_espece_mode())
    d <- cpue_espece_mode()
    ymax <- calc_ymax_boxplot(d$cpue_sortie_poids, d$NomVernaculaire)
    ggplot(cpue_espece_mode(), aes(x = NomVernaculaire, y = cpue_sortie_poids, fill = NomVernaculaire)) +
      geom_boxplot(outlier.alpha = 0.4) +
      scale_y_continuous(
        breaks = scales::pretty_breaks(n = 7),
        labels = scales::label_number(),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(ylim = c(0, ymax)) +
      
      labs(
        title = paste0("CPUE des 10 espèces les plus capturées - ", input$mode_peche_cpue),
        x = "Espèce", y = "CPUE (g/pêcheur/sortie)"
      ) +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        plot.title    = element_text(face = "bold", size = 16),
        axis.title.x = element_text(size = 15),
        axis.text.y   = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 12)
      )
  })
  
  # =====================================================================
  # ONGLET 3 : Tableaux résumés CPUE des 10 espèces les plus pêchées,
  # pour le mode de pêche choisi dans SON PROPRE menu déroulant
  # (indépendant de celui de l'onglet 2, mêmes 10 espèces cependant)
  # =====================================================================
  cpue_espece_mode_tableau <- reactive({
    req(data_avec_poids(), cpue_espece_top10(), input$mode_peche_cpue_tableau)
    
    data_avec_poids() %>%
      filter(Mode_peche == input$mode_peche_cpue_tableau) %>%
      semi_join(cpue_espece_top10(), by = "NomScientifique") %>%
      group_by(ID_sortie_corrige, NomScientifique, NomVernaculaire) %>%
      summarise(
        nb_captures    = sum(Nb_de_prises, na.rm = TRUE),
        poids_captures = sum(poids_final_total, na.rm = TRUE),
        nb_pecheurs    = first(Nb_de_pêcheurs),
        zone           = first(Nom_ZoneDePeche),
        date           = first(Heure_debut_session),
        mois           = first(mois),
        .groups = "drop"
      ) %>%
      mutate(
        cpue_sortie_nb    = nb_captures / nb_pecheurs,
        cpue_sortie_poids = poids_captures / nb_pecheurs,
        NomVernaculaire   = factor(NomVernaculaire, levels = ordre_especes_cpue())
      )
  })
  
  tableau_cpue_nb <- reactive({
    req(cpue_espece_mode_tableau())
    
    cpue_espece_mode_tableau() %>%
      group_by(NomVernaculaire, NomScientifique) %>%
      summarise(
        n_sorties = n(),
        moyenne   = mean(cpue_sortie_nb, na.rm = TRUE),
        mediane   = median(cpue_sortie_nb, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(NomVernaculaire = factor(NomVernaculaire, levels = ordre_especes_cpue())) %>%
      arrange(NomVernaculaire) %>%
      mutate(across(where(is.numeric), ~ round(.x, 2)))
  })
  
  tableau_cpue_poids <- reactive({
    req(cpue_espece_mode_tableau())
    
    cpue_espece_mode_tableau() %>%
      group_by(NomVernaculaire, NomScientifique) %>%
      summarise(
        n_sorties = n(),
        moyenne   = mean(cpue_sortie_poids, na.rm = TRUE),
        mediane   = median(cpue_sortie_poids, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(NomVernaculaire = factor(NomVernaculaire, levels = ordre_especes_cpue())) %>%
      arrange(NomVernaculaire) %>%
      mutate(across(where(is.numeric), ~ round(.x, 2)))
  })
  
  output$table_cpue_nb <- renderTable({
    req(tableau_cpue_nb())
    tableau_cpue_nb() %>%
      rename("Nom vernaculaire" = NomVernaculaire, "Nom scientifique" = NomScientifique,
             "Nb sorties" = n_sorties, "Moyenne" = moyenne, "Médiane" = mediane)
  })
  
  output$table_cpue_poids <- renderTable({
    req(tableau_cpue_poids())
    tableau_cpue_poids() %>%
      rename("Nom vernaculaire" = NomVernaculaire, "Nom scientifique" = NomScientifique,
             "Nb sorties" = n_sorties, "Moyenne" = moyenne, "Médiane" = mediane)
  })
  
  # =====================================================================
  # Biomasse totale / conservée / relâchée, sur poids brut ou corrigé,
  # filtrée par année si choisie
  # =====================================================================
  df_biomasse <- reactive({
    req(data_avec_poids())
    filtrer_annee(data_avec_poids(), input$annee_biomasse)
  })
  
  biomasse_totale <- reactive({
    req(df_biomasse())
    
    df_biomasse() %>%
      mutate(
        poids_conserve = if_else(Nb_de_prises > 0,
                                 poids_final_total * (Nb_conservées / Nb_de_prises), 0),
        poids_relache  = if_else(Nb_de_prises > 0,
                                 poids_final_total * (Nb_relâchées / Nb_de_prises), 0)
      ) %>%
      summarise(
        biomasse_totale_kg    = round(sum(poids_final_total, na.rm = TRUE) / 1000, 2),
        biomasse_conservee_kg = round(sum(poids_conserve, na.rm = TRUE) / 1000, 2),
        biomasse_relachee_kg  = round(sum(poids_relache, na.rm = TRUE) / 1000, 2)
      )
  })
  
  output$biomasse_texte <- renderUI({
    req(biomasse_totale())
    b <- biomasse_totale()
    
    tagList(
      p(strong("Totale : "), format(b$biomasse_totale_kg, big.mark = " "), " kg"),
      p(strong("Conservée : "), format(b$biomasse_conservee_kg, big.mark = " "), " kg"),
      p(strong("Relachée : "), format(b$biomasse_relachee_kg, big.mark = " "), " kg")
    )
  })
  
  
  # =====================================================================
  # Classe de taille
  # =====================================================================
  
  data_taille_espece <- reactive({
    req(df(), input$espece_taille)
    
    df() %>%
      filter(NomVernaculaire == input$espece_taille) %>%
      uncount(Nb_de_prises) %>%
      filter(!is.na(Taille_MM))
  })
  
  output$plot_taille_espece <- renderPlot({
    req(data_taille_espece())
    validate(need(nrow(data_taille_espece()) > 0, "Aucune donnée de taille pour cette espèce."))
    
    q95  <- quantile(data_taille_espece()$Taille_MM, probs = 0.998, na.rm = TRUE)
    xmax <- ceiling(q95 / 10) * 10
    
    ggplot(data_taille_espece(), aes(x = Taille_MM)) +
      geom_histogram(
        binwidth = 10,
        fill = "lightblue",
        color = "black",
        boundary = 0
      ) +
      scale_x_continuous(
        limits = c(0.5, xmax),
        breaks = scales::pretty_breaks(n = 10),
        expand = expansion(mult = c(0, 0.05))
      ) +
      labs(
        title = paste0("Distribution des tailles - ", input$espece_taille),
        x = "Taille (mm)", y = "Effectif"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 15),
        plot.title    = element_text(face = "bold", size = 16),
        axis.title.x = element_text(size = 15),
        axis.text.y   = element_text(size = 15),
        axis.title.y = element_text(size = 15)
      )
  })
  
 
  ##### RTP (FishBase uniquement, avec diagnostic par point) #####
  
  # -------------------------------------------------------------------
  # Référentiel FishBase calculé UNE SEULE FOIS par jeu de données chargé
  # (appel réseau via rfishbase::length_weight, pas à chaque changement
  # d'espèce dans le menu déroulant)
  # -------------------------------------------------------------------
  rtp_fishbase_ref <- reactive({
    req(df())
    
    tableau_a_b_especes %>%
      filter(NomScientifique %in% unique(df()$NomScientifique)) %>%
      dplyr::select(
        NomScientifique,
        a        = a_FB,
        b        = b_FB,
        origine  = origine_FB
      ) %>%
      distinct(NomScientifique, .keep_all = TRUE)
  })
  
  # -------------------------------------------------------------------
  # Données de l'espèce choisie : nuage de points diagnostiqués + courbe
  # théorique taille-poids (uniquement les prises à 1 seul individu, avec
  # tolérance de 20% autour des tailles/poids max connus)
  # -------------------------------------------------------------------
  rtp_espece_data <- reactive({
    req(df(), input$espece_rtp, rtp_fishbase_ref())
    
    data_coherente_FB <- df() %>%
      filter(NomVernaculaire == input$espece_rtp) %>%
      filter(Nb_de_prises == 1) %>%
      filter(
        Taille_MM <= Taille_max_mm * 1.2 | is.na(Taille_max_mm),
        Poids_g   <= Poids_max_g  * 1.2 | is.na(Poids_max_g)
      ) %>%
      mutate(Taille_MM = Taille_MM / 10) %>%
      rename(Taille_CM = Taille_MM) %>%
      mutate(outlier_valeur = Taille_CM == 0 | Poids_g == 0) %>%
      left_join(
        rtp_fishbase_ref() %>% dplyr::select(NomScientifique, a, b, origine),
        by = "NomScientifique"
      ) %>%
      mutate(
        Poids_RTP  = a * Taille_CM ^ b,
        Poids_diff = abs(Poids_g - Poids_RTP),
        outlier_RTP = Poids_diff / Poids_RTP > 1,
        outlier = outlier_valeur | outlier_RTP,
        diag = ifelse(outlier_valeur, "Val. manquante", "val. ok"),
        diag = ifelse(!outlier_valeur, paste0(diag, " / ",
                                              ifelse(is.na(outlier_RTP), "RTP inconnue",
                                                     ifelse(outlier_RTP, "RTP aberrante", "RTP ok"))), diag)
      )
    
    if (nrow(data_coherente_FB) == 0 || all(is.na(data_coherente_FB$a))) {
      return(list(disponible = FALSE))
    }
    
    courbe <- data_coherente_FB %>%
      filter(!outlier_valeur) %>%
      summarise(
        Taille_max_obs = max(Taille_CM, na.rm = TRUE),
        a = first(a), b = first(b), origine = first(origine)
      ) %>%
      reframe(Taille = seq(from = 0, to = Taille_max_obs, by = 0.5), a = a, b = b, origine = origine) %>%
      mutate(Poids = a * Taille ^ b)
    
    list(disponible = TRUE, points = data_coherente_FB, courbe = courbe)
  })
  
  output$plot_rtp_espece <- renderPlot({
    d <- rtp_espece_data()
    req(d)
    
    if (!isTRUE(d$disponible)) {
      ggplot() +
        annotate("text", x = 0, y = 0, label = "Pas de données disponibles pour cette espèce.", size = 6) +
        theme_void()
    } else {
      origine_sp <- unique(d$courbe$origine)
      subtitle_text <- paste0("origine des données : ", ifelse(is.na(origine_sp), "NA", origine_sp))
      
      ggplot() +
        geom_line(data = d$courbe, aes(x = Taille, y = Poids), color = "brown", linewidth = 1) +
        geom_point(data = d$points, aes(x = Taille_CM, y = Poids_g, color = diag)) +
        labs(
          title = input$espece_rtp,
          subtitle = subtitle_text,
          x = "Taille (cm)", y = "Poids (g)"
        ) +
        scale_color_manual(
          name = NULL,
          values = c(
            "val. ok / RTP ok"        = "green",
            "val. ok / RTP inconnue"  = "darkgreen",
            "val. ok / RTP aberrante" = "orange",
            "Val. manquante"          = "red"
          )
        ) +
        theme_bw() +
        theme(legend.position = "bottom",
              axis.text.x = element_text(size = 15),
              plot.title    = element_text(face = "bold", size = 16),
              axis.title.x = element_text(size = 15),
              axis.text.y   = element_text(size = 15),
              axis.title.y = element_text(size = 15), 
              legend.text = element_text(size = 14)
        )
      
    }
  })

  ############### CARTOGRAPHIE ###############
  
    # --- Détermine quel shapefile utiliser selon les zones présentes dans les données ---
    shp_actif <- reactive({
      req(df())
      
      zones_presentes <- unique(df()$ID_ZoneDePeche)
      
      if ( 4 %in% zones_presentes) {
        list(shp = calanque_shp , disponible = TRUE)
      } else if ( 7 %in% zones_presentes) {
        list(shp = banyuls_shp , disponible = TRUE)
      } else {
        list(shp = NULL, disponible = FALSE)
      }
    })
    
    # --- Mise à jour dynamique des choix d'espèces (uniquement noms vernaculaires) ---
    observe({
      especes_dispo <- sort(unique(df()$NomVernaculaire))
      
      updateSelectInput(
        session,
        "espece_cartographie",
        choices = especes_dispo,
        selected = especes_dispo[1]
      )
    })
    
    # ============================================================
    # CARTE 1 : par espèce sélectionnée
    # ============================================================
    
    shp_pression_espece <- reactive({
      req(input$espece_cartographie)
      req(shp_actif()$disponible)
      
      data_sp <- df() %>%
        filter(NomVernaculaire == input$espece_cartographie)
      
      captures_par_zone <- data_sp %>%
        group_by(Nom_SousZoneDePeche) %>%
        summarise(
          total_nb_captures = sum(Nb_de_prises, na.rm = TRUE),
          nb_sorties        = n(),
          .groups = "drop"
        )
      
      shp_actif()$shp %>%
        left_join(captures_par_zone, by = "Nom_SousZoneDePeche") %>%
        mutate(total_nb_captures = replace(total_nb_captures, is.na(total_nb_captures), 0))
    })
    
    output$cartographie_ui <- renderUI({
      if (isTRUE(shp_actif()$disponible)) {
        withSpinner(leafletOutput("cartographie", height = "450px"), type = 4, color = "steelblue")
      } else {
        div(
          class = "alert alert-warning",
          "Pas de sous-zone disponible dans les données, la cartographie n'est pas faisable."
        )
      }
    })
    
    output$cartographie <- renderLeaflet({
      req(shp_actif()$disponible)
      shp_data <- shp_pression_espece()
      
      pal <- colorNumeric(
        palette = "inferno",
        domain  = shp_data$total_nb_captures,
        reverse = TRUE
      )
      
      leaflet(shp_data) %>%
        addProviderTiles(providers$OpenStreetMap) %>%
        addPolygons(
          fillColor    = ~pal(total_nb_captures),
          fillOpacity  = 0.7,
          color        = "white",
          weight       = 1.5,
          smoothFactor = 0.5,
          popup = ~paste0(
            "<b>Zone : </b>", Nom_SousZoneDePeche, "<br>",
            "<b>Captures : </b>", total_nb_captures, " prises<br>",
            "<b>Sorties : </b>", nb_sorties
          ),
          highlight = highlightOptions(
            weight      = 3,
            color       = "#FFFFFF",
            fillOpacity = 0.9,
            bringToFront = TRUE
          ),
          label = ~paste0(Nom_SousZoneDePeche, " — ", total_nb_captures, " prises"),
          labelOptions = labelOptions(
            style     = list("font-weight" = "bold", padding = "4px 8px"),
            textsize  = "13px",
            direction = "auto"
          )
        ) %>%
        addLegend(
          position = "bottomright",
          pal      = pal,
          values   = ~total_nb_captures,
          title    = paste0("Nombre de prises<br><small>", input$espece_cartographie, "</small>"),
          opacity  = 0.8
        )
    })
    
    # ============================================================
    # CARTE 2 : toutes espèces confondues
    # ============================================================
    
    shp_pression_globale <- reactive({
      req(shp_actif()$disponible)
      
      captures_par_zone <- df() %>%
        group_by(Nom_SousZoneDePeche) %>%
        summarise(
          total_nb_captures = sum(Nb_de_prises, na.rm = TRUE),
          nb_sorties        = n(),
          .groups = "drop"
        )
      
      shp_actif()$shp %>%
        left_join(captures_par_zone, by = "Nom_SousZoneDePeche") %>%
        mutate(total_nb_captures = replace(total_nb_captures, is.na(total_nb_captures), 0))
    })
    
    output$cartographie_globale_ui <- renderUI({
      if (isTRUE(shp_actif()$disponible)) {
        withSpinner(leafletOutput("cartographie_globale", height = "450px"), type = 4, color = "steelblue")
      } else {
        div(
          class = "alert alert-warning",
          "Pas de sous-zone disponible dans les données, la cartographie n'est pas faisable."
        )
      }
    })
    
    output$cartographie_globale <- renderLeaflet({
      req(shp_actif()$disponible)
      shp_data <- shp_pression_globale()
      
      pal <- colorNumeric(
        palette = "inferno",
        domain  = shp_data$total_nb_captures,
        reverse = TRUE
      )
      
      leaflet(shp_data) %>%
        addProviderTiles(providers$OpenStreetMap) %>%
        addPolygons(
          fillColor    = ~pal(total_nb_captures),
          fillOpacity  = 0.7,
          color        = "white",
          weight       = 1.5,
          smoothFactor = 0.5,
          popup = ~paste0(
            "<b>Zone : </b>", Nom_SousZoneDePeche, "<br>",
            "<b>Captures : </b>", total_nb_captures, " prises<br>",
            "<b>Sorties : </b>", nb_sorties
          ),
          highlight = highlightOptions(
            weight      = 3,
            color       = "#FFFFFF",
            fillOpacity = 0.9,
            bringToFront = TRUE
          ),
          label = ~paste0(Nom_SousZoneDePeche, " — ", total_nb_captures, " prises"),
          labelOptions = labelOptions(
            style     = list("font-weight" = "bold", padding = "4px 8px"),
            textsize  = "13px",
            direction = "auto"
          )
        ) %>%
        addLegend(
          position = "bottomright",
          pal      = pal,
          values   = ~total_nb_captures,
          title    = "Nombre de prises<br><small>Toutes espèces</small>",
          opacity  = 0.8
        )
    })
  }

  

# Run the application
shinyApp(ui = ui, server = server)
