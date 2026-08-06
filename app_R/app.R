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
      Poids_theorique_SIH    = a_SIH * Taille_MM ^ b_SIH,                       # calcul des poids théoriques a partir des coefficients a et b du SIH
      Poids_theorique_FB     = a_FB * (Taille_MM / 10) ^ b_FB,                  # calcul des poids théoriques a partir des coefficients a et b de FishBase
      Poids_theorique_ObsBio = a_obs * (Taille_MM / 10) ^ b_obs,                # calcul des poids théoriques a partir des coefficients a et b de Obsbio (IFremer)
      
      nb_theo_non_na = rowSums(!is.na(across(c(Poids_theorique_SIH, Poids_theorique_FB, Poids_theorique_ObsBio)))),    # compte combien de poids théoriques sont disponibles (= combien de références)
      nb_med         = rowSums(across(c(origine_SIH, origine_FB, origine_obs)) == "med", na.rm = TRUE),                # Combien de ces coefficients sont calculés à partir de données de Méditerranée
      cas_equivalent_unique = nb_theo_non_na == 1 | (nb_theo_non_na > 1 & nb_med == 1),                                # Note si il y n'y a qu'une seule référence théorique ou une seule de Méditerranée
      
      # condition de validité ObsBio (robuste ET taille dans la gamme), appliquée partout de façon cohérente, y compris cas_equivalent_unique
      obs_valide = !is.na(Poids_theorique_ObsBio) &
        Taille_MM >= Taille_min_ratio & Taille_MM <= Taille_max_ratio & ratio_robuste == TRUE,  # code pour valider si ObsBio est valide = True -> utilisable
      
      # on selectionne un seul poids théorique par ligne . 
      poids_theo_ref = case_when(
        cas_equivalent_unique & origine_SIH == "med"                            ~ Poids_theorique_SIH,                         # poids théorique = SIH car seul dispo et de Med
        cas_equivalent_unique & origine_FB  == "med"                            ~ Poids_theorique_FB,                          # poids théorique = Fishbase car seul dispo et de Med
        cas_equivalent_unique & origine_obs == "med" & obs_valide               ~ Poids_theorique_ObsBio,                      # poids théorique = Obsbio car seul dispo et de Med
        cas_equivalent_unique & obs_valide                                      ~ Poids_theorique_ObsBio,                      # poids théorique = Obsbio car seul dispo 
        cas_equivalent_unique & !is.na(Poids_theorique_FB)                      ~ Poids_theorique_FB,                          # poids théorique = Fishbase car seul dispo
        cas_equivalent_unique & !is.na(Poids_theorique_SIH)                     ~ Poids_theorique_SIH,                         # poids théorique = SIH car seul dispo
        !cas_equivalent_unique & origine_obs == "med" & obs_valide              ~ Poids_theorique_ObsBio,                      # poids théorique = Obsbio car vient de Med 
        !cas_equivalent_unique & origine_FB == "med"                            ~ Poids_theorique_FB,                          # poids théorique = FB car vient de Med et Obsbio non/pas dispo
        !cas_equivalent_unique & origine_SIH == "med"                           ~ Poids_theorique_SIH,                         # poids théorique = SIH car car vient de Med et Obsbio et FB non/pas dispo
        !cas_equivalent_unique & !is.na(Poids_theorique_ObsBio) & obs_valide    ~ Poids_theorique_ObsBio,                      # poids théorique = Obsbio car ordre de priorité 
        !cas_equivalent_unique & !is.na(Poids_theorique_FB)                     ~ Poids_theorique_FB,                          # poids théorique = FB car ordre de priorité 
        !cas_equivalent_unique & !is.na(Poids_theorique_SIH)                    ~ Poids_theorique_SIH                          # poids théorique = SIH car ordre de priorité 
      ),
      
      ratio_si_total    = ifelse(Taille_MM < Taille_max_mm * 1.1, (Poids_g / Nb_de_prises) / poids_theo_ref, NA_real_),        # on regarde si c'est un poids total qui est renseigné
      ratio_si_moyenne  = ifelse(Taille_MM < Taille_max_mm * 1.1, Poids_g / poids_theo_ref, NA_real_),                         # si c'est une déclaration d'un individu moyen du groupe 
      ratio_si_total_kg = ifelse(Taille_MM < Taille_max_mm * 1.1, (Poids_g * 1000 / Nb_de_prises) / poids_theo_ref, NA_real_), # même chose mais on regarde s'il y a une erreur d'unité
      ratio_si_moy_kg   = ifelse(Taille_MM < Taille_max_mm * 1.1, (Poids_g * 1000) / poids_theo_ref, NA_real_),                # même chose mais on regarde s'il y a une erreur d'unité
      
      # diagnostique de poids
      diag_poids = case_when(
        is.na(poids_theo_ref)                                                   ~ "aucun poids théorique, poids obs gardé",
        Taille_MM > Taille_max_mm * 1.2 & !is.na(Taille_max_mm)                 ~ "taille aberrante car > à taille max +20%, poids obs gardé",
        Nb_de_prises == 0                                                       ~ "pas de capture",
        Taille_MM == 0                                                          ~ "pas de taille, poids obs gardé",
        Poids_g == 0                                                            ~ "poids absent",
        Nb_de_prises == 1 & between(ratio_si_moyenne, 0.5, 1.5)                 ~ "poids obs ok",
        Nb_de_prises == 1 & between(ratio_si_moy_kg, 0.5, 1.5)                  ~ "erreur unité : poids saisi en kg (individu)",
        Nb_de_prises == 1                                                       ~ "poids obs aberrant",
        Nb_de_prises > 1 & between(ratio_si_total, 0.5, 1.5)                    ~ "TOTAL GARDÉ",
        Nb_de_prises > 1 & between(ratio_si_moyenne, 0.5, 1.5)                  ~ "poids donné moyen d'un individu",
        Nb_de_prises > 1 & between(ratio_si_total_kg, 0.5, 1.5)                 ~ "TOTAL GARDÉ - erreur unité : poids saisi en kg",
        Nb_de_prises > 1 & between(ratio_si_moy_kg, 0.5, 1.5)                   ~ "erreur unité : poids saisi en kg (individu)",
        TRUE                                                                    ~ "poids aberrant"
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


mois_ordre <- c("janvier","février","mars","avril","mai","juin",
                "juillet","août","septembre","octobre","novembre","décembre")


# =============================================================================
# UI
# =============================================================================
# le UI c'est l'interface, ce que l'utilisateur.ice du Tableau de bord va manipuler. C'est donc les affichages, les menu déroulants pour selection de certaines données, thème, design... 

ui <- page_navbar(
  title = "Tableau de bord - CatchMachine",
  theme = bs_theme(version = 5),
  fillable=FALSE,
  
  header = tags$style(HTML("
    .dashboard-card { border: 1px solid #444; border-radius: 4px; margin-bottom: 15px; }
    .dashboard-card .card-header { font-weight: 600; background: #fff; border-bottom: 1px solid #444; }
    .card-header-flex { display:flex; justify-content:space-between; align-items:center; }
    .year-select { max-width: 160px; }
  ")),
  
  nav_panel("Informations générales",
            layout_columns(
              col_widths = c(12),
              div(
                class = "p-2 mb-2 d-flex justify-content-between align-items-center flex-wrap",
                style = "background-color:#C1CDCD; border-radius:12px;",
                textOutput("bandeau_texte_p1"),
                
              )
            ),
            
            layout_columns(
              col_widths = c(3, 5, 4), # largeur des colonnes pour la page 1, le total doit faire 12 
              fill = FALSE,
              
              # ================= COLONNE GAUCHE =================
              
              #chaque case est codée de la même façon : 
              
              # div(                              #dit que tout ce qui est dans ce div() est dans la même colonne                        
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
                       "Décoché : elles utilisent le poids déclaré brut tel quel."), 
                     
                     selectInput("annee_page1", "Année (concerne toute la page)", choices = "Toutes les années", width = "220px")
                     
                ),
                
                card(class = "dashboard-card",
                     card_header("Informations générales"),
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
                     card_header("Evolution temporelle annuelle des captures"),
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
                       div(
                         div(style = "font-weight: 600; margin-bottom: 8px;",
                             "Top 5 des espèces les plus capturées"
                         ),
                         checkboxGroupInput(
                           "mois_top10",
                           "Mois :",
                           choices = mois_ordre,
                           selected = mois_ordre,
                           inline = TRUE
                         )
                       )
                     ),
                     
                     withSpinner(plotOutput("plot_top10_prises", height = "280px"), type = 4, color = "steelblue"), # ggplot des 5 especes les plus capturées en nombre
                     withSpinner(plotOutput("plot_top10_biomasse", height = "280px"), type = 4, color = "steelblue") # ggplot des 5 especes les plus capturées en nombre qui montre la biomasse
                     
                ),
                
                card(class = "dashboard-card",
                     card_header("Sorties infructueuses"),
                     withSpinner(plotOutput("plot_infructueuses", height = "280px"), type = 4, color = "steelblue")
                ),
                
                card(class = "dashboard-card",
                     card_header("Biomasse"),
                     withSpinner(uiOutput("biomasse_texte"), type = 4, color = "steelblue") 
                )
              )
            )
  ),
  
  nav_panel("Informations complémentaires",                                     # page 2
            layout_columns(
              col_widths = c(12),
              div(
                class = "p-2 mb-2 d-flex justify-content-between align-items-center flex-wrap",
                style = "background-color:#C1CDCD; border-radius:12px;",
                textOutput("bandeau_texte_p2"),
                
              )
            ),
            layout_columns(
              col_widths = c(12),
              div(
                class = "dashboard-card",
                selectInput("espece_page2", "Choisir une espèce (toute la page)",
                            choices = NULL, selected = NULL, width = "300px")
                
              )
            ),
            
            
            layout_columns(
              col_widths = c(6,6),
              fill = FALSE,
              
              # ================= COLONNE GAUCHE =================
              div(
                class = "d-flex flex-column",
                
                card(class = "dashboard-card",
                     card_header("Distribution des tailles par espèce"),
                     withSpinner(plotOutput("plot_taille_espece", height = "450px"), type = 4, color = "steelblue")
                )
              ), 
              # ================= COLONNE CENTRALE =================
              div(
                class = "d-flex flex-column", 
                
                card(class = "dashboard-card",
                     card_header("Relation taille-poids (RTP) par espèce"), 
                     withSpinner(plotOutput("plot_rtp_espece", height = "450px"), type = 4, color = "steelblue")
                )
              )
              
              
              
            )
            
  ), 
  nav_panel("Cartographie", # onglet 3/page 3 : disponible pour les AMP avec des sous zones de pêche 
            
            layout_columns(
              col_widths = c(12),
              div(
                class = "p-2 mb-2",
                style = "background-color:#C1CDCD; border-radius:12px;",
                textOutput("bandeau_texte_p3")
              )
            ),
            
            
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
  # Gestion du poids utilisé dans les calculs (brut ou corrigé)
  # -------------------------------------------------------------------
  # Objectif : créer une colonne commune `poids_final_total` qui sera utilisée
  # ensuite pour calculer CPUE et biomasse, indépendamment du fait que l'utilisateur
  # active ou non la correction des poids via la case à cocher.
  # Cela évite de dupliquer le code de calcul selon l'état de la case.
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
  # Mise à jour dynamique des menus déroulants de l'interface
  # -------------------------------------------------------------------
  # Cette section remplit automatiquement les choix des menus déroulants
  # (années, modes de pêche, espèces) à partir des données chargées.
  # Cela rend l'application adaptable à n'importe quel jeu de données sans
  # modifier le code à chaque fois.
  # -------------------------------------------------------------------
  
  # `observeEvent(df(), { ... })` signifie : "à chaque fois que `df()` change,
  # exécute le code ci-dessous". C'est utile quand les données sont rechargées
  # ou mises à jour dynamiquement.
 
  # --- Années disponibles (pour la page 1) ---
  # On récupère toutes les années présentes dans la colonne `annee` du jeu de données,
  # on les trie, puis on ajoute une option "Toutes les années" pour permettre
  # de ne pas filtrer par année.
  
  observeEvent(df(), {
    annees <- sort(unique(df()$annee))                                           
    choix  <- c("Toutes les années", as.character(annees))
    
    # On met à jour le menu déroulant `annee_page1` avec ces choix.
    # `selected` définit l'option sélectionnée par défaut.
    
    updateSelectInput(session, "annee_page1", choices = choix, selected = "Toutes les années")
    
    
    # --- Modes de pêche disponibles ---
    # On récupère les modes de pêche uniques, on les trie, puis on retire le mode
    # "Indeterminable" qui n'est pas pertinent pour les analyses.
    
    modes <- sort(unique(df()$Mode_peche)) # sort les modes de pêche des données
    modes <- modes[modes != "Indeterminable"] # on enleve le mode indterminable
    
    # On met à jour deux menus différents qui utilisent les mêmes options de modes de pêche :
    # - `mode_peche_cpue` : utilisé dans l'onglet "CPUE espèces"
    # - `mode_peche_cpue_tableau` : utilisé dans l'onglet "Tableaux résumés"
    
    updateSelectInput(session, "mode_peche_cpue", choices = modes, selected = modes[1])
    updateSelectInput(session, "mode_peche_cpue_tableau", choices = modes, selected = modes[1])
    
    # --- Espèces disponibles pour la page 2 ---
    # On récupère tous les noms vernaculaires d'espèces présents dans les données,
    # on les trie, et on alimente le menu `espece_page2` (utilisé pour les graphiques
    # de distribution des tailles et RTP sur la page 2).
    especes_page2 <- sort(unique(df()$NomVernaculaire))
    updateSelectInput(session, "espece_page2", choices = especes_page2, selected = especes_page2[1])
    
    
  })
  
  # -------------------------------------------------------------------
  # Fonction utilitaire de filtrage par année
  # -------------------------------------------------------------------
  # Cette fonction est réutilisée à plusieurs endroits du code pour
  # filtrer un jeu de données selon l'année choisie par l'utilisateur.
  # Elle évite de répéter le même code de filtrage partout.
  # -------------------------------------------------------------------
  filtrer_annee <- function(data, annee_choisie) {
    if (is.null(annee_choisie) || annee_choisie == "Toutes les années") {       #Si aucune année n'est sélectionnée ou si l'option "Toutes les années" est choisie, on retourne le jeu de données tel quel (pas de filtrage).
      return(data)
    }
    
    dplyr::filter(data, as.character(annee) == annee_choisie)                   # Sinon, on filtre les lignes où l'année (convertie en caractère pour être sûr de bien comparer) correspond à l'année choisie.
  }
  
  
  # -------------------------------------------------------------------
  # Calcul de la limite haute d'un axe Y adapté à des boxplots
  # -------------------------------------------------------------------
  # Objectif : déterminer une valeur maximale pour l'axe Y d'un graphique
  # de type boxplot, de manière à ne pas couper les boîtes ni leurs moustaches.
  # Plutôt que d'utiliser un simple quantile des données brutes (ce qui peut
  # tronquer des boîtes dans le cas de données très étalées ou peu nombreuses),
  # on utilise les vraies moustaches supérieures calculées par `boxplot.stats()`,
  # groupe par groupe, puis on prend le maximum de ces moustaches, avec une
  # petite marge de sécurité.
  # -------------------------------------------------------------------
  calc_ymax_boxplot <- function(valeurs, groupes, marge = 1.05) {               # la marge est à changer ici si besoin
    d <- data.frame(valeur = valeurs, groupe = groupes) %>%
      filter(!is.na(valeur))
    
    if (nrow(d) == 0) return(NA_real_)
    
    moustaches_sup <- d %>%                                                     # Pour chaque groupe, on calcule la moustache supérieure du boxplot via `boxplot.stats()$stats[5]` (5ème élément = upper whisker).
      group_by(groupe) %>%
      summarise(sup = boxplot.stats(valeur)$stats[5], .groups = "drop") %>%
      pull(sup)                                                                 # On extrait le vecteur des moustaches supérieures
    
    max(moustaches_sup, na.rm = TRUE) * marge                                   #On prend la plus grande de ces moustaches et on lui applique une marge (par défaut 5 %) pour éviter que le haut du graphique ne colle exactementà la moustache et rende le rendu moins lisible.
  }
  
  
  # --------------------------------------------------------------------
  # bandeau de texte warning en haut de chaque page
  # --------------------------------------------------------------------
  
  output$bandeau_texte_p1 <- renderText({
    "ATTENTION : les résultats présentés ici ne sont qu'à titre d'information et doivent être interprétés avec nuance. Ils ne reflètent pas l'exacte réalité de la pêche récréative dans les zones étudiées."
  })
  
  output$bandeau_texte_p2 <- renderText({
    "ATTENTION : les résultats présentés ici ne sont qu'à titre d'information et doivent être interprétés avec nuance. Ils ne reflètent pas l'exacte réalité de la pêche récréative dans les zones étudiées."
  })
  
  output$bandeau_texte_p3 <- renderText({
    "ATTENTION : les résultats présentés ici ne sont qu'à titre d'information et doivent être interprétés avec nuance. Ils ne reflètent pas l'exacte réalité de la pêche récréative dans les zones étudiées."
  })
  
  
  # --------------------------------------------------------------------
  # Résumé des indicateurs clés (nb pêcheurs, sorties, captures)
  # --------------------------------------------------------------------
  # Ce bloc crée un jeu de données filtré selon l'année choisie sur la page 1
  # --------------------------------------------------------------------
  
  df_info <- reactive({
    req(df())
    filtrer_annee(df(), input$annee_page1)                                      # Application de la fonction utilitaire `filtrer_annee` définie plus haut :
                                                                                # - si "Toutes les années" est sélectionné, aucun filtrage n'est appliqué
                                                                                # - sinon, seules les lignes correspondant à l'année choisie sont conservées
  })
  
  output$resume_texte <- renderUI({
    req(df_info())
    
    nb_pecheurs     <- n_distinct(df_info()$Id_Abonné)                          # Nombre de pêcheurs distincts : on compte les identifiants uniques d'abonnés.
    nb_sorties      <- n_distinct(df_info()$ID_sortie_corrige)                  # Nombre de sorties de pêche distinctes : on compte les ID de sortie corrigés uniques.
    capture_totale  <- sum(df_info()$Nb_de_prises, na.rm = TRUE)                # Nombre total de captures : somme de la colonne `Nb_de_prises`, en ignorant les NA.

    
    tagList(
      p(strong(format(nb_pecheurs, big.mark = " ")), " pêcheurs déclarants sur la période sélectionnée."),
      p(strong(format(nb_sorties, big.mark = " ")), " sorties de pêche sur la période sélectionnée."),
      p(strong(format(capture_totale, big.mark = " ")), " captures déclarées sur la période sélectionnée.")
    )
  })
  
  # --------------------------------------------------------------------
  # Liste des espèces les plus pêchées (tableau dynamique, sans sélection)
  # --------------------------------------------------------------------
  # Ce bloc crée un tableau récapitulatif par espèce (nombre de captures,
  # biomasse totale, pourcentage du total) et l'affiche dans un widget DT
  # en mode "liste déroulante" (scrollable), sans permettre la sélection
  # de lignes : c'est purement un affichage informatif.
  # --------------------------------------------------------------------
  
  liste_especes <- reactive({
    req(data_avec_poids())                                                      # On s'assure que `data_avec_poids()` (le jeu de données avec la colonne `poids_final_total`) est disponible avant de continuer.
    data_avec_poids() %>%
      group_by(NomVernaculaire) %>%
      summarise(
        nb_captures = sum(Nb_de_prises, na.rm = TRUE),                              # - le nombre total de captures (somme des `Nb_de_prises`)
        biomasse_totale_kg = round(sum(poids_final_total, na.rm = TRUE) / 1000, 2), # - la biomasse totale en kg : somme des `poids_final_total` (en grammes) divisée par 1000 et arrondie à 2 décimales.
        .groups = "drop") %>%
      arrange(desc(nb_captures)) %>%
      mutate(pourcentage = round(100 * nb_captures / sum(nb_captures), 2)) %>%
      rename("Espèce" = NomVernaculaire, "Nb captures" = nb_captures, "% du total de capture" = pourcentage, "Biomasse totale (kg)"  = biomasse_totale_kg )  # Renommer les colonnes pour l'affichage final dans le tableau.
  })
  
  output$table_especes <- renderDT({
    datatable(
      liste_especes(),                                                          # Données à afficher
      rownames = FALSE,                                                         # Ne pas afficher les noms de lignes (index)
      selection = "none",                                                       # retire la possibilité de sélectionner des lignes
      options = list(
        dom = "t",
        scrollY = "553px",
        scrollCollapse = TRUE,
        paging = FALSE,
        ordering = FALSE
      )
    )
  })
  
  # --------------------------------------------------------------------
  # Top 5 des espèces pêchées, filtrable par année (page 1) et par mois
  # --------------------------------------------------------------------
  # Ce bloc crée :
  # - un jeu de données filtré par année et par mois sélectionnés
  # - un réactif qui calcule le top 5 des espèces (en nombre de prises)
  # - deux graphiques : l'un sur le nombre de prises, l'autre sur la biomasse
  # --------------------------------------------------------------------
  
  # ---------------------------------------------------------------------
  # 1. Données filtrées par année et par mois
  # ---------------------------------------------------------------------
  # `df_top10` est un réactif qui contient les données principales,
  # filtrées selon :
  # - l'année choisie sur la page 1 (`input$annee_page1`)
  # - les mois sélectionnés dans le widget `input$mois_top10`
  # (le nom `top10` est historique ; ici on affiche en réalité un top 5)
  # ---------------------------------------------------------------------
  df_top10 <- reactive({
    req(data_avec_poids())
    req(input$mois_top10)
    
    filtrer_annee(data_avec_poids(), input$annee_page1)%>%                      # Filtrage par année via la fonction utilitaire `filtrer_annee`.

      filter(mois %in% input$mois_top10)                                        # Filtrage par mois : on ne garde que les lignes dont la colonne `mois` correspond à l'un des mois sélectionnés dans `input$mois_top10`.
  })
  
  
  # ---------------------------------------------------------------------
  # 2. Calcul du top 5 des espèces (en nombre de prises)
  # ---------------------------------------------------------------------
  # `espece_top10` est un réactif qui :
  # - vérifie que les colonnes nécessaires sont présentes
  # - vérifie qu'il y a bien des données après filtrage
  # - agrège par espèce (nom scientifique) le nombre de prises et la biomasse
  # - sélectionne les 5 espèces les plus capturées
  # - rajoute les noms vernaculaires pour l'affichage
  # --------------------------------------------------------------------
  
  espece_top10 <- reactive({
    req(df_top10())
    
    validate(
      need(all(c("NomScientifique", "Nb_de_prises", "CodeFAO", "NomVernaculaire", "poids_final_total") %in% names(df_top10())),
           "Colonnes attendues manquantes après nettoyage (NomScientifique, Nb_de_prises, CodeFAO, NomVernaculaire, poids_final_total)"), 
      need(nrow(df_top10()) > 0, "Aucune donnée pour la sélection année/mois choisie")
    )
    
    top5 <- df_top10() %>%                                                      # calcul des informations de nombre de prises et biomasse pour les espèces 
      group_by(NomScientifique) %>%
      summarise(
        total = sum(Nb_de_prises, na.rm = TRUE),
        biomasse_totale_kg = round(sum(poids_final_total, na.rm = TRUE) / 1000, 2),
        .groups = "drop") %>%
      arrange(desc(total)) %>%
      slice_head(n = 5)                                                         # pour changer le nombre d'espèces affichées : changer le chiffre ici. c'est organisé dans l'ordre décroissant du nombre de prises, donc pour avoir les 10 sp les + capturées, mettre n = 10
    
    top5 %>%
      left_join(
        df_top10() %>% select(CodeFAO, NomVernaculaire, NomScientifique) %>% distinct(),
        by = "NomScientifique"
      ) %>%
      arrange(desc(total))
  })
  
  # ---------------------------------------------------------------------
  # 3. Graphique 1 : Top 5 des espèces - Nombre de prises
  # ---------------------------------------------------------------------
  # Barplot horizontal montrant le nombre total de prises par espèce
  # pour le top 5, avec les valeurs affichées à droite des barres.
  # ---------------------------------------------------------------------
  
  output$plot_top10_prises <- renderPlot({
    req(espece_top10())
    
    ggplot(espece_top10(), aes(x = total, y = reorder(NomVernaculaire, total))) +
      geom_col(fill = "skyblue1") +
      geom_text(aes(label = total), hjust = -0.2, size = 6) +
      labs(title = "Top 5 des espèces - Nombre de prises", x = "Nombre de prises", y = "Espèce") +
      theme_minimal() +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10),
                         expand = expansion(mult = c(0, 0.15))) +
      # Personnalisation de la taille et du style des textes.
      theme(
        text = element_text(size = 16),                                         # taille de base du texte
        plot.title = element_text(face = "bold", size = 16),                    # titre en gras
        axis.title = element_text(size = 14),                                   # titres d'axes
        axis.text = element_text(size = 13),                                    # textes des axes
        axis.text.y = element_text(size = 14)                                   # noms d'espèces (axe Y) un peu plus grands
      )
  })
  
  # ---------------------------------------------------------------------
  # 4. Graphique 2 : Top 5 des espèces - Biomasse (kg)
  # ---------------------------------------------------------------------
  # Barplot horizontal montrant la biomasse totale (en kg) par espèce
  # pour le même top 5, avec les valeurs affichées à droite des barres.
  # ---------------------------------------------------------------------
  
  
  output$plot_top10_biomasse <- renderPlot({
    req(espece_top10())
    
    ggplot(espece_top10(), aes(x = biomasse_totale_kg, y = reorder(NomVernaculaire, total))) +
      geom_col(fill = "seagreen3") +
      geom_text(aes(label = biomasse_totale_kg), hjust = -0.2, size = 6) +
      labs(title = "Top 5 des espèces - Biomasse", x = "Biomasse (kg)", y = "Espèce") +
      theme_minimal() +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10),
                         expand = expansion(mult = c(0, 0.20))) +
      theme(
        text = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        axis.text.y = element_text(size = 14)
      )
  })
  
  # ---------------------------------------------------------------------
  # f - Évolution temporelle des captures - filtrée par année (page 1)
  # ---------------------------------------------------------------------
  
  # Ce bloc crée :
  # - un jeu de données filtré par année
  # - un réactif qui agrège le nombre de prises par mois et par année
  # - un graphique en lignes montrant l'évolution mensuelle des captures
  #   pour chaque année, avec une zone grisée pour les mois futurs de
  #   l'année en cours (si présente).
  # ---------------------------------------------------------------------
  
  
  
  # ---------------------------------------------------------------------
  # 1. Données filtrées par année
  # ---------------------------------------------------------------------
  # `df_evolution` contient le jeu de données principal, filtré selon
  # l'année choisie sur la page 1 via `input$annee_page1`.
  # ---------------------------------------------------------------------
  
  
  df_evolution <- reactive({
    req(df())
    filtrer_annee(df(), input$annee_page1)
  })
  
  # ---------------------------------------------------------------------
  # 2. Agrégation du nombre de prises par mois et par année
  # ---------------------------------------------------------------------
  # `prise_mois` calcule, pour chaque année et chaque mois, le nombre total
  # de prises déclarées. Les mois sont ordonnés dans l'ordre calendaire
  # (janvier → décembre) pour un affichage correct sur le graphique.
  # ---------------------------------------------------------------------
  
  prise_mois <- reactive({
    req(df_evolution())
    
    ordre_mois <- c("janvier", "février", "mars", "avril", "mai", "juin",
                    "juillet", "août", "septembre", "octobre", "novembre", "décembre")
    
    df_evolution() %>%
      mutate(mois = factor(mois, levels = ordre_mois)) %>%
      group_by(annee, mois) %>%
      summarise(nb_prises = sum(Nb_de_prises, na.rm = TRUE), .groups = "drop") %>%      # Somme des prises (`Nb_de_prises`) pour chaque combinaison année/mois.
      arrange(annee, mois)
  })
  
  # ---------------------------------------------------------------------
  # 3. Graphique : évolution mensuelle des captures par année
  # ---------------------------------------------------------------------
  # Graphique en lignes avec points et étiquettes :
  # - une ligne par année
  # - le nombre de prises en ordonnée
  # - les mois en abscisse (janvier → décembre)
  # - une zone grisée à droite pour indiquer les mois "futurs" de l'année
  #   en cours (si cette année est présente dans les données).
  # ---------------------------------------------------------------------
  
  output$plot_evolution <- renderPlot({
    req(prise_mois())
    
    ordre_mois <- c("janvier", "février", "mars", "avril", "mai", "juin",
                    "juillet", "août", "septembre", "octobre", "novembre", "décembre")
    
    d <- prise_mois()
    annee_courante <- max(d$annee, na.rm = TRUE)
    d_max <- d %>%
      filter(annee == annee_courante) %>%
      filter(!is.na(mois)) %>%
      slice_max(mois, n = 1, with_ties = FALSE)
    
    mois_courant <- d_max$mois[1]
    
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
      
      labs(title = "Nombre de prises déclarées de pêche par mois",
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
  
  
  # ---------------------------------------------------------------------
  # Sorties infructueuses -> camembert, filtrable par année (page 1)
  # ---------------------------------------------------------------------
  # Ce bloc crée :
  # - un jeu de données filtré par année
  # - un réactif qui calcule, par zone de pêche, le nombre de sorties
  #   avec et sans captures (bredouilles)
  # - un réactif qui agrège ces counts au niveau global
  # - un graphique en camembert montrant la proportion de sorties
  #   infructueuses vs sorties avec captures.
  # ---------------------------------------------------------------------
  
  
  
  # ---------------------------------------------------------------------
  # 1. Données filtrées par année pour les sorties infructueuses
  # ---------------------------------------------------------------------
  # `df_infructueuses` est similaire à `df_evolution`, mais dédié à
  # l'analyse des sorties avec/sans captures.
  # ---------------------------------------------------------------------
  
  
  df_infructueuses <- reactive({
    req(df())
    filtrer_annee(df(), input$annee_page1)
  })
  
  
  # ---------------------------------------------------------------------
  # 2. Calcul du nombre de sorties avec/sans captures par zone
  # ---------------------------------------------------------------------
  # `sorties_bredouilles` :
  # - regroupe les données par zone de pêche et par sortie
  # - calcule le nombre total de captures par sortie
  # - classe chaque sortie comme "bredouille" (0 capture) ou non
  # - agrège par zone : nombre de sorties, nombre de bredouilles,
  #   et nombre de sorties avec captures.
  # ---------------------------------------------------------------------
  sorties_bredouilles <- reactive({
    req(df_infructueuses())
    
    df_infructueuses() %>%
      group_by(Nom_ZoneDePeche, ID_sortie_corrige) %>%
      summarise(nb_captures = sum(Nb_de_prises, na.rm = TRUE), .groups = "drop") %>%
      group_by(Nom_ZoneDePeche) %>%
      summarise(
        nb_sorties     = n_distinct(ID_sortie_corrige),                         # Nombre de sorties total
        nb_bredouilles = sum(nb_captures == 0),                                 # Nombre de sorties avec 0 capture (bredouilles/infructueuses).
        .groups = "drop"
      ) %>%
      mutate(nb_avec_captures = nb_sorties - nb_bredouilles)                    # Nombre de sorties avec >1 captures
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
  
  # ---------------------------------------------------------------------
  # 4. Graphique : camembert des sorties avec/sans captures
  # ---------------------------------------------------------------------
  # Camembert (pie chart) montrant la répartition entre :
  # - sorties avec captures
  # - sorties sans captures (bredouilles)
  # Les effectifs sont affichés au centre de chaque secteur.
  # ---------------------------------------------------------------------
  
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
          "nb_bredouilles" = "Nombre de sorties sans capture"
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
  
  # ---------------------------------------------------------------------
  # ONGLET 1 : CPUE générales (mode + zone)
  # ---------------------------------------------------------------------
  # Utilisation de data_avec_poids() (poids brut ou corrigé selon la case).
  # Échelle Y dynamique : pas de limite fixe, calculée via calc_ymax_boxplot().
  # ---------------------------------------------------------------------
  
  # ---------------------------------------------------------------------
  # Données de CPUE par sortie et par mode de pêche
  # ---------------------------------------------------------------------
  # À modifier si :
  # - changer la définition d'une "sortie" (group_by)
  # - ajouter/retirer des variables dans le résumé
  # ---------------------------------------------------------------------
  
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
        cpue_sortie_nb    = nb_captures / nb_pecheurs,                          # ind/pêcheur/sortie
        cpue_sortie_poids = poids_captures / nb_pecheurs                        # g/pêcheur/sortie
      ) %>%
      filter(Mode_peche != "Indeterminable")                                    # on exclue le mode indeterminable
  })
  
  # ---------------------------------------------------------------------
  # Graphique 1 : CPUE en nombre (ind/pêcheur/sortie) par mode de pêche
  # ---------------------------------------------------------------------
  # À modifier si :
  # -  changer le titre, les labels, la palette de couleurs, etc.
  # -  passer à un autre type de graphique
  # ---------------------------------------------------------------------
 
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
  
  # ---------------------------------------------------------------------
  # Graphique 2 : CPUE en poids (g/pêcheur/sortie) par mode de pêche
  # ---------------------------------------------------------------------
  # Même structure que le précédent, mais sur cpue_sortie_poids.
  # À modifier aux mêmes endroits (titres, type de géom, etc.).
  # ---------------------------------------------------------------------
  
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
  
  
  
  # ---------------------------------------------------------------------
  # ONGLET 2 : CPUE des 10 espèces les plus pêchées (par mode sélectionné)
  # ---------------------------------------------------------------------
  # Le top 10 est calculé une fois sur toutes les données, puis on filtre
  # par le mode de pêche choisi dans input$mode_peche_cpue.
  # Échelle Y dynamique via calc_ymax_boxplot().
  # ---------------------------------------------------------------------

  
  # ---------------------------------------------------------------------
  # Top 10 des espèces (tous modes confondus)
  # ---------------------------------------------------------------------
  # À modifier si :
  # - tu veux changer le nombre d'espèces (slice_head(n = ...))
  # - tu veux changer le critère de "top" (ex: biomasse au lieu de nb)
  # ---------------------------------------------------------------------
  
  cpue_espece_top10 <- reactive({
    req(data_avec_poids())
    
    data_avec_poids() %>%
      group_by(NomScientifique, NomVernaculaire) %>%
      summarise(total = sum(Nb_de_prises, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total)) %>%
      slice_head(n = 10)                                                        # changer ce chiffre pour un nombre d'espèces affichées différent
  })
  
  # ---------------------------------------------------------------------
  # Ordre des espèces (pour l'affichage et les facteurs)
  # ---------------------------------------------------------------------
  # À modifier pour changer l'ordre d'affichage (ex: par biomasse, par ordre
  #   alphabétique, etc.)
  # ---------------------------------------------------------------------
  
  ordre_especes_cpue <- reactive({
    req(cpue_espece_top10())
    cpue_espece_top10() %>%
      arrange(desc(total)) %>%
      pull(NomVernaculaire)
  })
  
  # ---------------------------------------------------------------------
  # CPUE par espèce pour le mode de pêche sélectionné (onglet 2)
  # ---------------------------------------------------------------------
  # À modifier si :
  # - changer le filtre de mode (ex: plusieurs modes)
  # - ajouter d'autres variables dans le résumé par sortie
  # ---------------------------------------------------------------------
  
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
  
  # ---------------------------------------------------------------------
  # Graphique 1 (onglet 2) : CPUE nb par espèce (mode choisi)
  # ---------------------------------------------------------------------
  # À modifier si :
  # - tu veux changer le titre, les labels, l'angle des étiquettes x, etc.
  # - tu veux afficher/masquer les outliers (outlier.alpha)
  # ---------------------------------------------------------------------
  
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
  
  # ---------------------------------------------------------------------
  # Graphique 2 (onglet 2) : CPUE poids par espèce (mode choisi)
  # ---------------------------------------------------------------------
  # Mêmes points de modification que pour le graphique en nombre.
  # ---------------------------------------------------------------------
  
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
  
  # ---------------------------------------------------------------------
  # ONGLET 3 : Tableaux résumés CPUE des 10 espèces (mode indépendant)
  # ---------------------------------------------------------------------
  # Utilise un menu de mode séparé : input$mode_peche_cpue_tableau.
  # Les 10 espèces sont les mêmes que pour l'onglet 2 (cpue_espece_top10).
  # ---------------------------------------------------------------------
  
  # ---------------------------------------------------------------------
  # CPUE par espèce pour le mode de pêche sélectionné (onglet 3)
  # ---------------------------------------------------------------------
  # À modifier si :
  # - tu veux changer le mode de filtrage (ex: plusieurs modes, zones, etc.)
  # - tu veux ajouter des colonnes dans le résumé par sortie
  # ---------------------------------------------------------------------
  
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
  
  # ---------------------------------------------------------------------
  # Tableau récapitulatif CPUE en nombre (moyenne, médiane, n_sorties)
  # ---------------------------------------------------------------------
  # À modifier si :
  # - changer les statistiques affichées (ex: ajouter sd, min, max)
  # - changer l'ordre ou le nom des colonnes
  # ---------------------------------------------------------------------
  
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
  
  # ---------------------------------------------------------------------
  # Tableau récapitulatif CPUE en poids (moyenne, médiane, n_sorties)
  # ---------------------------------------------------------------------
  # Même logique que pour le tableau en nombre, mais sur cpue_sortie_poids.
  # ---------------------------------------------------------------------
  
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
  
  # ---------------------------------------------------------------------
  # Affichage des tableaux dans l'UI (onglet 3)
  # ---------------------------------------------------------------------
  # À modifier si :
  # - changer les noms de colonnes affichés
  # - passer de renderTable à renderDT pour un tableau interactif
  # ---------------------------------------------------------------------
  
  output$table_cpue_nb <- renderTable({
    req(tableau_cpue_nb())
    tableau_cpue_nb() %>%
      rename("Nom vernaculaire" = NomVernaculaire,
             "Nom scientifique" = NomScientifique,
             "Nb sorties" = n_sorties, 
             "Moyenne" = moyenne, 
             "Médiane" = mediane)
  })
  
  output$table_cpue_poids <- renderTable({
    req(tableau_cpue_poids())
    tableau_cpue_poids() %>%
      rename("Nom vernaculaire" = NomVernaculaire, 
             "Nom scientifique" = NomScientifique,
             "Nb sorties" = n_sorties, 
             "Moyenne" = moyenne, 
             "Médiane" = mediane)
  })
  
  # ---------------------------------------------------------------------
  # Biomasse totale / conservée / relâchée
  # ---------------------------------------------------------------------
  # Utilisation de data_avec_poids() (poids brut ou corrigé selon la case).
  # Filtrage par année via input$annee_page1.
  # Affichage sous forme de texte dynamique (UI).
  # ---------------------------------------------------------------------
  
  # ---------------------------------------------------------------------
  # Données filtrées par année pour la biomasse
  # ---------------------------------------------------------------------
  # À modifier si :
  # - tu veux changer la période de filtrage (ex: ajouter un filtre par mois)
  # - tu veux utiliser un autre jeu de données de base
  # ---------------------------------------------------------------------
  
  df_biomasse <- reactive({
    req(data_avec_poids())
    filtrer_annee(data_avec_poids(), input$annee_page1)
  })
  
  # ---------------------------------------------------------------------
  # Calcul de la biomasse totale, conservée et relâchée (en kg)
  # ---------------------------------------------------------------------
  # À modifier si :
  # - changer la règle de répartition conservée/relâchée
  # - ajouter d'autres catégories (ex: biomasse par zone, par mode)
  # ---------------------------------------------------------------------
  biomasse_totale <- reactive({
    req(df_biomasse())
    
    df_biomasse() %>%                                                           # Estimation du poids conservé et relâché à partir des proportions Nb_conservées / Nb_de_prises et Nb_relâchées / Nb_de_prises
      mutate(
        poids_conserve = if_else(Nb_de_prises > 0,
                                 poids_final_total * (Nb_conservées / Nb_de_prises), 0),  #calcul du poids conservé
        poids_relache  = if_else(Nb_de_prises > 0,
                                 poids_final_total * (Nb_relâchées / Nb_de_prises), 0)    #calcul du poids relâché
      ) %>%
      summarise(
        biomasse_totale_kg    = round(sum(poids_final_total, na.rm = TRUE) / 1000, 2),
        biomasse_conservee_kg = round(sum(poids_conserve, na.rm = TRUE) / 1000, 2),
        biomasse_relachee_kg  = round(sum(poids_relache, na.rm = TRUE) / 1000, 2)
      )
  })
  
  # ---------------------------------------------------------------------
  # Affichage texte des biomasses dans l'UI
  # ---------------------------------------------------------------------
  # À modifier si :
  # -  changer les libellés, l'ordre, ajouter des % etc.
  # -  passer à un autre type d'affichage (tableau, etc.)
  # ---------------------------------------------------------------------
  
  output$biomasse_texte <- renderUI({
    req(biomasse_totale())
    b <- biomasse_totale()
    
    tagList(
      p(strong("Totale : "), format(b$biomasse_totale_kg, big.mark = " "), " kg"),
      p(strong("Conservée : "), format(b$biomasse_conservee_kg, big.mark = " "), " kg"),
      p(strong("Relâchée : "), format(b$biomasse_relachee_kg, big.mark = " "), " kg")
    )
  })
  
  # =====================================================================
  # PAGE 2
  # =====================================================================
  
  # ---------------------------------------------------------------------
  # Classe de taille (page 2 - espèce commune)
  # ---------------------------------------------------------------------
  # Histogramme des tailles pour l'espèce sélectionnée dans input$espece_page2.
  # Les lignes sont "désagrégées" (uncount) pour avoir une ligne par individu.
  # ---------------------------------------------------------------------
  
  # ---------------------------------------------------------------------
  # Données de taille pour l'espèce choisie
  # ---------------------------------------------------------------------
  # À modifier si :
  # - changer le filtre d'espèce (ex: plusieurs espèces)
  # - ajouter des filtres (taille min/max, zone, mode, etc.)
  # ---------------------------------------------------------------------
  
  data_taille_espece <- reactive({
    req(df(), input$espece_page2)
    
    df() %>%
      filter(NomVernaculaire == input$espece_page2) %>%
      uncount(Nb_de_prises) %>%                                                 # Transformation : une ligne par individu (Nb_de_prises fois la ligne)
      filter(!is.na(Taille_MM))                                                 # On ne garde que les lignes avec une taille valide
  })
  
  # ---------------------------------------------------------------------
  # Histogramme des tailles
  # ---------------------------------------------------------------------
  # À modifier pour :
  # - changer la largeur des bins (binwidth)
  # - changer la limite haute (actuellement ~ q95 arrondi à la dizaine)
  # - changer les couleurs, titres, etc.
  # ---------------------------------------------------------------------
  
  output$plot_taille_espece <- renderPlot({
    req(data_taille_espece())
    validate(need(nrow(data_taille_espece()) > 0, "Aucune donnée de taille pour cette espèce."))
    
    q95  <- quantile(data_taille_espece()$Taille_MM, probs = 0.995, na.rm = TRUE)  # Calcul d'une limite haute "raisonnable" : 99,5e percentile arrondi à la dizaine supérieure
    xmax <- ceiling(q95 / 10) * 10
    
    ggplot(data_taille_espece(), aes(x = Taille_MM)) +
      geom_histogram(
        binwidth = 10,                                                          # changer ici pour des classes plus fines/larges
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
        title = paste0("Distribution des tailles - ", input$espece_page2),
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
  
  
  # ---------------------------------------------------------------------
  # RTP (FishBase uniquement, avec diagnostic par point)
  # ---------------------------------------------------------------------
  # Construction d'un référentiel FishBase une fois par chargement de données,
  # puis affichage d'un nuage de points taille-poids avec courbe théorique
  # et diagnostic (RTP ok / aberrante / inconnue, valeurs manquantes).
  # ---------------------------------------------------------------------
  
  # ---------------------------------------------------------------------
  # Référentiel FishBase (calculé une seule fois)
  # ---------------------------------------------------------------------
  # À modifier si :
  # - tu veux changer les colonnes utilisées (ex: autres coefficients a/b)
  # - tu veux filtrer certaines espèces ou origines
  # ---------------------------------------------------------------------
  
  rtp_fishbase_ref <- reactive({
    req(df())
    
    tableau_a_b_especes %>%                                                     # On part du tableau a/b déjà chargé (tableau_a_b_especes)
      filter(NomScientifique %in% unique(df()$NomScientifique)) %>%
      dplyr::select(
        NomScientifique,
        a        = a_FB,
        b        = b_FB,
        origine  = origine_FB
      ) %>%
      distinct(NomScientifique, .keep_all = TRUE)
  })
  
  # ---------------------------------------------------------------------
  # Données de l'espèce choisie pour le graphique RTP
  # ---------------------------------------------------------------------
  # À modifier si :
  # - changer le critère de sélection des prises (ici Nb_de_prises == 1)
  # - changer la tolérance autour des tailles/poids max (ici * 1.2)
  # - changer les règles de détection d'outliers / diagnostic
  # ---------------------------------------------------------------------
  
  rtp_espece_data <- reactive({
    req(df(), input$espece_page2, rtp_fishbase_ref())
    
    data_coherente_FB <- df() %>%
      filter(NomVernaculaire == input$espece_page2) %>%
      filter(Nb_de_prises == 1) %>%                                             # On ne garde que les prises à 1 individu (pour éviter les moyennes implicites). a changer si necessaire 
      filter(
        Taille_MM <= Taille_max_mm * 1.2 | is.na(Taille_max_mm),                # changer la marge de tolérance ici. taille et poids dans ±20% des maxima connus
        Poids_g   <= Poids_max_g  * 1.2 | is.na(Poids_max_g)
      ) %>%
      # Conversion mm -> cm pour la courbe taille-poids
      mutate(Taille_MM = Taille_MM / 10) %>%
      rename(Taille_CM = Taille_MM) %>%
      # Détection des valeurs manifestement fausses (0)
      mutate(outlier_valeur = Taille_CM == 0 | Poids_g == 0) %>%
      # Jointure avec les paramètres a, b de FishBase
      left_join(
        rtp_fishbase_ref() %>% dplyr::select(NomScientifique, a, b, origine),
        by = "NomScientifique"
      ) %>%
      # Poids théorique selon la relation taille-poids (RTP)
      mutate(
        Poids_RTP  = a * Taille_CM ^ b,
        Poids_diff = abs(Poids_g - Poids_RTP),
        # Outlier RTP : écart relatif > 100% par rapport au poids théorique
        outlier_RTP = Poids_diff / Poids_RTP > 1,
        outlier = outlier_valeur | outlier_RTP,
        # Diagnostic textuel par point
        diag = ifelse(outlier_valeur, "Val. manquante", "val. ok"),
        diag = ifelse(!outlier_valeur, paste0(diag, " / ",
                                              ifelse(is.na(outlier_RTP), "RTP inconnue",
                                                     ifelse(outlier_RTP, "RTP aberrante", "RTP ok"))), diag)
      )
    
    # Si pas de données ou pas de paramètres a/b, on indique "non disponible"
    if (nrow(data_coherente_FB) == 0 || all(is.na(data_coherente_FB$a))) {
      return(list(disponible = FALSE))
    }
    
    # Construction de la courbe théorique (taille-poids) sur la plage observée
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
  
  # ---------------------------------------------------------------------
  # Graphique RTP : nuage de points + courbe théorique
  # ---------------------------------------------------------------------
  # À modifier si :
  # - changer les couleurs, le seuil d'outlier, les labels
  # - afficher/masquer certaines catégories de points
  # - changer la référence (ex: autre base que FishBase)
  # ---------------------------------------------------------------------
  
  output$plot_rtp_espece <- renderPlot({
    d <- rtp_espece_data()
    req(d)
    
    if (!isTRUE(d$disponible)) {
      ggplot() +
        annotate("text", x = 0, y = 0, label = "Pas de données disponibles pour cette espèce.", size = 6) +
        theme_void()
    } else {
      origine_sp <- unique(d$courbe$origine)
      subtitle_text <- paste0("référentiel FishBase, origine des données : ", ifelse(is.na(origine_sp), "NA", origine_sp))
      
      ggplot() +
        geom_line(data = d$courbe, aes(x = Taille, y = Poids), color = "brown", linewidth = 1) +
        geom_point(data = d$points, aes(x = Taille_CM, y = Poids_g, color = diag)) +
        labs(
          title = paste0("Courbe de ratio taille-poids - ", input$espece_page2),
          subtitle = subtitle_text,
          x = "Taille (cm)", y = "Poids (g)"
        ) +
        scale_color_manual(                                                     # Palette de couleurs pour les catégories de diagnostic
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
              plot.subtitle = element_text(size = 13),
              axis.text.y   = element_text(size = 15),
              axis.title.y = element_text(size = 15), 
              legend.text = element_text(size = 14)
        )
      
    }
  })
  
  # ---------------------------------------------------------------------
  # Graphique RTP : nuage de points + courbe théorique
  # ---------------------------------------------------------------------
  
  # ---------------------------------------------------------------------
  # Sélection du shapefile actif selon les zones présentes dans les données
  # ---------------------------------------------------------------------
  # À modifier si :
  # - ajout d'un nouveau site (ex: autre zone que 4 et 7)
  # - changement de la logique de sélection (ex: via un input utilisateur)
  # ---------------------------------------------------------------------
  
  
  # --- Détermine quel shapefile utiliser selon les zones présentes dans les données ---
  shp_actif <- reactive({
    req(df())
    
    zones_presentes <- unique(df()$ID_ZoneDePeche)
    # Choix du shapefile selon la zone détectée
    if (4 %in% zones_presentes) {
      # Site de la Calanque
      list(shp = calanque_shp, disponible = TRUE)
    } else if (7 %in% zones_presentes) {
      # Site de Banyuls
      list(shp = banyuls_shp, disponible = TRUE)                                # si nouveau shapefile est ajouté dans le jeu de données; le mettre dans global.R et ajouter ici 
    } else {                                                                    # rajouter ce code avec le bon ID de zone : 
                                                                                #} else if (ID_NOUVEAU %in% zones_presentes) {
                                                                                #list(shp = NOUVEAU_shp, disponible = TRUE)  
      # Aucun shapefile compatible
      list(shp = NULL, disponible = FALSE)
    }
  })
  
  # ---------------------------------------------------------------------
  # Mise à jour dynamique du menu des espèces (cartographie par espèce)
  # ---------------------------------------------------------------------
  # À modifier si :
  # - limiter la liste à certaines espèces (ex: top N, filtre par groupe)
  # - changer l'espèce sélectionnée par défaut
  # ---------------------------------------------------------------------
  
  observe({
    especes_dispo <- sort(unique(df()$NomVernaculaire))
    
    updateSelectInput(
      session,
      "espece_cartographie",
      choices = especes_dispo,
      selected = especes_dispo[1]                           # changer ici pour une autre espèce par défaut
    )
  })
  
  # ---------------------------------------------------------------------
  # CARTE 1 : pression de pêche par espèce sélectionnée
  # ---------------------------------------------------------------------
  # Cette carte affiche, par sous-zone, le nombre total de captures
  # pour l'espèce choisie dans input$espece_cartographie.
  # ---------------------------------------------------------------------
  
  # ---------------------------------------------------------------------
  # Préparation des données spatiales pour l'espèce choisie
  # ---------------------------------------------------------------------
  # À modifier si :
  # - changer la métrique (ex: biomasse au lieu de nb captures)
  # - ajouter des filtres (année, mode, etc.)
  # ---------------------------------------------------------------------
  
  shp_pression_espece <- reactive({
    req(input$espece_cartographie)
    req(shp_actif()$disponible)
    
    data_sp <- df() %>%
      filter(NomVernaculaire == input$espece_cartographie)
    
    # Agrégation par sous-zone : nombre de captures et nombre de sorties
    captures_par_zone <- data_sp %>%
      group_by(Nom_SousZoneDePeche) %>%
      summarise(
        total_nb_captures = sum(Nb_de_prises, na.rm = TRUE),
        nb_sorties        = n(),
        .groups = "drop"
      )
    
    # Jointure avec le shapefile (polygones des sous-zones)
    shp_actif()$shp %>%
      left_join(captures_par_zone, by = "Nom_SousZoneDePeche") %>%
      mutate(total_nb_captures = replace(total_nb_captures, is.na(total_nb_captures), 0)) # Remplace les NA (zones sans captures) par 0
  })
  
  # ---------------------------------------------------------------------
  # UI conditionnelle : affiche la carte ou un message d'avertissement
  # ---------------------------------------------------------------------
  # À modifier si :
  # - changement du message, la hauteur de la carte, etc.
  # ---------------------------------------------------------------------
  
  output$cartographie_ui <- renderUI({
    if (isTRUE(shp_actif()$disponible)) {
      withSpinner(leafletOutput("cartographie", height = "450px"), type = 4, color = "steelblue")     # hauteur de la carte ici :  height = "450px"
    } else {
      div(
        class = "alert alert-warning",
        "Pas de sous-zone disponible dans les données, la cartographie n'est pas faisable."
      )
    }
  })
  
  # ---------------------------------------------------------------------
  # Rendu de la carte Leaflet (par espèce)
  # ---------------------------------------------------------------------
  # À modifier pour :
  # - changer la palette de couleurs, le titre de la légende
  # - changer le contenu des popups / labels
  # - ajouter des couches (ex: points, fond différent, etc.)
  # ---------------------------------------------------------------------
  
  output$cartographie <- renderLeaflet({
    req(shp_actif()$disponible)
    shp_data <- shp_pression_espece()
    
    pal <- colorNumeric(
      palette = "inferno",
      domain  = shp_data$total_nb_captures,
      reverse = TRUE # <- inverser ou non la palette
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
  
  # ---------------------------------------------------------------------
  # CARTE 2 : pression de pêche toutes espèces confondues
  # ---------------------------------------------------------------------
  # Même logique que la carte 1, mais sans filtrer par espèce :
  # on agrège toutes les captures par sous-zone.
  # ---------------------------------------------------------------------
  
  # ---------------------------------------------------------------------
  # Préparation des données spatiales (toutes espèces)
  # ---------------------------------------------------------------------
  # À modifier si :
  # - tu veux changer la métrique (biomasse, CPUE, etc.)
  # - tu veux ajouter des filtres (année, mode, etc.)
  # ---------------------------------------------------------------------
  
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
  
  # ---------------------------------------------------------------------
  # UI conditionnelle pour la carte globale
  # ---------------------------------------------------------------------
  
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
  
  
  # ---------------------------------------------------------------------
  # Rendu de la carte Leaflet (toutes espèces)
  # ---------------------------------------------------------------------
  # À modifier aux mêmes endroits que pour la carte par espèce :
  # palette, popups, labels, légende, fond de carte, etc.
  # ---------------------------------------------------------------------
  
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