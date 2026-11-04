# ============================================================
# RHIZOSPHERE PRIMING MODEL - shiny app
# ============================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)

# ============================================================
# HELPERS: reference list for the About tab
# ============================================================
ref_item <- function(authors, source, doi, note = NULL) {
  tags$li(HTML(paste0(
    htmltools::htmlEscape(authors), " <em>", htmltools::htmlEscape(source), "</em> ",
    '<a href="https://doi.org/', doi, '" target="_blank">doi:', htmltools::htmlEscape(doi), '</a>',
    if (!is.null(note)) paste0(". ", htmltools::htmlEscape(note)) else "."
  )))
}

ref_group <- function(title, ...) {
  tagList(h5(strong(title)), tags$ul(...))
}

# ============================================================
# MODEL FUNCTIONS 
# ============================================================
microbial_cn <- function(cn_som = 22, cn_re) {
  cn_base <- 5 + (14 - 5) / (1 + exp(-0.49 * (cn_som - 13.23)))
  cn_base + 0.15 * (cn_re - 35)
}

model_step <- function(state, day, c_re_mean, cn_re_mean, priming_strength,
                       Km = 16, Km_plant = 50) {
  
  angle <- 2 * pi * (day - 1) / 60 + pi/2
  c_re <- pmax(0, c_re_mean + 5 * sin(angle))
  cn_re <- pmax(8, cn_re_mean + 20 * sin(angle))
  
  n_re <- c_re / cn_re
  cn_som <- 22
  cn_mo <- microbial_cn(cn_som, cn_re)
  cue <- 0.30
  
  growth_n_limited <- n_re * cn_mo
  c_consumed <- growth_n_limited / cue
  c_excess <- pmax(0, c_re - c_consumed)
  
  if (c_excess > 0) {
    growth_excess <- c_excess * cue
    n_mined_potential <- (growth_excess / cn_mo) * priming_strength
    mining_factor <- pmax(0.2, Km / (state$soil_N + Km))
    n_mined <- n_mined_potential * mining_factor
    r_nm_positive <- n_mined * cn_som
    n_to_soil <- n_mined * 0.5
  } else {
    n_mined <- r_nm_positive <- mining_factor <- n_to_soil <- 0
  }
  
  suppression <- if (cn_re < 35) ((35 - cn_re) / 25)^0.85 * 0.85 else 0
  suppression <- pmin(1, pmax(0, suppression))
  
  net_rpe <- r_nm_positive * (1 - suppression) * 1.52 -
    suppression * (c_re * 0.040) * 1.58
  
  total_growth <- growth_n_limited + ifelse(c_excess > 0, c_excess * cue, 0)
  
  grazing_consumption <- state$C_mo * 0.095
  microbial_respiration <- total_growth * (1 - cue)
  grazer_respiration <- grazing_consumption * 0.765
  
  mo_new <- max(4, state$C_mo + total_growth * 0.65 - grazing_consumption)
  grazer_new <- max(2, state$C_grazer + grazing_consumption * 0.235 - state$C_grazer * 0.015)
  
  new_plant_C <- state$plant_C + c_re * 0.25
  new_soil_N <- state$soil_N + n_to_soil + 0.005
  plant_N_uptake <- new_plant_C * 0.002 * (new_soil_N / (new_soil_N + Km_plant))
  new_soil_N <- pmax(0.1, new_soil_N - plant_N_uptake)
  
  cum_net_c_add <- c_re - microbial_respiration - grazer_respiration - net_rpe
  
  list(net_rpe = net_rpe, cum_net_c_add = cum_net_c_add,
       c_re = c_re, cn_re = cn_re,
       C_mo = mo_new, C_grazer = grazer_new, soil_N = new_soil_N,
       plant_N_uptake = plant_N_uptake, plant_C = new_plant_C,
       microbial_growth = total_growth,
       grazing_consumption = grazing_consumption,
       microbial_respiration = microbial_respiration,
       grazer_respiration = grazer_respiration)
}

# ============================================================
# HELPER: Run full simulation (returns sim df + metrics)
# ============================================================
run_full_simulation <- function(c_re_mean, cn_re_mean, priming_strength = 1.92,
                                Km = 16, Km_plant = 50,
                                init_C_mo = 8, init_C_grazer = 6,
                                init_soil_N = 15, init_plant_C = 20,
                                sim_days = 300) {
  
  state <- list(C_mo = init_C_mo, C_grazer = init_C_grazer,
                soil_N = init_soil_N, plant_C = init_plant_C)
  
  sim <- data.frame(day = 1:sim_days)
  sim$net_rpe <- sim$cum_net_c <- sim$c_re <- sim$cn_re <- 0
  sim$C_mo <- sim$C_grazer <- sim$soil_N <- sim$plant_N_uptake <- sim$plant_C <- 0
  sim$microbial_growth <- sim$grazing_consumption <-
    sim$microbial_respiration <- sim$grazer_respiration <- 0
  
  cum_net_c <- 0
  for (d in 1:sim_days) {
    out <- model_step(state, d, c_re_mean, cn_re_mean, priming_strength,
                      Km = Km, Km_plant = Km_plant)
    state <- list(C_mo = out$C_mo, C_grazer = out$C_grazer,
                  soil_N = out$soil_N, plant_C = out$plant_C)
    cum_net_c <- cum_net_c + out$cum_net_c_add
    
    sim$net_rpe[d] <- out$net_rpe
    sim$cum_net_c[d] <- cum_net_c
    sim$c_re[d] <- out$c_re
    sim$cn_re[d] <- out$cn_re
    sim$C_mo[d] <- out$C_mo
    sim$C_grazer[d] <- out$C_grazer
    sim$soil_N[d] <- out$soil_N
    sim$plant_N_uptake[d] <- out$plant_N_uptake
    sim$plant_C[d] <- out$plant_C
    sim$microbial_growth[d] <- out$microbial_growth
    sim$grazing_consumption[d] <- out$grazing_consumption
    sim$microbial_respiration[d] <- out$microbial_respiration
    sim$grazer_respiration[d] <- out$grazer_respiration
  }
  
  mean_rpe <- mean(sim$net_rpe, na.rm = TRUE)
  min_rpe  <- min(sim$net_rpe, na.rm = TRUE)
  max_rpe  <- max(sim$net_rpe, na.rm = TRUE)
  final_c  <- tail(sim$cum_net_c, 1)
  
  metrics <- data.frame(
    c_re_mean = c_re_mean,
    cn_re_mean = cn_re_mean,
    Km = Km,
    Km_plant = Km_plant,
    priming_strength = priming_strength,
    Mean_RPE = round(mean_rpe, 3),
    Min_RPE = round(min_rpe, 3),
    Max_RPE = round(max_rpe, 3),
    Final_Net_C_balance = round(final_c, 2),
    stringsAsFactors = FALSE
  )
  
  list(sim = sim, metrics = metrics)
}

# ============================================================
# PRESET SCENARIOS (with intended Km / Km_plant)
# ============================================================
presets <- list(
  "Mostly Positive" = list(
    name = "Mostly_positive",
    c_re_mean = 10, cn_re_mean = 43,
    Km = 6, Km_plant = 90, priming_strength = 1.92
  ),
  "Net RPE Balanced" = list(
    name = "NetRPE_balanced",
    c_re_mean = 19, cn_re_mean = 42,
    Km = 30, Km_plant = 50, priming_strength = 1.92
  ),
  "Symmetric Equal" = list(
    name = "Symmetric_equal",
    c_re_mean = 25, cn_re_mean = 40,
    Km = 35, Km_plant = 48, priming_strength = 1.92
  ),
  "Mostly Negative" = list(
    name = "Mostly_negative",
    c_re_mean = 15, cn_re_mean = 38,
    Km = 60, Km_plant = 45, priming_strength = 1.92
  )
)

# ============================================================
# UI
# ============================================================
ui <- navbarPage(
  "Rhizosphere Priming Effect (RPE) Model",
  
  # ---------- TAB 1: SINGLE SIMULATION ----------
  tabPanel("Single Simulation",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Exudate Inputs"),
        sliderInput("c_re_mean", "Mean daily C exudation (µg C g⁻¹ soil d⁻¹)", 
                    min = 5, max = 100, value = 19, step = 1),
        sliderInput("cn_re_mean", "Mean exudate C:N ratio", 
                    min = 5, max = 60, value = 42, step = 1),
        
        hr(),
        h4("Process Parameters"),
        sliderInput("priming_strength", "Priming strength multiplier", 
                    min = 0.5, max = 3.0, value = 1.92, step = 0.01),
        sliderInput("Km", "Km — microbial N-mining half-saturation", 
                    min = 1, max = 100, value = 30, step = 1),
        sliderInput("Km_plant", "Km_plant — plant N-uptake half-saturation", 
                    min = 1, max = 100, value = 50, step = 1),
        
        hr(),
        h4("Initial Conditions"),
        numericInput("init_C_mo", "Initial microbial biomass C", value = 8, min = 1, step = 1),
        numericInput("init_C_grazer", "Initial grazer biomass C", value = 6, min = 1, step = 1),
        numericInput("init_soil_N", "Initial soil mineral N", value = 15, min = 0.1, step = 0.5),
        numericInput("init_plant_C", "Initial plant C pool", value = 20, min = 1, step = 1),
        numericInput("sim_days", "Simulation length (days)", value = 300, min = 100, max = 1000, step = 50),
        
        hr(),
        actionButton("run_sim", "▶ Run Simulation", class = "btn-primary btn-lg", width = "100%"),
        
        hr(),
        h4("Quick Presets"),
        fluidRow(
          column(6, actionButton("preset1", "Mostly +", class = "btn-success btn-sm", width = "100%")),
          column(6, actionButton("preset2", "Balanced", class = "btn-info btn-sm", width = "100%"))
        ),
        br(),
        fluidRow(
          column(6, actionButton("preset3", "Symmetric", class = "btn-warning btn-sm", width = "100%")),
          column(6, actionButton("preset4", "Mostly −", class = "btn-danger btn-sm", width = "100%"))
        ),
        
        hr(),
        downloadButton("dl_sim_csv", "Download simulation data (.csv)", class = "btn-outline-primary btn-sm"),
        br(), br(),
        downloadButton("dl_rpe_pdf", "Download current RPE plot (.pdf)", class = "btn-outline-primary btn-sm")
      ),
      
      mainPanel(
        width = 9,
        tabsetPanel(
          tabPanel("Key Metrics",
            br(),
            fluidRow(
              column(3, wellPanel(h5("Mean Net RPE"), h3(textOutput("val_mean_rpe")), style = "text-align:center;")),
              column(3, wellPanel(h5("Min Net RPE"), h3(textOutput("val_min_rpe")), style = "text-align:center;")),
              column(3, wellPanel(h5("Max Net RPE"), h3(textOutput("val_max_rpe")), style = "text-align:center;")),
              column(3, wellPanel(h5("Final Net C Balance"), h3(textOutput("val_final_c")), style = "text-align:center;"))
            ),
            br(),
            h4("Simulation Summary"),
            tableOutput("metrics_table"),
            hr(),
            p(em("Net RPE > 0 indicates net positive priming (more soil C loss). Net RPE < 0 indicates net suppression or C stabilization."))
          ),
          
          tabPanel("RPE & Carbon Balance",
            plotOutput("plot_rpe", height = "320px"),
            plotOutput("plot_balance", height = "320px")
          ),
          
          tabPanel("Biomass & Fluxes",
            plotOutput("plot_biomass", height = "320px"),
            plotOutput("plot_fluxes", height = "320px")
          ),
          
          tabPanel("Soil N & Diagnostics",
            plotOutput("plot_soilN", height = "320px"),
            plotOutput("plot_diagnostics", height = "320px")
          ),
          
          tabPanel("Full 8-Panel View",
            plotOutput("plot_full", height = "950px")
          )
        )
      )
    )
  ),
  
  # ---------- TAB 2: SCENARIO COMPARISON ----------
  tabPanel("Scenario Comparison",
    fluidPage(
      h3("Compare the four targeted scenarios from the original study"),
      actionButton("run_compare", "Run / Refresh All Four Scenarios", class = "btn-primary"),
      br(), br(),
      plotOutput("plot_compare_2x2", height = "620px"),
      br(),
      h4("Scenario Metrics Summary"),
      tableOutput("compare_table")
    )
  ),
  
  # ---------- TAB 3: SENSITIVITY (Km × Km_plant) ----------
  tabPanel("Sensitivity Explorer (Km × Km_plant)",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Base Scenario"),
        selectInput("sens_preset", "Choose base preset", 
                    choices = names(presets), selected = "Net RPE Balanced"),
        helpText("The selected preset sets c_re_mean, cn_re_mean and priming_strength. Km and Km_plant are varied in the grid below."),
        
        hr(),
        h4("Grid Settings"),
        sliderInput("sens_km_min", "Km minimum", min = 1, max = 50, value = 1, step = 5),
        sliderInput("sens_km_max", "Km maximum", min = 10, max = 101, value = 101, step = 5),
        sliderInput("sens_kmp_min", "Km_plant minimum", min = 1, max = 50, value = 1, step = 5),
        sliderInput("sens_kmp_max", "Km_plant maximum", min = 10, max = 101, value = 101, step = 5),
        sliderInput("sens_step", "Step size (resolution)", min = 2, max = 10, value = 5, step = 1),
        
        hr(),
        actionButton("run_sens", "Run Sensitivity Grid", class = "btn-success btn-lg", width = "100%"),
        br(), br(),
        helpText("A 21×21 grid (step=5) runs in < 2 seconds. Larger grids are also fast."),
        downloadButton("dl_sens_csv", "Download sensitivity results (.csv)")
      ),
      
      mainPanel(
        width = 9,
        tabsetPanel(
          tabPanel("Mean RPE Heatmap",
            plotOutput("sens_heatmap_rpe", height = "650px")
          ),
          tabPanel("Final C Balance Heatmap",
            plotOutput("sens_heatmap_cbal", height = "650px")
          ),
          tabPanel("Results Table (top 20)",
            tableOutput("sens_table")
          )
        )
      )
    )
  ),
  
  # ---------- TAB 4: ABOUT ----------
  tabPanel("About / Instructions",
    fluidPage(
      h2("Rhizosphere Priming Effect (RPE) Interactive Explorer"),
      p("This Shiny app lets you interactively explore a mechanistic model of rhizosphere priming, microbial N-mining, plant N-uptake competition, and grazer dynamics."),
      h4("How to use"),
      tags$ol(
        tags$li(strong("Single Simulation tab:"), "Adjust sliders or click a preset, then press 'Run Simulation'. Explore the different plot tabs."),
        tags$li(strong("Scenario Comparison tab:"), "Run all four classic scenarios side-by-side (2×2 RPE time series + summary table)."),
        tags$li(strong("Sensitivity Explorer:"), "Pick a base scenario and vary microbial Km and plant Km_plant across a grid to see how they control mean RPE and net C balance (heatmaps)."),
        tags$li("All plots and data can be downloaded.")
      ),
      h4("Key Parameters"),
      tags$ul(
        tags$li(strong("c_re_mean / cn_re_mean:"), "Average daily root exudate C input and its C:N ratio (seasonal sine wave added internally)."),
        tags$li(strong("Km:"), "Controls how strongly microbes mine soil organic N when mineral N is low (lower Km = stronger mining response at low N)."),
        tags$li(strong("Km_plant:"), "Controls plant competitive strength for mineral N (lower Km_plant = stronger plant uptake)."),
        tags$li(strong("priming_strength:"), "Multiplier on potential N mining from excess C."),
        tags$li(strong("Net RPE:"), "Net rhizosphere priming effect on soil C (positive = net loss of soil C, negative = net gain/stabilization).")
      ),
      hr(),
      h4("Key references"),
      ref_group("Model basis",
        ref_item("Michel (2026)", "Biol Fertil Soils 62: 935–955.", "10.1007/s00374-026-02021-8"),
        ref_item("Chertov et al. (2022)", "Plants 11: 2605.", "10.3390/plants11192605")
      ),
      ref_group("Priming concept",
        ref_item("Kuzyakov, Friedel & Stahr (2000)", "Soil Biol Biochem 32: 1485–1498.", "10.1016/s0038-0717(00)00084-5"),
        ref_item("Cheng et al. (2014)", "New Phytol 201: 31–44.", "10.1111/nph.12440")
      ),
      ref_group("N mining, the Km term",
        ref_item("Craine, Morrow & Fierer (2007)", "Ecology 88: 2105–2113.", "10.1890/06-1847.1"),
        ref_item("Chen et al. (2014)", "Glob Change Biol 20: 2356–2367.", "10.1111/gcb.12475"),
        ref_item("Dijkstra et al. (2013)", "Front Microbiol 4.", "10.3389/fmicb.2013.00216")
      ),
      ref_group("Negative priming",
        ref_item("Blagodatskaya & Kuzyakov (2008)", "Biol Fertil Soils 45: 115–131.", "10.1007/s00374-008-0334-y"),
        ref_item("Michel et al. (2023)", "Biogeochemistry 162: 141–161.", "10.1007/s10533-022-00996-8")
      ),
      ref_group("Microbial CUE and C:N",
        ref_item("Manzoni et al. (2012)", "New Phytol 196: 79–91.", "10.1111/j.1469-8137.2012.04225.x"),
        ref_item("Mooshammer et al. (2014)", "Front Microbiol 5.", "10.3389/fmicb.2014.00022")
      ),
      ref_group("Plant–microbe N competition",
        ref_item("Kuzyakov & Xu (2013)", "New Phytol 198: 656–669.", "10.1111/nph.12235"),
        ref_item("Hodge, Robinson & Fitter (2000)", "Trends Plant Sci 5: 304–308.", "10.1016/S1360-1385(00)01656-3")
      ),
      ref_group("Grazers and the microbial loop",
        ref_item("Clarholm (1985)", "Soil Biol Biochem 17: 181–187.", "10.1016/0038-0717(85)90113-0"),
        ref_item("Bonkowski (2004)", "New Phytol 162: 617–631.", "10.1111/j.1469-8137.2004.01066.x")
      ),
      ref_group("Earlier priming models",
        ref_item("Blagodatsky et al. (2010)", "Soil Biol Biochem 42: 1275–1283.", "10.1016/j.soilbio.2010.04.005"),
        ref_item("Wutzler & Reichstein (2013)", "Biogeosciences 10: 2089–2103.", "10.5194/bg-10-2089-2013")
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {
  
  # Reactive value to store latest single-run results
  rv <- reactiveValues(
    sim = NULL,
    metrics = NULL,
    last_params = NULL
  )
  
  # ---- Run simulation when button pressed ----
  observeEvent(input$run_sim, {
    res <- run_full_simulation(
      c_re_mean = input$c_re_mean,
      cn_re_mean = input$cn_re_mean,
      priming_strength = input$priming_strength,
      Km = input$Km,
      Km_plant = input$Km_plant,
      init_C_mo = input$init_C_mo,
      init_C_grazer = input$init_C_grazer,
      init_soil_N = input$init_soil_N,
      init_plant_C = input$init_plant_C,
      sim_days = input$sim_days
    )
    rv$sim <- res$sim
    rv$metrics <- res$metrics
    rv$last_params <- list(
      c_re_mean = input$c_re_mean, cn_re_mean = input$cn_re_mean,
      Km = input$Km, Km_plant = input$Km_plant,
      priming_strength = input$priming_strength
    )
  })
  
  # ---- Preset buttons: update UI + run simulation directly ----
  observeEvent(input$preset1, {
    p <- presets[["Mostly Positive"]]
    updateSliderInput(session, "c_re_mean", value = p$c_re_mean)
    updateSliderInput(session, "cn_re_mean", value = p$cn_re_mean)
    updateSliderInput(session, "Km", value = p$Km)
    updateSliderInput(session, "Km_plant", value = p$Km_plant)
    updateSliderInput(session, "priming_strength", value = p$priming_strength)
    
    res <- run_full_simulation(p$c_re_mean, p$cn_re_mean, p$priming_strength,
                               p$Km, p$Km_plant)
    rv$sim <- res$sim
    rv$metrics <- res$metrics
    rv$last_params <- p
  })
  
  observeEvent(input$preset2, {
    p <- presets[["Net RPE Balanced"]]
    updateSliderInput(session, "c_re_mean", value = p$c_re_mean)
    updateSliderInput(session, "cn_re_mean", value = p$cn_re_mean)
    updateSliderInput(session, "Km", value = p$Km)
    updateSliderInput(session, "Km_plant", value = p$Km_plant)
    updateSliderInput(session, "priming_strength", value = p$priming_strength)
    
    res <- run_full_simulation(p$c_re_mean, p$cn_re_mean, p$priming_strength,
                               p$Km, p$Km_plant)
    rv$sim <- res$sim
    rv$metrics <- res$metrics
    rv$last_params <- p
  })
  
  observeEvent(input$preset3, {
    p <- presets[["Symmetric Equal"]]
    updateSliderInput(session, "c_re_mean", value = p$c_re_mean)
    updateSliderInput(session, "cn_re_mean", value = p$cn_re_mean)
    updateSliderInput(session, "Km", value = p$Km)
    updateSliderInput(session, "Km_plant", value = p$Km_plant)
    updateSliderInput(session, "priming_strength", value = p$priming_strength)
    
    res <- run_full_simulation(p$c_re_mean, p$cn_re_mean, p$priming_strength,
                               p$Km, p$Km_plant)
    rv$sim <- res$sim
    rv$metrics <- res$metrics
    rv$last_params <- p
  })
  
  observeEvent(input$preset4, {
    p <- presets[["Mostly Negative"]]
    updateSliderInput(session, "c_re_mean", value = p$c_re_mean)
    updateSliderInput(session, "cn_re_mean", value = p$cn_re_mean)
    updateSliderInput(session, "Km", value = p$Km)
    updateSliderInput(session, "Km_plant", value = p$Km_plant)
    updateSliderInput(session, "priming_strength", value = p$priming_strength)
    
    res <- run_full_simulation(p$c_re_mean, p$cn_re_mean, p$priming_strength,
                               p$Km, p$Km_plant)
    rv$sim <- res$sim
    rv$metrics <- res$metrics
    rv$last_params <- p
  })
  
  # ---- Reactive plot builders (used by multiple outputs) ----
  theme_custom <- theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_blank(),
          panel.grid.minor = element_blank())
  
  make_rpe_plot <- function(sim) {
    ggplot(sim, aes(day, net_rpe)) +
      geom_line(color = "#1f77b4", linewidth = 1.3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
      labs(x = "Day", y = "Net RPE") +
      theme_custom
  }
  
  make_balance_plot <- function(sim) {
    ggplot(sim, aes(day, cum_net_c)) +
      geom_line(color = "#2ca02c", linewidth = 1.3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
      labs(x = "Day", y = "Cumulative Net C Balance") +
      theme_custom
  }
  
  make_biomass_plot <- function(sim) {
    data_biomass <- sim %>%
      select(day, C_mo, C_grazer) %>%
      pivot_longer(-day, names_to = "compartment", values_to = "value")
    ggplot(data_biomass, aes(day, value, color = compartment)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("C_mo" = "#d62728", "C_grazer" = "#ff7f0e"),
                         labels = c("Microbial C", "Grazer C")) +
      labs(x = "Day", y = "Biomass C (µg g⁻¹)", color = NULL) +
      theme_custom
  }
  
  make_flux_plot <- function(sim) {
    data_flux <- sim %>%
      select(day, microbial_growth, grazing_consumption,
             microbial_respiration, grazer_respiration) %>%
      pivot_longer(-day, names_to = "flux", values_to = "value")
    ggplot(data_flux, aes(day, value, color = flux)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("microbial_growth" = "#1f77b4",
                                    "grazing_consumption" = "#9467bd",
                                    "microbial_respiration" = "#d62728",
                                    "grazer_respiration" = "#ff7f0e"),
                         labels = c("Microbial growth", "Grazing consumption",
                                    "Microbial respiration", "Grazer respiration")) +
      labs(x = "Day", y = "C flux (µg g⁻¹ d⁻¹)", color = NULL) +
      theme_custom
  }
  
  make_soilN_plot <- function(sim) {
    data_soil <- sim %>%
      select(day, soil_N, plant_N_uptake) %>%
      pivot_longer(-day, names_to = "variable", values_to = "value")
    ggplot(data_soil, aes(day, value, color = variable)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("soil_N" = "#2ca02c", "plant_N_uptake" = "#1f77b4"),
                         labels = c("Soil mineral N", "Plant N uptake")) +
      labs(x = "Day", y = "N (µg g⁻¹ or µg g⁻¹ d⁻¹)", color = NULL) +
      theme_custom
  }
  
  make_diagnostics_plot <- function(sim) {
    data_diag <- sim %>%
      mutate(cn_mo = microbial_cn(22, cn_re))
    p1 <- ggplot(data_diag, aes(day, cn_mo)) +
      geom_line(color = "#9467bd", linewidth = 1.2) +
      geom_hline(yintercept = 22, linetype = "dotted", color = "gray50") +
      labs(x = "Day", y = "Microbial C:N") +
      theme_custom
    p2 <- ggplot(data_diag, aes(day, cn_re)) +
      geom_line(color = "#d62728", linewidth = 1.2, linetype = "dashed") +
      geom_hline(yintercept = 35, linetype = "dotted", color = "gray50") +
      labs(x = "Day", y = "Exudate C:N") +
      theme_custom
    grid.arrange(p1, p2, ncol = 1)
  }
  
  # ---- Outputs for Single Simulation tab ----
  output$val_mean_rpe <- renderText({
    req(rv$metrics)
    rv$metrics$Mean_RPE
  })
  output$val_min_rpe <- renderText({
    req(rv$metrics)
    rv$metrics$Min_RPE
  })
  output$val_max_rpe <- renderText({
    req(rv$metrics)
    rv$metrics$Max_RPE
  })
  output$val_final_c <- renderText({
    req(rv$metrics)
    rv$metrics$Final_Net_C_balance
  })
  
  output$metrics_table <- renderTable({
    req(rv$metrics)
    rv$metrics
  }, rownames = FALSE, striped = TRUE, hover = TRUE)
  
  output$plot_rpe <- renderPlot({
    req(rv$sim)
    make_rpe_plot(rv$sim)
  })
  
  output$plot_balance <- renderPlot({
    req(rv$sim)
    make_balance_plot(rv$sim)
  })
  
  output$plot_biomass <- renderPlot({
    req(rv$sim)
    make_biomass_plot(rv$sim)
  })
  
  output$plot_fluxes <- renderPlot({
    req(rv$sim)
    make_flux_plot(rv$sim)
  })
  
  output$plot_soilN <- renderPlot({
    req(rv$sim)
    make_soilN_plot(rv$sim)
  })
  
  output$plot_diagnostics <- renderPlot({
    req(rv$sim)
    make_diagnostics_plot(rv$sim)
  })
  
  output$plot_full <- renderPlot({
    req(rv$sim)
    sim <- rv$sim
    
    data_rpe <- data.frame(day = sim$day, net_rpe = sim$net_rpe)
    data_biomass <- sim %>% select(day, C_mo, C_grazer) %>%
      pivot_longer(-day, names_to = "compartment", values_to = "value")
    data_flux <- sim %>% select(day, microbial_growth, grazing_consumption,
                                microbial_respiration, grazer_respiration) %>%
      pivot_longer(-day, names_to = "flux", values_to = "value")
    data_soilN <- sim %>% select(day, soil_N, plant_N_uptake) %>%
      pivot_longer(-day, names_to = "variable", values_to = "value")
    data_balance <- data.frame(day = sim$day, cum_net_c = sim$cum_net_c)
    data_exudates <- sim %>% select(day, c_re, cn_re) %>%
      mutate(n_from_exudates = c_re / cn_re)
    data_diagnostic <- sim %>% mutate(cn_mo = microbial_cn(22, cn_re))
    
    p1 <- ggplot(data_rpe, aes(day, net_rpe)) +
      geom_line(color = "#1f77b4", linewidth = 1.4) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      labs(y = "Net RPE", x = NULL) + theme_custom
    
    p2 <- ggplot(data_biomass, aes(day, value, color = compartment)) +
      geom_line(linewidth = 1.3) +
      scale_color_manual(values = c("C_mo" = "#d62728", "C_grazer" = "#ff7f0e")) +
      labs(y = "Biomass C", x = NULL) + theme_custom + theme(legend.title = element_blank())
    
    p3 <- ggplot(data_flux, aes(day, value, color = flux)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("microbial_growth"="#1f77b4","grazing_consumption"="#9467bd",
                                    "microbial_respiration"="#d62728","grazer_respiration"="#ff7f0e")) +
      labs(y = "C flux", x = NULL) + theme_custom + theme(legend.title = element_blank())
    
    p4 <- ggplot(data_soilN, aes(day, value, color = variable)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("soil_N"="#2ca02c","plant_N_uptake"="#1f77b4")) +
      labs(y = "Soil N / Plant uptake", x = NULL) + theme_custom + theme(legend.title = element_blank())
    
    p5 <- ggplot(data_balance, aes(day, cum_net_c)) +
      geom_line(color = "#2ca02c", linewidth = 1.4) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      labs(y = "Net C balance", x = NULL) + theme_custom
    
    p6 <- ggplot(data_exudates, aes(x = day)) +
      geom_line(aes(y = c_re), color = "#1f77b4", linewidth = 1.3) +
      geom_line(aes(y = n_from_exudates * 50), color = "#d62728", linewidth = 1.1, linetype = "dashed") +
      scale_y_continuous(name = "C input", sec.axis = sec_axis(~./50, name = "N input (scaled)")) +
      labs(x = NULL) + theme_custom
    
    p7 <- ggplot(data_diagnostic, aes(day, cn_mo)) +
      geom_line(color = "#9467bd", linewidth = 1.3) +
      geom_hline(yintercept = 22, color = "gray50", linetype = "dotted") +
      labs(y = "Microbial C:N", x = NULL) + theme_custom
    
    p8 <- ggplot(data_diagnostic, aes(day, cn_re)) +
      geom_line(color = "#d62728", linewidth = 1.2, linetype = "dashed") +
      geom_hline(yintercept = 35, color = "gray50", linetype = "dotted") +
      labs(y = "Exudate C:N", x = "Day") + theme_custom
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 2, heights = c(1,1,1,1))
  })
  
  # ---- Downloads for single run ----
  output$dl_sim_csv <- downloadHandler(
    filename = function() paste0("RPE_simulation_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$sim)
      write.csv(rv$sim, file, row.names = FALSE)
    }
  )
  
  output$dl_rpe_pdf <- downloadHandler(
    filename = function() paste0("RPE_plot_", Sys.Date(), ".pdf"),
    content = function(file) {
      req(rv$sim)
      pdf(file, width = 10, height = 6)
      print(make_rpe_plot(rv$sim) + 
              labs(title = paste0("Net RPE — c_re_mean=", input$c_re_mean,
                                  ", cn_re_mean=", input$cn_re_mean)))
      dev.off()
    }
  )
  
  # ============================================================
  # SCENARIO COMPARISON TAB
  # ============================================================
  compare_results <- reactiveVal(NULL)
  
  observeEvent(input$run_compare, {
    res_list <- list()
    metrics_list <- list()
    
    for (nm in names(presets)) {
      p <- presets[[nm]]
      res <- run_full_simulation(p$c_re_mean, p$cn_re_mean, p$priming_strength,
                                 p$Km, p$Km_plant, sim_days = 300)
      res$sim$Scenario <- p$name
      res_list[[nm]] <- res$sim
      metrics_list[[nm]] <- cbind(Scenario = p$name, res$metrics)
    }
    
    all_sims <- bind_rows(res_list)
    all_metrics <- bind_rows(metrics_list)
    
    compare_results(list(sims = all_sims, metrics = all_metrics))
  })
  
  output$plot_compare_2x2 <- renderPlot({
    req(compare_results())
    df <- compare_results()$sims
    
    ggplot(df, aes(day, net_rpe, color = Scenario)) +
      geom_line(linewidth = 1.1, alpha = 0.85) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      facet_wrap(~ Scenario, ncol = 2, scales = "free_y") +
      scale_color_manual(values = c("#2ca02c", "#1f77b4", "#ff7f0e", "#d62728")) +
      labs(title = "Net RPE Time Series — Four Targeted Scenarios",
           x = "Day", y = "Net RPE") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold", size = 12))
  })
  
  output$compare_table <- renderTable({
    req(compare_results())
    compare_results()$metrics
  }, striped = TRUE, hover = TRUE)
  
  # ============================================================
  # SENSITIVITY TAB (Km × Km_plant grid)
  # ============================================================
  sens_results <- reactiveVal(NULL)
  
  observeEvent(input$run_sens, {
    base <- presets[[input$sens_preset]]
    
    km_seq  <- seq(input$sens_km_min,  input$sens_km_max,  by = input$sens_step)
    kmp_seq <- seq(input$sens_kmp_min, input$sens_kmp_max, by = input$sens_step)
    
    grid <- expand.grid(Km = km_seq, Km_plant = kmp_seq)
    
    sens_df <- list()
    
    for (i in seq_len(nrow(grid))) {
      par <- grid[i, ]
      res <- run_full_simulation(
        c_re_mean = base$c_re_mean,
        cn_re_mean = base$cn_re_mean,
        priming_strength = base$priming_strength,
        Km = par$Km,
        Km_plant = par$Km_plant,
        sim_days = 300
      )
      sens_df[[i]] <- data.frame(
        Km = par$Km,
        Km_plant = par$Km_plant,
        mean_RPE = res$metrics$Mean_RPE,
        min_RPE = res$metrics$Min_RPE,
        max_RPE = res$metrics$Max_RPE,
        final_C_balance = res$metrics$Final_Net_C_balance
      )
    }
    
    sens_summary <- bind_rows(sens_df)
    sens_results(sens_summary)
  })
  
  output$sens_heatmap_rpe <- renderPlot({
    req(sens_results())
    df <- sens_results()
    
    ggplot(df, aes(Km, Km_plant, fill = mean_RPE)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = 0, name = "Mean Net RPE") +
      labs(title = paste("Mean Net RPE —", input$sens_preset),
           subtitle = "Effect of microbial Km and plant Km_plant",
           x = "Km (microbial N-mining)", y = "Km_plant (plant N uptake)") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "right")
  })
  
  output$sens_heatmap_cbal <- renderPlot({
    req(sens_results())
    df <- sens_results()
    
    ggplot(df, aes(Km, Km_plant, fill = final_C_balance)) +
      geom_tile() +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                           midpoint = 0, name = "Final Net C Balance") +
      labs(title = paste("Final Net C Balance —", input$sens_preset),
           subtitle = "Effect of microbial Km and plant Km_plant",
           x = "Km (microbial N-mining)", y = "Km_plant (plant N uptake)") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "right")
  })
  
  output$sens_table <- renderTable({
    req(sens_results())
    sens_results() %>%
      arrange(desc(abs(mean_RPE))) %>%
      head(20)
  }, striped = TRUE)
  
  output$dl_sens_csv <- downloadHandler(
    filename = function() paste0("RPE_sensitivity_", input$sens_preset, "_", Sys.Date(), ".csv"),
    content = function(file) {
      req(sens_results())
      write.csv(sens_results(), file, row.names = FALSE)
    }
  )
  
  # Optional: auto-run a default simulation on startup 
  observe({
    if (is.null(rv$sim)) {
      # Run a default balanced scenario on first load
      p <- presets[["Net RPE Balanced"]]
      res <- run_full_simulation(p$c_re_mean, p$cn_re_mean, p$priming_strength,
                                 p$Km, p$Km_plant)
      rv$sim <- res$sim
      rv$metrics <- res$metrics
      rv$last_params <- p
    }
  })
}

# ============================================================
# RUN THE APP
# ============================================================
shinyApp(ui = ui, server = server)

