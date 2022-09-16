library(shiny)
library(shinyBS)
library(deSolve)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(dplyr)

# Define UI for random distribution application
ui = fluidPage(

  # Application title
  titlePanel("Laboratoire 8 - Compétition interspécifique"),

  sidebarLayout(
    sidebarPanel(
      tags$div(
        tags$p("Entrez directement ou déplacez les curseurs pour modifier les valeurs d'entrée des paramètres.")),

    #  tags$div(tags$p("Paramètres de l'espèce 1")),

      numericInput(inputId = "N1", label = HTML("N<sub>1</sub>(0): Abondance initiale de l'espèce 1 - N1(0)"), min = 1, value = 200),

      numericInput("N2", HTML("N<sub>2</sub>(0): Abondance initiale de l'espèce 2 - N2(0)"), min = 1, value = 200),

      numericInput(inputId = "K1", label = HTML("K<sub>1</sub>: Capacité de support de l'espèce 1"), min = 1, value = 200),

      numericInput("K2", HTML("K<sub>2</sub>: Capacité de support de l'espèce 2"), min = 1, value = 200),

      sliderInput("r1", label = HTML("r<sub>1</sub>: Taux de croissance intrinsèque de l'espèce 1"),
                    value = 0.5,
                    min = 0,
                    max = 2,
                    step = 0.1),

      sliderInput("r2", HTML("r<sub>2</sub>: Taux de croissance intrinsèque de l'espèce 2"),
                    value = 0.5,
                    min = 0,
                    max = 2,
                    step = 0.1),

      sliderInput("alpha", HTML("&alpha;: Effet relatif de l'espèce 2 sur l'espèce 1"),
                  value = 5,
                  min = 0,
                  max = 10,
                  step = 0.05),

      #tags$div(tags$p("Paramètres de l'espèce 2")),
      sliderInput("beta", HTML("&beta;: Effet relatif de l'espèce 1 sur l'espèce 2"),
                  value = 5,
                  min = 0,
                  max = 10,
                  step = 0.05),

      numericInput("time", "Durée de la simulation", min = 1, max=1000, value = 50),
      #width=5 #largeur du slider
    ),

    mainPanel(
      withMathJax(),
      tabsetPanel(type = "pills",
                          tabPanel("Graphiques",
                          h4("Le premier graphique montre le nombre d'individus en fonction du temps."),
                          plotOutput("plot_numbers"),
                          h4("Le deuxième graphique montre l'espace des phases du modèle de compétition de Lotka-Volterra à deux espèces."),
                          plotOutput("phase_portrait"),
                          #h4("Trajectory of the overall population growth rate")
                          ),
                  tabPanel("Tableau",
                           br(),
                           tableOutput("df_table"),
                           )
                         ),
    tags$div(
      tags$p("Les graphiques ci-haut montrent le nombre d'individus de l'espèce 1 et de l'espèce 2 prédit par le modèle de compétition de Lotka-Volterra:
             $$\\begin{align}
             \\frac{\\textrm{dN1}}{\\textrm{dt}} &= r1 N1 (1-\\frac{N1 + \\alpha\\ N2}{K1}) \\\\
             \\frac{\\textrm{dN2}}{\\textrm{dt}} &= r2 N2 (1-\\frac{N2 + \\beta\\ N1}{K2}) \\\\
             \\end{align}$$
             Où \\(alpha\\) est l'effet relatif de la présence de l'espèce 2 sur l'espèce 1, \\(beta\\) l'effet relatif de la présence de l'espèce 1 sur l'espèce 2, \\(K1\\) la capacité de support (ou capacité biotique) de l'espèce 1, \\(K2\\) la capacité de support (ou capacité biotique) de l'espèce 2 et  \\(r1\\) et \\(r2\\) le taux de croissance de l'espèce 1 et 2, respectivement.
             ")
    )
    )
  )
)

#------ Function for Lotka-Volterra system ------------------------------------
lvmodel = function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    # Differential equations
    dN1 <- r1*N1*(1 - (N1 + alpha*N2)/K1)
    dN2 <- r2*N2*(1 - (N2 + beta*N1)/K2)
    list(c(dN1, dN2))
  })
}

#----- Function to calculate solution to Lotka-Volterra system ----------------
# Out for the table in Shiny App
lotka_volterra = function(N1_in, N2_in, K1_in, K2_in, r1_in, r2_in, alpha_in, beta_in, time_limit){
  times = seq(0, time_limit, by = 1)
  parameters = c(K1 = K1_in, K2 = K2_in, r1 = r1_in, r2 = r2_in, alpha = alpha_in, beta=beta_in)
  state = c(N1 = N1_in, N2 = N2_in)
  results = data.frame(ode(y = state, times = times, func = lvmodel,
                           parms = parameters))
  return(results)
}

#----- Function to calculate solution to Lotka-Volterra system ----------------
# Needed to return the value of alpha, beta, K1 and K2 to compute isoclines
lotka_volterra_all = function(N1_in, N2_in, K1_in, K2_in, r1_in, r2_in, alpha_in, beta_in, time_limit){
  times = seq(0, time_limit, by = 1)
  parameters = c(K1 = K1_in, K2 = K2_in, r1 = r1_in, r2 = r2_in, alpha = alpha_in, beta=beta_in)
  state = c(N1 = N1_in, N2 = N2_in)
  results = data.frame(ode(y = state, times = times, func = lvmodel,
                           parms = parameters),K1_in,alpha_in,K2_in,beta_in)
  return(results)
}

server = function(input, output, session) {
  # ------------- get data from slider inputs ---------------------------------
  # Out for the table in Shiny App
  data_table = reactive({
    lotka_volterra(N1_in = input$N1,
                   N2_in = input$N2,
                   K1_in = input$K1,
                   K2_in = input$K2,
                   r1_in = input$r1,
                   r2_in = input$r2,
                   alpha_in=input$alpha,
                   beta_in=input$beta,
                   time_limit = input$time)
  })

  data = reactive({
    lotka_volterra_all(N1_in = input$N1,
                   N2_in = input$N2,
                   K1_in = input$K1,
                   K2_in = input$K2,
                   r1_in = input$r1,
                   r2_in = input$r2,
                   alpha_in=input$alpha,
                   beta_in=input$beta,
                   time_limit = input$time)
  })

  # -------- Generate plot of numbers of predator/prey -----------------------
output$plot_numbers <- renderPlot(width = 650,height = 350,{
    out_long = gather(data(), species, num_anim, N1:N2)
    color_pal <- c("#E41A1C", "#377EB8", "#4DAF4A")
    ggplot(out_long, aes(x = time, y = num_anim, colour = species)) +
      geom_line(size = 2) +
      labs(x = "Temps", y = "Abondance de l'espèce 1 et 2")+
          # title = "Nombre d'individus en fonction du temps")
      #scale_colour_discrete(name = "Espèces", labels = c("Espèce 1", "Espèce 2")) +
      scale_color_manual(values = color_pal,
                         labels = c("1", "2")) +
      labs(colour="Espèces",linetype="Espèces",shape="Espèces")+
      theme(text = element_text(size = 15))
  })

output$df_table <- renderTable({data_table()})
  # -------- Generate phase portrait ------------------------------------------
  output$phase_portrait <- renderPlot(width = 650,height = 350,{
      color_pal <- c("#E41A1C", "#377EB8", "#4DAF4A")
    ggplot(data(), aes(x = N1, y = N2)) +
      geom_point() +
      labs(x = "Nombre de l'espèce 1", y = "Nombre de l'espèce 2")+
          # ,title = "Phase Plot of competition Model")
      scale_color_manual(values = color_pal,
                         labels = c("1", "2")) +
     #geom_path(aes(x = N1, y = N2), size = 1) +
     #isoclines
     geom_segment(data(), mapping=aes(x = 0, y = K1_in/alpha_in, xend = K1_in, yend = 0, color = "Species 1"), size = 2) +
     geom_segment(data(), mapping=aes(x = 0, y = K2_in, xend = K2_in/beta_in, yend = 0, color = "Species 2"), size = 2) +
     geom_point(x = first(data()$N1), y = first(data()$N2),
                pch = 21, size = 3, fill = "black") +
     geom_point(x = last(data()$N1), y = last(data()$N2),
                pch = 21, size = 3, fill = "white", stroke = 1) +
     labs(colour="Espèces",linetype="Espèces",shape="Espèces")+
      theme(text = element_text(size = 15))
  })

}
shinyApp(ui = ui, server = server)
