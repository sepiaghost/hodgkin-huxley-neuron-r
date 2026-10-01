# Interactive Hodgkin-Huxley neuron
# Run with: shiny::runApp()   (from this folder)
#
# Move the stimulus current slider and watch the neuron cross its threshold and
# fire. This is a rebuild, in R, of an Excel model with a stimulus-current bar.

library(shiny)
library(ggplot2)
source("hodgkin_huxley.R")

stim_start <- 10   # ms, when the current step begins

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: Georgia, 'Times New Roman', serif; color: #1E2A2A; }
    h2 { font-weight: 600; margin-bottom: 4px; }
    .lead { color: #4A5858; margin-bottom: 20px; }
    .readout { font-size: 15px; line-height: 1.7; margin-top: 18px; }
    .readout b { color: #0B6E6E; }
  "))),

  titlePanel("Hodgkin-Huxley neuron"),
  p(class = "lead",
    "Change the stimulus current and see when the neuron fires. ",
    "Below a threshold nothing happens; above it, a full spike."),

  sidebarLayout(
    sidebarPanel(
      sliderInput("I0", "Stimulus current (\u00b5A/cm\u00b2)",
                  min = 0, max = 20, value = 2, step = 0.1),
      sliderInput("duration", "Stimulus duration (ms)",
                  min = 1, max = 50, value = 30, step = 1),
      div(class = "readout", htmlOutput("readout")),
      width = 3
    ),
    mainPanel(
      plotOutput("voltage_plot", height = "320px"),
      plotOutput("gate_plot", height = "240px"),
      width = 9
    )
  )
)

server <- function(input, output, session) {

  sim <- reactive({
    simulate_hh(I0 = input$I0, t_on = stim_start,
                t_off = stim_start + input$duration, t_end = 80)
  })

  output$voltage_plot <- renderPlot({
    s <- sim()
    ggplot(s, aes(time, V)) +
      annotate("rect", xmin = stim_start, xmax = stim_start + input$duration,
               ymin = -Inf, ymax = Inf, fill = "#F2C572", alpha = 0.25) +
      geom_hline(yintercept = -65, colour = "grey60", linetype = "dashed") +
      geom_line(colour = "#0B6E6E", linewidth = 0.9) +
      coord_cartesian(ylim = c(-85, 50)) +
      labs(x = NULL, y = "Membrane potential (mV)",
           title = "Membrane potential",
           subtitle = "Shaded band = stimulus. Dashed line = resting potential.") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })

  output$gate_plot <- renderPlot({
    s <- sim()
    long <- data.frame(
      time = rep(s$time, 3),
      value = c(s$m, s$h, s$n),
      gate = rep(c("m  (Na+ activation)", "h  (Na+ inactivation)",
                   "n  (K+ activation)"), each = nrow(s))
    )
    ggplot(long, aes(time, value, colour = gate)) +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = c("#0B6E6E", "#B5462F", "#3B5BA5")) +
      labs(x = "Time (ms)", y = "Open probability", colour = NULL,
           title = "Ion channel gates") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  })

  output$readout <- renderUI({
    s <- sim()
    n_spikes <- count_spikes(s)
    HTML(sprintf(
      "Spikes: <b>%d</b><br>Peak potential: <b>%.1f mV</b>",
      n_spikes, max(s$V)
    ))
  })
}

shinyApp(ui, server)
