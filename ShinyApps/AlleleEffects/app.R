library(shiny)
library(ggplot2)



# --- Shiny UI ---

ui <- fluidPage(
  titlePanel("Quantitative Genetics: Allele Substitution & Variance"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("p", "Allele Frequency (p):", min = 0.01, max = 0.99, value = 0.5, step = 0.01),
      sliderInput("a", "Additive Effect (a):", min = 0, max = 5, value = 1, step = 0.1),
      sliderInput("d", "Dominance Effect (d):", min = -5, max = 5, value = 0, step = 0.1),
      sliderInput("sigma", "Environmental Standard Deviation (sigma):", min = 0.05, max = 2, value = 0.25, step = 0.05),
      sliderInput("N", "Population Size (N):", min = 100, max = 5000, value = 1000, step = 100)
    ),
    
    mainPanel(
      plotOutput("distPlot"),
      hr(),
      plotOutput("effectPlot")
    )
  )
)

# --- Shiny Server ---

server <- function(input, output) {
  # --- Define the Functions ---
  
  plot_allele_contribution <- function(p, a, d, sigma = 0.25, N = 5000){
    
    freq_0    <- (1 - p)^2
    freq_half <- 2 * p * (1 - p)
    freq_1    <- p^2
    
    n_0    <- round(N * freq_0)
    n_half <- round(N * freq_half)
    n_1    <- N - n_0 - n_half 
    
    sub_effect = a + d*((1-p)-p)
    
    genotypes <- c(rep("0", n_0), rep("1/2", n_half), rep("1", n_1))
    mean_values <- c(rep(-a, n_0), rep(d, n_half), rep(a, n_1))
    phenotypes <- rnorm(N, mean = mean_values, sd = sigma)
    
    df <- data.frame(Genotype = factor(genotypes, levels = c("0", "1/2", "1")), 
                     Phenotype = phenotypes)

    Va = 2*p*(1-p)*sub_effect^2
    Vd = (2*p*(1-p)*d)^2 
    Ve = var(df$Phenotype) - Va - Vd
    
    dist_0    <- function(x) dnorm(x, mean = -a, sd = sigma) * freq_0
    dist_half <- function(x) dnorm(x, mean = d,  sd = sigma) * freq_half
    dist_1    <- function(x) dnorm(x, mean = a,  sd = sigma) * freq_1
    
    up = max(c(2, a, d))
    lw = min(c(-a, d, -2))
    
    ggplot(df, aes(x = Phenotype)) +
      geom_histogram(aes(y = after_stat(density)), bins = 40, 
                     fill = "grey85", color = "grey60", alpha = 0.7) +
      stat_function(fun = dist_0, aes(color = "0"), linewidth = 1.2, alpha = 0.7) +
      stat_function(fun = dist_half, aes(color = "1/2"), linewidth = 1.2, alpha = 0.7) +
      stat_function(fun = dist_1, aes(color = "1"), linewidth = 1.2, alpha = 0.7) +
      scale_x_continuous(limits = c(lw - 4*sigma, up + 4*sigma)) + 
      scale_color_manual(values = c("0" = "#E41A1C", "1/2" = "#4DAF4A", "1" = "#377EB8"), name = "Genotype") +
      labs(
        title = "Locus Contribution to Phenotypic Variance",
        subtitle = bquote("Parameters:" ~ p == .(p) * "," ~ a == .(a) * "," ~ d == .(d) * "," ~ sigma == .(sigma) * "," ~ Va == .(Va |> round(3)) * "," ~ Vd == .(Vd |> round(3)) * "," ~ Ve == .(Ve |> round(3))),
        x = "Phenotypic Value",
        y = "Probability Density"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", size = 14), axis.title = element_text(face = "bold"))
  }
  
  plot_average_effect <- function(p, a, d, sigma = 0.25, N = 1000){
    freq_0    <- (1 - p)^2
    freq_half <- 2 * p * (1 - p)
    freq_1    <- p^2
    
    n_0    <- round(N * freq_0)
    n_half <- round(N * freq_half)
    n_1    <- N - n_0 - n_half 
    
    sub_effect = a + d*((1-p)-p)
    
    genotypes <- c(rep(0, n_0), rep(0.5, n_half), rep(1, n_1))
    mean_values <- c(rep(-a, n_0), rep(d, n_half), rep(a, n_1))
    phenotypes <- rnorm(N, mean = mean_values, sd = sigma)
    
    df <- data.frame(Genotype = genotypes, Phenotype = phenotypes)
    
    up = max(a, d)
    lw = min(-a, d)
    
    ggplot(df, aes(x = Genotype, y = Phenotype)) +
      geom_jitter(width = 0.02, alpha = 0.4, color = "steelblue") +
      stat_summary(fun = mean, geom = "point", color = "darkred", size = 4, shape = 18) +
      geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 1) +
      scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "1/2", "1")) +
      scale_y_continuous(limits = c(lw - 3*sigma, up + 3*sigma)) + 
      labs(
        title = "Average Effect of Allele Substitution",
        subtitle = bquote("Parameters:" ~ p == .(p) * "," ~ a == .(a) * "," ~ d == .(d) * "," ~ sigma == .(sigma) * "  |" ~ alpha == .(sub_effect)),
        x = "Within-Individual Allele Frequency",
        y = "Phenotype Value"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", size = 14), axis.title = element_text(face = "bold"))
  }
  
  output$distPlot <- renderPlot({
    plot_allele_contribution(p = input$p, a = input$a, d = input$d, sigma = input$sigma, N = input$N)
  })
  
  output$effectPlot <- renderPlot({
    plot_average_effect(p = input$p, a = input$a, d = input$d, sigma = input$sigma, N = input$N)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)