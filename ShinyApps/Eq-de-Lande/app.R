#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(MASS)

if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Visualizing the Lande Equation"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("corr",
                        "Genetic correlation:",
                        min = -1,
                        max = 1,
                        step = 0.01,
                        value = 0.5)
            ,
            sliderInput("beta_norm",
                        "Strength of selection:",
                        min = 0,
                        max = 2,
                        step = 0.01,
                        value = 0)
            ,
            sliderInput("beta_angle",
                        "Direction of selection:",
                        min = 0,
                        max = 360,
                        step = 1,
                        value = 60),
            checkboxInput("show_response", "Show response vector (delta Z)", FALSE),
            checkboxInput("next_gen", "Show next generation.", FALSE)
            
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        G = matrix(c(1, input$corr, input$corr, 1), 2, 2)
        set.seed(input$corr*100)
        x <- mvrnorm(n = 100, Sigma = G, mu = c(10, 10))
        df = data.frame(x = x[,1], y = x[,2])
        # draw the histogram with the specified number of bins
        p = ggplot(df, aes(x, y)) + geom_point() + theme_cowplot() + labs(x = "Trait x", y = "Trait y")  + stat_ellipse() + coord_fixed()
        if(input$beta_norm > 0){
            beta_x = cos(input$beta_angle / 180 * pi) * input$beta_norm
            beta_y = sin(input$beta_angle / 180 * pi) * input$beta_norm
            beta = c(beta_x, beta_y)
            dZ = G %*% beta
            nx = mvrnorm(n = 100, Sigma = G, mu = c(10, 10) + dZ)
            n_df = data.frame(x = nx[,1], y = nx[,2])
            
            p = p + geom_segment(aes(x = 10, y = 10, xend = 10 + beta_x, yend = 10 + beta_y), size = 1, color = "tomato3", arrow = arrow(length = unit(0.5, "cm")))
 
            if(input$show_response){
                p = p + geom_segment(aes(x = 10, y = 10, xend = 10 + dZ[1], yend = 10 + dZ[2]), 
                                     size = 1, color = "blue", arrow = arrow(length = unit(0.5, "cm")))
            }
            if(input$next_gen){
                p = p + geom_point(data = n_df, color = "orange") + stat_ellipse(data = n_df, color = "orange")
            }
        }
        p
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
