library(shiny)
library(ggplot2)
library(cowplot)
library(wesanderson)
library(mvtnorm)
library(matrixStats)
library(MASS)
library(bslib) 
library(shinycssloaders)
library(ContourFunctions)

Norm = function(x) sqrt(sum(x^2))

Normalize = function(x) x / Norm(x)

gplotW_bar_trajectory <-
  function(run, space_size = 6, xlimits = c(-space_size, 
                                            space_size), 
           ylimits = c(-space_size, space_size), resolution = 0.2,
           mypalette = colorRampPalette(c("white", 
                                          wes_palette(10, 
                                                      name = "Zissou1", 
                                                      type = "continuous"), 
                                          "darkred")),
           log = FALSE, main = "", ...){
    
    W_bar = W_bar_factory(run$theta)
    x <- seq(xlimits[1], xlimits[2], resolution)
    y <- seq(ylimits[1], ylimits[2], resolution)
    X <- as.matrix(expand.grid(x, y))
    Z <- vector()
    
    for(i in 1:nrow(X)){
      Z[i] <- W_bar(c(X[i,1], X[i,2]))
    }
    
    if(log) { Z = Z - logSumExp(Z)
    } else Z = exp(Z - logSumExp(Z))
    
    # 1. Format peak data
    theta_df <- data.frame(run$theta)
    colnames(theta_df) <- c("X1", "X2")
    
    # 2. Format trajectory data into segments (z_t to z_{t+1})
    traj_raw <- data.frame(run$trajectory)
    n_gens <- nrow(traj_raw)
    
    # Calculate the start and absolute end of each step
    x_start <- traj_raw[1:(n_gens-1), 1]
    y_start <- traj_raw[1:(n_gens-1), 2]
    x_end_full <- traj_raw[2:n_gens, 1]
    y_end_full <- traj_raw[2:n_gens, 2]
    
    # 3. Create a gap by scaling back the end coordinates 
    # (0.85 means the arrow covers 85% of the distance, leaving a 15% gap)
    arrow_scale <- 0.85 
    
    seg_df <- data.frame(
      x = x_start,
      y = y_start,
      xend = x_start + (x_end_full - x_start) * arrow_scale,
      yend = y_start + (y_end_full - y_start) * arrow_scale
    )
    
    gcf_grid(x, y, Z, xlim = xlimits, ylim = ylimits, 
             color.palette = mypalette, mainminmax = FALSE, 
             mainminmax_minmax = FALSE, ...) +
      
      # Dashed origin axes
      geom_segment(aes(x = 0, xend = 0, y = ylimits[1], yend = ylimits[2]), 
                   color = "gray50", linetype = "dashed", alpha = 0.6) + 
      geom_segment(aes(y = 0, yend = 0, x = xlimits[1], xend = xlimits[2]), 
                   color = "gray50", linetype = "dashed", alpha = 0.6) +
      
      # Adaptive Peaks
      geom_point(data = theta_df, aes(x = X1, y = X2), 
                 shape = 24, fill = "#F1C40F", color = "black", size = 3, stroke = 1) +
      
      # Generational Arrows with Gaps
      geom_segment(data = seg_df, aes(x = x, y = y, xend = xend, yend = yend),
                   color = "#2C3E50", linewidth = 0.6, 
                   arrow = arrow(type = "closed", length = unit(0.06, "inches")),
                   linejoin = "mitre") +
      
      # Starting position dot
      geom_point(data = traj_raw[1, , drop = FALSE], aes(x = X1, y = X2), 
                 color = "#2C3E50", size = 2.5) +
      
      ggtitle(main) + 
      coord_fixed() + 
      theme_void() + 
      theme(legend.position = "none")
  }

diff_cut_off = 1e-4
max_gens = 10000
max_stand_still = 100
space_size = 6

mypalette = colorRampPalette(c(wes_palette(10, name = "Zissou1", type = "continuous"), "darkred"))(50)

vector_cor = function(x, y) abs(x %*% y/(Norm(x)*Norm(y)))

W_bar_factory = function(theta_matrix, w_cov = diag(dim(theta_matrix)[2])) {
  function(x) logSumExp(apply(theta_matrix, 1, function(theta) dmvnorm(x, mean = theta, w_cov, log = T)))
}

W_bar_gradient_factory = function(theta_matrix, w_cov = NULL){
  if(is.null(w_cov)){
    function(x) rowSums(apply(theta_matrix, 1, function(theta) - dmvnorm(x, mean = theta) * t(x - theta)))/exp(W_bar_factory(theta_matrix)(x))
  } else{
    function(x) rowSums(apply(theta_matrix, 1, function(theta) - dmvnorm(x, mean = theta, w_cov) * solve(w_cov, x - theta)))/exp(W_bar_factory(theta_matrix, w_cov)(x))
  }
}

randomPeaks = function(n = n_peaks, p = n_traits, x = rep(1, p), intervals = 1, prop = 1, dz_limits,
                       max_uniform = n * 100, sigma_init = 2, sigma_step = 0.01, verbose = FALSE){
  steps = length(intervals)
  counter = vector("numeric", steps)
  n_per = ceiling(n * prop)
  peaks = matrix(0, n, p)
  k = 1
  attempts = 1
  while(k <= n & attempts < max_uniform){
    attempts = attempts + 1
    rpeak = Normalize(rnorm(p))
    corr = vector_cor(x, rpeak)
    for(i in 1:steps) {
      if(corr < intervals[i]){
        if(counter[i] < n_per[i]){
          counter[i] = counter[i] + 1
          if(verbose) print(counter)
          peaks[k,] = rpeak * runif(1, dz_limits[1], dz_limits[2])
          k = k + 1
        }
        break
      }
    }
  }
  if(k < n){
    mask = which(counter != n_per)
    mask = c(mask[1]-1, mask)
    target_intervals = intervals[mask]
    sigma = sigma_init
    while(k <= n){
      rpeak = Normalize(x + rnorm(p, 0, sigma))
      corr = vector_cor(x, rpeak)
      if(corr < target_intervals[1]){
        sigma = sigma - sigma_step
      } else if(corr > target_intervals[length(target_intervals)]) {
        sigma = sigma + sigma_step
      } else { for(i in 1:steps) {
        if(corr < intervals[i]){
          if(counter[i] < n_per[i]){
            counter[i] = counter[i] + 1
            if(verbose) print(counter)
            peaks[k,] = rpeak * runif(1, dz_limits[1], dz_limits[2])
            k = k + 1
            mask = which(counter != n_per)
            mask = c(mask[1]-1, mask)
            target_intervals = intervals[mask]
          }
          break
        }
      }
      }
    }
  }
  peaks[sample(1:n, n),]
}

calculateTrajectory <- function (start_position, G, W_bar, W_bar_grad, scale = 2) {
  p = dim(G)[1]
  trajectory = matrix(NA, max_gens, p)
  betas = matrix(NA, max_gens, p)
  current_position = start_position
  stand_still_counter = 0
  net_beta = rep(0, p)
  gen = 1
  while(gen <= max_gens){
    trajectory[gen,] = current_position
    beta = W_bar_grad(as.vector(current_position))
    betas[gen,] = beta
    net_beta = net_beta + beta
    next_position = current_position + (G/scale)%*%beta
    if(Norm(next_position) > space_size*2) stop("Out of bounds")
    if(Norm(next_position - current_position) < diff_cut_off){
      stand_still_counter = stand_still_counter + 1
    }
    if(stand_still_counter > max_stand_still){
      break
    }
    current_position = next_position
    gen = gen+1
  }
  trajectory = unique(trajectory[!is.na(trajectory[,1]),])
  betas = betas[!is.na(betas[,1]),]
  net_dz = trajectory[dim(trajectory)[1],] - start_position
  return(list(start_position = start_position,
              trajectory = trajectory,
              betas = betas,
              net_beta = net_beta,
              net_dz = net_dz))
}

runSimulation = function(G_type = c("Diagonal", "Integrated"), G = NULL,
                         n_peaks = 1, p, rho = 0.7, scale = 6, peakPool = NULL, theta = NULL){
  G_type = match.arg(G_type)
  if(is.null(G)){
    if(G_type == "Diagonal"){
      G = G_factory(p, 0.1)
      gmax = eigen(G)$vectors[,1]
    } else if(G_type == "Integrated"){
      G = G_factory(p, rho)
      gmax = eigen(G)$vectors[,1]
    } else stop("Unknown G type")
  } else
    gmax = eigen(G)$vectors[,1]
  if(n_peaks == 1){
    Surface_type = "Single"
  } else
    Surface_type = "Multiple"
  if(is.null(theta)) theta = matrix(peakPool[sample(1:nrow(peakPool), n_peaks),], n_peaks, p)
  W_bar = W_bar_factory(theta)
  W_bar_grad = W_bar_gradient_factory(theta)
  trajectory = calculateTrajectory(rep(0, p), G, W_bar, W_bar_grad, scale = scale)
  trajectory$G_type = G_type
  trajectory$G = G
  trajectory$gmax = gmax
  trajectory$theta = theta
  trajectory$z = trajectory$trajectory[dim(trajectory$trajectory)[1],]
  trajectory$W_bar = W_bar
  trajectory$W_bar_grad = W_bar_grad
  trajectory$Surface_type = Surface_type
  trajectory$normZ = Norm(trajectory$z)
  return(trajectory)
}

G_factory = function(p, rho, sigma = 0.1){
  while(TRUE){
    G = matrix(rnorm(p*p, rho, sigma), p, p)
    G = (G + t(G))/2
    diag(G) = rnorm(p, 1, sigma)
    tryCatch({chol(G); break}, error = function(x) FALSE)
  }
  G
}

space_size <- 6
diff_cut_off <- 1e-4
max_gens <- 10000
max_stand_still <- 100

ui <- page_sidebar(
  title = "Evolutionary Trajectories on Adaptive Surfaces",
  theme = bs_theme(version = 5, bootswatch = "flatly"), # Modern, clean theme
  
  sidebar = sidebar(
    sliderInput("seed", "Simulation random seed:", min = 0, max = 100, step = 1, value = 42),
    sliderInput("n_peaks", "Number of adaptive peaks:", min = 1, max = 10, step = 1, value = 1),
    hr(),
    sliderInput("corr_1", "Genetic correlation (Pop 1):", min = -1, max = 1, step = 0.05, value = 0.8),
    sliderInput("corr_2", "Genetic correlation (Pop 2):", min = -1, max = 1, step = 0.05, value = 0.0)
  ),
  
  # Main panel 
  card(
    full_screen = TRUE,
    plotOutput("distPlot", height = "600px") |> shinycssloaders::withSpinner() # Optional: add shinycssloaders::withSpinner()
  )
)

server <- function(input, output, session) {
  
  # 1. REACTIVE: Generate the landscape (Peaks and Theta). 
  # This ONLY runs when 'seed' or 'n_peaks' changes.
  landscape_data <- reactive({
    set.seed(input$seed * input$n_peaks)
    p <- 2
    peakPool <- randomPeaks(100, p = p, dz_limits = c(3, space_size), 
                            intervals = c(1), prop = c(1))
    theta <- matrix(peakPool[sample(1:nrow(peakPool), input$n_peaks),], input$n_peaks, p)
    
    return(theta)
  })
  
  # 2. REACTIVE: Simulate Population 1. 
  # ONLY runs when the landscape changes or 'corr_1' changes.
  pop1_sim <- reactive({
    theta <- landscape_data()
    set.seed(input$seed * input$n_peaks)
    
    # Ellipse Data
    G1 <- matrix(c(1, input$corr_1, input$corr_1, 1), 2, 2)
    x1 <- mvrnorm(n = 100, Sigma = G1, mu = c(10, 10))
    df1 <- data.frame(x = x1[,1], y = x1[,2])
    p_xG <- ggplot(df1, aes(x, y)) + geom_point(alpha=0.5) + stat_ellipse(color="blue") + 
      coord_fixed() + theme_void()
    
    # Trajectory
    sim <- runSimulation("Integrated", rho = input$corr_1, n_peaks = input$n_peaks, 
                         p = 2, scale = 4, theta = theta)
    
    # Generate plot (assuming gplotW_bar_trajectory is defined)
    p_main <- gplotW_bar_trajectory(sim, space_size = 8)
    
    # Combine
    ggdraw(p_main) + draw_plot(p_xG, .6, .8, .2, .2)
  })
  
  # 3. REACTIVE: Simulate Population 2. 
  # ONLY runs when the landscape changes or 'corr_2' changes.
  pop2_sim <- reactive({
    theta <- landscape_data()
    set.seed(input$seed * input$n_peaks)
    
    # Ellipse Data
    G2 <- matrix(c(1, input$corr_2, input$corr_2, 1), 2, 2)
    x2 <- mvrnorm(n = 100, Sigma = G2, mu = c(10, 10))
    df2 <- data.frame(x = x2[,1], y = x2[,2])
    p_yG <- ggplot(df2, aes(x, y)) + geom_point(alpha=0.5) + stat_ellipse(color="red") + 
      coord_fixed() + theme_void()
    
    # Trajectory
    sim <- runSimulation("Integrated", rho = input$corr_2, n_peaks = input$n_peaks, 
                         p = 2, scale = 4, theta = theta)
    
    # Generate plot
    p_main <- gplotW_bar_trajectory(sim, space_size = 8)
    
    # Combine
    ggdraw(p_main) + draw_plot(p_yG, .6, .8, .2, .2)
  })
  
  # 4. RENDER: Combine the two plots
  output$distPlot <- renderPlot({
    plot_grid(pop1_sim(), pop2_sim(), labels = c("Population 1", "Population 2"))
  })
}

shinyApp(ui, server)
