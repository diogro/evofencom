#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
if(!require(ContourFunctions)){install.packages("ContourFunctions"); library(ContourFunctions)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
theme_set(theme_cowplot())
if(!require(wesanderson)) { install.packages("wesanderson"); library(wesanderson) }
#if(!require(grDevices)) { install.packages("grDevices"); library(grDevices) }

if(!require(mvtnorm)) { install.packages("mvtnorm"); library(mvtnorm) }
if(!require(matrixStats)){install.packages("matrixStats"); library(matrixStats)}
if(!require(MASS)){install.packages("MASS"); library(MASS)}

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
    gcf_grid(x, y, Z, xlim = xlimits, ylim = ylimits, 
             color.palette = mypalette, mainminmax = FALSE, 
             mainminmax_minmax = FALSE, ...) +
        geom_point(data=data.frame(run$theta), aes(X1, X2), shape = 17) +
        geom_point(data=data.frame(run$trajectory), aes(X1, X2), shape = 19) +
        ggtitle(main) + coord_fixed() + theme_void() + 
        theme(legend.position = "none") +
        geom_segment(aes(x = 0, xend = 0, y = ylimits[1], yend = ylimits[2])) + 
        geom_segment(aes(y = 0, yend = 0, x = xlimits[1], xend = xlimits[2]))
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


# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    output$distPlot <- renderPlot({
        p = 2
        set.seed(input$seed * input$n_peaks)
        peakPool = randomPeaks(100, p = 2, dz_limits = c(3, space_size),
                               intervals = c(1), prop = c(1))
        theta = matrix(peakPool[sample(1:nrow(peakPool), input$n_peaks),], input$n_peaks, p)

        set.seed(input$seed * input$n_peaks)
        G1 = matrix(c(1, input$corr_1, input$corr_1, 1), 2, 2)
        x1 <- mvrnorm(n = 100, Sigma = G1, mu = c(10, 10))
        df1 = data.frame(x = x1[,1], y = x1[,2])
        p_xG = ggplot(df1, aes(x, y)) + geom_point() + stat_ellipse() + coord_fixed() + theme_nothing()

        set.seed(input$seed * input$n_peaks)
        G2 = matrix(c(1, input$corr_2, input$corr_2, 1), 2, 2)
        x2 <- mvrnorm(n = 100, Sigma = G2, mu = c(10, 10))
        df2 = data.frame(x = x2[,1], y = x2[,2])
        p_yG = ggplot(df2, aes(x, y)) + geom_point() + stat_ellipse() + coord_fixed() + theme_nothing()

        set.seed(input$seed * input$n_peaks)
        x = runSimulation("Integrated", rho = input$corr_1, n_peaks = input$n_peaks, p = 2, scale = 4, theta = theta)
        p_x = gplotW_bar_trajectory(x, 8)
        p_x = ggdraw(p_x) +  draw_plot(p_xG, .6, .8, .2, .2)
        y = runSimulation("Integrated", rho = input$corr_2, n_peaks = input$n_peaks, p = 2, scale = 4, theta = theta)
        p_y = gplotW_bar_trajectory(y, 8)
        p_y = ggdraw(p_y) +  draw_plot(p_yG, .6, .8, .2, .2)
        p_xy = plot_grid(p_x, p_y, labels = c("Population 1", "Population 2"))
        p_xy
    })

})

