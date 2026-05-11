if(!require(rsconnect)){install.packages("rsconnect"); library(rsconnect)}
rsconnect::setAccountInfo(name='damelo', 
                          token='95346A5518D36EC95824A1AC86E965D1', 
                          secret='SecretKey')

# Deploy
if(!require(here)){install.packages("here"); library(here)}
rsconnect::deployApp(here("ShinyApps/Eq-de-Lande"))
rsconnect::deployApp(here("ShinyApps/Superficie-Adaptativa"))
rsconnect::deployApp(here("ShinyApps/AlleleEffects"))

# Run local
if(!require(shiny)){install.packages("shiny"); library(shiny)}
runGitHub("diogro", "evofencom", subdir = "ShinyApps/Eq-de-Lande")
runGitHub("diogro", "evofencom", subdir = "ShinyApps/Superficie-Adaptativa")

