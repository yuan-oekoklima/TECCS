
library(shiny)
library(markdown)
library(data.table)
library(stringr)
library(DT)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(rdwd)
rdwd::updateRdwd() # to update the latest data if rdwd package doesn't work properly


df.station <- fread("Station_DWD.csv")
df.station$startyear <- as.numeric(substr(df.station$von_datum, 1, 4))
df.station$endyear <- as.numeric(substr(df.station$bis_datum, 1, 4))
df.station$choice <- paste0(df.station$Stationsname, " (", df.station$Stations_id, ")", sep = "")
choice.station <- df.station$choice[df.station$Bundesland == "Bayern" &
                                      df.station$startyear <= 2010 &
                                      df.station$endyear >= 2019]


ui <- navbarPage("Twig Experiment Climate Change Simulator TECCS",
        tabPanel("Plot",
          sidebarLayout(
            sidebarPanel(
              fileInput("FileUpload", 
                        "Upload your experimental data (.csv / .txt)", 
                        accept = c(".csv", ".txt")),
              hr(),
              
              h2("Bud Development Model"),
              numericInput("Tinside", 
                           label = "Estimate the mean temperature (°C) inside your home",
                           value = 20),
              selectizeInput("ClimateStation",
                             label = "Select your nearest available climate station",
                             choices = choice.station),
              numericInput("ThresholdChilling", 
                          label = "Choose the chilling temperature threshold (°C)",
                          value = 5),
              numericInput("ThresholdForcing", 
                          label = "Choose the forcing temperature threshold (°C)",
                          value = 5),
              hr(),
              
              h2("Climate Change Simulator"),
              numericInput("SimuDate", 
                        label = "Choose the base year for your climate change simulation (i.e. the year of the corresponding autumn / winter before the simulated spring)",
                        value = 2015),
              sliderInput("WinterScenario", 
                          label = "Choose the winter warming scenario",
                          min = -1, max = 5,
                          value = 0, step = 0.5, post = " °C", 
                          animate = animationOptions(interval = 2000, loop = TRUE)),
              sliderInput("SpringScenario", 
                          label = "Choose the spring warming scenario",
                          min = -1, max = 5,
                          value = 0, step = 0.5, post = " °C", 
                          animate = animationOptions(interval = 2000, loop = TRUE)),
              width = 3
            ),
            mainPanel(
              plotOutput("ChillForcePlot")
            )
          )
        ),
        
        tabPanel("Literature-based model",
          sidebarLayout(
            sidebarPanel(
              helpText("Hazel - Menzel et al. (2020)"),
              helpText("BBCH: 7"),
              helpText("Starting date of experiment: 2015-11-21"),
              helpText("Temperature: recording from own instrument"),
              helpText("Climate station: Landshut-Reithof (13710)"),
              helpText("Chilling/Forcing temperature threshold: 5 °C"),
              numericInput("SimuDateModel", 
                           label = "Choose the base year for your climate change simulation (i.e. the year of the corresponding autumn / winter before the simulated spring)",
                           value = 2015),
              sliderInput("WinterScenarioModel", 
                          label = "Choose the winter warming scenario",
                          min = -1, max = 5,
                          value = 0, step = 0.5, post = " °C", 
                          animate = animationOptions(interval = 2000, loop = TRUE)),
              sliderInput("SpringScenarioModel", 
                          label = "Choose the spring warming scenario",
                          min = -1, max = 5,
                          value = 0, step = 0.5, post = " °C", 
                          animate = animationOptions(interval = 2000, loop = TRUE)),
              width = 3
            ),
            mainPanel(
              plotOutput("ChillForceModelPlot")
            )
          )
        ),
        
        tabPanel("About", 
                 includeHTML("About.html")
        )
)



server <- function(input, output) {
  
  output$ChillForcePlot <- renderPlot({
    
    # file input
    file <- input$FileUpload
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext %in% c("csv", "txt"), "Please check your data format!"))
    df.exp <- fread(file$datapath, skip = 2)
    
    # subtract information from input data
    info.exp <- str_extract_all(readLines(file$datapath, n = 2), "\\(?[0-9,.]+\\)?")
    year.start <- as.numeric(info.exp[[1]][1])
    month.start <- as.numeric(info.exp[[1]][2])
    day.start <- as.numeric(info.exp[[1]][3])
    bbch.exp <- as.numeric(info.exp[[2]][1])
    
    # temperature dataset imported from DWD
    id.station <- df.station$Stations_id[df.station$choice == input$ClimateStation]
    link <- selectDWD(id = id.station, res="daily", var="kl", per="h") # only historical data set
    file <- dataDWD(link, read=FALSE)
    clim <- readDWD(file, varnames=TRUE)
    df.T <- clim %>%
      mutate(MESS_DATUM = as.Date(MESS_DATUM)) %>%
      filter(MESS_DATUM >= as.Date(paste0(year.start, "-", month.start, "-", day.start, sep = ""))) %>%
      select(MESS_DATUM, TMK.Lufttemperatur)
    
    # calculate chilling and forcing units outside
    df.T$chilling[df.T$TMK.Lufttemperatur < input$ThresholdChilling] <- 1
    df.T$chilling[is.na(df.T$chilling)] <- 0
    df.T$chilling.sum <- cumsum(df.T$chilling)
    
    df.T$forcing <- df.T$TMK.Lufttemperatur
    df.T$forcing[df.T$forcing < input$ThresholdForcing] <- input$ThresholdForcing
    df.T$forcing.dd <- df.T$forcing - input$ThresholdForcing
    df.T$forcing.dd.sum <- cumsum(df.T$forcing.dd)
    
    df.exp$chilling.outside <- NA
    df.exp$forcing.outside <- NA
    for (i in 1:nrow(df.exp)){
      if (df.exp$Days_outside[i] == 0){
        df.exp$chilling.outside[i] <- 0
        df.exp$forcing.outside[i] <- 0
      } else{
        df.exp$chilling.outside[i] <- df.T$chilling.sum[df.exp$Days_outside[i]]
        df.exp$forcing.outside[i] <- df.T$forcing.dd.sum[df.exp$Days_outside[i]]
      }
    }
    
    # calculate forcing units inside
    if (input$Tinside <= input$ThresholdForcing){
      de.exp$forcing.inside <- 0
    } else{
      df.exp$forcing.inside <- df.exp$Days_inside * (input$Tinside - input$ThresholdForcing)
    }
    
    df.exp$chilling <- df.exp$chilling.outside
    df.exp$forcing <- df.exp$forcing.inside + df.exp$forcing.outside
    
    # fit simple model
    # y = a + blnx
    model.exp <- lm(forcing ~ log(chilling), df.exp[df.exp$chilling != 0,])
    cf.exp <- round(coef(model.exp), 3)
    eq.exp <- paste0("Forcing = ", cf.exp[1],
                     ifelse(sign(cf.exp[2])==1, " + ", " - "), abs(cf.exp[2]), " ln(Chilling), Adjusted R-squared: ", 
                     round(summary(model.exp)$adj.r.squared, 3))
    predict.exp <- data.frame(chilling = seq(1, 200, 1))
    predict.exp <- cbind(predict.exp, predict(model.exp, predict.exp, interval = "confidence"))

    
    # climate change simulation
    df.T.simu <- clim %>%
      mutate(MESS_DATUM = as.Date(MESS_DATUM)) %>%
      filter(MESS_DATUM >= as.Date(paste0(input$SimuDate, "-", month.start, "-", day.start, sep = ""))) %>%
      select(MESS_DATUM, TMK.Lufttemperatur)
    
    # adjust T scenario
    df.T.simu$TMK.Lufttemperatur[month(df.T.simu$MESS_DATUM) %in% c(12, 1, 2)] <- df.T.simu$TMK.Lufttemperatur[month(df.T.simu$MESS_DATUM) %in% c(12, 1, 2)] + input$WinterScenario
    df.T.simu$TMK.Lufttemperatur[month(df.T.simu$MESS_DATUM) %in% c(3, 4, 5)] <- df.T.simu$TMK.Lufttemperatur[month(df.T.simu$MESS_DATUM) %in% c(3, 4, 5)] + input$SpringScenario
    
    # calculate chilling and forcing units for simulation
    df.T.simu$chilling[df.T.simu$TMK.Lufttemperatur < input$ThresholdChilling] <- 1
    df.T.simu$chilling[is.na(df.T.simu$chilling)] <- 0
    df.T.simu$chilling.sum <- cumsum(df.T.simu$chilling)
    
    df.T.simu$forcing <- df.T.simu$TMK.Lufttemperatur
    df.T.simu$forcing[df.T.simu$forcing < input$ThresholdForcing] <- input$ThresholdForcing
    df.T.simu$forcing.dd <- df.T.simu$forcing - input$ThresholdForcing
    df.T.simu$forcing.dd.sum <- cumsum(df.T.simu$forcing.dd)
    
    # calculate the intersection
    df.T.simu$forcing.exp <- cf.exp[1] + cf.exp[2] * log(df.T.simu$chilling.sum)
    nrow.intersect <- head(which(df.T.simu$forcing.dd.sum > df.T.simu$forcing.exp), 1)
    
    df.T.simu$yeartime <- year(df.T.simu$MESS_DATUM)
    df.T.simu$weektime <- week(df.T.simu$MESS_DATUM)
    df.T.simu$timeindex <- paste0(df.T.simu$yeartime, "-W", df.T.simu$weektime, sep = "")
    df.T.simu.unique <- df.T.simu[!duplicated(df.T.simu$timeindex),] %>%
      filter(MESS_DATUM <= df.T.simu$MESS_DATUM[nrow.intersect])
    df.T.simu.unique.text <- df.T.simu %>%
      filter(MESS_DATUM %in% c(df.T.simu$MESS_DATUM[1], df.T.simu$MESS_DATUM[nrow.intersect]))
    
    df.ggplot1 <- data.frame(chilling = predict.exp$chilling,
                             forcing = predict.exp$fit,
                             lwr = predict.exp$lwr,
                             upr = predict.exp$upr,
                             group = "Experiment")
    df.ggplot2 <- data.frame(chilling = df.T.simu[1:nrow.intersect,]$chilling.sum,
                             forcing = df.T.simu[1:nrow.intersect,]$forcing.dd.sum,
                             lwr = NA,
                             upr = NA,
                             group = "Temperature")
    df.ggplot <- rbind(df.ggplot1, df.ggplot2)
    df.ggplot$group <- factor(df.ggplot$group, levels = c("Experiment", "Temperature"))
    
    annotations <- data.frame(
      xpos = c(-Inf,-Inf,Inf,Inf),
      ypos =  c(-Inf, Inf,-Inf,Inf),
      annotateText = c("","",paste0("BBCH: ", bbch.exp, sep = ""),eq.exp),
      hjustvar = c(0,0,1,1) ,
      vjustvar = c(0,1,-1,2)) #<- adjust
    
    ggplot(df.ggplot)+
      geom_point(data = df.exp[df.exp$chilling != 0,], aes(chilling ,forcing), pch = 21, col = "#377eb8", alpha = .5, size = 2)+
      geom_line(aes(x = chilling, y = forcing, col = group), size = 2)+
      scale_color_manual(name = "", values = c("#377eb8", "#e41a1c"),
                         label = c("Model fit based on experimental data (points)",  
                                   "Chilling-Forcing conditions based on temperature data from\nthe climate station and the scenario chosen in the simulator;\nLabels indicate the start date and the modeled date\nof reaching the studied BBCH stage"))+
      geom_label_repel(data = df.T.simu.unique.text, aes(x=chilling.sum, y=forcing.dd.sum, label=as.character(MESS_DATUM)),
                       fontface = 'bold', size = 6,
                       box.padding = unit(0.25, "lines"), point.padding = unit(0.5, "lines"))+
      geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), 
                fontface = 'bold', size = 6)+
      ylim(0, 2000)+
      xlab(paste0("Chilling [days with Tmean < ", input$ThresholdChilling, " °C]", sep=""))+
      ylab(paste0("Forcing [degree days with Tmean > ", input$ThresholdForcing, " °C]", sep = ""))+
      theme_bw()+
      theme(legend.position = "bottom", text = element_text(size=20))
    
  }, height = 800, width = 1000)
  
  output$ChillForceModelPlot <- renderPlot({
    # Landshut-Reithof (13710)
    link <- selectDWD(id = 13710, res="daily", var="kl", per="h")
    file <- dataDWD(link, read=FALSE)
    clim <- readDWD(file, varnames=TRUE)
    df.T <- clim %>%
      mutate(MESS_DATUM = as.Date(MESS_DATUM)) %>%
      filter(as.Date(MESS_DATUM) >= as.Date("2015-11-21")) %>%
      select(MESS_DATUM, TMK.Lufttemperatur)
    
    # calculate the coefficients a and b from paper exponential fit
    predict.model <- data.frame(chilling = seq(1, 200, 1))
    predict.model$forcing <- -330.301 + 1296.723 * exp(-0.01376 * predict.model$chilling)
    
    # climate change simulation
    df.T.simu <- clim %>%
      mutate(MESS_DATUM = as.Date(MESS_DATUM)) %>%
      filter(MESS_DATUM >= as.Date(paste0(input$SimuDateModel, "-11-21", sep = ""))) %>%
      select(MESS_DATUM, TMK.Lufttemperatur)
    
    # adjust T scenario
    df.T.simu$TMK.Lufttemperatur[month(df.T.simu$MESS_DATUM) %in% c(12, 1, 2)] <- df.T.simu$TMK.Lufttemperatur[month(df.T.simu$MESS_DATUM) %in% c(12, 1, 2)] + input$WinterScenarioModel
    df.T.simu$TMK.Lufttemperatur[month(df.T.simu$MESS_DATUM) %in% c(3, 4, 5)] <- df.T.simu$TMK.Lufttemperatur[month(df.T.simu$MESS_DATUM) %in% c(3, 4, 5)] + input$SpringScenarioModel
    
    # calculate chilling and forcing units for simulation
    df.T.simu$chilling[df.T.simu$TMK.Lufttemperatur < 5] <- 1
    df.T.simu$chilling[is.na(df.T.simu$chilling)] <- 0
    df.T.simu$chilling.sum <- cumsum(df.T.simu$chilling)
    
    df.T.simu$forcing <- df.T.simu$TMK.Lufttemperatur
    df.T.simu$forcing[df.T.simu$forcing < 5] <- 5
    df.T.simu$forcing.dd <- df.T.simu$forcing - 5
    df.T.simu$forcing.dd.sum <- cumsum(df.T.simu$forcing.dd)
    
    # calculate the intersection
    df.T.simu$forcing.exp <- -330.301 + 1296.723 * exp(-0.01376 * df.T.simu$chilling.sum)
    nrow.intersect <- head(which(df.T.simu$forcing.dd.sum > df.T.simu$forcing.exp), 1)
    
    df.T.simu$yeartime <- year(df.T.simu$MESS_DATUM)
    df.T.simu$weektime <- week(df.T.simu$MESS_DATUM)
    df.T.simu$timeindex <- paste0(df.T.simu$yeartime, "-W", df.T.simu$weektime, sep = "")
    df.T.simu.unique <- df.T.simu[!duplicated(df.T.simu$timeindex),] %>%
      filter(MESS_DATUM <= df.T.simu$MESS_DATUM[nrow.intersect])
    df.T.simu.unique.text <- df.T.simu %>%
      filter(MESS_DATUM %in% c(df.T.simu$MESS_DATUM[1], df.T.simu$MESS_DATUM[nrow.intersect]))
    
    df.ggplot1 <- data.frame(chilling = predict.model$chilling,
                             forcing = predict.model$forcing,
                             group = "Model")
    df.ggplot2 <- data.frame(chilling = df.T.simu[1:nrow.intersect,]$chilling.sum,
                             forcing = df.T.simu[1:nrow.intersect,]$forcing.dd.sum,
                             group = "Temperature")
    df.ggplot <- rbind(df.ggplot1, df.ggplot2)
    df.ggplot$group <- factor(df.ggplot$group, levels = c("Model", "Temperature"))
    
    eq.model <- "Forcing = -330.301 + 1296.723 * exp(-0.01376 * Chilling)"
    annotations <- data.frame(
      xpos = c(-Inf,-Inf,Inf,Inf),
      ypos =  c(-Inf, Inf,-Inf,Inf),
      annotateText = c("","","",eq.model),
      hjustvar = c(0,0,1,1) ,
      vjustvar = c(0,1,-1,2)) #<- adjust
    
    ggplot(df.ggplot)+
      geom_line(aes(x = chilling, y = forcing, col = group), size = 2)+
      scale_color_manual(name = "", values = c("#377eb8", "#e41a1c"),
                         label = c("Model fit from literature for comparison",  
                                   "Chilling-Forcing conditions based on temperature data from\nthe climate station and the scenario chosen in the simulator;\nLabels indicate the start date and the modeled date\nof reaching the studied BBCH stage"))+
      geom_label_repel(data = df.T.simu.unique.text, aes(x=chilling.sum, y=forcing.dd.sum, label=as.character(MESS_DATUM)),
                       fontface = 'bold', size = 6,
                       box.padding = unit(0.25, "lines"), point.padding = unit(0.5, "lines"))+
      geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),
                fontface = 'bold', size = 6)+
      ylim(0, 2000)+
      xlab(paste0("Chilling [days with Tmean < 5 °C]", sep=""))+
      ylab(paste0("Forcing [degree days with Tmean > 5 °C]", sep = ""))+
      theme_bw()+
      theme(legend.position = "bottom", text = element_text(size=20))
    
  }, height = 800, width = 1000)
  
}

# Run the application 
shinyApp(ui = ui, server = server)
