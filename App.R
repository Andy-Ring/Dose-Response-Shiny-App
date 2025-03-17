# Load required packages
library(shiny)
library(bslib)          # For custom Bootstrap theming
library(ggplot2)
library(drc)
library(dplyr)
library(broom)          # For tidying ANOVA results
library(thematic)

thematic_shiny()

# Helper function to clean group names
cleanGroupName <- function(name) {
  if (grepl(":", name)) {
    parts <- strsplit(name, ":")[[1]]
    if (length(parts) >= 2) return(trimws(parts[2]))
    else return(name)
  }
  trimws(name)
}

ui <- fluidPage(
  theme = bs_theme(
    bootswatch = "darkly",
    secondary = "#BA0C2F",
    "table-bg" = "primary"
  ),
  titlePanel("Dose–Response Analysis"),
  sidebarLayout(
    sidebarPanel(
      HTML('<img src="logo.png" width="100%" height="auto">'),
      br(), br(),
      # Removed the downloadTemplate button
      fileInput("dataFile", "Upload Your Data (CSV)",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      # UI for selecting which columns to use
      uiOutput("columnSelectionUI"),
      uiOutput("recommendationUI"),
      tags$hr(),
      radioButtons("modelChoice", "Select Curve Model:",
                   choices = list("Recommended" = "recommended",
                                  "4-parameter logistic (LL.4)" = "LL.4",
                                  "3-parameter logistic (LL.3)" = "LL.3",
                                  "2-parameter logistic (LL.2)" = "LL.2",
                                  "Linear Model" = "linear"),
                   selected = "recommended"),
      tags$hr(),
      h4("Plot Customization"),
      textInput("plotTitle", "Plot Title:", value = "Dose–Response Curve"),
      textInput("xLabel", "X Axis Label:", value = "Dose"),
      textInput("yLabel", "Y Axis Label:", value = "Response"),
      checkboxInput("logScale", "Logarithmic Dose Axis", value = FALSE),
      checkboxInput("sciIC50", "Display IC50 in scientific notation", value = TRUE),
      tags$hr(),
      actionButton("runAnalysis", "Run Analysis"),
      tags$hr(),
      tags$p("Created by Andy Ring"),
      tags$p("Version 1.1.0 | March 17th, 2025")
    ),
    mainPanel(
      layout_columns(
        card(card_header("Dose-Response Plot"), plotOutput("doseResponsePlot")),
        card(card_header("IC50 Values"), tableOutput("ic50Table")),
        card(card_header("ANOVA Results"), tableOutput("anovaTable")),
        downloadButton("downloadPlot", "Download Plot"),
        downloadButton("downloadData", "Download Analysis Summary"),
        col_widths = c(12, 4, 8, 4, 4, -4)
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Read the raw CSV data
  rawData <- reactive({
    req(input$dataFile)
    read.csv(input$dataFile$datapath)
  })
  
  # Create a UI to allow the user to choose the columns for dose, response, and group
  output$columnSelectionUI <- renderUI({
    req(rawData())
    data <- rawData()
    tagList(
      selectInput("doseColumn", "Select Dose Column:", choices = names(data)),
      selectInput("responseColumn", "Select Response Column:", choices = names(data)),
      selectInput("groupColumn", "Select Group Column (if applicable):", 
                  choices = c("None", names(data)))
    )
  })
  
  # Process the data by renaming the selected columns to standard names
  processedData <- reactive({
    req(rawData(), input$doseColumn, input$responseColumn)
    data <- rawData()
    data$Dose <- data[[input$doseColumn]]
    data$Response <- data[[input$responseColumn]]
    if (!is.null(input$groupColumn) && input$groupColumn != "None") {
      data$Group <- data[[input$groupColumn]]
    }
    data
  })
  
  #### 2. Pre-Analysis Recommendations ####
  recommendations <- reactive({
    req(processedData())
    data <- processedData()
    rec <- list()
    
    if("Dose" %in% names(data)) {
      doseRange <- range(data$Dose, na.rm = TRUE)
      rec$model <- if (diff(doseRange) > 50) "LL.4" else "linear"
    } else {
      rec$model <- "LL.4"
    }
    
    if("Group" %in% names(data) && length(unique(na.omit(data$Group))) > 1) {
      rec$test <- "Two-Way ANOVA (Dose as factor and Group as factor) will be run automatically."
    } else {
      rec$test <- "No grouping variable detected (two-way ANOVA not applicable)."
    }
    rec
  })
  
  output$recommendationUI <- renderUI({
    req(processedData())
    rec <- recommendations()
    tagList(
      h4("Recommendations (Based on Uploaded Data):"),
      p(strong("Curve Model:"), rec$model),
      p(strong("Statistical Test:"), rec$test)
    )
  })
  
  #### 3. Run Analysis ####
  analysisResults <- eventReactive(input$runAnalysis, {
    req(processedData())
    data <- processedData()
    
    selectedModel <- input$modelChoice
    if(selectedModel == "recommended") {
      selectedModel <- recommendations()$model
    }
    
    result <- list(data = data, selectedModel = selectedModel)
    
    if("Group" %in% names(data) && length(unique(na.omit(data$Group))) > 1) {
      if(selectedModel %in% c("LL.4", "LL.3", "LL.2")) {
        fit <- try(drm(Response ~ Dose, curveid = Group, data = data,
                       fct = switch(selectedModel,
                                    "LL.4" = LL.4(),
                                    "LL.3" = LL.3(),
                                    "LL.2" = LL.2())),
                   silent = TRUE)
        if(inherits(fit, "try-error")) {
          result$fit <- NULL
          result$ic50 <- NA
        } else {
          result$fit <- fit
          # Use type = "relative" to calculate IC50 as the midpoint between asymptotes
          ed <- try(ED(fit, 50, type = "relative"), silent = TRUE)
          result$ic50 <- if (inherits(ed, "try-error")) NA else ed[,1]
        }
      } else if(selectedModel == "linear") {
        lmFit <- try(lm(Response ~ Dose * as.factor(Group), data = data), silent = TRUE)
        if(inherits(lmFit, "try-error")) {
          result$fit <- NULL
          result$ic50 <- NA
        } else {
          result$fit <- lmFit
          coefs <- coef(lmFit)
          groups_levels <- levels(as.factor(data$Group))
          ic50s <- list()
          for(g in groups_levels) {
            if(g == groups_levels[1]) {
              int_val <- coefs["(Intercept)"]
              slope_val <- coefs["Dose"]
            } else {
              int_val <- coefs["(Intercept)"] + coefs[paste0("as.factor(Group)", g)]
              slope_val <- coefs["Dose"] + coefs[paste0("Dose:as.factor(Group)", g)]
            }
            subdata <- data[data$Group == g, ]
            target <- (max(subdata$Response, na.rm = TRUE) + min(subdata$Response, na.rm = TRUE)) / 2
            ic50s[[g]] <- ifelse(slope_val == 0, NA, (target - int_val) / slope_val)
          }
          result$ic50 <- ic50s
        }
      }
      anovaFit <- aov(Response ~ as.factor(Dose) * as.factor(Group), data = data)
      result$anovaModel <- anovaFit
      
    } else {
      ### Single Group Analysis ###
      if(selectedModel %in% c("LL.4", "LL.3", "LL.2")) {
        fit <- try(drm(Response ~ Dose, data = data,
                       fct = switch(selectedModel,
                                    "LL.4" = LL.4(),
                                    "LL.3" = LL.3(),
                                    "LL.2" = LL.2())),
                   silent = TRUE)
        if(inherits(fit, "try-error")) {
          result$fit <- NULL
          result$ic50 <- NA
        } else {
          result$fit <- fit
          # Use type = "relative" for a relative IC50 calculation
          ed <- try(ED(fit, 50, type = "relative"), silent = TRUE)
          result$ic50 <- if (inherits(ed, "try-error")) NA else ed[1,1]
        }
      } else if(selectedModel == "linear") {
        fit <- try(lm(Response ~ Dose, data = data), silent = TRUE)
        if(inherits(fit, "try-error")) {
          result$fit <- NULL
          result$ic50 <- NA
        } else {
          result$fit <- fit
          coefs <- coef(fit)
          target <- (max(data$Response, na.rm = TRUE) + min(data$Response, na.rm = TRUE)) / 2
          result$ic50 <- ifelse(coefs["Dose"] == 0, NA, (target - coefs["(Intercept)"]) / coefs["Dose"])
        }
      }
      result$anovaModel <- NULL
    }
    result
  })
  
  #### 4. Create Plot Object ####
  plotObject <- reactive({
    req(analysisResults())
    res <- analysisResults()
    data <- res$data
    
    if("Group" %in% names(data) && length(unique(na.omit(data$Group))) > 1) {
      summaryData <- data %>%
        group_by(Dose, Group) %>%
        summarise(meanResponse = mean(Response, na.rm = TRUE),
                  seResponse = sd(Response, na.rm = TRUE) / sqrt(n()),
                  .groups = "drop")
      p <- ggplot(summaryData, aes(x = Dose, y = meanResponse, color = as.factor(Group))) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = meanResponse - seResponse, ymax = meanResponse + seResponse),
                      width = 0.1) +
        labs(title = input$plotTitle, x = input$xLabel, y = input$yLabel, color = "Group")
      if(!is.null(res$fit)) {
        dose_seq <- seq(min(data$Dose, na.rm = TRUE), max(data$Dose, na.rm = TRUE),
                        length.out = 100)
        for(g in unique(data$Group)) {
          newData <- data.frame(Dose = dose_seq, Group = g)
          pred <- predict(res$fit, newdata = newData)
          pred_df <- data.frame(Dose = dose_seq, pred = pred, Group = g)
          p <- p + geom_line(data = pred_df, aes(x = Dose, y = pred, color = as.factor(Group)),
                             size = 1)
        }
      }
    } else {
      summaryData <- data %>%
        group_by(Dose) %>%
        summarise(meanResponse = mean(Response, na.rm = TRUE),
                  seResponse = sd(Response, na.rm = TRUE) / sqrt(n()),
                  .groups = "drop")
      p <- ggplot(summaryData, aes(x = Dose, y = meanResponse)) +
        geom_point(size = 3, color = "blue") +
        geom_errorbar(aes(ymin = meanResponse - seResponse, ymax = meanResponse + seResponse),
                      width = 0.1, color = "blue") +
        labs(title = input$plotTitle, x = input$xLabel, y = input$yLabel)
      if(!is.null(res$fit)) {
        dose_seq <- seq(min(data$Dose, na.rm = TRUE), max(data$Dose, na.rm = TRUE),
                        length.out = 100)
        pred <- predict(res$fit, newdata = data.frame(Dose = dose_seq))
        pred_df <- data.frame(Dose = dose_seq, pred = pred)
        p <- p + geom_line(data = pred_df, aes(x = Dose, y = pred),
                           color = "red", size = 1)
      }
    }
    
    if(input$logScale) {
      p <- p + scale_x_log10()
    }
    p
  })
  
  output$doseResponsePlot <- renderPlot({
    plotObject()
  })
  
  #### 5. Render IC50 Table ####
  output$ic50Table <- renderTable({
    req(analysisResults())
    res <- analysisResults()
    sci <- input$sciIC50
    
    if(length(res$ic50) > 1) {
      groups <- names(res$ic50)
      if(is.null(groups)) {
        groups <- paste0("Group ", seq_along(res$ic50))
      } else {
        groups <- sapply(groups, cleanGroupName)
      }
      df <- data.frame(
        Group = groups,
        IC50 = sapply(res$ic50, function(x) {
          if(is.na(x)) "NA" else format(x, scientific = sci)
        }),
        stringsAsFactors = FALSE
      )
      df
    } else {
      groupName <- if(is.null(names(res$ic50))) "All" else cleanGroupName(names(res$ic50))
      df <- data.frame(
        Group = groupName,
        IC50 = if(is.na(res$ic50)) "NA" else format(res$ic50, scientific = sci),
        stringsAsFactors = FALSE
      )
      df
    }
  }, striped = TRUE, hover = TRUE)
  
  #### 6. Render Cleaned-Up ANOVA Results ####
  output$anovaTable <- renderTable({
    req(analysisResults())
    res <- analysisResults()
    if(!is.null(res$anovaModel)) {
      tidy_anova <- tidy(res$anovaModel)
      tidy_anova$term <- gsub("as.factor\\((.*?)\\)", "\\1", tidy_anova$term)
      tidy_anova$p.value <- sapply(tidy_anova$p.value, function(x) format(x, scientific = TRUE))
      tidy_anova
    } else {
      data.frame(Message = "No two-way ANOVA performed (single group data).")
    }
  }, striped = TRUE, hover = TRUE)
  
  #### 7. Download Handlers ####
  output$downloadPlot <- downloadHandler(
    filename = function() { "dose_response_plot.png" },
    content = function(file) {
      req(analysisResults())
      res <- analysisResults()
      data <- res$data
      p <- NULL
      
      if("Group" %in% names(data) && length(unique(na.omit(data$Group))) > 1) {
        summaryData <- data %>%
          group_by(Dose, Group) %>%
          summarise(
            meanResponse = mean(Response, na.rm = TRUE),
            seResponse = sd(Response, na.rm = TRUE) / sqrt(n()),
            .groups = "drop"
          )
        p <- ggplot(summaryData, aes(x = Dose, y = meanResponse, color = as.factor(Group))) +
          geom_point(size = 3) +
          geom_errorbar(aes(ymin = meanResponse - seResponse, ymax = meanResponse + seResponse),
                        width = 0.1) +
          labs(title = input$plotTitle, x = input$xLabel, y = input$yLabel, color = "Group") +
          theme_light()
        if (!is.null(res$fit)) {
          dose_seq <- seq(min(data$Dose, na.rm = TRUE), max(data$Dose, na.rm = TRUE), length.out = 100)
          for (g in unique(data$Group)) {
            newData <- data.frame(Dose = dose_seq, Group = g)
            pred <- predict(res$fit, newdata = newData)
            pred_df <- data.frame(Dose = dose_seq, pred = pred, Group = g)
            p <- p + geom_line(data = pred_df, aes(x = Dose, y = pred, color = as.factor(Group)),
                               size = 1)
          }
        }
      } else {
        summaryData <- data %>%
          group_by(Dose) %>%
          summarise(
            meanResponse = mean(Response, na.rm = TRUE),
            seResponse = sd(Response, na.rm = TRUE) / sqrt(n()),
            .groups = "drop"
          )
        p <- ggplot(summaryData, aes(x = Dose, y = meanResponse)) +
          geom_point(size = 3, color = "blue") +
          geom_errorbar(aes(ymin = meanResponse - seResponse, ymax = meanResponse + seResponse),
                        width = 0.1, color = "blue") +
          labs(title = input$plotTitle, x = input$xLabel, y = input$yLabel) +
          theme_light()
        if (!is.null(res$fit)) {
          dose_seq <- seq(min(data$Dose, na.rm = TRUE), max(data$Dose, na.rm = TRUE), length.out = 100)
          pred <- predict(res$fit, newdata = data.frame(Dose = dose_seq))
          pred_df <- data.frame(Dose = dose_seq, pred = pred)
          p <- p + geom_line(data = pred_df, aes(x = Dose, y = pred),
                             color = "red", size = 1)
        }
      }
      
      if (input$logScale) {
        p <- p + scale_x_log10()
      }
      
      ggsave(file, plot = p, device = "png", width = 12, height = 8, dpi = 800)
    }
  )
  
  output$downloadData <- downloadHandler(
    filename = function() { "analysis_summary.txt" },
    content = function(file) {
      res <- analysisResults()
      sink(file)
      cat("=== Analysis Summary ===\n")
      cat("Curve Model Used:\n", res$selectedModel, "\n\n")
      cat("IC50 (ED50) value(s):\n")
      print(res$ic50)
      cat("\n")
      if(!is.null(res$anovaModel)) {
        cat("ANOVA Results:\n")
        tidy_anova <- tidy(res$anovaModel)
        tidy_anova$term <- gsub("as.factor\\((.*?)\\)", "\\1", tidy_anova$term)
        tidy_anova$p.value <- sapply(tidy_anova$p.value, function(x) format(x, scientific = TRUE))
        print(tidy_anova)
      } else {
        cat("No two-way ANOVA performed (single group data).\n")
      }
      sink()
    }
  )
}

shinyApp(ui = ui, server = server)
