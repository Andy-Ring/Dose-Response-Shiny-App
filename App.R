# Load required packages
library(shiny)
library(ggplot2)
library(drc)
library(dplyr)

# Define the UI
ui <- fluidPage(
  titlePanel("Dose–Response Curve Fitting, IC50 Calculation, and Two-Way ANOVA"),
  sidebarLayout(
    sidebarPanel(
      # Download CSV template button
      downloadButton("downloadTemplate", "Download Data Template"),
      tags$hr(),
      
      # File input for data
      fileInput("dataFile", "Upload Your Data (CSV)",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      
      # Display recommendations (based on uploaded data)
      uiOutput("recommendationUI"),
      tags$hr(),
      
      # Choose the curve model
      radioButtons("modelChoice", "Select Curve Model:",
                   choices = list("Recommended" = "recommended",
                                  "4-parameter logistic (LL.4)" = "LL.4",
                                  "3-parameter logistic (LL.3)" = "LL.3",
                                  "2-parameter logistic (LL.2)" = "LL.2",
                                  "Linear Model" = "linear"),
                   selected = "recommended"),
      tags$hr(),
      
      # Plot customization options
      h4("Plot Customization"),
      textInput("plotTitle", "Plot Title:", value = "Dose–Response Curve"),
      textInput("xLabel", "X Axis Label:", value = "Dose"),
      textInput("yLabel", "Y Axis Label:", value = "Response"),
      checkboxInput("logScale", "Logarithmic Dose Axis", value = FALSE),
      tags$hr(),
      
      # Run analysis button
      actionButton("runAnalysis", "Run Analysis"),
      tags$hr(),
      
      # Download buttons for the plot and summary
      downloadButton("downloadPlot", "Download High-Res Plot"),
      br(), br(),
      downloadButton("downloadData", "Download Analysis Summary")
    ),
    mainPanel(
      plotOutput("doseResponsePlot"),
      verbatimTextOutput("analysisSummary")
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  #### 1. CSV Template Download ####
  output$downloadTemplate <- downloadHandler(
    filename = function() { "data_template.csv" },
    content = function(file) {
      # Template columns: Dose (numeric), Response (numeric), and Group (optional)
      template <- data.frame(
        Dose = c(0.1, 0.5, 1, 5, 10),
        Response = rep(NA, 5),
        Group = rep(NA, 5)
      )
      write.csv(template, file, row.names = FALSE)
    }
  )
  
  #### 2. Pre-Analysis Recommendations ####
  recommendations <- reactive({
    req(input$dataFile)
    data <- read.csv(input$dataFile$datapath)
    rec <- list()
    
    # Simple model recommendation: if the dose range is large (>50), use a logistic model; otherwise, linear.
    if("Dose" %in% names(data)) {
      doseRange <- range(data$Dose, na.rm = TRUE)
      if(diff(doseRange) > 50) {
        rec$model <- "LL.4"
      } else {
        rec$model <- "linear"
      }
    } else {
      rec$model <- "LL.4"
    }
    
    # Recommend two-way ANOVA if a grouping variable exists with more than one level.
    if("Group" %in% names(data) && length(unique(na.omit(data$Group))) > 1) {
      rec$test <- "Two-Way ANOVA (Dose as factor and Group as factor) will be run automatically."
    } else {
      rec$test <- "No grouping variable detected (two-way ANOVA not applicable)."
    }
    rec
  })
  
  output$recommendationUI <- renderUI({
    req(input$dataFile)
    rec <- recommendations()
    tagList(
      h4("Recommendations (Based on Uploaded Data):"),
      p(strong("Curve Model:"), rec$model),
      p(strong("Statistical Test:"), rec$test)
    )
  })
  
  #### 3. Run Analysis ####
  analysisResults <- eventReactive(input$runAnalysis, {
    req(input$dataFile)
    data <- read.csv(input$dataFile$datapath)
    
    # Determine which curve model to use.
    selectedModel <- input$modelChoice
    if(selectedModel == "recommended") {
      selectedModel <- recommendations()$model
    }
    
    result <- list(data = data, selectedModel = selectedModel)
    
    # Multiple groups: fit curves per group and run two-way ANOVA.
    if("Group" %in% names(data) && length(unique(na.omit(data$Group))) > 1) {
      groups <- unique(na.omit(data$Group))
      if(selectedModel %in% c("LL.4", "LL.3", "LL.2")) {
        fit <- try(
          drm(Response ~ Dose, curveid = Group, data = data,
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
          ed <- try(ED(fit, 50, type = "absolute"), silent = TRUE)
          if(inherits(ed, "try-error")){
            result$ic50 <- NA
          } else {
            result$ic50 <- ed[,1]  # ED returns a matrix with rows for each group
          }
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
            ic50_val <- ifelse(slope_val == 0, NA, (target - int_val) / slope_val)
            ic50s[[g]] <- ic50_val
          }
          result$ic50 <- ic50s
        }
      }
      # Run two-way ANOVA (treat Dose as a factor and Group as a factor)
      anovaFit <- aov(Response ~ as.factor(Dose) * as.factor(Group), data = data)
      result$anovaResult <- summary(anovaFit)
      result$testResult <- paste(capture.output(result$anovaResult), collapse = "\n")
      
    } else {
      ### Single Group Analysis ###
      if(selectedModel %in% c("LL.4", "LL.3", "LL.2")) {
        fit <- try(
          drm(Response ~ Dose, data = data,
              fct = switch(selectedModel,
                           "LL.4" = LL.4(),
                           "LL.3" = LL.3(),
                           "LL.2" = LL.2())),
          silent = TRUE)
        if(inherits(fit, "try-error")){
          result$fit <- NULL
          result$ic50 <- NA
        } else {
          result$fit <- fit
          ed <- try(ED(fit, 50, type = "absolute"), silent = TRUE)
          if(inherits(ed, "try-error")){
            result$ic50 <- NA
          } else {
            result$ic50 <- ed[1,1]
          }
        }
      } else if(selectedModel == "linear") {
        fit <- try(lm(Response ~ Dose, data = data), silent = TRUE)
        if(inherits(fit, "try-error")){
          result$fit <- NULL
          result$ic50 <- NA
        } else {
          result$fit <- fit
          coefs <- coef(fit)
          target <- (max(data$Response, na.rm = TRUE) + min(data$Response, na.rm = TRUE)) / 2
          ic50_val <- ifelse(coefs["Dose"] == 0, NA, (target - coefs["(Intercept)"]) / coefs["Dose"])
          result$ic50 <- ic50_val
        }
      }
      result$testResult <- "No two-way ANOVA performed (single group data)."
    }
    return(result)
  })
  
  #### 4. Create Plot Object with Mean ± SE and Overlay Fitted Curve(s) ####
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
        geom_errorbar(aes(ymin = meanResponse - seResponse, ymax = meanResponse + seResponse), width = 0.1) +
        labs(title = input$plotTitle, x = input$xLabel, y = input$yLabel, color = "Group") +
        theme_light()
      
      # Overlay fitted curves if available
      if(!is.null(res$fit)) {
        dose_seq <- seq(min(data$Dose, na.rm = TRUE), max(data$Dose, na.rm = TRUE), length.out = 100)
        for(g in unique(data$Group)) {
          newData <- data.frame(Dose = dose_seq, Group = g)
          pred <- predict(res$fit, newdata = newData)
          pred_df <- data.frame(Dose = dose_seq, pred = pred, Group = g)
          p <- p + geom_line(data = pred_df, aes(x = Dose, y = pred, color = as.factor(Group)), size = 1)
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
        geom_errorbar(aes(ymin = meanResponse - seResponse, ymax = meanResponse + seResponse), width = 0.1, color = "blue") +
        labs(title = input$plotTitle, x = input$xLabel, y = input$yLabel) +
        theme_light()
      
      # Overlay fitted curve if available
      if(!is.null(res$fit)) {
        dose_seq <- seq(min(data$Dose, na.rm = TRUE), max(data$Dose, na.rm = TRUE), length.out = 100)
        pred <- predict(res$fit, newdata = data.frame(Dose = dose_seq))
        pred_df <- data.frame(Dose = dose_seq, pred = pred)
        p <- p + geom_line(data = pred_df, aes(x = Dose, y = pred), color = "red", size = 1)
      }
    }
    
    # Optionally apply logarithmic scale to the Dose axis.
    if(input$logScale) {
      p <- p + scale_x_log10()
    }
    return(p)
  })
  
  #### 5. Render Plot and Analysis Summary ####
  output$doseResponsePlot <- renderPlot({
    plotObject()
  })
  
  output$analysisSummary <- renderPrint({
    req(analysisResults())
    res <- analysisResults()
    cat("=== Analysis Summary ===\n")
    cat("Curve Model Used:\n")
    cat("  ", res$selectedModel, "\n\n")
    cat("IC50 (ED50) value(s):\n")
    print(res$ic50)
    cat("\nTwo-Way ANOVA Results (if applicable):\n")
    if(!is.null(res$anovaResult)) {
      cat(res$testResult, "\n")
    } else {
      cat("No two-way ANOVA performed (single group data).\n")
    }
  })
  
  #### 6. Download Handlers ####
  # Download high-resolution plot as PNG
  output$downloadPlot <- downloadHandler(
    filename = function() { "dose_response_plot.png" },
    content = function(file) {
      ggsave(file, plot = plotObject(), device = "png", width = 12, height = 8, dpi = 800)
    }
  )
  
  # Download analysis summary as a text file
  output$downloadData <- downloadHandler(
    filename = function() { "analysis_summary.txt" },
    content = function(file) {
      res <- analysisResults()
      sink(file)
      cat("=== Analysis Summary ===\n")
      cat("Curve Model Used:\n")
      cat("  ", res$selectedModel, "\n\n")
      cat("IC50 (ED50) value(s):\n")
      print(res$ic50)
      cat("\nTwo-Way ANOVA Results (if applicable):\n")
      if(!is.null(res$anovaResult)) {
        cat(res$testResult, "\n")
      } else {
        cat("No two-way ANOVA performed (single group data).\n")
      }
      sink()
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)









