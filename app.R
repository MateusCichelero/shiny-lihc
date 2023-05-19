################################################################################
#EXERCICIO FINAL - CIENCIA DE DADOS EM SAUDE
#ALUNO: MATEUS CICHELERO DA SILVA
#MAIO 2023
################################################################################
library(SummarizedExperiment)
library(survival)
library(survminer)
library(circlize)
library(RColorBrewer)
library(palmerpenguins)
library(ComplexHeatmap)
library(dplyr)
################################################################################
### Loading anda processing data
################################################################################
setwd("~/Documents/cursos/ai-specialization/ds_saude/exercicio3")
load(file = "./tcgaLIHCdata_preprocessed.RData")

lihc_df <- colData(tcgaLIHCdata) #used in survival analysis
lihc_df$Age_Category = case_when(lihc_df$Age > 69  ~ '69y - 90y',
                                 lihc_df$Age > 61  & lihc_df$Age <= 69 ~ '61y - 69y',
                                 lihc_df$Age > 52  & lihc_df$Age <= 61 ~ '52y - 61y',
                                 lihc_df$Age > 0  & lihc_df$Age <= 52 ~ '0y - 52y') # end function

#--- Extract data matrix and metadata
gexp <- assay(tcgaLIHCdata)
rowAnnotation <- rowData(tcgaLIHCdata)
colAnnotation <- colData(tcgaLIHCdata)
rownames(gexp) <- rowAnnotation$SYMBOL

### Remove genes with low counts
idx <- rowSums(gexp!=0)/ncol(gexp)
gexp <- gexp[idx>0.3,]


################################################################################
### User interface
################################################################################

ui <- fluidPage(
  
  navbarPage(
    "TGCA-LIHC Data Analysis",
    
    tabPanel("ComplexHeatmap Clustering",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("selection_type_heatmap", "Feature Selection Strategy:",
                              c("Abundance Coefficient"="ac", "Spearman Correlation"="sc")
                 )
               ),
               mainPanel(
                 plotOutput("plot_heatmap")
               )
             )
    ),
    
    tabPanel("Survival Analysis",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("stratify_survival", "Stratify by:",
                              c("Tumor Stage"="ts", "Gender"="gender", "Age Group"="age")
                 )
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Kaplan-Meier Plot", 
                            plotOutput("plot_km")),
                   tabPanel("Cox Regression Forest Plot", 
                            plotOutput("plot_cox"))
                 )
               
    
    )
    
    
             )
  
)))

################################################################################
### Server
################################################################################
server <- function(input, output){

  
  output$plot_heatmap <- renderPlot({
    
    if (input$selection_type_heatmap == "ac") {
      
      ### Filter the 'gexp' matrix using feature abundance
      idx <- sort.list(apply(gexp, 1, mean), decreasing = T)[1:100]
      gexp_filt <- gexp[idx,]
      ### Filter the 'colAnnotation', removing NAs
      #--- Get sample annotations, and remove NAs
      colAnnotation_filt <- colAnnotation[,c("Tumor_Stage"), drop=F]
      colAnnotation_filt <- colAnnotation_filt[complete.cases(colAnnotation_filt),, drop=F]
      gexp_filt <- gexp_filt[ ,rownames(colAnnotation_filt)]
      ### Re-scale the data using a columnwise rank transformation
      x <- gexp_filt
      x <- t(apply(x, 1, rank));x <- x/max(x)
      x <- t(scale(t(x), center = TRUE, scale = F)) #z-score
      ### Run semi- or unsupervised clustering analysis 
      #--- Set col annotations
      colAnnotation_filt$Tumor_Stage <- as.factor(colAnnotation_filt$Tumor_Stage)
      pal1 <- brewer.pal(4,"Set1")
      names(pal1) <- levels(colAnnotation_filt$Tumor_Stage)
      top_annotation <- columnAnnotation(df=colAnnotation_filt, 
                                         col=list('Tumor_Stage'=pal1))
      #--- Set a color scheme
      pal2 <- rev(brewer.pal(7,"RdYlBu"))
      bks <- quantile(as.numeric(x), probs = seq(0,1, length.out = length(pal2)))
      colors <- colorRamp2(breaks = bks, colors = pal2)
      
      #--- Run clustering analysis and plot a large heatmap with ComplexHeatmap
      p <- Heatmap(x, col = colors, name = "RNA-seq", 
              column_split = colAnnotation_filt$Tumor_Stage,
              show_row_names = F, show_column_names = F, 
              top_annotation=top_annotation,
              clustering_method_rows = "ward.D2", 
              clustering_distance_rows="spearman",
              clustering_method_columns = "ward.D2", 
              clustering_distance_columns = "spearman")
    }
    
    if (input$selection_type_heatmap == "sc") {
      
      idx <- cor(t(gexp), colAnnotation$Tumor_Stage, method = "spearman",
                 use="complete.obs")
      idx <- sort.list(abs(idx), decreasing = T)[1:100]
      gexp_filt <- gexp[idx,]
      ### Filter the 'colAnnotation', removing NAs
      #--- Get sample annotations, and remove NAs
      colAnnotation_filt <- colAnnotation[,c("Tumor_Stage"), drop=F]
      colAnnotation_filt <- colAnnotation_filt[complete.cases(colAnnotation_filt),, drop=F]
      gexp_filt <- gexp_filt[ ,rownames(colAnnotation_filt)]
      ### Re-scale the data using a columnwise rank transformation
      x <- gexp_filt
      x <- t(apply(x, 1, rank));x <- x/max(x)
      x <- t(scale(t(x), center = TRUE, scale = F)) #z-score
      ### Run semi- or unsupervised clustering analysis 
      #--- Set col annotations
      colAnnotation_filt$Tumor_Stage <- as.factor(colAnnotation_filt$Tumor_Stage)
      pal1 <- brewer.pal(4,"Set1")
      names(pal1) <- levels(colAnnotation_filt$Tumor_Stage)
      top_annotation <- columnAnnotation(df=colAnnotation_filt, 
                                         col=list('Tumor_Stage'=pal1))
      #--- Set a color scheme
      pal2 <- rev(brewer.pal(7,"RdYlBu"))
      bks <- quantile(as.numeric(x), probs = seq(0,1, length.out = length(pal2)))
      colors <- colorRamp2(breaks = bks, colors = pal2)
      
      #--- Run clustering analysis and plot a large heatmap with ComplexHeatmap
      p <- Heatmap(x, col = colors, name = "RNA-seq", 
                   column_split = colAnnotation_filt$Tumor_Stage,
                   show_row_names = F, show_column_names = F, 
                   top_annotation=top_annotation,
                   clustering_method_rows = "ward.D2", 
                   clustering_distance_rows="spearman",
                   clustering_method_columns = "ward.D2", 
                   clustering_distance_columns = "spearman")
      
    }
    p
    
    
  })
  
  output$plot_km <- renderPlot({
    if (input$stratify_survival == "ts") {
      # Fit survival curves
      fit <- survfit(Surv(OS.time.months, OS) ~ Tumor_Stage, data = lihc_df) 

      
      # Customized survival curves
      p <- ggsurvplot(fit, data = lihc_df,
                 surv.median.line = "hv", # Add medians survival
                 
                 # Change legends: title & labels
                 legend.title = "Tumor Stage",
                 legend.labs = c("1", "2", "3", "4"),
                 xlab = "Months",
                 # Add p-value and tervals
                 pval = TRUE,
                 
                 conf.int = TRUE,
                 # Add risk table
                 risk.table = TRUE,
                 tables.height = 0.2,
                 tables.theme = theme_cleantable(),
                 
                 # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
                 # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
                 #palette = c("#E7B800", "#2E9FDF"),
                 ggtheme = theme_bw() # Change ggplot2 theme
      )
      
    }
    
    if (input$stratify_survival == "gender") {
      # Fit survival curves
      fit <- survfit(Surv(OS.time.months, OS) ~ gender, data = lihc_df) 
      
      # Customized survival curves
      p <- ggsurvplot(fit, data = lihc_df,
                 surv.median.line = "hv", # Add medians survival
                 
                 # Change legends: title & labels
                 legend.title = "Gender",
                 legend.labs = c("Female", "Male"),
                 xlab = "Months",
                 # Add p-value and intervals
                 pval = TRUE,
                 
                 conf.int = TRUE,
                 # Add risk table
                 risk.table = TRUE,
                 tables.height = 0.2,
                 tables.theme = theme_cleantable(),
                 
                 # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
                 # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
                 #palette = c("#E7B800", "#2E9FDF"),
                 ggtheme = theme_bw() # Change ggplot2 theme
      )
      
    }
    
    if (input$stratify_survival == "age") {
      # Fit survival curves
      fit <- survfit(Surv(OS.time.months, OS) ~ Age_Category, data = lihc_df) 

      # Customized survival curves
      p <- ggsurvplot(fit, data = lihc_df,
                 surv.median.line = "hv", # Add medians survival
                 
                 # Change legends: title & labels
                 legend.title = "Age",
                 #legend.labs = c("Female", "Male"),
                 xlab = "Months",
                 # Add p-value and intervals
                 pval = TRUE,
                 
                 conf.int = TRUE,
                 # Add risk table
                 risk.table = TRUE,
                 tables.height = 0.2,
                 tables.theme = theme_cleantable(),
                 
                 # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
                 # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
                 #palette = c("#E7B800", "#2E9FDF"),
                 ggtheme = theme_bw() # Change ggplot2 theme
      )
      
    }
    
    p
  })
  
  output$plot_cox <- renderPlot({
    if (input$stratify_survival == "ts") {
      fit <- coxph(Surv(OS.time.months, OS) ~ Stage, data = lihc_df)
      p <- ggforest(fit)
    }
    
    if (input$stratify_survival == "gender") {
      fit <- coxph(Surv(OS.time.months, OS) ~ gender, data = lihc_df)
      p <- ggforest(fit)
      
    }
    
    if (input$stratify_survival == "age") {
      fit <- coxph(Surv(OS.time.months, OS) ~ Age_Category, data = lihc_df)
      p <- ggforest(fit)
      
    }
    
    p

  })

  
}

shinyApp(ui, server)