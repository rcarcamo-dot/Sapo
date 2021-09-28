library(shiny)

shinyUI(fluidPage(
    # App title ----
    titlePanel("Sapoapp"),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            radioButtons("radioButtons", 
                               h3("Select input format"), 
                               choices = list("Upload FASTA file" = 1, 
                                              "Paste sequence" = 2),
                               selected = 0),
            
            conditionalPanel(
                condition ="input.radioButtons=='1'",
                fileInput("file1", h3("Choose Fasta File"),
                          multiple = TRUE,
                          accept = c("text/fasta",
                                     ".fa",
                                     ".fasta"))
                ),
            
            conditionalPanel(
                condition ="input.radioButtons=='2'",
                #selectizeInput("text", h3("Input words here"), 
                               #choices = NULL, 
                               #multiple = FALSE, 
                               #options = list(create = TRUE)),
                textAreaInput('txt', h3('Input sequence:'), value = "", 
                              placeholder = "", width = "600px", height="200px"),
                tags$hr(),
                
                
                
                actionButton(inputId = "searchButton", label = "Run")
            )
            
            
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            tabsetPanel(
                tabPanel("Genotype", tableOutput("selected_var") ), 
                tabPanel("Phylogenetic analysis", verbatimTextOutput("summary"))
            )
        
    )
)))

