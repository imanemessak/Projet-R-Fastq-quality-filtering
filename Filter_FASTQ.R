library('shiny')
library('Biostrings')
library('ShortRead')
library('FastqCleaner')
library('ggplot2')
library("Rqc")
filfas<-function(fq){
  
  filter1 <- nFilter(threshold=1) # Keeps only reads without Ns
  filter2 <- polynFilter(threshold=20, nuc=c("A","T","G","C")) # Removes reads with >=20 of one nucleotide
  filter <- compose(filter1, filter2)
  
  return (fq[filter(fq)])
  
}

#filterFastq("f1.fq", "Test1z.fq", filter=filfas)

filtrim3<-function(fq,ffilter){
  # create ShortReadQ object
  my_read <- readFastq(fq)
  
  # apply the filter 
  filtered <- trim3q_filter(my_read, rm.3qual = 28)
  filterfile<-writeFastq(filtered, ffilter, mode="w", full=FALSE, compress=TRUE)
  return (filterfile)
  
}

filcom<-function(fq,ffilter){
  my_read <- readFastq(fq)
  # apply the filter
  filtered <- complex_filter(my_read)
  filterfile<-writeFastq(filtered, ffilter, mode="w", full=FALSE, compress=TRUE)
  return (filterfile)
  
}

filadpt<-function(fq,ffilter){
  my_read <- readFastq(fq)
  # trim adapter
  filtered <- adapter_filter(my_read, rc.R = TRUE)
  filterfile<-writeFastq(filtered, ffilter, mode="w", full=FALSE, compress=TRUE)
  return (filterfile)
}



filn<-function(fq,ffilter){
  my_read <- readFastq(fq)
  # trim adapter
  filtered <- n_filter(my_read, rm.N = 3)
  
  filterfile<-writeFastq(filtered, ffilter, mode="w", full=FALSE, compress=TRUE)
  return (filterfile)
  
}


filfix<-function(fq,ffilter){
  my_read <- readFastq(fq)
  
  # apply the filter 
  filtered3 <- fixed_filter(my_read, trim5 = 5)
  filterfile<-writeFastq(filtered3, ffilter, mode="w", full=FALSE, compress=TRUE)
  return (filterfile)
}



# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("FASTQ Filter"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    
    # Sidebar panel for inputs ----
    sidebarPanel(
     
    
      fileInput("file1", "Choose FASTQ File",
                multiple = FALSE,
                accept = c("text/fq",
                           ".fq")),
      # Horizontal line ----
      tags$hr(),
      radioButtons("readFile", "readFastq",
                   choices = c(None="none",readLines="readlines",Quality = "quality",
                               Reads = "reads",
                              Id = "id",Cycle_GC="cycle",Cycle_Base_CallsLine="cyclebase")
      ),
      
      
      # Horizontal line ----
      tags$hr(),
      
      textInput("textt", "Destination File ")  ,
      # Input: Select number of rows to display ----
      tags$hr(),
      h4("Methods"),
        actionButton("filter_fastqq", "filter_fastq"),
        actionButton("trim3q_filterr", "trim3q_filter"),
        actionButton("adapter_filterr", "adapter_filter"),
        actionButton("n_filterr", "n_filter"),
        actionButton("fixed_filterr", "fixed_filter"),
        actionButton("complex_filterr", "complex_filter")
      
    


),

    
      mainPanel(
     
      
      # Output: Data file ----
      tableOutput("contents")
    )
    
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  output$contents <- renderTable({
   
  req(input$file1)
    
    # when reading semicolon separated files,
    # having a comma separator causes `read.csv` to error
    tryCatch(
      {
        df <- readLines(input$file1$datapath)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    observeEvent( input$filter_fastqq,{ filterFastq(input$file1$datapath,input$textt, filter=filfas)
              #print(plotQualityProfile(input$file1$datapath))
    
      })
    
    observeEvent( input$trim3q_filterr,{ filtrim3(input$file1$datapath,input$textt) 
      #print(plotQualityProfile(input$textt))
     
      
      })
    observeEvent( input$adapter_filterr, { filadpt(input$file1$datapath,input$textt)
      #print(plotQualityProfile(input$textt))
      })
    observeEvent( input$n_filterr,{ filn(input$file1$datapath,input$textt) 
      #print(plotQualityProfile(input$textt))
      })
    observeEvent( input$fixed_filterr,{ filfix(input$file1$datapath,input$textt) 
      #print(plotQualityProfile(input$textt))
      })
    observeEvent( input$complex_filterr,{ filcom(input$file1$datapath,input$textt) 
      #print(plotQualityProfile(input$textt))
      })
    
    
    if(input$readFile == "none") {
      
      return()
    }
   else if(input$readFile == "reads") {
      
     df <- readFastq(input$file1$datapath)
     qs<-sread(df)[1:20]
     
    return(qs)
   
   }else if(input$readFile == "quality") {
     df <- readFastq(input$file1$datapath)
     qs<-as(quality(df), "matrix") [1:4,1:12]
     print(plotQualityProfile(input$file1$datapath))
  
     return(qs)
   }
    else if(input$readFile == "id") {
      df <- readFastq(input$file1$datapath)
      qs<-id(df)[1:20]
      return(qs)
    
    } 
    else if (input$readFile == "cycle"){
      qa <- rqcQA(input$file1$datapath)
   
      print(rqcCycleGCPlot(qa))
     
    }
    else if (input$readFile == "cyclebase"){
      qa <- rqcQA(input$file1$datapath)
      print(rqcCycleBaseCallsLinePlot(qa))
      
      
    }
    
 
    
    else {
      
      return(df)
    }
    
    
    
      
  })
}

shinyApp(ui, server)   
        
        
        
      
      
                                        