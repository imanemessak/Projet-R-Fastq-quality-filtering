library("shiny")
library("rBLAST")
library("seqinr")
library("shinyjs")
library("readr")
library("annotate")
library("R4X")
library("httr") 
library("jsonlite")
library("RCurl")
library("stringr")
library("lubridate")
library("Biostrings")
library("shinydashboard")
library('shiny')
library('Biostrings')
library('ShortRead')
library('FastqCleaner')
library('ggplot2')
library("Rqc")

if (interactive()) {
  
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



  ui<-fluidPage(
    dashboardPage(
      dashboardHeader(title =  "Medical Genomic Analysis",
                      titleWidth = 400,
                      dropdownMenu(
                        type = "notifications", 
                        icon = icon("question-circle"),
                        badgeStatus = NULL,
                        headerText = HTML("<h3> Search keyword for sequences in  the NCBI Database:</h3>
                                          <table style=`width:100%`>
                                          <tr>
                                          <th>Argument</th>
                                          <th>Example</th>
                                          <th>Restricts your search to sequences</th></tr>
                                          <tr>
                                          <td>  “AC=” </td>
                                          <td>“AC=NC_001477”</td>
                                          <td>With a particular accession number</td></tr>
                                          <tr>
                                          <td>“SP=”</td>
                                          <td>“SP=Chlamydia”</td>
                                          <td>From a particular organism or taxon</td></tr>
                                          <tr>
                                          <td>“M=”</td>
                                          <td>“M=mRNA”</td>
                                          <td>Of a specific type (eg. mRNA)</td></tr>
                                          <tr>
                                          <td> “J=”</td>
                                          <td>“J=Nature”</td>
                                          <td>Described in a paper published in a particular journal</td></tr>
                                          <tr>
                                          <td> “R=”</td>
                                          <td>“R=Nature/460/352”</td>
                                          <td>Described in a paper in a particular journal, volume and start-page</td></tr>
                                          <tr>
                                          <td> “AU=”</td>
                                          <td>“AU=Smith”</td>
                                          <td>Described in a paper, or submitted to NCBI, by a particular author </td></tr>
                                          </table>
                                          
                                          <table style=`width:100%`  >
                                          <tr>
                                          <th> Input keyword</th>
                                          <th>Searches for sequences</th>
                                          </tr>
                                          <tr>
                                          <td>  “AC=NC_001477”</td>
                                          <td>With accession number NC_001477</td></tr>
                                          <tr>
                                          <td>  “R=Nature/460/352”</td>
                                          <td>Published in Nature 460:352-358</td></tr>
                                          <tr>
                                          <td>  “SP=Chlamydia trachomatis”</td>
                                          <td>From the bacterium Chlamydia trachomatis</td></tr>
                                          <tr>
                                          <td>   “AU=Berriman”</td>
                                          <td>Published in a paper, or submitted to NCBI, by someone called Berriman</td></tr>
                                          <tr>
                                          <td>   “K=flagellin OR K=fibrinogen”</td>
                                          <td>Which have the keyword ‘flagellin’ or ‘fibrinogen’</td></tr>
                                          <tr>
                                          <td>   “SP=Mycobacterium leprae AND K=dnaA”</td>
                                          <td>Which are from M. leprae, and have the keyword “dnaA”</td></tr>
                                          <tr>
                                          <td>   “SP=Homo sapiens AND K=colon cancer”</td>
                                          <td>Which are from human, and have the keyword “colon cancer”</td></tr>
                                          <tr>
                                          <td>   “SP=Homo sapiens AND K=malaria”</td>
                                          <td>Which are from human, and have the keyword “malaria”</td></tr>
                                          <tr>
                                          <td>  “SP=Homo sapiens AND M=mrna”</td>
                                          <td>Which are mRNA sequences from human</td></tr>
                                          <tr>
                                          <td>  “SP=Bacteria”</td>
                                          <td>Which are sequences from Bacteria</td>
                                          </tr>
                                          </table>") )  
                        ),
      dashboardSidebar(width=400,
                       sidebarMenu(
                         menuItem("About", icon = icon("info"), tabName = "about"),
                         menuItem("Quality Controle/Filitring", icon = icon("filter"), tabName = "filtre"
                                  
                         ),
                         menuItem("Search/Alignement BLAST", icon = icon("list"), tabName = "keyalign",
                                  
                         
                         menuSubItem("KeyWord", icon = icon("search"), tabName = "keyword"
                                     
                         ),
                         menuSubItem("DNA sequence", icon = icon("dna"), tabName = "dna"),
                         menuSubItem("File", icon = icon("file-import"), tabName = "file-import"
                                     
                         ),
                         menuSubItem("API", icon = icon("exchange-alt"), tabName = "API"
                                     
                         )
                         )
                       )),
      dashboardBody(
        #CSS 
        tags$head(tags$style(HTML('
                                  .main-header .logo {
                                  font-family: "Georgia", Times, "Times New Roman", serif;
                                  font-weight: bold;
                                  font-size: 24px;
                                  
                                  }
                                  th {
                                  border-bottom	: 1px solid black;
                                  border-collapse: collapse;
                                  }
                                  td, th {
                                  padding: 3px;
                                  }
                                  .navbar-nav>.messages-menu>.dropdown-menu, .navbar-nav>.notifications-menu>.dropdown-menu, .navbar-nav>.tasks-menu>.dropdown-menu {
                                  width: 690px;
                                  
                                  }
                                  .set2{
                                  font-family: "Georgia", Times, "Times New Roman", serif;
                                  font-weight: italic;
                                  font-size: 14px;
                                  }
                                  .set2 form.well { 
                                  background-color: #fff7e6;
                                  border: 1px inset #b3b3ff;
                                  }
                                  .set1{
                                  font-family: "Georgia", Times, "Times New Roman", serif;
                                  font-weight: italic;
                                  font-size: 15px;
                                  }
                                  .h3, h3 {
                                  font-family: "Georgia", Times, "Times New Roman", serif;
                                  font-weight: bold;
                                  font-size: 15px;
                                  }
                                  .set1 form.well { 
                                  background-color: #fff7e6;
                                  border: 1px inset #b3b3ff;
                                  }
                                  .skin-blue .main-header .navbar .sidebar-toggle:hover{
                                  background-color: #ffccb3 ;
                                  }
                                  
                                  .skin-blue .sidebar-menu>li.active>a, .skin-blue .sidebar-menu>li:hover>a {
                                  color: #ff884d;
                                  background: #fff2e6;
                                  border-left-color: #ffccb3;
                                  }
                                  .skin-blue .sidebar a {
                                  color: #ffaa80;
                                  font-weight: bold;
                                  }
                                  .skin-blue .main-header .navbar {
                                  background-color: #ffddcc;
                                  }
                                  .container-fluid{
                                  background-color: #fff7e6 !important;
                                  }
                                  
                                  .content-wrapper, .right-side {
                                  
                                  background-color: #fff7e6 !important;
                                  }
                                  .main-sidebar, .skin-blue .wrapper {
                                  background-color: #fff7e6 !important;
                                  }
                                  
                                  .skin-blue .main-header .logo {
                                  background-color: #ffccb3;
                                  }
                                  .skin-blue .main-sidebar {
                                  background-color: #fff7e6!important;
                                  }
                                  .skin-blue .treeview-menu>li>a {
                                  color: #ff884d;
                                  background: #fff2e6;
                                  border-left-color: #ffccb3;
                                  }
                                  .skin-blue .sidebar-menu>li>.treeview-menu {
                                  border-left-color: #ffccb3;
                                  background: #fff7e6 !important;
                                  }
                                  .skin-blue .main-header .logo:hover {
                                  background-color: #ffddcc;
                                  }
                                  body{
                                  background-color: #fff7e6 !important;
                                  }
                                  .box-header {
                                  color: #ffaa80;
                                  }
                                  .box {
                                  background: #fff;
                                  border-top: 3px solid #ff9966;
                                  
                                  }
                                  
                                  '))),
        useShinyjs(),
        
        div(class="set1",
            
            tabItems(
              tabItem(tabName = "about",
                      box(width = 800,title = "About",
                          h4("Today, the sequence is the first known information concerning a protein, even before it has been purified and its role has been known. Genomic analysis of this sequence allows us to better understand the diversity of living things, to build phylogenetic trees or to identify genes associated with diseases.")
                         , h4(
                            
                            
                            
                          )
                      )
                      
                      
              ),
              
                            tabItem(tabName = "API",
                      box(width = 800,title = "Web API-BLAST",
                          fileInput("infile", label = h3("File input")),
                          textAreaInput("outfile", label = h3("File output "), placeholder =  "Enter output file `namefile.out`..."),
                          actionButton("blast"," Blast+",icon = icon("check-double")),
                          actionButton("reset", "Reset ALL",icon = icon("broom")),
                          actionButton("resettab", "Reset Main Panel",icon = icon("trash-alt"))
                      )
                      
                      
              ),
              
              tabItem(tabName = "dna",
                      box(width = 800,title = "Sequence",
                          textAreaInput("text", label = h3("Sequence input"), placeholder =  "Enter sequence..."),
                          actionButton("GO"," Blast",icon = icon("check")),
                          actionButton("resetSequence", "Reset Sequence",icon = icon("broom")),
                          actionButton("resettab1", "Reset Main Panel",icon = icon("trash-alt"))
                      )),
              
              tabItem(tabName = "keyword",
                      box(width = 800,title = "Search By Keyword",
                          textAreaInput("mot", label = h3("Search"), placeholder =  "Enter keyword..."),
                          selectInput("select", label = h3("ACNUC sub-database"),
                                      choices = list("genbank","embl","emblwgs","swissprot","ensembl","hogenom7dna","hogenom7","hogenom","hogenomdna","hovergendna",
                                                     "hovergen","hogenom5","hogenom5dna","hogenom4","hogenom4dna","homolens","homolensdna","hobacnucl","hobacprot",
                                                     "phever2","phever2dna","refseq","refseq16s","greviews","bacterial","archaeal","protozoan","ensprotists","ensfungi",
                                                     "ensmetazoa","ensplants","ensemblbacteria","mito","polymorphix","emglib","refseqViruses","ribodb","taxodb"),selected ="genbank"),
                          textAreaInput("name", label = h3("File name"), placeholder =  "Enter file name `namefile.fasta`..."),
                          actionButton("rec","Search",icon = icon("search")),
                          actionButton("resetsearch", "Reset All",icon = icon("broom")),
                          actionButton("resettab2", "Reset Main Panel",icon = icon("trash-alt"))
                      )
                      
              ),
              
              
              tabItem(tabName = "file-import",
                      box(width = 800,title = "Upload FASTA File", 
                          fileInput("file", label = h3("File input")),
                          actionButton("go"," Blast",icon = icon("check")),
                          actionButton("resetfile", "Reset File",icon = icon("broom")),
                          actionButton("resettab3", "Reset Main Panel",icon = icon("trash-alt"))
                      )),
              tabItem(tabName = "filtre",
                      box(width = 800,title = "FASTQ Filter",
                          fileInput("file1", "Choose FASTQ File", multiple = FALSE,accept = c("text/fq",".fq")),
                          radioButtons("readFile", "Read Fastq",
                                       choices = c(None="none",ReadLines="readlines",Quality = "quality",
                                                   Reads = "reads",
                                                   Id = "id",Cycle_GC="cycle",Cycle_Base_CallsLine="cyclebase")
                          ),
                          
                          
                          # Horizontal line ----
                          tags$hr(),
                          textAreaInput("textt", label = h3("Destination File name"), placeholder =  "Enter file name `namefile.fq`..."),
                          
                          # Input: Select number of rows to display ----
                          
                          h3("Filtring Methods"),
                          actionButton("filter_fastqq", "Filter Fastq",icon = icon("filter")),
                          actionButton("trim3q_filterr", "Trim3q Filter",icon = icon("filter")),
                          actionButton("adapter_filterr", "Adapter Filter",icon = icon("filter")),
                          actionButton("n_filterr", "N Filter",icon = icon("filter")),
                          actionButton("fixed_filterr", "Fixed Filter",icon = icon("filter")),
                          actionButton("complex_filterr", "Complex Filter",icon = icon("filter")),
                          tags$hr(),
                          actionButton("resetsearch", "Reset All",icon = icon("broom")),
                          actionButton("resettab2", "Reset Main Panel",icon = icon("trash-alt"))
                      )
                      
              )
              
              
              
              
              #end tabItems
            )
        ),#end div
        div(class="set2",
            
            mainPanel(width = 800,
                      tableOutput("contents"),
                      verbatimTextOutput("search"),
                      dataTableOutput("table"),
                      plotOutput("graphe2"),
                      plotOutput("graphe")
            ) 
        )#end div
        )#enddashboardbody
        )#enddashboard
        )#endfluidpage
  
  
  server <- function(input, output) {
    file_RNA<-function(x){
      Sys.which("blastn")
      library(rBLAST)
      bl <- blast(db="./16SMicrobialDB/16SMicrobial")
      chemin<-getwd()
      nam<-paste(chemin,"/",x$name,sep = "")
      seq <- readRNAStringSet(nam)
      names(seq) <-  sapply(strsplit(names(seq), " "), "[", 1)
      cl <- predict(bl, seq[1,])
      return(cl)
    }
    sequence<- function(y) {
      Sys.which("blastn")
      library(rBLAST)
      bl <- blast(db="./16SMicrobialDB/16SMicrobial")
      seqe<-RNAStringSet(y)
      cl<- predict(bl, seqe)
      return(cl)
    }
    
    observeEvent(input$go,{
      file<-input$file
      cl<-file_RNA(file)
      output$table<-renderDataTable({ 
        cl})
      output$graphe<-renderPlot({
        barplot(cl$Perc.Ident,main = "SubjectID by Perc.Ident",names.arg = cl$SubjectID,col = "#ffddcc")})
      output$graphe2<-renderPlot({
        barplot(cl[1:10,3],main = "10 first SubjectID by Perc.Ident",names.arg = cl[1:10,2],col = "#ffddcc")})
    })
    observeEvent(input$GO,{
      fill<-input$text
      cl<-sequence(fill)
      output$table<-renderDataTable({ 
        cl
      })
      output$graphe<-renderPlot({
        barplot(cl$Perc.Ident,main = "SubjectID by Perc.Ident",names.arg = cl$SubjectID,col = "#ffddcc")})
      output$graphe2<-renderPlot({
        barplot(cl[1:10,3],main = "10 first SubjectID by Perc.Ident",names.arg = cl[1:10,2],col = "#ffddcc")})
      
    })
    search_base<-function(){
      library(seqinr)
      BANK<-input$select 
      query<-input$mot
      choosebank(BANK)
      humtRNAs<-query("humtRNAs", query)
      myseqs <- getSequence(humtRNAs)         
      mynames <- getName(humtRNAs)
      namefile<-input$name
      write.fasta(myseqs, mynames, file.out=namefile)
      closebank()
      filefasta<-read.fasta(namefile)
      return(filefasta)
    }
    observeEvent(input$rec,{
      mot_cle<-search_base()
      
      output$search<-renderPrint({ 
        mot_cle
      })
    })
    observeEvent(input$resetfile, {
      shinyjs::reset("file")
    })
    observeEvent(input$resetsearch, {
      shinyjs::reset("mot")
      shinyjs:: reset("name")
      shinyjs::reset("select")
    })
    observeEvent(input$resetSequence, {
      shinyjs:: reset("text")
    })
    cmdCreate <- function(inf, outf){
      chemin<-getwd()
      name<-paste(chemin,"/",inf,sep = "")
      names<-paste(chemin,"/",outf,sep = "")
      paste("blastn -db nr -query ", name, " -remote -out ",names, sep = "")
    }
    observeEvent(input$blast,{
      cm<-cmdCreate(input$infile,input$outfile)
      t<- system(cm)
      output$search<-renderPrint({
        t
        ff<- read_file(input$outfile)
        test<-str_split(ff,"\n")
        test
      })
      
    })
    
    observeEvent(input$reset, {
      shinyjs::reset("infile")
      shinyjs:: reset("outfile")
    })
    observeEvent(input$resettab, {
      output$table<-renderDataTable({data=NULL})
      output$search<-renderPrint({data=NULL})
      
    }) 
    observeEvent(input$resettab1, {
      output$table<-renderDataTable({data=NULL})
      output$search<-renderPrint({data=NULL})
      
    }) 
    observeEvent(input$resettab2, {
      output$table<-renderDataTable({data=NULL})
      output$search<-renderPrint({data=NULL})
      
    }) 
    observeEvent(input$resettab3, {
      output$table<-renderDataTable({data=NULL})
      output$search<-renderPrint({data=NULL})
      
    }) 
    
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
  }

