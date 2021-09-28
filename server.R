#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(Biostrings)
library(crisscrosslinker)
library(ape)
library(phangorn)
library(seqinr)
library(msa)
library(tidyverse)
library(ggplot2)
library(ggtree)
library(ggrepel)
library(treeio)

shinyServer(function(input,output,session){
  
  
  output$selected_var <- renderTable({
    #validate(need((input$txt | input$file1) != '', 'Enter one type of input!'))
    if (input$txt != ""){
      query <- input$txt
      myquerypath <- file.path("../fastas/query.fa")
      if (startsWith(query, ">")){
        writeLines(query, myquerypath)
      } else {
        writeLines(paste0(">Query\n",query), myquerypath)
      }
      myquery <- readLines(myquerypath)
    } else if (input$txt == "" && input$file1 != "") {
    observe({
      file1 = input$file1
      if (is.null(file1)) {
        return(NULL)
      }
      data1 = read.fasta(file1$datapath)
      dnabin_to_fasta <- lapply(data1, function(x) as.character(x[1:length(x)]))
      cat(file="../fastas/query.fa", paste(paste0(">",names(dnabin_to_fasta)),
                                 sapply(dnabin_to_fasta, paste, collapse=""), sep="\n"), sep="\n");
    })
      myquerypath <- file.path("../fastas/query.fa")
      myquery <- readLines(myquerypath)
    } else {
      stop ("Enter one kind of input")
    }
    
    #this makes sure the fasta is formatted properly
    
    results_folder <- file.path("../results")
    myquerypath <- file.path("../fastas/query.fa")
    myquery <- readLines(myquerypath)
    header <- gsub(">", "", myquery[grep("^>", myquery)])
    query_length <- getLength(myquery)
    #querystring <- toupper(c2s(query[[1]]))  #extract query sequence and capitalize it
    
    #Load templates fasta files
    gi_1_templates <- seqinr::read.fasta (file = '../fastas/GI_1.fa', seqtype = "DNA", as.string = TRUE)
    gi_2_templates <- seqinr::read.fasta (file = '../fastas/GI_2.fa', seqtype = "DNA", as.string = TRUE)
    gi_3_templates <- seqinr::read.fasta (file = '../fastas/GI_3.fa', seqtype = "DNA", as.string = TRUE)
    gi_4_templates <- seqinr::read.fasta (file = '../fastas/GI_4.fa', seqtype = "DNA", as.string = TRUE)
    gi_5_templates <- seqinr::read.fasta (file = '../fastas/GI_5.fa', seqtype = "DNA", as.string = TRUE)
    gi_6_templates <- seqinr::read.fasta (file = '../fastas/GI_6.fa', seqtype = "DNA", as.string = TRUE)
    gi_7_templates <- seqinr::read.fasta (file = '../fastas/GI_7.fa', seqtype = "DNA", as.string = TRUE)
    gii_1_templates <- seqinr::read.fasta (file = '../fastas/GII_1.fa', seqtype = "DNA", as.string = TRUE)
    gii_2_templates <- seqinr::read.fasta (file = '../fastas/GII_2.fa', seqtype = "DNA", as.string = TRUE)
    gii_3_templates <- seqinr::read.fasta (file = '../fastas/GII_3.fa', seqtype = "DNA", as.string = TRUE)
    gii_4_templates <- seqinr::read.fasta (file = '../fastas/GII_4.fa', seqtype = "DNA", as.string = TRUE)
    gii_5_templates <- seqinr::read.fasta (file = '../fastas/GII_5.fa', seqtype = "DNA", as.string = TRUE)
    gii_6_templates <- seqinr::read.fasta (file = '../fastas/GII_6.fa', seqtype = "DNA", as.string = TRUE)
    gii_7_templates <- seqinr::read.fasta (file = '../fastas/GII_7.fa', seqtype = "DNA", as.string = TRUE)
    gii_8_templates <- seqinr::read.fasta (file = '../fastas/GII_8.fa', seqtype = "DNA", as.string = TRUE)
    giv_i_templates <- seqinr::read.fasta (file = '../fastas/GIV_1.fa', seqtype = "DNA", as.string = TRUE)
    gv_i_templates <- seqinr::read.fasta (file = '../fastas/GV_1.fa', seqtype = "DNA", as.string = TRUE)
    
    #List of the files that contain templates for each genotype
    list_files <- c(gi_1_templates, gi_2_templates, gi_3_templates, gi_4_templates, gi_5_templates,
                    gi_6_templates, gi_7_templates, gii_1_templates, gii_2_templates, gii_3_templates,
                    gii_4_templates, gii_5_templates, gii_6_templates, gii_7_templates, gii_8_templates,
                    giv_i_templates, gv_i_templates)
    
    headers <- attr(list_files, "names")
    list_genogroups <- str_extract(headers, "G[IV]++")
    list_genotypes <- gsub('.{1}$', '', str_extract(headers, "G.+/"))
    #length(genotypes)
    #length(genogroups)
    

    #This iterates 7 times over each file, but we have to manage how to print the messages without this loop
    
    for(i in 1:length(list_files)) { 
      list_templates <- list()
      nam <- paste("seq", i, sep = "")
      assign(nam, toupper(c2s(list_files[[i]])))
      header <- gsub(">", "", myquery[grep("^>", myquery)])
      write.fasta(sequences= toupper(c2s(list_files[[i]])), names=nam, file.out="../fastas/template.fa")
      fasta_file_names <- c(myquerypath, '../fastas//template.fa')
      combined_fasta <- fasta.combine(fasta_file_names = fasta_file_names)
      write.fasta(sequences=combined_fasta,names=names(combined_fasta),file.out="../fastas/query_template.fa")
      query_template <- file.path("../fastas/query_template.fa")
      query_template = readLines("../fastas/query_template.fa")
      query_template = toupper(query_template)
      writeLines(query_template, "../fastas/query_template.fa")
      mysequences = readDNAStringSet("../fastas/query_template.fa")
      myClustalWAlignment <- msa(mysequences, "ClustalW", type = "dna")
      myClustalWAlignment <- msaConvert(myClustalWAlignment, "seqinr::alignment")
      d <- dist.alignment(myClustalWAlignment, matrix = "identity")
      #print(d)
      if (d <= 0.488){
        genogroup <- list_genogroups[i]
        print(genogroup)
        if (d <= 0.169){
          genotype <- list_genotypes[i]
          print(genotype)
          break
        }
      }
    }
    summary_table <- data.frame(
    Name    = header,
    Length  = query_length,
    Genogroup = genogroup,
    Genotype = genotype
  )
  summary_table
})
    
    

    
})# sanitize.text.function = function(x) x)
  
#})

