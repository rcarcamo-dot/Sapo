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
    query <- input$txt
    myquerypath <- file.path("../fastas/query.fa")
    if (startsWith(query, ">")){
      writeLines(query, myquerypath)
    } else {
      writeLines(paste0(">Query\n",query), myquerypath)
    }
    
    #this makes sure the fasta is formatted properly
    
    results_folder <- file.path("../results")
    myquerypath <- file.path("../fastas/query.fa")
    myquery <- readLines(myquerypath)
    header <- gsub(">", "", myquery[grep("^>", myquery)])
    query_length <- getLength(query)
    querystring <- toupper(c2s(query[[1]]))  #extract query sequence and capitalize it
    #Load templates fasta files
    gi_1_templates <- seqinr::read.fasta (file = '../fastas/GI_1.fa', seqtype = "DNA", as.string = TRUE)
    gi_2_templates <- seqinr::read.fasta (file = '../fastas/GI_2.fa', seqtype = "DNA", as.string = TRUE)
    gii_1_templates <- seqinr::read.fasta (file = '../fastas/GII_1.fa', seqtype = "DNA", as.string = TRUE)
    gii_3_templates <- seqinr::read.fasta (file = '../fastas/GII_3.fa', seqtype = "DNA", as.string = TRUE)
    gii_4_templates <- seqinr::read.fasta (file = '../fastas/GII_4.fa', seqtype = "DNA", as.string = TRUE)
    gii_5_templates <- seqinr::read.fasta (file = '../fastas/GII_5.fa', seqtype = "DNA", as.string = TRUE)
    gv_i_templates <- seqinr::read.fasta (file = '../fastas/GV_1.fa', seqtype = "DNA", as.string = TRUE)
    
    headers <- attr(list_files, "names")
    genogroups <- str_extract(headers, "G[IV]++")
    genogroups
    genotypes <- str_extract(headers, "G.+/")
    genotypes
    
    
    #List of the files that contain templates for each genotype
    list_files <- c(gi_1_templates, gi_2_templates, gii_1_templates, gii_3_templates, gii_4_templates, gii_5_templates, gv_i_templates)
    genotypes <- attr(gi_1_templates, "names")
    for (i in 1:length(genotypes)){
      #print(genotypes[i])
      genotype <- str_extract(genotypes[i], ".G.")
      print(genotype)
    }
      
    attr(gi_1_templates, "names")
    genotype <- gsub(">", "", list_files[1][grep("^/G*.*/", list_files[1])])
    list_genotypes <- c('../fastas/GI_1.fa', '../fastas/GI_2.fa', '../fastas/GII_1.fa', '../fastas/GII_3.fa', '../fastas/GII_4.fa', '../fastas/GII_5.fa', '../fastas/GV_1.fa')
    #This iterates 7 times over each file, but we have to manage how to print the messages without this loop
    for (r in 1:length(list_genotypes)){ #Iterate through all the template files
      list_templates <- list()
      my_positive_string <- paste("This sequence belongs to", substr(list_genotypes[[r]], start = 11, stop = 15), "genotype.", sep = " ")
      my_negative_string <- paste("This sequence does not belong to", substr(list_genotypes[[r]], start = 11, stop = 15), "genotype.", sep = " ")
      for(i in 1:length(list_files)) { 
        nam <- paste("seq", i, sep = "")
        assign(nam, toupper(c2s(list_files[[i]])))
        attr(list_files[[i]], "names")
        #print(attr(list_files[[i]], "names"))
        #print (list_files[[i]])
        print(nam)
        genotype <- gsub(">", "", list_files[i][grep("^>*/", list_files[i])])
        list_templates[[nam]] <- (c2s(list_files[[i]]))
        write.fasta(sequences= toupper(c2s(list_files[[i]])), names=nam, file.out="../fastas/template.fa")
        fasta_file_names <- c(myquerypath, '../fastas//template.fa')
        combined_fasta <- fasta.combine(fasta_file_names = fasta_file_names)
        write.fasta(sequences=combined_fasta,names=names(combined_fasta),file.out="../fastas/query_template.fa")
        query_template <- file.path("../fastas/query_template.fa")
        query_template = readLines("../fastas/query_template.fa")
        query_template = toupper(query_template)
        writeLines(query_template, "../fastas/query_template.fa")
        mysequences = readDNAStringSet("../fastas/query_template.fa")
        #print(toupper(mysequences))
        myClustalWAlignment <- msa(mysequences, "ClustalW", type = "dna")
        myClustalWAlignment <- msaConvert(myClustalWAlignment, "seqinr::alignment")
        d <- dist.alignment(myClustalWAlignment, matrix = "identity")
        #print(d)
        if (d <= 0.169){
          genotype <- my_positive_string
          break
        } else {
          print (my_negative_string, quote = FALSE)
          next
        }
      }
    }
    summary_table <- data.frame(
      Name    = header,
      Length  = query_length,
      Genotype = genotype
    )
    summary_table
  })
    
    

    
  })# sanitize.text.function = function(x) x)
  
#})

