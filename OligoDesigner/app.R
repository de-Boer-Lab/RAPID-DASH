library(shiny)
library(DT)
library(Biostrings)
library(shinythemes)


# Primer prefix templates
FWD_TEMPLATE <- "ATAAGGATCCGGTCTCA%sGTAAAACGACGGCCAGT"
REV_TEMPLATE <- "ATAATGTACAGGTCTCT%sAGGAAACAGCTATGACCATG"

# Terminal and intermediate overhangs
terminal_overhangs <- c("GGTA", "CGTA")
intermediate_overhangs <- c(
  "TAGA", "CTCC", "ATCA", "CTGA", "AGCG", "AAGG", "CATC",
  "ACCT", "GCGA", "ACTG", "ATAC", "GAAA", "AAAT", "GAGC",
  "CGGC", "CCCA", "CCAC", "CAAG", "AATA")

# Reverse complement helper
rev_comp <- function(seq) {
  as.character(reverseComplement(DNAString(seq)))
}

# Function to generate all gRNA oligos and primer pairs
generate_gRNA_array <- function(spacers, spacer_names, num_units) {
  spacers <- toupper(trimws(spacers))
  spacers <- spacers[spacers != ""]
  
  if (length(spacers) != num_units) {
    return(data.frame(Error = "⚠️ Number of spacers doesn't match gRNA unit count ⚠️"))
  }
  
  if (length(spacer_names) != num_units) {
    return(data.frame(Error = "⚠️ Number of gRNA names doesn't match gRNA unit count ⚠️"))
  }
  if (length(spacers) != length(spacer_names)) {
    return(data.frame(Error = "⚠️ Number of spacers doesn't match gRNA names ⚠️"))
  }
  
  # Validate all spacers
  valid <- sapply(spacers, function(s) nchar(s) == 20 && grepl("^[ACGTN]+$", s))
  if (!all(valid)) {
    return(data.frame(Error = "⚠️ All spacers must be 20bp A/C/G/T/N only ⚠️"))
  }
  
  # Assign overhangs
  if (num_units - 2 > (length(intermediate_overhangs))) {
    return(data.frame(Error = "⚠️ Not enough unique intermediate overhangs. You can only design up to 20 gRNA units  ⚠️"))
  }
  
  # Hard limit to 20 gRNAs
  if (num_units > 20) {
    return(data.frame(Error = "⚠️ You can only design up to 20 gRNA units ⚠️"))
  }
  
  overhangs <- c(
    terminal_overhangs[1],
    sample(intermediate_overhangs, num_units - 1),
    terminal_overhangs[2]
  )
  fwd_ovhgs <- overhangs[1:(length(overhangs) - 1)]
  rev_ovhgs <- sapply(overhangs[2:length(overhangs)], rev_comp)

  rows <- vector("list", length = num_units)
  
  for (i in seq_len(length(fwd_ovhgs))) {
    spacer <- spacers[i]
    oligo <- paste0("GTGGAAAGGACGAAACACCg", spacer, "gttttagagctaGAAAtag")
    
    fwd_ovhg <- fwd_ovhgs[i]
    rev_ovhg <- rev_ovhgs[i]
    
    fwd_primer <- sprintf(FWD_TEMPLATE, fwd_ovhg)
    rev_primer <- sprintf(REV_TEMPLATE, rev_ovhg)
    
    rows[[i]] <- list(
      gRNA = spacer_names[i],
      Spacer = spacer,
      Oligo = oligo,
      Fwd_Primer = fwd_primer,
      Rev_Primer = rev_primer,
      Fwd_Overhang = fwd_ovhg,
      Rev_Overhang = rev_ovhg
    )
  }
  
  results <- do.call(rbind.data.frame, rows)
  return(results)
}

generate_idt_ready_table <- function(results) {
  if (!ncol(results) == 7) {
    return(data.frame(""))
  }
  clean_names <- gsub("\\s+", "_", results$gRNA)
  
  data.frame(
    Name = c(paste0(clean_names, "_spacer"),
             paste0(clean_names, "_fwd"),
             paste0(clean_names, "_rev")),
    Sequence = c(results$Oligo, results$Fwd_Primer, results$Rev_Primer),
    stringsAsFactors = FALSE
  )
}


ui <- fluidPage(
  
  theme = shinytheme("cosmo"),
  
  titlePanel("🧬 RAPID-DASH Oligo Designer"),
  
  sidebarLayout(
    sidebarPanel(
      class = "sidebar",
      
      h4("🔧 Guide Array Configuration"),
      
      numericInput("gRNAunits", "Number of gRNAs in array:", value = NA, min = 1, max = 20),
      textInput("spacer_names", "gRNA names (comma-separated):", placeholder = "Name of the gRNAs (comma separated)"),
      textInput("spacers", "gRNA spacer sequences:", placeholder = "gRNA sequences (5'-3') (comma separated)"),
      
      actionButton("generate", "🚀 Generate Oligos"),
      
      br(), br(),
      p("This tool implements the oligo design for gRNA array assembly described in:"),
      tags$em("RAPID-DASH: Single-Day Assembly of Guide RNA Arrays for Multiplexed CRISPR-Cas9 Applications"),
      br(),
      HTML('<a href="https://www.biorxiv.org/content/10.1101/2025.04.09.648054v1" target="_blank">🔗 View Preprint</a>'),
      br(), br(),
      HTML('<div style="text-align: center;"><img src="rapidash.jpg" width="300"></div>')
    ),
    
    mainPanel(
      class = "mainpanel",
      DTOutput("oligo_table"),
      uiOutput("download_ui"),
    )
  )
)

# Server
server <- function(input, output) {
  
  observeEvent(input$generate, {
    
    if (is.na(input$gRNAunits) || input$spacers == "" || input$spacer_names == "") {
      return(NULL)
    }
    
    spacers <- unlist(strsplit(input$spacers, ","))
    spacer_names <- unlist(strsplit(input$spacer_names, ","))
    results <- generate_gRNA_array(spacers, spacer_names, input$gRNAunits)
    
  output$oligo_table <- renderDT({
      datatable(results, options = list(pageLength = 30), rownames = FALSE)
    })
  
  output$download_table <- downloadHandler(
      "RAPID-DASH-Oligos.tsv",
      content = function(file){write.table(results, file, sep = "\t", quote = F)}
    )
  
  output$download_ui <- renderUI({
    tagList(
      downloadButton("download_table", "Download Table"),
      br(), br(),
      h4("Copy-paste table to order from IDT:"),
      
      tags$div(
        style = "position: relative;",
        tags$pre(id = "clipboard_text", style = "white-space: pre-wrap; background: #f9f9f9; padding: 10px; border: 1px solid #ccc; border-radius: 5px;",
                 textOutput("idt_clipboard"))
      )
    )
  })
  
  idt_df <- generate_idt_ready_table(results)
  
  output$idt_clipboard <- renderText({
    df <- generate_idt_ready_table(results)
    paste(apply(df, 1, function(x) paste(x, collapse = ",")), collapse = "\n")
  })
  
  })
}

# Run App
shinyApp(ui = ui, server = server)
