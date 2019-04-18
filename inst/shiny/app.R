library(magrittr)
library(shiny)

# to do: valiate GSE accession, 'GSE\\d+'

#" optional
#" in mainPanel
#       h3('current settings'),
#       uiOutput('current_settings', T),
#" in observeEvent(input$accession)
#       output$current_settings = renderUI(
#       	if (input$task == '1') {
#       		p('output dir:', code(parseDirPath(roots, input$output_dir)))
#       	} else {
#       		tagList(
#       			p('matrix file:', code(parseFilePaths(roots, input$matrix_file)$datapath)),
#       			p('SOFT file:', code(parseFilePaths(roots, input$soft_file)$datapath)),
#       			p('output dir:', code(parseDirPath(roots, input$output_dir)))
#       		)
#       	}
#       )

# Define UI for application that draws a histogram
ui <- fluidPage(
	tags$style('h2 {text-align: center;}'),
	tags$style('style + .row, .container-fluid .col-sm-6 {display: flex; justify-content: center;}'),
	tags$style('div.col-sm-12{max-width:600px;}'),
	tags$style('button a:hover {text-decoration: none}'),

	shinyjs::useShinyjs(),

	sidebarLayout(
	    sidebarPanel(
	    	width = 12,
	    	h3('prepare GSEA input files from microarray in GSE'),
	    	br(),
	    	selectInput(
	    	 	'task', h4('What do you want to do?'), selected = 1,
                choices = list(
                	'A. run App using GSE accession' = 1,
                	'B-1. download GSE matrix and SOFT file' = 2,
                	'B-2. run App using GSE matrix and SOFT file' = 3)
	    	),
			fluidRow(
				column(4, h4('GSE accession'), style = 'text-align:center'),
				column(6, textInput('accession', NULL, placeholder = 'GSE19161'))
			),
	    	div(id = 'fill_accession'),
			conditionalPanel(
				'input.task == "1"',
				helpText('All right, we will download the matrix and SOFT file for you. But we suggest you run', code('B-1'), 'and', code('B-2'), 'i.e., download the files yourself and use them as input.')
			),
			conditionalPanel(
				'input.task == "2"',
				helpText('Click the following buttons to download files, then run', code('B-2'), '.'),
				# 'If you can\'t download files, check whether you have filled correct', code('GSE accession'), '.'
				div(
					id = 'task2',
					fluidRow(
						column(6, actionButton('get_matrix', a('Download matrix file'))),
						column(6, actionButton('get_soft', a('Download SOFT file')))
					)
				)
			),
			conditionalPanel(
				'input.task == "3"',
				helpText('If you are confused, run', code('B-1'), 'first.'),
				fluidRow(
					column(6, shinyFiles::shinyFilesButton(
					    'matrix_file', 'Choose matrix file', 'Choose series matrix file', T)),
					column(6, shinyFiles::shinyFilesButton(
					    'soft_file', 'Choose SOFT file', 'Choose platform annotation (SOFT) file', T))
				),
			br()
			),
			conditionalPanel(
				'input.task == "1" || input.task == "3"',
				fluidRow(
					column(6, shinyFiles::shinyDirButton(
					    'output_dir', 'Chose output directory', 'where the result files would be put')),
					column(6, actionButton('run', 'Run', icon('play-circle')))
				),
				br(),
				div(id = 'instruction'),
				textInput('gene', h5('interested gene ', span('(optional)', style = 'font-size:80%')), placeholder = 'EIF4G2'),
		    	div(id = 'correct_gene'),
				helpText('If you want to identify gene sets correlated with a gene of interest (gene neighbors), you can fill its offical HUGO gene symbol in the aboving blank (multiple value should be separated by "', code(', '), '").', br(), 'Such as ', code('ANKLE2'), ' or ', code('BRCA1, BRCA2'), '.')
			)
		),
	    mainPanel(width = 0)
	)
)

set_element <- function(id, inner_html = '') {
	shinyjs::runjs(paste0("document.getElementById('", id, "').innerHTML = '", inner_html, "'", collapse = ''))
}

set_download_link <- function(parent_id, href) {
	shinyjs::runjs(paste0("a = document.querySelector('#", parent_id," > a'); a.href='", href, "'; a.download = ''", collapse = ''))
}

#' @param gene string.
# gene <- 'BRCA1, BRCA2'
# gene <- 'BRCB1, BRCE2'
validate_symbol <- function(gene) {
	non_symbol <- stringr::str_split(gene, ', ')[[1]] %>% setdiff(hgnc::hugo_symbol);
	if (length(non_symbol) > 0)
		paste0('<code>', non_symbol, '</code>', collapse = ', ') %>%
			paste0(., ' is not valid HUGO gene symbol!') %>% HTML()
	else
		NULL
}

validate_input_files <- function(matrix_file, soft_file, accession) {
	if (length(matrix_file) == 0) return('please select matrix file!')

	expect_matrix_basename = paste0(accession, '_series_matrix.txt.gz')
	if (basename(matrix_file) != expect_matrix_basename)
		return(paste0('matrix file name should be <code>', expect_matrix_basename, '</code>. please select another file!', collapse = ''))


	if (length(soft_file) == 0) return('please select SOFT file!')

	expect_soft_basename = paste0(accession, '_family.soft.gz')
	if (basename(soft_file) != expect_soft_basename)
		return(paste0('SOFT file name should be <code>', expect_soft_basename, '</code>. please select another file!', collapse = ''))

	return(NULL)
}



# Define server logic required to draw a histogramp
server <- function(input, output, session) {
	roots <- c(home = normalizePath('~')) # Downloads = normalizePath('~/Downloads'), root = '/root'
	shinyFiles::shinyFileChoose(input, 'soft_file', roots = roots, filetypes = 'gz')
	shinyFiles::shinyFileChoose(input, 'matrix_file', roots = roots, filetypes = 'gz')
	shinyFiles::shinyDirChoose(input, 'output_dir', roots = roots)

	shinyjs::runjs("document.getElementById('get_matrix').onclick = function() {document.querySelector('#get_matrix > a').click()}")
	shinyjs::runjs("document.getElementById('get_soft').onclick = function() {document.querySelector('#get_soft > a').click()}")


	observeEvent(input$accession, {
		if (input$accession == '')
			set_element('fill_accession', 'please fill <code>GSE accession</code>!')
		else
			set_element('fill_accession')

		set_download_link('get_matrix', rGEO::gse_matrix_ftp(input$accession))
		set_download_link('get_soft', rGEO::gse_soft_ftp(input$accession))
	})

	observeEvent(input$run, {
		#updateTextInput(session, 'accession', value = 'GSE19161')
		run_error  = F;
		gene_error = F;
		#" validate gene symbol
		gene <- input$gene;
		if (gene != '') {
			gene_message <- validate_symbol(gene);
			if (is.null(gene_message))
				set_element('correct_gene')
			else {gene_error = T; set_element('correct_gene', gene_message)}
		}
		#" validate output directory
		output_dir = shinyFiles::parseDirPath(roots, input$output_dir)
		if (length(output_dir) == 0) {run_error = T; set_element('instruction', 'please select output directory!')}
		#" validate input files
		if (input$task == '3') {
			matrix_file = shinyFiles::parseFilePaths(roots, input$matrix_file)$datapath %>% as.character()
			soft_file = shinyFiles::parseFilePaths(roots, input$soft_file)$datapath %>% as.character()
			file_messages <- validate_input_files(matrix_file, soft_file, input$accession)
			if (!is.null(file_messages)) {run_error = T; set_element('instruction', file_messages)}
		}
		if (!run_error) set_element('instruction');

		if (!run_error && input$accession != '' && !gene_error) {
			if (input$task == '1') {
				matrix_url <- rGEO::gse_matrix_ftp(input$accession);
				matrix_file <- file.path(tempdir(), basename(matrix_url));
				download.file(matrix_url, matrix_file);

				soft_url   <- rGEO::gse_soft_ftp(input$accession);
				soft_file <- file.path(tempdir(), basename(soft_url));
				download.file(soft_url, soft_file);
			}

			if (input$task %in% c('1', '3'))
				qGSEA::make_gsea_input(matrix_file, soft_file, output_dir, gene)
		}
	})
}


# Run the application
shinyApp(ui = ui, server = server)

# java -Xmx512m -cp /media/computer/opt/GSEA/gsea-3.0.jar xtools.gsea.Gsea -gmx /media/computer/opt/GSEA/msigdb_v6.1/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt -res GSE19161.txt -cls GSE19161.cls#EIF4G2 -chip GSE19161.chip -out . -collapse true -nperm 10 -permute phenotype -rpt_label GSE19161 -metric Pearson -gui false

