# .txt --------------

#' @title write GSEA expression dataset input file
#'
#' @param matrix tibble. usually the return value of [rGEO::read_gse_matrix()], `NULL` is okay.
#' @param path string. path to output file
#'
#' @return logical scalar
#'
#' @seealso [read_txt]
#'
#' @examples NULL
#'
#' @export

write_txt <- function(matrix, path) {
	if (is.null(matrix)) return(F)
	if (nrow(matrix) == 0L) {warning('empty matrix'); return(F)}
	if (ncol(matrix) < 3L) {warning('at least two samples is needed for GSEA');return(F)}

	matrix %>%
		tibble::add_column(Description = 'NA', .after = 1) %>%
		readr::write_tsv(path, na = '');

	return(T)
}



#' @title read .txt file info a `matrix`
#'
#' @return tibble. identical to [rGEO::read_gse_matrix()]
#'
#' @seealso [write_txt]
#'
#' @export

# file <- 'tests/testthat/output/GSE51280.txt'
read_txt <- function(file) {
	if (!file.exists(file)) return(NULL);

	readr::read_tsv(
		file,
		col_types = readr::cols(
			.default = readr::col_double(),
			ID_REF = readr::col_character(),
			Description = readr::col_skip()
		)
	) %>% {attr(., 'spec') <- NULL; .}
}






# .chip --------------

#' @title write GSEA microarray chip annotation input file
#'
#' @section to do: fetch description from hgnc later
#'
#' @param chip tibble. usually the return value of [rGEO::read_gse_soft()], `NULL` is okay.
#' @param path string. path to output file
#'
#' @return logical scalar
#'
#' @seealso [read_chip]
#'
#' @examples `NULL`
#'
#' @export

write_chip <- function(chip, path) {
	if (is.null(chip)) return(F)
	if (nrow(chip) == 0L) {warning('empty chip'); return(F)}

	chip %>%
		dplyr::select('Probe Set ID' = 1, 'Gene Symbol' = 2) %>%
		tibble::add_column('Gene Title' = 'NA') %>%
		readr::write_tsv(path);

	return(T)
}



#' @title read .chip file info a `chip`
#'
#' @return tibble. identical to [read_gse_soft]
#'
#' @seealso [write_chip]
#'
#' @export

# file <- 'tests/testthat/output/GSE51280.chip'
read_chip <- function(file) {
	if (!file.exists(file)) return(NULL);

	readr::read_tsv(file, T, 'cc-') %>% dplyr::select(ID_REF = 1, symbol = 2) %>%
		{attr(., 'spec') <- NULL; .}
}






# .cls --------------

#' @title write GSEA phenotype data input file
#'
#' @description only support continuous phenotype for now
#'
#' @details Although GSEA doesn't require every variable in a .cls file to be the same length, I do. Since they should be used with the same expression data (such as GSE008.txt + GSE008.chip + GSE008.cls, on the other hand, let GSE006.txt and GSE007.txt use different variable in the same .cls file is definitely insane)
#'
#' @param phenotype tibble. usually the return value of [make_phenotype]
#' @param path string. path to output file
#'
#' @return logical scalar
#'
#' @examples `NULL`
#'
#' @export

# phenotype <- tibble::tibble(a=rnorm(5), b = rpois(5, 2))
# path <- 'tests/testthat/output/test.cls'
write_cls <- function(phenotype, path) {
	if (is.null(phenotype)) return(F)
	if (ncol(phenotype) == 0L || nrow(phenotype) == 0L) {warning('empty phenotype'); return(F)}

	na2nan  <- . %>% {.[is.na(.)] = NaN; .}
	Inf2nan <- . %>% {.[is.infinite(.)] = NaN; .}
	phenotype %<>% dplyr::mutate_all(na2nan) %>% dplyr::mutate_all(Inf2nan)

	colnames(phenotype) %>%
		sapply(. %>% {paste0('#', ., '\n', paste0(phenotype[[.]], collapse = ' '), '\n')}) %>%
		c('#numeric\n', .) %>% readr::write_lines(path)

	return(T)
}



#' @title make phenotype tibble for [write_cls]
#'
#' @param matrix tibble. usually the return value of [rGEO::read_gse_matrix()]
#' @param chip tibble. usually the return value of [rGEO::read_gse_soft()]
#' @param gene character. HUGO gene symbol of your interested gene
#'
#' @return tibble. each column represents a _phenotype_, its name is the phenotype's name (hugo gene symbol), its value contains expression value of all samples in the .txt file
#' @export
#' 
#' @section internal note: `max(c(NA,NA), na.rm = TRUE)` produces `-Inf`
#'
#' @section to do: reject too many NA
#' \code{
#' 		if (sum(is.infinite(value)) / length(value) <= 0.25) {}
#'      testthat::test_that('make_gsea_input()', {
#' 	        dataset <- gse51280 %>% {.[[2]][1:5] = 'ANKLE2'; .}
#' 	        gene <- 'ANKLE2'
#' 	        output_file <- 'tests/testthat/output/GSE51280'
#' 	        testthat::expect_true(make_gsea_input(dplyr::mutate_at(dataset, 3:4, . %>% {NA}), gene, output_file));
#' 	       testthat::expect_false(make_gsea_input(dplyr::mutate_at(dataset, 3:10, . %>% {NA}), gene, '/tmp/foobar'));
#'     });
#' }
#'
#' @examples NULL

# matrix <- rGEO::read_gse_matrix(system.file('extdata/GSE51280_series_matrix.txt.gz', package = 'rGEO'))
# chip <- rGEO::read_gse_soft(system.file('extdata/GSE51280_family.soft.gz', package = 'rGEO'))
# gene <- c('CCNB1', 'NEUROG1')
make_phenotype <- function(matrix, chip, gene) {
	maximum <- . %>% {suppressWarnings(max(., na.rm = T))}

	if (is.null(matrix) || is.null(chip)) return(NULL)
	dataset <- dplyr::inner_join(chip, matrix, 'ID_REF') %>% dplyr::select(-1)
	if (nrow(dataset) == 0L) return(NULL)

	dataset %>%
		dplyr::filter(.data$symbol %in% gene) %>%
		dplyr::group_by_at('symbol') %>%
		dplyr::summarise_if(is.numeric, suppressWarnings(dplyr::funs(maximum))) %>%
		plyr::dlply('symbol', . %>% {
			name = .[1, 1, drop = T];
			value = unlist(.[1, -1, drop = T]);
			tibble::tibble(!!name := value)
		}) %>% dplyr::bind_cols()
}






# others --------------



#' @title prepare GSEA input files from GSE raw data
#' 
#' @description You can use first example if you already have a [phenotype labels file](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Phenotype_Labels), or use second example if you are interested in gene sets co-expressing with certain gene(s).
#' 
#' @details both `NA` and `Inf` are not allowed in input files (`.txt` & `.cls`). 
#' 
#' @param matrix_file string. path to GSE series matrix file
#' @param soft_file string. path to GSE SOFT file
#' @param output_dir string, path to output directory
#' @param gene character. HUGO gene symbol(s) of your interested gene(s)
#' 
#' @return `NULL`
#' 
#' @examples
#' \dontrun{
#'   make_gsea_input(
#'     system.file('extdata/GSE51280_series_matrix.txt.gz', package = 'rGEO'),
#'     system.file('extdata/GSE51280_family.soft.gz', package = 'rGEO')
#'   )
#' 
#'   make_gsea_input(
#'     system.file('extdata/GSE19161_series_matrix.txt.gz', package = 'rGEO'),
#'     system.file('extdata/GSE19161_family.soft.gz', package = 'rGEO'),
#'     '.', 'EIF4G2'
#'   )
#' }
#' 
#' @export


# to do:  and testthat
# matrix_file = system.file('extdata/GSE19161_series_matrix.txt.gz', package = 'rGEO');
# soft_file = system.file('extdata/GSE19161_family.soft.gz', package = 'rGEO');
# output_dir = 'tests/testthat/output';
# gene = c('EIF4G2', 'PARK7');
make_gsea_input <- function(matrix_file, soft_file, output_dir = '.', gene = NULL) {
    accession <- stringr::str_extract(matrix_file, 'GSE\\d+')
	output_file <- output_dir %>% stringr::str_remove('/$') %>% paste0('/', accession)

    matrix <- rGEO::read_gse_matrix(matrix_file)
    qGSEA::write_txt(matrix, paste0(output_file, '.txt'))

    chip <- rGEO::read_gse_soft(soft_file)
    qGSEA::write_chip(chip, paste0(output_file, '.chip'))

    if (!is.null(gene)) {
        phenotype <- qGSEA::make_phenotype(matrix, chip, gene)
        qGSEA::write_cls(phenotype, paste0(output_file, '.cls'))
    }
    
    NULL
}















