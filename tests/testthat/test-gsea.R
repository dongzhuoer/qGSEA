testthat::context('Testing gsea.R')
setwd(here::here(''))  # workspace is reset per file



# .txt --------
testthat::test_that('write_txt()', {
	testthat::expect_false(write_txt(NULL, 'doubi'))

	# reject dataset with only one sample
	testthat::expect_warning(
		tibble::tibble(ID_REF = character(), GSM1 = numeric()) %>% write_txt('doubi')
	) %>% testthat::expect_false()

	testthat::expect_true(
	    system.file('extdata/GSE51280_series_matrix.txt.gz', package = 'rGEO') %>% 
	        rGEO::read_gse_matrix() %>% write_txt('tests/testthat/output/GSE51280.txt')
	)

	testthat::expect_true(
		system.file('extdata/GSE19161_series_matrix.txt.gz', package = 'rGEO') %>% 
	        rGEO::read_gse_matrix() %>% write_txt('tests/testthat/output/GSE19161.txt')
	)
});



testthat::test_that('read_txt()', {
	testthat::expect_identical(read_txt('non/existent/file'), NULL)

	testthat::expect_equal(
	    system.file('extdata/GSE51280_series_matrix.txt.gz', package = 'rGEO') %>% 
	        rGEO::read_gse_matrix(),
	    read_txt('tests/testthat/output/GSE51280.txt')
	)

	testthat::expect_equal(
	    system.file('extdata/GSE19161_series_matrix.txt.gz', package = 'rGEO') %>% 
	        rGEO::read_gse_matrix(),
	    read_txt('tests/testthat/output/GSE19161.txt')
	)
});






# .chip --------
testthat::test_that('write_chip()', {
	testthat::expect_false(write_chip(NULL, 'doubi'))

	testthat::expect_warning(
		tibble::tibble(ID_REF = character(), symbol = character()) %>% write_chip('doubi')
	) %>% testthat::expect_false()

	testthat::expect_true(
	    system.file('extdata/GSE51280_family.soft.gz', package = 'rGEO') %>% 
	        rGEO::read_gse_soft() %>% write_chip('tests/testthat/output/GSE51280.chip')
	)

	testthat::expect_true(
		system.file('extdata/GSE19161_family.soft.gz', package = 'rGEO') %>% 
	        rGEO::read_gse_soft() %>% write_chip('tests/testthat/output/GSE19161.chip')
	)
});



testthat::test_that('read_chip()', {
	testthat::expect_identical(read_chip('non/existent/file'), NULL)

	testthat::expect_equal(
	    system.file('extdata/GSE51280_family.soft.gz', package = 'rGEO') %>% 
	        rGEO::read_gse_soft(),
	    read_chip('tests/testthat/output/GSE51280.chip')
	)

	testthat::expect_equal(
	    system.file('extdata/GSE19161_family.soft.gz', package = 'rGEO') %>% 
	        rGEO::read_gse_soft(),
	    read_chip('tests/testthat/output/GSE19161.chip')
	)

});






# .cls --------
testthat::test_that('write_cls()', {
	testthat::expect_false(write_cls(NULL, 'doubi'))
	testthat::expect_warning(write_cls(tibble::tibble()))
	testthat::expect_warning(write_cls(tibble::tibble(x = numeric())))

 	testthat::expect_true(
		tibble::tibble(a = c(1, NA, 2, 3), b = c(1, 2, Inf, 4)) %>%
			write_cls('tests/testthat/output/test.cls')
	)
 	testthat::expect_identical(
		readr::read_file('tests/testthat/output/test.cls'),
		'#numeric\n\n#a\n1 NaN 2 3\n\n#b\n1 2 NaN 4\n\n'
	)
});



testthat::test_that('make_phenotype()', {
	testthat::expect_identical(make_phenotype(NULL, tibble::tibble(), 'doubi'), NULL)
	testthat::expect_identical(make_phenotype(tibble::tibble(), NULL, 'doubi'), NULL)

	testthat::expect_identical(
		make_phenotype(tibble::tibble(ID_REF = 1:2), tibble::tibble(ID_REF = 3:4), 'doubi'),
		NULL
	)

	testthat::expect_identical(
		make_phenotype(tibble::tibble(ID_REF = 1:4, c(4, 7, 5, 4), c(8, 5, 11, 6)),
					   tibble::tibble(ID_REF = 1:4, symbol = c('a', 'b', 'a', 'b')),
					   c('a', 'b')
		),
		tibble::tibble(a = c(5, 11), b = c(7, 6))
	)

	matrix <- read_txt('tests/testthat/output/GSE51280.txt')
	chip <- read_chip('tests/testthat/output/GSE51280.chip')
	gene <- c('CCNB1', 'NEUROG1')
	phenotype <- make_phenotype(matrix, chip, gene)

	testthat::expect_true(tibble::is_tibble(phenotype))
	testthat::expect_identical(colnames(phenotype), gene)
 	testthat::expect_identical(dim(phenotype), c(24L, 2L))
 	testthat::expect_identical(sapply(phenotype, class) %>% {names(.) <- NULL; .}, rep('numeric', 2))

 	testthat::expect_warning(make_phenotype(matrix, chip, 'doubi') %>% write_cls(tempfile()))
});




# GSEA ------


testthat::test_that('make_gsea_input()', {
    testthat::expect_identical(
        make_gsea_input(
            system.file('extdata/GSE51280_series_matrix.txt.gz', package = 'rGEO'),
            system.file('extdata/GSE51280_family.soft.gz', package = 'rGEO'),
            'tests/testthat/output/'
        ),
        NULL
    )
    
    testthat::expect_identical(
        make_gsea_input(
            system.file('extdata/GSE19161_series_matrix.txt.gz', package = 'rGEO'),
            system.file('extdata/GSE19161_family.soft.gz', package = 'rGEO'),
            'tests/testthat/output', 
            c('EIF4G2', 'PARK7')
        ),
        NULL
    )
})


testthat::test_that('run GSEA using real data GSE19161', {
	testthat::skip('skip running GSEA')

    gsea_command <- . %>% paste0('cd tests/testthat/output; java -Xmx512m -cp /media/computer/opt/GSEA/gsea-3.0.jar xtools.gsea.Gsea -gmx /media/computer/opt/GSEA/msigdb_v6.1/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt -res ', ., '.txt -cls ', ., '.cls#ANKLE2 -chip ', ., '.chip -out . -collapse true -nperm 10 -permute phenotype -rpt_label ', ., ' -metric Pearson -gui false')

	'GSE19161' %>% gsea_command %>% stringr::str_replace('ANKLE2', 'EIF4G2')  %>% system
});


# GSE51280: none of the gene sets passed size thresholds










