context('knit_print working properly')

if(!exists('cgr_exp')) skip('Corrgrapher did not create the object')

test_that('knit_print works with numerical data',{
         expect_s3_class(knit_print.corrgrapher(cgr_exp), 'knit_asis')
         expect_s3_class(knit_print.corrgrapher(cgr_df), 'knit_asis')
})

test_that('knit_print works with mixed data',{
          expect_s3_class(knit_print.corrgrapher(tit_cgr_exp), 'knit_asis')
          expect_s3_class(knit_print.corrgrapher(cgr_df_mixed), 'knit_asis')
})
