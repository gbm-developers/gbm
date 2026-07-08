# pmml (2.6.0)

* GitHub: <https://github.com/mhahsler/r-pmml>
* Email: <mailto:mhahsler@lyle.smu.edu>
* GitHub mirror: <https://github.com/cran/pmml>

Run `revdepcheck::revdep_details(, "pmml")` for more info

## Newly broken

*   checking tests ...
     ```
       Running 'testthat.R'
      ERROR
     Running the tests in 'tests/testthat.R' failed.
     Last 13 lines of output:
         'test_pmml_integration_transformations.R:445:3'
       • skip until export issue is resolved (1): 'test_pmml.nnet.R:66:3'
       • {neighbr} is not installed (6): 'test_pmml.neighbr.R:6:3',
         'test_pmml.neighbr.R:29:3', 'test_pmml.neighbr.R:57:3',
         'test_pmml.neighbr.R:87:3', 'test_pmml.neighbr.R:116:3',
         'test_pmml.neighbr.R:147:3'
       
       ══ Failed tests ════════════════════════════════════════════════════════════════
       ── Failure ('test_pmml.gbm.R:11:3'): pmml.gbm final Segment contains modelName attribute ──
       Expected `... <- NULL` to produce warnings.
       
       [ FAIL 1 | WARN 25 | SKIP 58 | PASS 363 ]
       Error:
       ! Test failures.
       Execution halted
     ```

