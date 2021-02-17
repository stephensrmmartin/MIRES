
## Resubmission

This is a resubmission. In this version, I addressed the notes from the incoming CRAN check:

* Moved all user-facing fn tests to one test_that call, requiring only one fit to be estimated.
* Changed adapt_delta to .8; removed print tests
* Changed cores = 2, chains = 2 to 1, 1, respectively. I believe Windows is misreporting processing time when parallelized.
* Wrapped example in mires() in `dontrun`, because it runs for longer than CRAN permits.
* Added tests via testthat for user-facing functions in lieu of the dontrun example.
