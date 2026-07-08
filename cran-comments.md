# gbm 2.3.0 (resubmission)

## Resubmission note

The previous submission of gbm 2.3.0 failed the incoming pretest on
r-devel-linux-arm64 with two test failures in `test-least-squares.R`
("gaussian predictions recover the signal"), while passing on Windows and
Debian x86_64.

We traced the failure to a genuine, long-standing bug in the package's C++
code that the arm64 check environment exposed. The categorical
split-direction codes that `predict.gbm()` uses were stored in
`std::vector<char>` holding the values -1 (left branch) and +1 (right
branch). Plain `char` is signed on x86_64 and on Apple's arm64 ABI, but
unsigned on the Linux aarch64 ABI, so on linux-arm64 the stored -1 became
255. As a result, `predict()` routed every left-branch observation of every
categorical split to the missing-value branch, producing systematically
biased predictions. Model fitting itself was unaffected (the training code
uses a different, correctly typed representation), which is why only the
prediction-accuracy test failed and only on that platform.

The fix is to use `std::vector<signed char>` (src/node.h, src/gbm.h). We
verified on a native Linux aarch64 runner (Ubuntu 24.04, gcc) that the
previously failing checks now pass, and that `predict()` now agrees with the
internal training fit to within floating-point tolerance on all platforms. A
regression test guarding the serialized-prediction path has been added.

## Reverse dependencies

`revdepcheck` flags one regression in `pmml`: `test_pmml.gbm.R` asserts that
`gbm(..., distribution = "multinomial")` emits a warning
("...ill-advised... currently broken... use at your own risk."). This
submission intentionally removes that warning, since multinomial support has
been fixed and is now tested (see NEWS.md "Behavior changes"). This warning
is present in gbm 2.2.3, the version currently on CRAN, so `pmml`'s test
passes today and will only start failing once gbm 2.3.0 is published. We have
notified the pmml maintainer (cc'd) ahead of this submission so the test's
`expect_warning()` wrapper can be removed to match gbm's new behavior.

## Test environments

* Windows (win-builder, R-devel)
* Debian Linux x86_64 (R-devel)
* Ubuntu 24.04 Linux aarch64 (GitHub Actions, R release)
* macOS arm64 (GitHub Actions, R release)

## R CMD check results

0 errors | 0 warnings | 0 notes
