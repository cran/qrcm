qrcm 3.0
=============

Changes with respect to version 2.2
------------------

* bug with survival 2.41 fixed in a better way (timefix = FALSE).
* added new function: iqrL and associated auxiliary functions
* renamed internal functions
* defined a test.fit method. Allow R = 0.
* improved starting points, better estimation of Jacobian and outer product
* handling of potentially discrete data in start.iqr and trans
* fixed warnings with PDF < 0, and wrong definition of covar.ok
* fixed predict (contrasts)
* fixed plf with scalar input
* new dependency pch >= 2.0
* functions for quantile crossing
* added obj.function for censored and truncated data