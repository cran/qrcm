qrcm 2.0
=============

Important changes with respect to version 1.0
------------------
* the main function **iqr** has been expanded to allow censored and truncated data via **Surv()**.
* the function **test.fit** has been expanded accordingly.
* the chi-squared goodness-of-fit test provided by **summary.iqr** appeared unreliable and was removed.
* **plot.iqr** has been improved

Changes in the package structure
----------------
* new internal functions have been created to handle truncation and censoring.
* the algorithm has been improved to be faster and more stable.
* the package depends on **pch (>= 1.2)**.