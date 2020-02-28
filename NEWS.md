qrcm 2.2
=============

Changes with respect to version 2.1
------------------
* bug fixed. Replaced class(obj) == "try-error" with inherits(obj, "try-error")
* Replaced pch:::fun with fun <- getFromNamespace(fun, ns = "pch"), added getFromNamespace to imports
* require pch >= 1.4
* updated my e-mail address
* updated reference to Frumento and Bottai (2017) with issue and page numbers