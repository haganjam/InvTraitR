## R CMD check results

* Namespace in Imports field not imported from: ‘doParallel’
* * bdc 1.1.4 depends on it without declaring it; tests fail without having it installed - fix will come with 1.1.5
* Regarding `no visible binding for global varible X` for column names:
* * Should we really just define a local variable with the same name to silence the linter?
That seems odd.

## Downstream dependencies

none yet