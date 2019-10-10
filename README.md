# unit-gamma-tests
Computer code for performing improved testing inferences in unit gamma regressions

ARTICLE: Modified likelihood ratio tests for unit gamma regressions

AUTHORS: Ana C. Guedes, Francisco Cribari-Neto, Patrícia L. Espinheira

JOURNAL: Journal of Applied Statistics

FILES: 

README.txt : this file, contains instructions

data_application-1.mat : data 

application-1.ox : Ox code for reproducing some of the resuts in Section 5.1 

application-1.out : output

application-2.ox : Ox code to test a restriction on the mean submodel parameters

application-2.out : output

PROGRAMMING LANGUAGE: Ox - see https://www.doornik.com (OBS. The console version of Ox is freely distributed for academic use. It is currently available for Windows, Linux, OS-X (Apple), and several Unix platforms.) 

EXECUTION: 

[Command line] oxl application-1.ox > application-1.out 
[Command line] oxl application-2.ox > application-2.out 

(On some systems, one must use oxl64 instead of oxl.) 

[IDE, via OxEdit] Open the file application.ox in OxEdit (Ox IDE) and select (top menu): Run -> Ox 

PROGRAM DESCRIPTION application-1.ox: The program can be used to test the null hypothesis of constant precision.

PROGRAM DESCRIPTION application-2.ox: The program can be used to test the null hypothesis \beta_4=0 (in GammaModel--II). This result is not informed in the articlle. We provide this computer code to exemplify a testing inference on the parameters that index the mean submodel. 

