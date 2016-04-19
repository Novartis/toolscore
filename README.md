# Research code: tool score calculation 
## The compound target tool score
The use of potent and selective chemical tools with well-defined targets can help elucidate biological processes driving phenotypes in phenotypic screens. However, identification of selective compounds en masse to create targeted screening sets is nontrivial. A systematic approach is needed to prioritize probes, which prevents the repeated use of published but unselective compounds. Here we performed a metaanalysis of integrated large-scale, heterogeneous bioactivity data in order to create an evidence-based, quantitative metric to systematically rank tool compounds for targets. Our tool score (TS) was then tested on hundreds of compounds by assessing their activity profiles in a panel of 41 cell-based pathway assays. We demonstrate that high-TS tools show more reliably selective phenotypic profiles than lower-TS compounds. Additionally, we highlight frequently-tested compounds that are non-selective tools and distinguish target-family polypharmacology from cross-family promiscuity. TS can therefore be used to prioritize compounds from heterogeneous databases for phenotypic screening.

## USAGE
```
Rscript example_toolscore_calculation_script.R
```

NOTE: The R script assumes the input data "example_inputdata_toolscore_calculation.csv" is under the same directory.

## CONTACT
yuan-2.wang@novartis.com, jeremy.jenkins@novartis.com
