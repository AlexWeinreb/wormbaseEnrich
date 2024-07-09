

Wrap the Tissue Enrichment Analysis from Wormbase.

The web version is accessible [on Wormbase](http://www.wormbase.org/tools/enrichment/tea/tea.cgi). The Python package is hosted on [PyPI](https://pypi.org/project/tissue_enrichment_analysis/), source code [on Github](https://github.com/dangeles/TissueEnrichmentAnalysis).


# Installation

The R package can be installed with:
```r
devtools::install_github("AlexWeinreb/wormbaseEnrich")
```

The Python package then needs to be installed with:
```r
wormbaseEnrich::install_tea()
```

note that you can select a Python virtual environment beforehand if you wish the package to be installed there.




# Usage

Assuming we want to use a virtualenv named `'r-tea'`, we can first perform the installation:

```r
# reticulate::virtualenv_remove('r-tea')
# reticulate::virtualenv_create('r-tea')
reticulate::use_virtualenv('r-tea')
wormbaseEnrich::install_tea()
```

If the installation was done previously, in a new session we can call
```r
reticulate::use_virtualenv('r-tea')
library(wormbaseEnrich)

tea <- reticulate::import("tissue_enrichment_analysis", delay_load = TRUE, convert = FALSE)
```

We can fetch one or several dictionaries from Wormbase:

```r
my_tissue_dict <- fetch_dictionary("tissue")
```

And them use it:
```r


df <- enrichment_analysis(xx, my_tissue_dict)
```




