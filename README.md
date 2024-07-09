

Reimplements the Tissue Enrichment Analysis from Wormbase.

The web version is accessible [on Wormbase](http://www.wormbase.org/tools/enrichment/tea/tea.cgi). The Python package is hosted on [PyPI](https://pypi.org/project/tissue_enrichment_analysis/), source code [on Github](https://github.com/dangeles/TissueEnrichmentAnalysis). This package is written entirely in R, and gives the same results except for the FDR correction.


# Installation

This R package can be installed with:
```r
devtools::install_github("AlexWeinreb/wormbaseEnrich")
```


# Usage

We first need to download a dictionary.

```r
library(wormbaseEnrich)

dict <- fetch_dictionary("tissue")
```

We can visualize its contents, it has one column for gene IDs, and then one column per term (tissue in this case), containing 0 or 1.

```r
dim(dict)
#> [1] 20845   567

dict[1:4,1:4]
#>             wbid Z3 WBbt:0004575 ALA WBbt:0003955 ABplpaap WBbt:0006077
#> 1 WBGene00022250               0                0                     0
#> 2 WBGene00012193               0                0                     0
#> 3 WBGene00011886               0                0                     0
#> 4 WBGene00001043               0                0                     0

table(dict$`Z3 WBbt:0004575`)
#> 
#>     0     1 
#> 20786    59 
```

We can then use this dictionary to test lists of genes. Note, we have to give wormbase gene IDs (I use [wbData](https://github.com/AlexWeinreb/wbData) for conversions):

```r
gene_list <- c("WBGene00007913","WBGene00006890","WBGene00007802","WBGene00000449",
               "WBGene00009386","WBGene00021950","WBGene00012670","WBGene00004301")

res <- enrichment_analysis(gene_list, dict)

res
#> # A tibble: 10 × 7
#>    term_name            term_id     expected observed enrichment_fc p_value     FDR
#>    <chr>                <chr>          <dbl>    <int>         <dbl>   <dbl>   <dbl>
#>  1 ABplpapp             WBbt:00064…   0.0191        1          52.4 1.73e-4 0.0128 
#>  2 ABprpapaa            WBbt:00064…   0.0196        1          51.0 1.84e-4 0.0128 
#>  3 ABaraapa             WBbt:00058…   0.0196        1          51.0 1.84e-4 0.0128 
#>  4 ABaraapp             WBbt:00061…   0.0202        1          49.6 1.94e-4 0.0128 
#>  5 ABprpappp            WBbt:00058…   0.0202        1          49.6 1.94e-4 0.0128 
#>  6 ABplpapa             WBbt:00060…   0.0207        1          48.3 2.05e-4 0.0128 
#>  7 ABprpappa            WBbt:00060…   0.0207        1          48.3 2.05e-4 0.0128 
#>  8 ABalpapa             WBbt:00065…   0.0218        1          45.9 2.27e-4 0.0128 
#>  9 anterior arcade cell WBbt:00057…   0.385         5          13.0 2.91e-6 0.00165
#> 10 IL socket cell       WBbt:00084…   0.289         3          10.4 2.20e-4 0.0128 
```

Finally we can plot these results.

```r
plot_enrichment_results(res)
```


This is a standard `{ggplot2}` plot, it can be modified the usual way and saved with:

```r
my_gg <- plot_enrichment_results(res)

ggsave("filename.png", my_gg, width = 10, height = 10, units = "cm")
```


In addition, one should provide a background gene list:
```r
enrichment_analysis(gene_list, dict, background_genes = background_genes)
```



# Comparison with Python package

We can run the official Python package using `{reticulate}`. First, we create a virtual environment in which to install the Python packages:

```r
# install.packages("reticulate")
# reticulate::virtualenv_remove('r-tea')
# reticulate::virtualenv_create('r-tea')

reticulate::use_virtualenv('r-tea')
reticulate::py_install("tissue_enrichment_analysis")
```

If the installation was done previously, in a new session we can call
```r
reticulate::use_virtualenv('r-tea')


tea <- reticulate::import("tissue_enrichment_analysis", delay_load = TRUE, convert = FALSE)
```

We can fetch one or several dictionaries from Wormbase:

```r
my_tissue_dict_py <- tea$fetch_dictionary("tissue")
```

Note this is a Python object, we can visualize it, not directly manipulate it:
```r
class(my_tissue_dict_py)
#> [1] "pandas.core.frame.DataFrame"        "pandas.core.generic.NDFrame"       
#> [3] "pandas.core.base.PandasObject"      "pandas.core.accessor.DirNamesMixin"
#> [5] "pandas.core.indexing.IndexingMixin" "pandas.core.arraylike.OpsMixin"    
#> [7] "python.builtin.object"  

my_tissue_dict_py
#>                  wbid  ...  ABplappap WBbt:0006067
#> 0      WBGene00022250  ...                     0.0
#> 1      WBGene00012193  ...                     0.0
#> 2      WBGene00011886  ...                     0.0
#> 3      WBGene00001043  ...                     0.0
#> 4      WBGene00050906  ...                     0.0
#> ...               ...  ...                     ...
#> 20840  WBGene00014249  ...                     0.0
#> 20841  WBGene00006244  ...                     0.0
#> 20842  WBGene00011054  ...                     0.0
#> 20843  WBGene00017264  ...                     0.0
#> 20844  WBGene00004744  ...                     0.0
#> 
#> [20845 rows x 567 columns]
```


We can then perform the enrichment test of a gene list against the dictionary:
```r

gene_list <- c("WBGene00007913","WBGene00006890","WBGene00007802","WBGene00000449","WBGene00009386","WBGene00021950","WBGene00012670","WBGene00004301")

res_py <- tea$enrichment_analysis(
    gene_list = gene_list,
    tissue_df = my_tissue_dict_py
  )

res_py
#>                                   Term  Expected  ...   P value   Q value
#> 39   anterior arcade cell WBbt:0005794  0.384664  ...  0.000003  0.001646
#> 3          IL socket cell WBbt:0008418  0.288770  ...  0.000220  0.049068
#> 63               ABplpapp WBbt:0006420  0.019070  ...  0.000173  0.049068
#> 18              ABprpapaa WBbt:0006446  0.019615  ...  0.000184  0.049068
#> 72              ABprpappp WBbt:0005847  0.020159  ...  0.000194  0.049068
#> 73               ABalpapa WBbt:0006573  0.021794  ...  0.000227  0.049068
#> 79               ABplpapa WBbt:0006087  0.020704  ...  0.000205  0.049068
#> 71               ABaraapp WBbt:0006153  0.020159  ...  0.000194  0.049068
#> 101             ABprpappa WBbt:0006088  0.020704  ...  0.000205  0.049068
#> 104              ABaraapa WBbt:0005853  0.019615  ...  0.000184  0.049068
#> 
#> [10 rows x 6 columns]
```

Note how the Q value differs but other values identical.


