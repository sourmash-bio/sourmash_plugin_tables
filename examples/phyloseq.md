
# Sourmash phyloseq object creation

## Gather

The `sourmash gather` command was used to create 5 gather csv files
found in `examples/data/gather`. The parameters used to generate these
samples can be found in the file names.

The `podar-ref.zip` database is a small toy database for easy testing.

``` bash
for sig in (awk -F, 'NR>1 {print $21}' examples/data/tiny_example.csv); do
    sourmash gather examples/data/sigs/${sig}.sig examples/data/podar-ref.zip -k 31 --scaled 1000 --dna -o examples/data/gather/${sig}.k31-s1000-t3000.csv --threshold-bp 3000
done
```

## Tax

The taxonomy file is located in `examples/data/` and is
`podar-red.tax.csv`.

## Metadata

The metadata for the 5 example signatures used in this vingette is from
the [CuratedMetagenomicData repository on
GitHub](https://github.com/waldronlab/curatedMetagenomicData)

An example metadata file was generated with:

``` bash
head -1 data/big_metadata.csv > tiny_example.csv
tail -n +2 data/big_metadata.csv | shuf | head -n 5 >> tiny_example.csv
```

## Sourmash gather files to table

The phyloseq object expects an n by m matrix. To achieve this across a
large number of samples, I wrote the `tables` plugin for sourmash.

This plugin collects all the gather information and creates an OTU
object from our sourmash gather data.

An example of the tables plugin operation:

``` bash
sourmash scripts gather_tables examples/data/gather/*  --format dense -o examples/data/dense-gather.k31-s1000-t3000.csv --column f_unique_weighted --filter 0
```

| match_name                                                                                       | ERR1136663 | ERR4087400 | ERR4562131 | SRR059888 | SRR9217435 |
|:-------------------------------------------------------------------------------------------------|-----------:|-----------:|-----------:|----------:|-----------:|
| CP000139.1 Bacteroides vulgatus ATCC 8482, complete genome                                       |  0.0083015 |  0.0055121 |  0.0058230 | 0.0000336 |  0.0000579 |
| AE015928.1 Bacteroides thetaiotaomicron VPI-5482, complete genome                                |  0.0020327 |  0.0003969 |  0.0004230 | 0.0000000 |  0.0000000 |
| NZ_DS996397.1 Desulfovibrio piger ATCC 29098 Scfld442, whole genome shotgun sequence             |  0.0010764 |  0.0000000 |  0.0001549 | 0.0000000 |  0.0000000 |
| CP001071.1 Akkermansia muciniphila ATCC BAA-835, complete genome                                 |  0.0005676 |  0.0041199 |  0.0000430 | 0.0000098 |  0.0000000 |
| NZ_KE136524.1 Enterococcus faecalis V583 acyDH-supercont2.1, whole genome shotgun sequence       |  0.0000236 |  0.0000693 |  0.0000000 | 0.0000000 |  0.0000000 |
| CP001472.1 Acidobacterium capsulatum ATCC 51196, complete genome                                 |  0.0001166 |  0.0000000 |  0.0000000 | 0.0000000 |  0.0000000 |
| AE009951.2 Fusobacterium nucleatum subsp. nucleatum ATCC 25586, complete genome                  |  0.0002815 |  0.0000000 |  0.0000000 | 0.0032173 |  0.0000565 |
| CP001013.1 Leptothrix cholodnii SP-6, complete genome                                            |  0.0000094 |  0.0000000 |  0.0000000 | 0.0001588 |  0.0001924 |
| NZ_JGWU01000001.1 Bordetella bronchiseptica D989 ctg7180000008197, whole genome shotgun sequence |  0.0000000 |  0.0000378 |  0.0000000 | 0.0002794 |  0.0003983 |
| NZ_KQ961402.1 Zymomonas mobilis strain ATCC 31823 Scaffold1, whole genome shotgun sequence       |  0.0000000 |  0.0000000 |  0.0003728 | 0.0000000 |  0.0000293 |
| AP009380.1 Porphyromonas gingivalis ATCC 33277 DNA, complete genome                              |  0.0000000 |  0.0000000 |  0.0000574 | 0.0002765 |  0.0000086 |
| NC_007951.1 Burkholderia xenovorans LB400 chromosome 1, complete sequence                        |  0.0000000 |  0.0000000 |  0.0000258 | 0.0000000 |  0.0000043 |
| AE017226.1 Treponema denticola ATCC 35405, complete genome                                       |  0.0000000 |  0.0000000 |  0.0000000 | 0.0000726 |  0.0000000 |
| CP000850.1 Salinispora arenicola CNS-205, complete genome                                        |  0.0000000 |  0.0000000 |  0.0000000 | 0.0003465 |  0.0000000 |
| AE006470.1 Chlorobium tepidum TLS, complete genome                                               |  0.0000000 |  0.0000000 |  0.0000000 | 0.0000754 |  0.0000000 |
| AE017221.1 Thermus thermophilus HB27, complete genome                                            |  0.0000000 |  0.0000000 |  0.0000000 | 0.0000034 |  0.0000000 |
| NC_011663.1 Shewanella baltica OS223, complete genome                                            |  0.0000000 |  0.0000000 |  0.0000000 | 0.0000000 |  0.0001623 |

## R script to create sourmash phyloseq object

### Load necessary libraries

``` r
library(phyloseq)
library(stringr)
```

### Load sourmash OTU and tax data

``` r
df <- read.csv('data/dense-gather.k31-s1000-t3000.csv')
tax <- read.csv('data/podar-ref.tax.csv')
```

### Manipulate the taxonomy strings in both dataframes to have matching row names

``` r
str(df)
```

    ## 'data.frame':    17 obs. of  6 variables:
    ##  $ match_name: chr  "CP000139.1 Bacteroides vulgatus ATCC 8482, complete genome" "AE015928.1 Bacteroides thetaiotaomicron VPI-5482, complete genome" "NZ_DS996397.1 Desulfovibrio piger ATCC 29098 Scfld442, whole genome shotgun sequence" "CP001071.1 Akkermansia muciniphila ATCC BAA-835, complete genome" ...
    ##  $ ERR1136663: num  8.30e-03 2.03e-03 1.08e-03 5.68e-04 2.36e-05 ...
    ##  $ ERR4087400: num  5.51e-03 3.97e-04 0.00 4.12e-03 6.93e-05 ...
    ##  $ ERR4562131: num  0.005823 0.000423 0.000155 0.000043 0 ...
    ##  $ SRR059888 : num  3.36e-05 0.00 0.00 9.81e-06 0.00 ...
    ##  $ SRR9217435: num  5.79e-05 0.00 0.00 0.00 0.00 ...

``` r
df[,1]
```

    ##  [1] "CP000139.1 Bacteroides vulgatus ATCC 8482, complete genome"                                      
    ##  [2] "AE015928.1 Bacteroides thetaiotaomicron VPI-5482, complete genome"                               
    ##  [3] "NZ_DS996397.1 Desulfovibrio piger ATCC 29098 Scfld442, whole genome shotgun sequence"            
    ##  [4] "CP001071.1 Akkermansia muciniphila ATCC BAA-835, complete genome"                                
    ##  [5] "NZ_KE136524.1 Enterococcus faecalis V583 acyDH-supercont2.1, whole genome shotgun sequence"      
    ##  [6] "CP001472.1 Acidobacterium capsulatum ATCC 51196, complete genome"                                
    ##  [7] "AE009951.2 Fusobacterium nucleatum subsp. nucleatum ATCC 25586, complete genome"                 
    ##  [8] "CP001013.1 Leptothrix cholodnii SP-6, complete genome"                                           
    ##  [9] "NZ_JGWU01000001.1 Bordetella bronchiseptica D989 ctg7180000008197, whole genome shotgun sequence"
    ## [10] "NZ_KQ961402.1 Zymomonas mobilis strain ATCC 31823 Scaffold1, whole genome shotgun sequence"      
    ## [11] "AP009380.1 Porphyromonas gingivalis ATCC 33277 DNA, complete genome"                             
    ## [12] "NC_007951.1 Burkholderia xenovorans LB400 chromosome 1, complete sequence"                       
    ## [13] "AE017226.1 Treponema denticola ATCC 35405, complete genome"                                      
    ## [14] "CP000850.1 Salinispora arenicola CNS-205, complete genome"                                       
    ## [15] "AE006470.1 Chlorobium tepidum TLS, complete genome"                                              
    ## [16] "AE017221.1 Thermus thermophilus HB27, complete genome"                                           
    ## [17] "NC_011663.1 Shewanella baltica OS223, complete genome"

``` r
str_view(head(df[,1]), "\\.\\d+")
```

    ## [1] │ CP000139<.1> Bacteroides vulgatus ATCC 8482, complete genome
    ## [2] │ AE015928<.1> Bacteroides thetaiotaomicron VPI-5482, complete genome
    ## [3] │ NZ_DS996397<.1> Desulfovibrio piger ATCC 29098 Scfld442, whole genome shotgun sequence
    ## [4] │ CP001071<.1> Akkermansia muciniphila ATCC BAA-835, complete genome
    ## [5] │ NZ_KE136524<.1> Enterococcus faecalis V583 acyDH-supercont2<.1>, whole genome shotgun sequence
    ## [6] │ CP001472<.1> Acidobacterium capsulatum ATCC 51196, complete genome

``` r
word(df[, 1], 1, 3)
```

    ##  [1] "CP000139.1 Bacteroides vulgatus"            
    ##  [2] "AE015928.1 Bacteroides thetaiotaomicron"    
    ##  [3] "NZ_DS996397.1 Desulfovibrio piger"          
    ##  [4] "CP001071.1 Akkermansia muciniphila"         
    ##  [5] "NZ_KE136524.1 Enterococcus faecalis"        
    ##  [6] "CP001472.1 Acidobacterium capsulatum"       
    ##  [7] "AE009951.2 Fusobacterium nucleatum"         
    ##  [8] "CP001013.1 Leptothrix cholodnii"            
    ##  [9] "NZ_JGWU01000001.1 Bordetella bronchiseptica"
    ## [10] "NZ_KQ961402.1 Zymomonas mobilis"            
    ## [11] "AP009380.1 Porphyromonas gingivalis"        
    ## [12] "NC_007951.1 Burkholderia xenovorans"        
    ## [13] "AE017226.1 Treponema denticola"             
    ## [14] "CP000850.1 Salinispora arenicola"           
    ## [15] "AE006470.1 Chlorobium tepidum"              
    ## [16] "AE017221.1 Thermus thermophilus"            
    ## [17] "NC_011663.1 Shewanella baltica"

``` r
df[,1] <- word(df[,1], 1, 3)
df[,1]
```

    ##  [1] "CP000139.1 Bacteroides vulgatus"            
    ##  [2] "AE015928.1 Bacteroides thetaiotaomicron"    
    ##  [3] "NZ_DS996397.1 Desulfovibrio piger"          
    ##  [4] "CP001071.1 Akkermansia muciniphila"         
    ##  [5] "NZ_KE136524.1 Enterococcus faecalis"        
    ##  [6] "CP001472.1 Acidobacterium capsulatum"       
    ##  [7] "AE009951.2 Fusobacterium nucleatum"         
    ##  [8] "CP001013.1 Leptothrix cholodnii"            
    ##  [9] "NZ_JGWU01000001.1 Bordetella bronchiseptica"
    ## [10] "NZ_KQ961402.1 Zymomonas mobilis"            
    ## [11] "AP009380.1 Porphyromonas gingivalis"        
    ## [12] "NC_007951.1 Burkholderia xenovorans"        
    ## [13] "AE017226.1 Treponema denticola"             
    ## [14] "CP000850.1 Salinispora arenicola"           
    ## [15] "AE006470.1 Chlorobium tepidum"              
    ## [16] "AE017221.1 Thermus thermophilus"            
    ## [17] "NC_011663.1 Shewanella baltica"

``` r
df[,1] <- str_remove_all(df[,1], "\\.\\d+")
df[,1]
```

    ##  [1] "CP000139 Bacteroides vulgatus"            
    ##  [2] "AE015928 Bacteroides thetaiotaomicron"    
    ##  [3] "NZ_DS996397 Desulfovibrio piger"          
    ##  [4] "CP001071 Akkermansia muciniphila"         
    ##  [5] "NZ_KE136524 Enterococcus faecalis"        
    ##  [6] "CP001472 Acidobacterium capsulatum"       
    ##  [7] "AE009951 Fusobacterium nucleatum"         
    ##  [8] "CP001013 Leptothrix cholodnii"            
    ##  [9] "NZ_JGWU01000001 Bordetella bronchiseptica"
    ## [10] "NZ_KQ961402 Zymomonas mobilis"            
    ## [11] "AP009380 Porphyromonas gingivalis"        
    ## [12] "NC_007951 Burkholderia xenovorans"        
    ## [13] "AE017226 Treponema denticola"             
    ## [14] "CP000850 Salinispora arenicola"           
    ## [15] "AE006470 Chlorobium tepidum"              
    ## [16] "AE017221 Thermus thermophilus"            
    ## [17] "NC_011663 Shewanella baltica"

``` r
df
```

    ##                                   match_name   ERR1136663   ERR4087400
    ## 1              CP000139 Bacteroides vulgatus 8.301477e-03 5.512089e-03
    ## 2      AE015928 Bacteroides thetaiotaomicron 2.032678e-03 3.968704e-04
    ## 3            NZ_DS996397 Desulfovibrio piger 1.076401e-03 0.000000e+00
    ## 4           CP001071 Akkermansia muciniphila 5.676425e-04 4.119893e-03
    ## 5          NZ_KE136524 Enterococcus faecalis 2.355363e-05 6.929483e-05
    ## 6         CP001472 Acidobacterium capsulatum 1.165905e-04 0.000000e+00
    ## 7           AE009951 Fusobacterium nucleatum 2.814659e-04 0.000000e+00
    ## 8              CP001013 Leptothrix cholodnii 9.421452e-06 0.000000e+00
    ## 9  NZ_JGWU01000001 Bordetella bronchiseptica 0.000000e+00 3.779718e-05
    ## 10             NZ_KQ961402 Zymomonas mobilis 0.000000e+00 0.000000e+00
    ## 11         AP009380 Porphyromonas gingivalis 0.000000e+00 0.000000e+00
    ## 12         NC_007951 Burkholderia xenovorans 0.000000e+00 0.000000e+00
    ## 13              AE017226 Treponema denticola 0.000000e+00 0.000000e+00
    ## 14            CP000850 Salinispora arenicola 0.000000e+00 0.000000e+00
    ## 15               AE006470 Chlorobium tepidum 0.000000e+00 0.000000e+00
    ## 16             AE017221 Thermus thermophilus 0.000000e+00 0.000000e+00
    ## 17              NC_011663 Shewanella baltica 0.000000e+00 0.000000e+00
    ##      ERR4562131    SRR059888   SRR9217435
    ## 1  5.823024e-03 3.356448e-05 5.792533e-05
    ## 2  4.229973e-04 0.000000e+00 0.000000e+00
    ## 3  1.548600e-04 0.000000e+00 0.000000e+00
    ## 4  4.301667e-05 9.811156e-06 0.000000e+00
    ## 5  0.000000e+00 0.000000e+00 0.000000e+00
    ## 6  0.000000e+00 0.000000e+00 0.000000e+00
    ## 7  0.000000e+00 3.217285e-03 5.649507e-05
    ## 8  0.000000e+00 1.587858e-04 1.923693e-04
    ## 9  0.000000e+00 2.793598e-04 3.983260e-04
    ## 10 3.728112e-04 0.000000e+00 2.932023e-05
    ## 11 5.735556e-05 2.765197e-04 8.581530e-06
    ## 12 2.581000e-05 0.000000e+00 4.290765e-06
    ## 13 0.000000e+00 7.255092e-05 0.000000e+00
    ## 14 0.000000e+00 3.464887e-04 0.000000e+00
    ## 15 0.000000e+00 7.539099e-05 0.000000e+00
    ## 16 0.000000e+00 3.356448e-06 0.000000e+00
    ## 17 0.000000e+00 0.000000e+00 1.623339e-04

``` r
rownames(df) <- df[,1]
rownames(df)
```

    ##  [1] "CP000139 Bacteroides vulgatus"            
    ##  [2] "AE015928 Bacteroides thetaiotaomicron"    
    ##  [3] "NZ_DS996397 Desulfovibrio piger"          
    ##  [4] "CP001071 Akkermansia muciniphila"         
    ##  [5] "NZ_KE136524 Enterococcus faecalis"        
    ##  [6] "CP001472 Acidobacterium capsulatum"       
    ##  [7] "AE009951 Fusobacterium nucleatum"         
    ##  [8] "CP001013 Leptothrix cholodnii"            
    ##  [9] "NZ_JGWU01000001 Bordetella bronchiseptica"
    ## [10] "NZ_KQ961402 Zymomonas mobilis"            
    ## [11] "AP009380 Porphyromonas gingivalis"        
    ## [12] "NC_007951 Burkholderia xenovorans"        
    ## [13] "AE017226 Treponema denticola"             
    ## [14] "CP000850 Salinispora arenicola"           
    ## [15] "AE006470 Chlorobium tepidum"              
    ## [16] "AE017221 Thermus thermophilus"            
    ## [17] "NC_011663 Shewanella baltica"

``` r
df <- df[,-1]

DT::datatable(df)
```

<div class="datatables html-widget html-fill-item" id="htmlwidget-4853653c9b64d3a4362a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-4853653c9b64d3a4362a">{"x":{"filter":"none","vertical":false,"data":[["CP000139 Bacteroides vulgatus","AE015928 Bacteroides thetaiotaomicron","NZ_DS996397 Desulfovibrio piger","CP001071 Akkermansia muciniphila","NZ_KE136524 Enterococcus faecalis","CP001472 Acidobacterium capsulatum","AE009951 Fusobacterium nucleatum","CP001013 Leptothrix cholodnii","NZ_JGWU01000001 Bordetella bronchiseptica","NZ_KQ961402 Zymomonas mobilis","AP009380 Porphyromonas gingivalis","NC_007951 Burkholderia xenovorans","AE017226 Treponema denticola","CP000850 Salinispora arenicola","AE006470 Chlorobium tepidum","AE017221 Thermus thermophilus","NC_011663 Shewanella baltica"],[0.008301477048164819,0.00203267830687083,0.001076400911054425,0.0005676424935757473,2.355363043882769e-05,0.0001165904706721971,0.0002814658837439909,9.421452175531076e-06,0,0,0,0,0,0,0,0,0],[0.005512088798175656,0.0003968703934686472,0,0.004119892656007862,6.929483060563682e-05,0,0,0,3.779718033034736e-05,0,0,0,0,0,0,0,0],[0.005823023670641407,0.0004229972870818062,0.0001548600237452036,4.301667326255657e-05,0,0,0,0,0,0.0003728111682754902,5.735556435007542e-05,2.581000395753394e-05,0,0,0,0,0],[3.356448072391872e-05,0,0,9.811155903914703e-06,0,0,0.00321728457154424,0.0001587858126554616,0.0002793597549483081,0,0.0002765196835024381,0,7.255091602631662e-05,0.0003464887163961456,7.539098747218667e-05,3.356448072391872e-06,0],[5.792532924471091e-05,0,0,0,0,0,5.649507420163164e-05,0.0001923693032941634,0.00039832602949758,2.932022838312528e-05,8.581530258475691e-06,4.290765129237846e-06,0,0,0,0,0.0001623339473894985]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>ERR1136663<\/th>\n      <th>ERR4087400<\/th>\n      <th>ERR4562131<\/th>\n      <th>SRR059888<\/th>\n      <th>SRR9217435<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"ERR1136663","targets":1},{"name":"ERR4087400","targets":2},{"name":"ERR4562131","targets":3},{"name":"SRR059888","targets":4},{"name":"SRR9217435","targets":5}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

``` r
str(tax)
```

    ## 'data.frame':    64 obs. of  10 variables:
    ##  $ accession   : chr  "AE000782" "NC_000909" "NC_003272" "AE009441" ...
    ##  $ taxid       : int  224325 243232 103690 178306 186497 190304 188937 190192 246200 194439 ...
    ##  $ superkingdom: chr  "Archaea" "Archaea" "Bacteria" "Archaea" ...
    ##  $ phylum      : chr  "Euryarchaeota" "Euryarchaeota" "Cyanobacteria" "Crenarchaeota" ...
    ##  $ class       : chr  "Archaeoglobi" "Methanococci" "" "Thermoprotei" ...
    ##  $ order       : chr  "Archaeoglobales" "Methanococcales" "Nostocales" "Thermoproteales" ...
    ##  $ family      : chr  "Archaeoglobaceae" "Methanocaldococcaceae" "Nostocaceae" "Thermoproteaceae" ...
    ##  $ genus       : chr  "Archaeoglobus" "Methanocaldococcus" "Nostoc" "Pyrobaculum" ...
    ##  $ species     : chr  "Archaeoglobus fulgidus" "Methanocaldococcus jannaschii" "Nostoc sp. PCC 7120" "Pyrobaculum aerophilum" ...
    ##  $ strain      : chr  "Archaeoglobus fulgidus DSM 4304" "Methanocaldococcus jannaschii DSM 2661" "" "Pyrobaculum aerophilum str. IM2" ...

``` r
tax[,1]
```

    ##  [1] "AE000782"        "NC_000909"       "NC_003272"       "AE009441"       
    ##  [5] "AE009950"        "AE009951"        "AE010299"        "AE009439"       
    ##  [9] "NC_003911"       "AE006470"        "AE015928"        "AL954747"       
    ## [13] "BX119912"        "BX571656"        "AE017180"        "AE017226"       
    ## [17] "BX950229"        "AE017221"        "BA000001"        "BA000023"       
    ## [21] "NC_007951"       "CP000492"        "NC_008751"       "CP000568"       
    ## [25] "CP000561"        "CP000609"        "CP000607"        "CP000660"       
    ## [29] "CP000667"        "CP000679"        "CP000702"        "CP000139"       
    ## [33] "NC_009665"       "CP000816"        "CP000850"        "CP000909"       
    ## [37] "CP000924"        "CP000969"        "CP001013"        "CP001071"       
    ## [41] "AP009380"        "NC_010730"       "CP001097"        "CP001110"       
    ## [45] "CP001130"        "NZ_CH959311"     "NZ_CH959317"     "CP001251"       
    ## [49] "NC_011663"       "CP000916"        "NZ_DS996397"     "CP001230"       
    ## [53] "CP001472"        "AP009153"        "CP001941"        "NC_013968"      
    ## [57] "NZ_KE136524"     "NZ_KQ961402"     "NZ_CP015081"     "NZ_ABZS01000228"
    ## [61] "NZ_JGWU01000001" "NZ_FWDH01000003" "NC_009972"       "NC_005213"

``` r
tax[,9]
```

    ##  [1] "Archaeoglobus fulgidus"              
    ##  [2] "Methanocaldococcus jannaschii"       
    ##  [3] "Nostoc sp. PCC 7120"                 
    ##  [4] "Pyrobaculum aerophilum"              
    ##  [5] "Pyrococcus furiosus"                 
    ##  [6] "Fusobacterium nucleatum"             
    ##  [7] "Methanosarcina acetivorans"          
    ##  [8] "Methanopyrus kandleri"               
    ##  [9] "Ruegeria pomeroyi"                   
    ## [10] "Chlorobaculum tepidum"               
    ## [11] "Bacteroides thetaiotaomicron"        
    ## [12] "Nitrosomonas europaea"               
    ## [13] "Rhodopirellula baltica"              
    ## [14] "Wolinella succinogenes"              
    ## [15] "Geobacter sulfurreducens"            
    ## [16] "Treponema denticola"                 
    ## [17] "Methanococcus maripaludis"           
    ## [18] "Thermus thermophilus"                
    ## [19] "Pyrococcus horikoshii"               
    ## [20] "Sulfolobus tokodaii"                 
    ## [21] "Paraburkholderia xenovorans"         
    ## [22] "Chlorobium phaeobacteroides"         
    ## [23] "Desulfovibrio vulgaris"              
    ## [24] "Ruminiclostridium thermocellum"      
    ## [25] "Pyrobaculum calidifontis"            
    ## [26] "Methanococcus maripaludis"           
    ## [27] "Chlorobium phaeovibrioides"          
    ## [28] "Pyrobaculum arsenaticum"             
    ## [29] "Salinispora tropica"                 
    ## [30] "Caldicellulosiruptor saccharolyticus"
    ## [31] "Thermotoga petrophila"               
    ## [32] "Bacteroides vulgatus"                
    ## [33] "Shewanella baltica"                  
    ## [34] "Ignicoccus hospitalis"               
    ## [35] "Salinispora arenicola"               
    ## [36] "Chloroflexus aurantiacus"            
    ## [37] "Thermoanaerobacter pseudethanolicus" 
    ## [38] "Thermotoga sp. RQ2"                  
    ## [39] "Leptothrix cholodnii"                
    ## [40] "Akkermansia muciniphila"             
    ## [41] "Porphyromonas gingivalis"            
    ## [42] "Sulfurihydrogenibium sp. YO3AOP1"    
    ## [43] "Chlorobium limicola"                 
    ## [44] "Pelodictyon phaeoclathratiforme"     
    ## [45] "Hydrogenobaculum sp. Y04AAS1"        
    ## [46] "Sulfitobacter sp. EE-36"             
    ## [47] "Sulfitobacter sp. NAS-14.1"          
    ## [48] "Dictyoglomus turgidum"               
    ## [49] "Shewanella baltica"                  
    ## [50] "Thermotoga neapolitana"              
    ## [51] "Desulfovibrio piger"                 
    ## [52] "Persephonella marina"                
    ## [53] "Acidobacterium capsulatum"           
    ## [54] "Gemmatimonas aurantiaca"             
    ## [55] "Aciduliprofundum boonei"             
    ## [56] "Haloferax volcanii"                  
    ## [57] "Enterococcus faecalis"               
    ## [58] "Zymomonas mobilis"                   
    ## [59] "Deinococcus radiodurans"             
    ## [60] "Sulfurihydrogenibium yellowstonense" 
    ## [61] "Bordetella bronchiseptica"           
    ## [62] "Caldicellulosiruptor bescii"         
    ## [63] "Herpetosiphon aurantiacus"           
    ## [64] "Nanoarchaeum equitans"

``` r
tax <- tax[,-2]

str_view(head(tax$accession), ".*")
```

    ## [1] │ <AE000782><>
    ## [2] │ <NC_000909><>
    ## [3] │ <NC_003272><>
    ## [4] │ <AE009441><>
    ## [5] │ <AE009950><>
    ## [6] │ <AE009951><>

``` r
#tax$ident <- str_remove_all(tax$ident, "\\..*")
rownames(tax) <- str_glue_data(tax, "{accession} {species}")
tax <- tax[,-1]

OTU = otu_table(df, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax))

physeq = phyloseq(OTU, TAX)
physeq
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 15 taxa and 5 samples ]
    ## tax_table()   Taxonomy Table:    [ 15 taxa by 8 taxonomic ranks ]

## Create a metadata dataframe to use with the sourmash phyloseq object

Read in the csv file that contains your sample metadata

``` r
all_metadata <- read.csv("data/tiny_example.csv")

metadata <- sample_data(all_metadata)

str(metadata)
```

    ## 'data.frame':    5 obs. of  142 variables:
    ## Formal class 'sample_data' [package "phyloseq"] with 4 slots
    ##   ..@ .Data    :List of 142
    ##   .. ..$ : chr  "ZeeviD_2015" "MetaCardis_2020_a" "MetaCardis_2020_a" "GhensiP_2019" ...
    ##   .. ..$ : chr  "PNP_Main_553" "M0x20MCx1812" "M0x20MCx1360" "SP_101SPI_T016" ...
    ##   .. ..$ : chr  "PNP_Main_553" "M0x20MCx1812" "M0x20MCx1360" "sub_101" ...
    ##   .. ..$ : chr  "stool" "stool" "stool" "oralcavity" ...
    ##   .. ..$ : chr  "no" "no" "yes" NA ...
    ##   .. ..$ : chr  "control" "T2D" "IGT" "control" ...
    ##   .. ..$ : chr  "healthy" "T2D" "IGT;MS" "healthy" ...
    ##   .. ..$ : int  51 63 51 NA 27
    ##   .. ..$ : chr  "adult" "adult" "adult" "adult" ...
    ##   .. ..$ : chr  "female" "male" "female" NA ...
    ##   .. ..$ : chr  "ISR" "DEU" "DEU" "ITA" ...
    ##   .. ..$ : chr  NA "Leipzig" "Leipzig" NA ...
    ##   .. ..$ : chr  "no" "no" "no" "no" ...
    ##   .. ..$ : chr  "IlluminaHiSeq" "IonProton" "IonProton" "IlluminaHiSeq" ...
    ##   .. ..$ : chr  NA NA NA NA ...
    ##   .. ..$ : int  26590418 34880489 34880489 33127901 22699609
    ##   .. ..$ : int  28756422 15374062 20888458 20061461 94436101
    ##   .. ..$ : num  2.52e+09 2.33e+09 3.20e+09 1.99e+09 8.71e+09
    ##   .. ..$ : int  20 NA NA 75 60
    ##   .. ..$ : int  80 NA NA 100 101
    ##   .. ..$ : chr  "ERR1136663" "ERR4562131" "ERR4087400" "SRR9217435" ...
    ##   .. ..$ : num  NA NA 6.78 NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : chr  "no" "antidiab;metformin;dppiv" "antihta;ca2_cbl" NA ...
    ##   .. ..$ : chr  "Jacob_Wirbel;Paolo_Manghi" "Paolo_Manghi" "Paolo_Manghi" "Paolo_Manghi" ...
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : chr  "no" NA NA NA ...
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : num  23.7 27.7 50.1 NA 23
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : chr  NA "yes" "no" "no" ...
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : chr  "rectal_swab" NA NA "subgingival_plaque" ...
    ##   .. ..$ : chr  "no_immuno_suppressive;no_T2D;no_T1D;no_related_treatments;no_psychiatric_diseases;no_gastro_intestinal_disorder;non_celiac" NA NA NA ...
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : num  NA NA 210 NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : num  NA NA 240 NA NA
    ##   .. ..$ : num  NA NA 5.73 NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : chr  NA NA NA "yes" ...
    ##   .. ..$ : chr  NA NA NA "implant" ...
    ##   .. ..$ : chr  NA NA NA "yes" ...
    ##   .. ..$ : int  NA NA NA 4 NA
    ##   .. ..$ : int  NA NA NA 3 NA
    ##   .. ..$ : int  NA NA NA 4 NA
    ##   .. ..$ : int  NA NA NA 3 NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. ..$ : logi  NA NA NA NA NA
    ##   .. .. [list output truncated]
    ##   ..@ names    : chr  "study_name" "sample_id" "subject_id" "body_site" ...
    ##   ..@ row.names: chr  "sa1" "sa2" "sa3" "sa4" ...
    ##   ..@ .S3Class : chr "data.frame"

This file contained many columns we do not need. Lets remove those.

``` r
metadata <- metadata[, colSums(!is.na(metadata)) > 0]
str(metadata)
```

    ## 'data.frame':    5 obs. of  41 variables:
    ## Formal class 'sample_data' [package "phyloseq"] with 4 slots
    ##   ..@ .Data    :List of 41
    ##   .. ..$ : chr  "ZeeviD_2015" "MetaCardis_2020_a" "MetaCardis_2020_a" "GhensiP_2019" ...
    ##   .. ..$ : chr  "PNP_Main_553" "M0x20MCx1812" "M0x20MCx1360" "SP_101SPI_T016" ...
    ##   .. ..$ : chr  "PNP_Main_553" "M0x20MCx1812" "M0x20MCx1360" "sub_101" ...
    ##   .. ..$ : chr  "stool" "stool" "stool" "oralcavity" ...
    ##   .. ..$ : chr  "no" "no" "yes" NA ...
    ##   .. ..$ : chr  "control" "T2D" "IGT" "control" ...
    ##   .. ..$ : chr  "healthy" "T2D" "IGT;MS" "healthy" ...
    ##   .. ..$ : int  51 63 51 NA 27
    ##   .. ..$ : chr  "adult" "adult" "adult" "adult" ...
    ##   .. ..$ : chr  "female" "male" "female" NA ...
    ##   .. ..$ : chr  "ISR" "DEU" "DEU" "ITA" ...
    ##   .. ..$ : chr  NA "Leipzig" "Leipzig" NA ...
    ##   .. ..$ : chr  "no" "no" "no" "no" ...
    ##   .. ..$ : chr  "IlluminaHiSeq" "IonProton" "IonProton" "IlluminaHiSeq" ...
    ##   .. ..$ : chr  NA NA NA NA ...
    ##   .. ..$ : int  26590418 34880489 34880489 33127901 22699609
    ##   .. ..$ : int  28756422 15374062 20888458 20061461 94436101
    ##   .. ..$ : num  2.52e+09 2.33e+09 3.20e+09 1.99e+09 8.71e+09
    ##   .. ..$ : int  20 NA NA 75 60
    ##   .. ..$ : int  80 NA NA 100 101
    ##   .. ..$ : chr  "ERR1136663" "ERR4562131" "ERR4087400" "SRR9217435" ...
    ##   .. ..$ : num  NA NA 6.78 NA NA
    ##   .. ..$ : chr  "no" "antidiab;metformin;dppiv" "antihta;ca2_cbl" NA ...
    ##   .. ..$ : chr  "Jacob_Wirbel;Paolo_Manghi" "Paolo_Manghi" "Paolo_Manghi" "Paolo_Manghi" ...
    ##   .. ..$ : chr  "no" NA NA NA ...
    ##   .. ..$ : num  23.7 27.7 50.1 NA 23
    ##   .. ..$ : chr  NA "yes" "no" "no" ...
    ##   .. ..$ : chr  "rectal_swab" NA NA "subgingival_plaque" ...
    ##   .. ..$ : chr  "no_immuno_suppressive;no_T2D;no_T1D;no_related_treatments;no_psychiatric_diseases;no_gastro_intestinal_disorder;non_celiac" NA NA NA ...
    ##   .. ..$ : num  NA NA 210 NA NA
    ##   .. ..$ : num  NA NA 240 NA NA
    ##   .. ..$ : num  NA NA 5.73 NA NA
    ##   .. ..$ : chr  NA NA NA "yes" ...
    ##   .. ..$ : chr  NA NA NA "implant" ...
    ##   .. ..$ : chr  NA NA NA "yes" ...
    ##   .. ..$ : int  NA NA NA 4 NA
    ##   .. ..$ : int  NA NA NA 3 NA
    ##   .. ..$ : int  NA NA NA 4 NA
    ##   .. ..$ : int  NA NA NA 3 NA
    ##   .. ..$ : int  NA NA 4 NA NA
    ##   .. ..$ : chr  "PRJEB11532" "PRJEB38742" "PRJEB37249" "PRJNA547717" ...
    ##   ..@ names    : chr  "study_name" "sample_id" "subject_id" "body_site" ...
    ##   ..@ row.names: chr  "sa1" "sa2" "sa3" "sa4" ...
    ##   ..@ .S3Class : chr "data.frame"

Format the dataframe to work well with the sourmash phyloseq object

``` r
names(metadata) <- tolower(names(metadata))
rownames(metadata) <- metadata$ncbi_accession
DT::datatable(metadata)
```

<div class="datatables html-widget html-fill-item" id="htmlwidget-b18601fa64fe400df915" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-b18601fa64fe400df915">{"x":{"filter":"none","vertical":false,"data":[["ERR1136663","ERR4562131","ERR4087400","SRR9217435","SRR059888"],["ZeeviD_2015","MetaCardis_2020_a","MetaCardis_2020_a","GhensiP_2019","HMP_2012"],["PNP_Main_553","M0x20MCx1812","M0x20MCx1360","SP_101SPI_T016","SRS013723"],["PNP_Main_553","M0x20MCx1812","M0x20MCx1360","sub_101","HMP_2012_159268001"],["stool","stool","stool","oralcavity","oralcavity"],["no","no","yes",null,null],["control","T2D","IGT","control","control"],["healthy","T2D","IGT;MS","healthy","healthy"],[51,63,51,null,27],["adult","adult","adult","adult","adult"],["female","male","female",null,"male"],["ISR","DEU","DEU","ITA","USA"],[null,"Leipzig","Leipzig",null,null],["no","no","no","no","no"],["IlluminaHiSeq","IonProton","IonProton","IlluminaHiSeq","IlluminaHiSeq"],[null,null,null,null,"Qiagen"],[26590418,34880489,34880489,33127901,22699609],[28756422,15374062,20888458,20061461,94436101],[2515847316,2326962418,3200734644,1987501438,8712278548],[20,null,null,75,60],[80,null,null,100,101],["ERR1136663","ERR4562131","ERR4087400","SRR9217435","SRR059888"],[null,null,6.78,null,null],["no","antidiab;metformin;dppiv","antihta;ca2_cbl",null,null],["Jacob_Wirbel;Paolo_Manghi","Paolo_Manghi","Paolo_Manghi","Paolo_Manghi","Paolo_Manghi"],["no",null,null,null,null],[23.7,27.6795859109,50.1440220699,null,23],[null,"yes","no","no",null],["rectal_swab",null,null,"subgingival_plaque","supragingival_plaque"],["no_immuno_suppressive;no_T2D;no_T1D;no_related_treatments;no_psychiatric_diseases;no_gastro_intestinal_disorder;non_celiac",null,null,null,null],[null,null,210.3648,null,null],[null,null,240.0247,null,null],[null,null,5.73,null,null],[null,null,null,"yes",null],[null,null,null,"implant",null],[null,null,null,"yes",null],[null,null,null,4,null],[null,null,null,3,null],[null,null,null,4,null],[null,null,null,3,null],[null,null,4,null,null],["PRJEB11532","PRJEB38742","PRJEB37249","PRJNA547717","PRJNA48479"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>study_name<\/th>\n      <th>sample_id<\/th>\n      <th>subject_id<\/th>\n      <th>body_site<\/th>\n      <th>antibiotics_current_use<\/th>\n      <th>study_condition<\/th>\n      <th>disease<\/th>\n      <th>age<\/th>\n      <th>age_category<\/th>\n      <th>gender<\/th>\n      <th>country<\/th>\n      <th>location<\/th>\n      <th>non_westernized<\/th>\n      <th>sequencing_platform<\/th>\n      <th>dna_extraction_kit<\/th>\n      <th>pmid<\/th>\n      <th>number_reads<\/th>\n      <th>number_bases<\/th>\n      <th>minimum_read_length<\/th>\n      <th>median_read_length<\/th>\n      <th>ncbi_accession<\/th>\n      <th>hscrp<\/th>\n      <th>treatment<\/th>\n      <th>curator<\/th>\n      <th>pregnant<\/th>\n      <th>bmi<\/th>\n      <th>smoker<\/th>\n      <th>body_subsite<\/th>\n      <th>uncurated_metadata<\/th>\n      <th>ldl<\/th>\n      <th>triglycerides<\/th>\n      <th>hba1c<\/th>\n      <th>ever_smoker<\/th>\n      <th>dental_sample_type<\/th>\n      <th>history_of_periodontitis<\/th>\n      <th>ppd_m<\/th>\n      <th>ppd_b<\/th>\n      <th>ppd_d<\/th>\n      <th>ppd_l<\/th>\n      <th>bristol_score<\/th>\n      <th>ncbi_bioproject<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[8,16,17,18,19,20,22,26,30,31,32,36,37,38,39,40]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"study_name","targets":1},{"name":"sample_id","targets":2},{"name":"subject_id","targets":3},{"name":"body_site","targets":4},{"name":"antibiotics_current_use","targets":5},{"name":"study_condition","targets":6},{"name":"disease","targets":7},{"name":"age","targets":8},{"name":"age_category","targets":9},{"name":"gender","targets":10},{"name":"country","targets":11},{"name":"location","targets":12},{"name":"non_westernized","targets":13},{"name":"sequencing_platform","targets":14},{"name":"dna_extraction_kit","targets":15},{"name":"pmid","targets":16},{"name":"number_reads","targets":17},{"name":"number_bases","targets":18},{"name":"minimum_read_length","targets":19},{"name":"median_read_length","targets":20},{"name":"ncbi_accession","targets":21},{"name":"hscrp","targets":22},{"name":"treatment","targets":23},{"name":"curator","targets":24},{"name":"pregnant","targets":25},{"name":"bmi","targets":26},{"name":"smoker","targets":27},{"name":"body_subsite","targets":28},{"name":"uncurated_metadata","targets":29},{"name":"ldl","targets":30},{"name":"triglycerides","targets":31},{"name":"hba1c","targets":32},{"name":"ever_smoker","targets":33},{"name":"dental_sample_type","targets":34},{"name":"history_of_periodontitis","targets":35},{"name":"ppd_m","targets":36},{"name":"ppd_b","targets":37},{"name":"ppd_d","targets":38},{"name":"ppd_l","targets":39},{"name":"bristol_score","targets":40},{"name":"ncbi_bioproject","targets":41}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

Finally, merge the metadata with the sourmash phyloseq object

``` r
pseq_final <- merge_phyloseq(physeq, metadata)
pseq_final
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 15 taxa and 5 samples ]
    ## sample_data() Sample Data:       [ 5 samples by 41 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 15 taxa by 8 taxonomic ranks ]
