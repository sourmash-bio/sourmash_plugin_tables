
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

    ##                                                                                          match_name
    ## 1                                        CP000139.1 Bacteroides vulgatus ATCC 8482, complete genome
    ## 2                                 AE015928.1 Bacteroides thetaiotaomicron VPI-5482, complete genome
    ## 3              NZ_DS996397.1 Desulfovibrio piger ATCC 29098 Scfld442, whole genome shotgun sequence
    ## 4                                  CP001071.1 Akkermansia muciniphila ATCC BAA-835, complete genome
    ## 5        NZ_KE136524.1 Enterococcus faecalis V583 acyDH-supercont2.1, whole genome shotgun sequence
    ## 6                                  CP001472.1 Acidobacterium capsulatum ATCC 51196, complete genome
    ## 7                   AE009951.2 Fusobacterium nucleatum subsp. nucleatum ATCC 25586, complete genome
    ## 8                                             CP001013.1 Leptothrix cholodnii SP-6, complete genome
    ## 9  NZ_JGWU01000001.1 Bordetella bronchiseptica D989 ctg7180000008197, whole genome shotgun sequence
    ## 10       NZ_KQ961402.1 Zymomonas mobilis strain ATCC 31823 Scaffold1, whole genome shotgun sequence
    ## 11                              AP009380.1 Porphyromonas gingivalis ATCC 33277 DNA, complete genome
    ## 12                        NC_007951.1 Burkholderia xenovorans LB400 chromosome 1, complete sequence
    ## 13                                       AE017226.1 Treponema denticola ATCC 35405, complete genome
    ## 14                                        CP000850.1 Salinispora arenicola CNS-205, complete genome
    ## 15                                               AE006470.1 Chlorobium tepidum TLS, complete genome
    ## 16                                            AE017221.1 Thermus thermophilus HB27, complete genome
    ## 17                                            NC_011663.1 Shewanella baltica OS223, complete genome
    ##      ERR1136663   ERR4087400   ERR4562131    SRR059888   SRR9217435
    ## 1  8.301477e-03 5.512089e-03 5.823024e-03 3.356448e-05 5.792533e-05
    ## 2  2.032678e-03 3.968704e-04 4.229973e-04 0.000000e+00 0.000000e+00
    ## 3  1.076401e-03 0.000000e+00 1.548600e-04 0.000000e+00 0.000000e+00
    ## 4  5.676425e-04 4.119893e-03 4.301667e-05 9.811156e-06 0.000000e+00
    ## 5  2.355363e-05 6.929483e-05 0.000000e+00 0.000000e+00 0.000000e+00
    ## 6  1.165905e-04 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
    ## 7  2.814659e-04 0.000000e+00 0.000000e+00 3.217285e-03 5.649507e-05
    ## 8  9.421452e-06 0.000000e+00 0.000000e+00 1.587858e-04 1.923693e-04
    ## 9  0.000000e+00 3.779718e-05 0.000000e+00 2.793598e-04 3.983260e-04
    ## 10 0.000000e+00 0.000000e+00 3.728112e-04 0.000000e+00 2.932023e-05
    ## 11 0.000000e+00 0.000000e+00 5.735556e-05 2.765197e-04 8.581530e-06
    ## 12 0.000000e+00 0.000000e+00 2.581000e-05 0.000000e+00 4.290765e-06
    ## 13 0.000000e+00 0.000000e+00 0.000000e+00 7.255092e-05 0.000000e+00
    ## 14 0.000000e+00 0.000000e+00 0.000000e+00 3.464887e-04 0.000000e+00
    ## 15 0.000000e+00 0.000000e+00 0.000000e+00 7.539099e-05 0.000000e+00
    ## 16 0.000000e+00 0.000000e+00 0.000000e+00 3.356448e-06 0.000000e+00
    ## 17 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.623339e-04

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
df
```

    ##                                             ERR1136663   ERR4087400
    ## CP000139 Bacteroides vulgatus             8.301477e-03 5.512089e-03
    ## AE015928 Bacteroides thetaiotaomicron     2.032678e-03 3.968704e-04
    ## NZ_DS996397 Desulfovibrio piger           1.076401e-03 0.000000e+00
    ## CP001071 Akkermansia muciniphila          5.676425e-04 4.119893e-03
    ## NZ_KE136524 Enterococcus faecalis         2.355363e-05 6.929483e-05
    ## CP001472 Acidobacterium capsulatum        1.165905e-04 0.000000e+00
    ## AE009951 Fusobacterium nucleatum          2.814659e-04 0.000000e+00
    ## CP001013 Leptothrix cholodnii             9.421452e-06 0.000000e+00
    ## NZ_JGWU01000001 Bordetella bronchiseptica 0.000000e+00 3.779718e-05
    ## NZ_KQ961402 Zymomonas mobilis             0.000000e+00 0.000000e+00
    ## AP009380 Porphyromonas gingivalis         0.000000e+00 0.000000e+00
    ## NC_007951 Burkholderia xenovorans         0.000000e+00 0.000000e+00
    ## AE017226 Treponema denticola              0.000000e+00 0.000000e+00
    ## CP000850 Salinispora arenicola            0.000000e+00 0.000000e+00
    ## AE006470 Chlorobium tepidum               0.000000e+00 0.000000e+00
    ## AE017221 Thermus thermophilus             0.000000e+00 0.000000e+00
    ## NC_011663 Shewanella baltica              0.000000e+00 0.000000e+00
    ##                                             ERR4562131    SRR059888
    ## CP000139 Bacteroides vulgatus             5.823024e-03 3.356448e-05
    ## AE015928 Bacteroides thetaiotaomicron     4.229973e-04 0.000000e+00
    ## NZ_DS996397 Desulfovibrio piger           1.548600e-04 0.000000e+00
    ## CP001071 Akkermansia muciniphila          4.301667e-05 9.811156e-06
    ## NZ_KE136524 Enterococcus faecalis         0.000000e+00 0.000000e+00
    ## CP001472 Acidobacterium capsulatum        0.000000e+00 0.000000e+00
    ## AE009951 Fusobacterium nucleatum          0.000000e+00 3.217285e-03
    ## CP001013 Leptothrix cholodnii             0.000000e+00 1.587858e-04
    ## NZ_JGWU01000001 Bordetella bronchiseptica 0.000000e+00 2.793598e-04
    ## NZ_KQ961402 Zymomonas mobilis             3.728112e-04 0.000000e+00
    ## AP009380 Porphyromonas gingivalis         5.735556e-05 2.765197e-04
    ## NC_007951 Burkholderia xenovorans         2.581000e-05 0.000000e+00
    ## AE017226 Treponema denticola              0.000000e+00 7.255092e-05
    ## CP000850 Salinispora arenicola            0.000000e+00 3.464887e-04
    ## AE006470 Chlorobium tepidum               0.000000e+00 7.539099e-05
    ## AE017221 Thermus thermophilus             0.000000e+00 3.356448e-06
    ## NC_011663 Shewanella baltica              0.000000e+00 0.000000e+00
    ##                                             SRR9217435
    ## CP000139 Bacteroides vulgatus             5.792533e-05
    ## AE015928 Bacteroides thetaiotaomicron     0.000000e+00
    ## NZ_DS996397 Desulfovibrio piger           0.000000e+00
    ## CP001071 Akkermansia muciniphila          0.000000e+00
    ## NZ_KE136524 Enterococcus faecalis         0.000000e+00
    ## CP001472 Acidobacterium capsulatum        0.000000e+00
    ## AE009951 Fusobacterium nucleatum          5.649507e-05
    ## CP001013 Leptothrix cholodnii             1.923693e-04
    ## NZ_JGWU01000001 Bordetella bronchiseptica 3.983260e-04
    ## NZ_KQ961402 Zymomonas mobilis             2.932023e-05
    ## AP009380 Porphyromonas gingivalis         8.581530e-06
    ## NC_007951 Burkholderia xenovorans         4.290765e-06
    ## AE017226 Treponema denticola              0.000000e+00
    ## CP000850 Salinispora arenicola            0.000000e+00
    ## AE006470 Chlorobium tepidum               0.000000e+00
    ## AE017221 Thermus thermophilus             0.000000e+00
    ## NC_011663 Shewanella baltica              1.623339e-04

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
tax
```

    ##                                                     superkingdom
    ## AE000782 Archaeoglobus fulgidus                          Archaea
    ## NC_000909 Methanocaldococcus jannaschii                  Archaea
    ## NC_003272 Nostoc sp. PCC 7120                           Bacteria
    ## AE009441 Pyrobaculum aerophilum                          Archaea
    ## AE009950 Pyrococcus furiosus                             Archaea
    ## AE009951 Fusobacterium nucleatum                        Bacteria
    ## AE010299 Methanosarcina acetivorans                      Archaea
    ## AE009439 Methanopyrus kandleri                           Archaea
    ## NC_003911 Ruegeria pomeroyi                             Bacteria
    ## AE006470 Chlorobaculum tepidum                          Bacteria
    ## AE015928 Bacteroides thetaiotaomicron                   Bacteria
    ## AL954747 Nitrosomonas europaea                          Bacteria
    ## BX119912 Rhodopirellula baltica                         Bacteria
    ## BX571656 Wolinella succinogenes                         Bacteria
    ## AE017180 Geobacter sulfurreducens                       Bacteria
    ## AE017226 Treponema denticola                            Bacteria
    ## BX950229 Methanococcus maripaludis                       Archaea
    ## AE017221 Thermus thermophilus                           Bacteria
    ## BA000001 Pyrococcus horikoshii                           Archaea
    ## BA000023 Sulfolobus tokodaii                             Archaea
    ## NC_007951 Paraburkholderia xenovorans                   Bacteria
    ## CP000492 Chlorobium phaeobacteroides                    Bacteria
    ## NC_008751 Desulfovibrio vulgaris                        Bacteria
    ## CP000568 Ruminiclostridium thermocellum                 Bacteria
    ## CP000561 Pyrobaculum calidifontis                        Archaea
    ## CP000609 Methanococcus maripaludis                       Archaea
    ## CP000607 Chlorobium phaeovibrioides                     Bacteria
    ## CP000660 Pyrobaculum arsenaticum                         Archaea
    ## CP000667 Salinispora tropica                            Bacteria
    ## CP000679 Caldicellulosiruptor saccharolyticus           Bacteria
    ## CP000702 Thermotoga petrophila                          Bacteria
    ## CP000139 Bacteroides vulgatus                           Bacteria
    ## NC_009665 Shewanella baltica                            Bacteria
    ## CP000816 Ignicoccus hospitalis                           Archaea
    ## CP000850 Salinispora arenicola                          Bacteria
    ## CP000909 Chloroflexus aurantiacus                       Bacteria
    ## CP000924 Thermoanaerobacter pseudethanolicus            Bacteria
    ## CP000969 Thermotoga sp. RQ2                             Bacteria
    ## CP001013 Leptothrix cholodnii                           Bacteria
    ## CP001071 Akkermansia muciniphila                        Bacteria
    ## AP009380 Porphyromonas gingivalis                       Bacteria
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1              Bacteria
    ## CP001097 Chlorobium limicola                            Bacteria
    ## CP001110 Pelodictyon phaeoclathratiforme                Bacteria
    ## CP001130 Hydrogenobaculum sp. Y04AAS1                   Bacteria
    ## NZ_CH959311 Sulfitobacter sp. EE-36                     Bacteria
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1                  Bacteria
    ## CP001251 Dictyoglomus turgidum                          Bacteria
    ## NC_011663 Shewanella baltica                            Bacteria
    ## CP000916 Thermotoga neapolitana                         Bacteria
    ## NZ_DS996397 Desulfovibrio piger                         Bacteria
    ## CP001230 Persephonella marina                           Bacteria
    ## CP001472 Acidobacterium capsulatum                      Bacteria
    ## AP009153 Gemmatimonas aurantiaca                        Bacteria
    ## CP001941 Aciduliprofundum boonei                         Archaea
    ## NC_013968 Haloferax volcanii                             Archaea
    ## NZ_KE136524 Enterococcus faecalis                       Bacteria
    ## NZ_KQ961402 Zymomonas mobilis                           Bacteria
    ## NZ_CP015081 Deinococcus radiodurans                     Bacteria
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense     Bacteria
    ## NZ_JGWU01000001 Bordetella bronchiseptica               Bacteria
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii             Bacteria
    ## NC_009972 Herpetosiphon aurantiacus                     Bacteria
    ## NC_005213 Nanoarchaeum equitans                          Archaea
    ##                                                                  phylum
    ## AE000782 Archaeoglobus fulgidus                           Euryarchaeota
    ## NC_000909 Methanocaldococcus jannaschii                   Euryarchaeota
    ## NC_003272 Nostoc sp. PCC 7120                             Cyanobacteria
    ## AE009441 Pyrobaculum aerophilum                           Crenarchaeota
    ## AE009950 Pyrococcus furiosus                              Euryarchaeota
    ## AE009951 Fusobacterium nucleatum                           Fusobacteria
    ## AE010299 Methanosarcina acetivorans                       Euryarchaeota
    ## AE009439 Methanopyrus kandleri                            Euryarchaeota
    ## NC_003911 Ruegeria pomeroyi                              Proteobacteria
    ## AE006470 Chlorobaculum tepidum                                 Chlorobi
    ## AE015928 Bacteroides thetaiotaomicron                     Bacteroidetes
    ## AL954747 Nitrosomonas europaea                           Proteobacteria
    ## BX119912 Rhodopirellula baltica                          Planctomycetes
    ## BX571656 Wolinella succinogenes                          Proteobacteria
    ## AE017180 Geobacter sulfurreducens                        Proteobacteria
    ## AE017226 Treponema denticola                               Spirochaetes
    ## BX950229 Methanococcus maripaludis                        Euryarchaeota
    ## AE017221 Thermus thermophilus                       Deinococcus-Thermus
    ## BA000001 Pyrococcus horikoshii                            Euryarchaeota
    ## BA000023 Sulfolobus tokodaii                              Crenarchaeota
    ## NC_007951 Paraburkholderia xenovorans                    Proteobacteria
    ## CP000492 Chlorobium phaeobacteroides                           Chlorobi
    ## NC_008751 Desulfovibrio vulgaris                         Proteobacteria
    ## CP000568 Ruminiclostridium thermocellum                      Firmicutes
    ## CP000561 Pyrobaculum calidifontis                         Crenarchaeota
    ## CP000609 Methanococcus maripaludis                        Euryarchaeota
    ## CP000607 Chlorobium phaeovibrioides                            Chlorobi
    ## CP000660 Pyrobaculum arsenaticum                          Crenarchaeota
    ## CP000667 Salinispora tropica                             Actinobacteria
    ## CP000679 Caldicellulosiruptor saccharolyticus                Firmicutes
    ## CP000702 Thermotoga petrophila                              Thermotogae
    ## CP000139 Bacteroides vulgatus                             Bacteroidetes
    ## NC_009665 Shewanella baltica                             Proteobacteria
    ## CP000816 Ignicoccus hospitalis                            Crenarchaeota
    ## CP000850 Salinispora arenicola                           Actinobacteria
    ## CP000909 Chloroflexus aurantiacus                           Chloroflexi
    ## CP000924 Thermoanaerobacter pseudethanolicus                 Firmicutes
    ## CP000969 Thermotoga sp. RQ2                                 Thermotogae
    ## CP001013 Leptothrix cholodnii                            Proteobacteria
    ## CP001071 Akkermansia muciniphila                        Verrucomicrobia
    ## AP009380 Porphyromonas gingivalis                         Bacteroidetes
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1                    Aquificae
    ## CP001097 Chlorobium limicola                                   Chlorobi
    ## CP001110 Pelodictyon phaeoclathratiforme                       Chlorobi
    ## CP001130 Hydrogenobaculum sp. Y04AAS1                         Aquificae
    ## NZ_CH959311 Sulfitobacter sp. EE-36                      Proteobacteria
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1                   Proteobacteria
    ## CP001251 Dictyoglomus turgidum                              Dictyoglomi
    ## NC_011663 Shewanella baltica                             Proteobacteria
    ## CP000916 Thermotoga neapolitana                             Thermotogae
    ## NZ_DS996397 Desulfovibrio piger                          Proteobacteria
    ## CP001230 Persephonella marina                                 Aquificae
    ## CP001472 Acidobacterium capsulatum                        Acidobacteria
    ## AP009153 Gemmatimonas aurantiaca                       Gemmatimonadetes
    ## CP001941 Aciduliprofundum boonei                          Euryarchaeota
    ## NC_013968 Haloferax volcanii                              Euryarchaeota
    ## NZ_KE136524 Enterococcus faecalis                            Firmicutes
    ## NZ_KQ961402 Zymomonas mobilis                            Proteobacteria
    ## NZ_CP015081 Deinococcus radiodurans                 Deinococcus-Thermus
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense           Aquificae
    ## NZ_JGWU01000001 Bordetella bronchiseptica                Proteobacteria
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii                  Firmicutes
    ## NC_009972 Herpetosiphon aurantiacus                         Chloroflexi
    ## NC_005213 Nanoarchaeum equitans                           Nanoarchaeota
    ##                                                                     class
    ## AE000782 Archaeoglobus fulgidus                              Archaeoglobi
    ## NC_000909 Methanocaldococcus jannaschii                      Methanococci
    ## NC_003272 Nostoc sp. PCC 7120                                            
    ## AE009441 Pyrobaculum aerophilum                              Thermoprotei
    ## AE009950 Pyrococcus furiosus                                  Thermococci
    ## AE009951 Fusobacterium nucleatum                            Fusobacteriia
    ## AE010299 Methanosarcina acetivorans                       Methanomicrobia
    ## AE009439 Methanopyrus kandleri                                Methanopyri
    ## NC_003911 Ruegeria pomeroyi                           Alphaproteobacteria
    ## AE006470 Chlorobaculum tepidum                                  Chlorobia
    ## AE015928 Bacteroides thetaiotaomicron                         Bacteroidia
    ## AL954747 Nitrosomonas europaea                         Betaproteobacteria
    ## BX119912 Rhodopirellula baltica                            Planctomycetia
    ## BX571656 Wolinella succinogenes                     Epsilonproteobacteria
    ## AE017180 Geobacter sulfurreducens                     Deltaproteobacteria
    ## AE017226 Treponema denticola                                 Spirochaetia
    ## BX950229 Methanococcus maripaludis                           Methanococci
    ## AE017221 Thermus thermophilus                                  Deinococci
    ## BA000001 Pyrococcus horikoshii                                Thermococci
    ## BA000023 Sulfolobus tokodaii                                 Thermoprotei
    ## NC_007951 Paraburkholderia xenovorans                  Betaproteobacteria
    ## CP000492 Chlorobium phaeobacteroides                            Chlorobia
    ## NC_008751 Desulfovibrio vulgaris                      Deltaproteobacteria
    ## CP000568 Ruminiclostridium thermocellum                        Clostridia
    ## CP000561 Pyrobaculum calidifontis                            Thermoprotei
    ## CP000609 Methanococcus maripaludis                           Methanococci
    ## CP000607 Chlorobium phaeovibrioides                             Chlorobia
    ## CP000660 Pyrobaculum arsenaticum                             Thermoprotei
    ## CP000667 Salinispora tropica                               Actinobacteria
    ## CP000679 Caldicellulosiruptor saccharolyticus                  Clostridia
    ## CP000702 Thermotoga petrophila                                Thermotogae
    ## CP000139 Bacteroides vulgatus                                 Bacteroidia
    ## NC_009665 Shewanella baltica                          Gammaproteobacteria
    ## CP000816 Ignicoccus hospitalis                               Thermoprotei
    ## CP000850 Salinispora arenicola                             Actinobacteria
    ## CP000909 Chloroflexus aurantiacus                            Chloroflexia
    ## CP000924 Thermoanaerobacter pseudethanolicus                   Clostridia
    ## CP000969 Thermotoga sp. RQ2                                   Thermotogae
    ## CP001013 Leptothrix cholodnii                          Betaproteobacteria
    ## CP001071 Akkermansia muciniphila                         Verrucomicrobiae
    ## AP009380 Porphyromonas gingivalis                             Bacteroidia
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1                      Aquificae
    ## CP001097 Chlorobium limicola                                    Chlorobia
    ## CP001110 Pelodictyon phaeoclathratiforme                        Chlorobia
    ## CP001130 Hydrogenobaculum sp. Y04AAS1                           Aquificae
    ## NZ_CH959311 Sulfitobacter sp. EE-36                   Alphaproteobacteria
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1                Alphaproteobacteria
    ## CP001251 Dictyoglomus turgidum                               Dictyoglomia
    ## NC_011663 Shewanella baltica                          Gammaproteobacteria
    ## CP000916 Thermotoga neapolitana                               Thermotogae
    ## NZ_DS996397 Desulfovibrio piger                       Deltaproteobacteria
    ## CP001230 Persephonella marina                                   Aquificae
    ## CP001472 Acidobacterium capsulatum                         Acidobacteriia
    ## AP009153 Gemmatimonas aurantiaca                         Gemmatimonadetes
    ## CP001941 Aciduliprofundum boonei                                         
    ## NC_013968 Haloferax volcanii                                 Halobacteria
    ## NZ_KE136524 Enterococcus faecalis                                 Bacilli
    ## NZ_KQ961402 Zymomonas mobilis                         Alphaproteobacteria
    ## NZ_CP015081 Deinococcus radiodurans                            Deinococci
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense             Aquificae
    ## NZ_JGWU01000001 Bordetella bronchiseptica              Betaproteobacteria
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii                    Clostridia
    ## NC_009972 Herpetosiphon aurantiacus                          Chloroflexia
    ## NC_005213 Nanoarchaeum equitans                                          
    ##                                                                      order
    ## AE000782 Archaeoglobus fulgidus                            Archaeoglobales
    ## NC_000909 Methanocaldococcus jannaschii                    Methanococcales
    ## NC_003272 Nostoc sp. PCC 7120                                   Nostocales
    ## AE009441 Pyrobaculum aerophilum                            Thermoproteales
    ## AE009950 Pyrococcus furiosus                                Thermococcales
    ## AE009951 Fusobacterium nucleatum                           Fusobacteriales
    ## AE010299 Methanosarcina acetivorans                      Methanosarcinales
    ## AE009439 Methanopyrus kandleri                              Methanopyrales
    ## NC_003911 Ruegeria pomeroyi                                Rhodobacterales
    ## AE006470 Chlorobaculum tepidum                                Chlorobiales
    ## AE015928 Bacteroides thetaiotaomicron                        Bacteroidales
    ## AL954747 Nitrosomonas europaea                            Nitrosomonadales
    ## BX119912 Rhodopirellula baltica                           Planctomycetales
    ## BX571656 Wolinella succinogenes                          Campylobacterales
    ## AE017180 Geobacter sulfurreducens                       Desulfuromonadales
    ## AE017226 Treponema denticola                                Spirochaetales
    ## BX950229 Methanococcus maripaludis                         Methanococcales
    ## AE017221 Thermus thermophilus                                    Thermales
    ## BA000001 Pyrococcus horikoshii                              Thermococcales
    ## BA000023 Sulfolobus tokodaii                                  Sulfolobales
    ## NC_007951 Paraburkholderia xenovorans                      Burkholderiales
    ## CP000492 Chlorobium phaeobacteroides                          Chlorobiales
    ## NC_008751 Desulfovibrio vulgaris                        Desulfovibrionales
    ## CP000568 Ruminiclostridium thermocellum                      Clostridiales
    ## CP000561 Pyrobaculum calidifontis                          Thermoproteales
    ## CP000609 Methanococcus maripaludis                         Methanococcales
    ## CP000607 Chlorobium phaeovibrioides                           Chlorobiales
    ## CP000660 Pyrobaculum arsenaticum                           Thermoproteales
    ## CP000667 Salinispora tropica                             Micromonosporales
    ## CP000679 Caldicellulosiruptor saccharolyticus       Thermoanaerobacterales
    ## CP000702 Thermotoga petrophila                               Thermotogales
    ## CP000139 Bacteroides vulgatus                                Bacteroidales
    ## NC_009665 Shewanella baltica                               Alteromonadales
    ## CP000816 Ignicoccus hospitalis                           Desulfurococcales
    ## CP000850 Salinispora arenicola                           Micromonosporales
    ## CP000909 Chloroflexus aurantiacus                           Chloroflexales
    ## CP000924 Thermoanaerobacter pseudethanolicus        Thermoanaerobacterales
    ## CP000969 Thermotoga sp. RQ2                                  Thermotogales
    ## CP001013 Leptothrix cholodnii                              Burkholderiales
    ## CP001071 Akkermansia muciniphila                        Verrucomicrobiales
    ## AP009380 Porphyromonas gingivalis                            Bacteroidales
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1                     Aquificales
    ## CP001097 Chlorobium limicola                                  Chlorobiales
    ## CP001110 Pelodictyon phaeoclathratiforme                      Chlorobiales
    ## CP001130 Hydrogenobaculum sp. Y04AAS1                          Aquificales
    ## NZ_CH959311 Sulfitobacter sp. EE-36                        Rhodobacterales
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1                     Rhodobacterales
    ## CP001251 Dictyoglomus turgidum                              Dictyoglomales
    ## NC_011663 Shewanella baltica                               Alteromonadales
    ## CP000916 Thermotoga neapolitana                              Thermotogales
    ## NZ_DS996397 Desulfovibrio piger                         Desulfovibrionales
    ## CP001230 Persephonella marina                                  Aquificales
    ## CP001472 Acidobacterium capsulatum                        Acidobacteriales
    ## AP009153 Gemmatimonas aurantiaca                          Gemmatimonadales
    ## CP001941 Aciduliprofundum boonei                                          
    ## NC_013968 Haloferax volcanii                                 Haloferacales
    ## NZ_KE136524 Enterococcus faecalis                          Lactobacillales
    ## NZ_KQ961402 Zymomonas mobilis                             Sphingomonadales
    ## NZ_CP015081 Deinococcus radiodurans                          Deinococcales
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense            Aquificales
    ## NZ_JGWU01000001 Bordetella bronchiseptica                  Burkholderiales
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii         Thermoanaerobacterales
    ## NC_009972 Herpetosiphon aurantiacus                      Herpetosiphonales
    ## NC_005213 Nanoarchaeum equitans                             Nanoarchaeales
    ##                                                                                                family
    ## AE000782 Archaeoglobus fulgidus                                                      Archaeoglobaceae
    ## NC_000909 Methanocaldococcus jannaschii                                         Methanocaldococcaceae
    ## NC_003272 Nostoc sp. PCC 7120                                                             Nostocaceae
    ## AE009441 Pyrobaculum aerophilum                                                      Thermoproteaceae
    ## AE009950 Pyrococcus furiosus                                                          Thermococcaceae
    ## AE009951 Fusobacterium nucleatum                                                     Fusobacteriaceae
    ## AE010299 Methanosarcina acetivorans                                                Methanosarcinaceae
    ## AE009439 Methanopyrus kandleri                                                        Methanopyraceae
    ## NC_003911 Ruegeria pomeroyi                                                          Rhodobacteraceae
    ## AE006470 Chlorobaculum tepidum                                                          Chlorobiaceae
    ## AE015928 Bacteroides thetaiotaomicron                                                  Bacteroidaceae
    ## AL954747 Nitrosomonas europaea                                                      Nitrosomonadaceae
    ## BX119912 Rhodopirellula baltica                                                     Planctomycetaceae
    ## BX571656 Wolinella succinogenes                                                     Helicobacteraceae
    ## AE017180 Geobacter sulfurreducens                                                      Geobacteraceae
    ## AE017226 Treponema denticola                                                          Spirochaetaceae
    ## BX950229 Methanococcus maripaludis                                                   Methanococcaceae
    ## AE017221 Thermus thermophilus                                                              Thermaceae
    ## BA000001 Pyrococcus horikoshii                                                        Thermococcaceae
    ## BA000023 Sulfolobus tokodaii                                                            Sulfolobaceae
    ## NC_007951 Paraburkholderia xenovorans                                                Burkholderiaceae
    ## CP000492 Chlorobium phaeobacteroides                                                    Chlorobiaceae
    ## NC_008751 Desulfovibrio vulgaris                                                  Desulfovibrionaceae
    ## CP000568 Ruminiclostridium thermocellum                                               Ruminococcaceae
    ## CP000561 Pyrobaculum calidifontis                                                    Thermoproteaceae
    ## CP000609 Methanococcus maripaludis                                                   Methanococcaceae
    ## CP000607 Chlorobium phaeovibrioides                                                     Chlorobiaceae
    ## CP000660 Pyrobaculum arsenaticum                                                     Thermoproteaceae
    ## CP000667 Salinispora tropica                                                       Micromonosporaceae
    ## CP000679 Caldicellulosiruptor saccharolyticus       Thermoanaerobacterales Family III. Incertae Sedis
    ## CP000702 Thermotoga petrophila                                                         Thermotogaceae
    ## CP000139 Bacteroides vulgatus                                                          Bacteroidaceae
    ## NC_009665 Shewanella baltica                                                           Shewanellaceae
    ## CP000816 Ignicoccus hospitalis                                                     Desulfurococcaceae
    ## CP000850 Salinispora arenicola                                                     Micromonosporaceae
    ## CP000909 Chloroflexus aurantiacus                                                     Chloroflexaceae
    ## CP000924 Thermoanaerobacter pseudethanolicus                                  Thermoanaerobacteraceae
    ## CP000969 Thermotoga sp. RQ2                                                            Thermotogaceae
    ## CP001013 Leptothrix cholodnii                                                                        
    ## CP001071 Akkermansia muciniphila                                                      Akkermansiaceae
    ## AP009380 Porphyromonas gingivalis                                                  Porphyromonadaceae
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1                                        Hydrogenothermaceae
    ## CP001097 Chlorobium limicola                                                            Chlorobiaceae
    ## CP001110 Pelodictyon phaeoclathratiforme                                                Chlorobiaceae
    ## CP001130 Hydrogenobaculum sp. Y04AAS1                                                    Aquificaceae
    ## NZ_CH959311 Sulfitobacter sp. EE-36                                                  Rhodobacteraceae
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1                                               Rhodobacteraceae
    ## CP001251 Dictyoglomus turgidum                                                        Dictyoglomaceae
    ## NC_011663 Shewanella baltica                                                           Shewanellaceae
    ## CP000916 Thermotoga neapolitana                                                        Thermotogaceae
    ## NZ_DS996397 Desulfovibrio piger                                                   Desulfovibrionaceae
    ## CP001230 Persephonella marina                                                     Hydrogenothermaceae
    ## CP001472 Acidobacterium capsulatum                                                  Acidobacteriaceae
    ## AP009153 Gemmatimonas aurantiaca                                                    Gemmatimonadaceae
    ## CP001941 Aciduliprofundum boonei                                                                     
    ## NC_013968 Haloferax volcanii                                                           Haloferacaceae
    ## NZ_KE136524 Enterococcus faecalis                                                     Enterococcaceae
    ## NZ_KQ961402 Zymomonas mobilis                                                       Sphingomonadaceae
    ## NZ_CP015081 Deinococcus radiodurans                                                    Deinococcaceae
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense                               Hydrogenothermaceae
    ## NZ_JGWU01000001 Bordetella bronchiseptica                                              Alcaligenaceae
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii         Thermoanaerobacterales Family III. Incertae Sedis
    ## NC_009972 Herpetosiphon aurantiacus                                                Herpetosiphonaceae
    ## NC_005213 Nanoarchaeum equitans                                                       Nanoarchaeaceae
    ##                                                                    genus
    ## AE000782 Archaeoglobus fulgidus                            Archaeoglobus
    ## NC_000909 Methanocaldococcus jannaschii               Methanocaldococcus
    ## NC_003272 Nostoc sp. PCC 7120                                     Nostoc
    ## AE009441 Pyrobaculum aerophilum                              Pyrobaculum
    ## AE009950 Pyrococcus furiosus                                  Pyrococcus
    ## AE009951 Fusobacterium nucleatum                           Fusobacterium
    ## AE010299 Methanosarcina acetivorans                       Methanosarcina
    ## AE009439 Methanopyrus kandleri                              Methanopyrus
    ## NC_003911 Ruegeria pomeroyi                                     Ruegeria
    ## AE006470 Chlorobaculum tepidum                             Chlorobaculum
    ## AE015928 Bacteroides thetaiotaomicron                        Bacteroides
    ## AL954747 Nitrosomonas europaea                              Nitrosomonas
    ## BX119912 Rhodopirellula baltica                           Rhodopirellula
    ## BX571656 Wolinella succinogenes                                Wolinella
    ## AE017180 Geobacter sulfurreducens                              Geobacter
    ## AE017226 Treponema denticola                                   Treponema
    ## BX950229 Methanococcus maripaludis                         Methanococcus
    ## AE017221 Thermus thermophilus                                    Thermus
    ## BA000001 Pyrococcus horikoshii                                Pyrococcus
    ## BA000023 Sulfolobus tokodaii                                  Sulfolobus
    ## NC_007951 Paraburkholderia xenovorans                   Paraburkholderia
    ## CP000492 Chlorobium phaeobacteroides                          Chlorobium
    ## NC_008751 Desulfovibrio vulgaris                           Desulfovibrio
    ## CP000568 Ruminiclostridium thermocellum                Ruminiclostridium
    ## CP000561 Pyrobaculum calidifontis                            Pyrobaculum
    ## CP000609 Methanococcus maripaludis                         Methanococcus
    ## CP000607 Chlorobium phaeovibrioides                           Chlorobium
    ## CP000660 Pyrobaculum arsenaticum                             Pyrobaculum
    ## CP000667 Salinispora tropica                                 Salinispora
    ## CP000679 Caldicellulosiruptor saccharolyticus       Caldicellulosiruptor
    ## CP000702 Thermotoga petrophila                                Thermotoga
    ## CP000139 Bacteroides vulgatus                                Bacteroides
    ## NC_009665 Shewanella baltica                                  Shewanella
    ## CP000816 Ignicoccus hospitalis                                Ignicoccus
    ## CP000850 Salinispora arenicola                               Salinispora
    ## CP000909 Chloroflexus aurantiacus                           Chloroflexus
    ## CP000924 Thermoanaerobacter pseudethanolicus          Thermoanaerobacter
    ## CP000969 Thermotoga sp. RQ2                                   Thermotoga
    ## CP001013 Leptothrix cholodnii                                 Leptothrix
    ## CP001071 Akkermansia muciniphila                             Akkermansia
    ## AP009380 Porphyromonas gingivalis                          Porphyromonas
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1          Sulfurihydrogenibium
    ## CP001097 Chlorobium limicola                                  Chlorobium
    ## CP001110 Pelodictyon phaeoclathratiforme                     Pelodictyon
    ## CP001130 Hydrogenobaculum sp. Y04AAS1                   Hydrogenobaculum
    ## NZ_CH959311 Sulfitobacter sp. EE-36                        Sulfitobacter
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1                     Sulfitobacter
    ## CP001251 Dictyoglomus turgidum                              Dictyoglomus
    ## NC_011663 Shewanella baltica                                  Shewanella
    ## CP000916 Thermotoga neapolitana                               Thermotoga
    ## NZ_DS996397 Desulfovibrio piger                            Desulfovibrio
    ## CP001230 Persephonella marina                              Persephonella
    ## CP001472 Acidobacterium capsulatum                        Acidobacterium
    ## AP009153 Gemmatimonas aurantiaca                            Gemmatimonas
    ## CP001941 Aciduliprofundum boonei                        Aciduliprofundum
    ## NC_013968 Haloferax volcanii                                   Haloferax
    ## NZ_KE136524 Enterococcus faecalis                           Enterococcus
    ## NZ_KQ961402 Zymomonas mobilis                                  Zymomonas
    ## NZ_CP015081 Deinococcus radiodurans                          Deinococcus
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense Sulfurihydrogenibium
    ## NZ_JGWU01000001 Bordetella bronchiseptica                     Bordetella
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii         Caldicellulosiruptor
    ## NC_009972 Herpetosiphon aurantiacus                        Herpetosiphon
    ## NC_005213 Nanoarchaeum equitans                             Nanoarchaeum
    ##                                                                                  species
    ## AE000782 Archaeoglobus fulgidus                                   Archaeoglobus fulgidus
    ## NC_000909 Methanocaldococcus jannaschii                    Methanocaldococcus jannaschii
    ## NC_003272 Nostoc sp. PCC 7120                                        Nostoc sp. PCC 7120
    ## AE009441 Pyrobaculum aerophilum                                   Pyrobaculum aerophilum
    ## AE009950 Pyrococcus furiosus                                         Pyrococcus furiosus
    ## AE009951 Fusobacterium nucleatum                                 Fusobacterium nucleatum
    ## AE010299 Methanosarcina acetivorans                           Methanosarcina acetivorans
    ## AE009439 Methanopyrus kandleri                                     Methanopyrus kandleri
    ## NC_003911 Ruegeria pomeroyi                                            Ruegeria pomeroyi
    ## AE006470 Chlorobaculum tepidum                                     Chlorobaculum tepidum
    ## AE015928 Bacteroides thetaiotaomicron                       Bacteroides thetaiotaomicron
    ## AL954747 Nitrosomonas europaea                                     Nitrosomonas europaea
    ## BX119912 Rhodopirellula baltica                                   Rhodopirellula baltica
    ## BX571656 Wolinella succinogenes                                   Wolinella succinogenes
    ## AE017180 Geobacter sulfurreducens                               Geobacter sulfurreducens
    ## AE017226 Treponema denticola                                         Treponema denticola
    ## BX950229 Methanococcus maripaludis                             Methanococcus maripaludis
    ## AE017221 Thermus thermophilus                                       Thermus thermophilus
    ## BA000001 Pyrococcus horikoshii                                     Pyrococcus horikoshii
    ## BA000023 Sulfolobus tokodaii                                         Sulfolobus tokodaii
    ## NC_007951 Paraburkholderia xenovorans                        Paraburkholderia xenovorans
    ## CP000492 Chlorobium phaeobacteroides                         Chlorobium phaeobacteroides
    ## NC_008751 Desulfovibrio vulgaris                                  Desulfovibrio vulgaris
    ## CP000568 Ruminiclostridium thermocellum                   Ruminiclostridium thermocellum
    ## CP000561 Pyrobaculum calidifontis                               Pyrobaculum calidifontis
    ## CP000609 Methanococcus maripaludis                             Methanococcus maripaludis
    ## CP000607 Chlorobium phaeovibrioides                           Chlorobium phaeovibrioides
    ## CP000660 Pyrobaculum arsenaticum                                 Pyrobaculum arsenaticum
    ## CP000667 Salinispora tropica                                         Salinispora tropica
    ## CP000679 Caldicellulosiruptor saccharolyticus       Caldicellulosiruptor saccharolyticus
    ## CP000702 Thermotoga petrophila                                     Thermotoga petrophila
    ## CP000139 Bacteroides vulgatus                                       Bacteroides vulgatus
    ## NC_009665 Shewanella baltica                                          Shewanella baltica
    ## CP000816 Ignicoccus hospitalis                                     Ignicoccus hospitalis
    ## CP000850 Salinispora arenicola                                     Salinispora arenicola
    ## CP000909 Chloroflexus aurantiacus                               Chloroflexus aurantiacus
    ## CP000924 Thermoanaerobacter pseudethanolicus         Thermoanaerobacter pseudethanolicus
    ## CP000969 Thermotoga sp. RQ2                                           Thermotoga sp. RQ2
    ## CP001013 Leptothrix cholodnii                                       Leptothrix cholodnii
    ## CP001071 Akkermansia muciniphila                                 Akkermansia muciniphila
    ## AP009380 Porphyromonas gingivalis                               Porphyromonas gingivalis
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1              Sulfurihydrogenibium sp. YO3AOP1
    ## CP001097 Chlorobium limicola                                         Chlorobium limicola
    ## CP001110 Pelodictyon phaeoclathratiforme                 Pelodictyon phaeoclathratiforme
    ## CP001130 Hydrogenobaculum sp. Y04AAS1                       Hydrogenobaculum sp. Y04AAS1
    ## NZ_CH959311 Sulfitobacter sp. EE-36                              Sulfitobacter sp. EE-36
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1                        Sulfitobacter sp. NAS-14.1
    ## CP001251 Dictyoglomus turgidum                                     Dictyoglomus turgidum
    ## NC_011663 Shewanella baltica                                          Shewanella baltica
    ## CP000916 Thermotoga neapolitana                                   Thermotoga neapolitana
    ## NZ_DS996397 Desulfovibrio piger                                      Desulfovibrio piger
    ## CP001230 Persephonella marina                                       Persephonella marina
    ## CP001472 Acidobacterium capsulatum                             Acidobacterium capsulatum
    ## AP009153 Gemmatimonas aurantiaca                                 Gemmatimonas aurantiaca
    ## CP001941 Aciduliprofundum boonei                                 Aciduliprofundum boonei
    ## NC_013968 Haloferax volcanii                                          Haloferax volcanii
    ## NZ_KE136524 Enterococcus faecalis                                  Enterococcus faecalis
    ## NZ_KQ961402 Zymomonas mobilis                                          Zymomonas mobilis
    ## NZ_CP015081 Deinococcus radiodurans                              Deinococcus radiodurans
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense  Sulfurihydrogenibium yellowstonense
    ## NZ_JGWU01000001 Bordetella bronchiseptica                      Bordetella bronchiseptica
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii                  Caldicellulosiruptor bescii
    ## NC_009972 Herpetosiphon aurantiacus                            Herpetosiphon aurantiacus
    ## NC_005213 Nanoarchaeum equitans                                    Nanoarchaeum equitans
    ##                                                                                             strain
    ## AE000782 Archaeoglobus fulgidus                                    Archaeoglobus fulgidus DSM 4304
    ## NC_000909 Methanocaldococcus jannaschii                     Methanocaldococcus jannaschii DSM 2661
    ## NC_003272 Nostoc sp. PCC 7120                                                                     
    ## AE009441 Pyrobaculum aerophilum                                    Pyrobaculum aerophilum str. IM2
    ## AE009950 Pyrococcus furiosus                                          Pyrococcus furiosus DSM 3638
    ## AE009951 Fusobacterium nucleatum                                                                  
    ## AE010299 Methanosarcina acetivorans                                 Methanosarcina acetivorans C2A
    ## AE009439 Methanopyrus kandleri                                          Methanopyrus kandleri AV19
    ## NC_003911 Ruegeria pomeroyi                                                Ruegeria pomeroyi DSS-3
    ## AE006470 Chlorobaculum tepidum                                           Chlorobaculum tepidum TLS
    ## AE015928 Bacteroides thetaiotaomicron                        Bacteroides thetaiotaomicron VPI-5482
    ## AL954747 Nitrosomonas europaea                                    Nitrosomonas europaea ATCC 19718
    ## BX119912 Rhodopirellula baltica                                        Rhodopirellula baltica SH 1
    ## BX571656 Wolinella succinogenes                                    Wolinella succinogenes DSM 1740
    ## AE017180 Geobacter sulfurreducens                                     Geobacter sulfurreducens PCA
    ## AE017226 Treponema denticola                                        Treponema denticola ATCC 35405
    ## BX950229 Methanococcus maripaludis                                    Methanococcus maripaludis S2
    ## AE017221 Thermus thermophilus                                            Thermus thermophilus HB27
    ## BA000001 Pyrococcus horikoshii                                           Pyrococcus horikoshii OT3
    ## BA000023 Sulfolobus tokodaii                                            Sulfolobus tokodaii str. 7
    ## NC_007951 Paraburkholderia xenovorans                            Paraburkholderia xenovorans LB400
    ## CP000492 Chlorobium phaeobacteroides                           Chlorobium phaeobacteroides DSM 266
    ## NC_008751 Desulfovibrio vulgaris                                        Desulfovibrio vulgaris DP4
    ## CP000568 Ruminiclostridium thermocellum                  Ruminiclostridium thermocellum ATCC 27405
    ## CP000561 Pyrobaculum calidifontis                               Pyrobaculum calidifontis JCM 11548
    ## CP000609 Methanococcus maripaludis                                    Methanococcus maripaludis C5
    ## CP000607 Chlorobium phaeovibrioides                             Chlorobium phaeovibrioides DSM 265
    ## CP000660 Pyrobaculum arsenaticum                                 Pyrobaculum arsenaticum DSM 13514
    ## CP000667 Salinispora tropica                                           Salinispora tropica CNB-440
    ## CP000679 Caldicellulosiruptor saccharolyticus        Caldicellulosiruptor saccharolyticus DSM 8903
    ## CP000702 Thermotoga petrophila                                         Thermotoga petrophila RKU-1
    ## CP000139 Bacteroides vulgatus                                       Bacteroides vulgatus ATCC 8482
    ## NC_009665 Shewanella baltica                                              Shewanella baltica OS185
    ## CP000816 Ignicoccus hospitalis                                        Ignicoccus hospitalis KIN4/I
    ## CP000850 Salinispora arenicola                                       Salinispora arenicola CNS-205
    ## CP000909 Chloroflexus aurantiacus                                 Chloroflexus aurantiacus J-10-fl
    ## CP000924 Thermoanaerobacter pseudethanolicus        Thermoanaerobacter pseudethanolicus ATCC 33223
    ## CP000969 Thermotoga sp. RQ2                                                                       
    ## CP001013 Leptothrix cholodnii                                            Leptothrix cholodnii SP-6
    ## CP001071 Akkermansia muciniphila                              Akkermansia muciniphila ATCC BAA-835
    ## AP009380 Porphyromonas gingivalis                              Porphyromonas gingivalis ATCC 33277
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1                                                        
    ## CP001097 Chlorobium limicola                                           Chlorobium limicola DSM 245
    ## CP001110 Pelodictyon phaeoclathratiforme                      Pelodictyon phaeoclathratiforme BU-1
    ## CP001130 Hydrogenobaculum sp. Y04AAS1                                                             
    ## NZ_CH959311 Sulfitobacter sp. EE-36                                                               
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1                                                            
    ## CP001251 Dictyoglomus turgidum                                      Dictyoglomus turgidum DSM 6724
    ## NC_011663 Shewanella baltica                                              Shewanella baltica OS223
    ## CP000916 Thermotoga neapolitana                                    Thermotoga neapolitana DSM 4359
    ## NZ_DS996397 Desulfovibrio piger                                     Desulfovibrio piger ATCC 29098
    ## CP001230 Persephonella marina                                           Persephonella marina EX-H1
    ## CP001472 Acidobacterium capsulatum                            Acidobacterium capsulatum ATCC 51196
    ## AP009153 Gemmatimonas aurantiaca                                      Gemmatimonas aurantiaca T-27
    ## CP001941 Aciduliprofundum boonei                                      Aciduliprofundum boonei T469
    ## NC_013968 Haloferax volcanii                                                Haloferax volcanii DS2
    ## NZ_KE136524 Enterococcus faecalis                                       Enterococcus faecalis V583
    ## NZ_KQ961402 Zymomonas mobilis                                                                     
    ## NZ_CP015081 Deinococcus radiodurans                                     Deinococcus radiodurans R1
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense       Sulfurihydrogenibium yellowstonense SS-5
    ## NZ_JGWU01000001 Bordetella bronchiseptica                           Bordetella bronchiseptica D989
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii                                                       
    ## NC_009972 Herpetosiphon aurantiacus                              Herpetosiphon aurantiacus DSM 785
    ## NC_005213 Nanoarchaeum equitans                                       Nanoarchaeum equitans Kin4-M

``` r
OTU = otu_table(df, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax))

OTU
```

    ## OTU Table:          [17 taxa and 5 samples]
    ##                      taxa are rows
    ##                                             ERR1136663   ERR4087400
    ## CP000139 Bacteroides vulgatus             8.301477e-03 5.512089e-03
    ## AE015928 Bacteroides thetaiotaomicron     2.032678e-03 3.968704e-04
    ## NZ_DS996397 Desulfovibrio piger           1.076401e-03 0.000000e+00
    ## CP001071 Akkermansia muciniphila          5.676425e-04 4.119893e-03
    ## NZ_KE136524 Enterococcus faecalis         2.355363e-05 6.929483e-05
    ## CP001472 Acidobacterium capsulatum        1.165905e-04 0.000000e+00
    ## AE009951 Fusobacterium nucleatum          2.814659e-04 0.000000e+00
    ## CP001013 Leptothrix cholodnii             9.421452e-06 0.000000e+00
    ## NZ_JGWU01000001 Bordetella bronchiseptica 0.000000e+00 3.779718e-05
    ## NZ_KQ961402 Zymomonas mobilis             0.000000e+00 0.000000e+00
    ## AP009380 Porphyromonas gingivalis         0.000000e+00 0.000000e+00
    ## NC_007951 Burkholderia xenovorans         0.000000e+00 0.000000e+00
    ## AE017226 Treponema denticola              0.000000e+00 0.000000e+00
    ## CP000850 Salinispora arenicola            0.000000e+00 0.000000e+00
    ## AE006470 Chlorobium tepidum               0.000000e+00 0.000000e+00
    ## AE017221 Thermus thermophilus             0.000000e+00 0.000000e+00
    ## NC_011663 Shewanella baltica              0.000000e+00 0.000000e+00
    ##                                             ERR4562131    SRR059888
    ## CP000139 Bacteroides vulgatus             5.823024e-03 3.356448e-05
    ## AE015928 Bacteroides thetaiotaomicron     4.229973e-04 0.000000e+00
    ## NZ_DS996397 Desulfovibrio piger           1.548600e-04 0.000000e+00
    ## CP001071 Akkermansia muciniphila          4.301667e-05 9.811156e-06
    ## NZ_KE136524 Enterococcus faecalis         0.000000e+00 0.000000e+00
    ## CP001472 Acidobacterium capsulatum        0.000000e+00 0.000000e+00
    ## AE009951 Fusobacterium nucleatum          0.000000e+00 3.217285e-03
    ## CP001013 Leptothrix cholodnii             0.000000e+00 1.587858e-04
    ## NZ_JGWU01000001 Bordetella bronchiseptica 0.000000e+00 2.793598e-04
    ## NZ_KQ961402 Zymomonas mobilis             3.728112e-04 0.000000e+00
    ## AP009380 Porphyromonas gingivalis         5.735556e-05 2.765197e-04
    ## NC_007951 Burkholderia xenovorans         2.581000e-05 0.000000e+00
    ## AE017226 Treponema denticola              0.000000e+00 7.255092e-05
    ## CP000850 Salinispora arenicola            0.000000e+00 3.464887e-04
    ## AE006470 Chlorobium tepidum               0.000000e+00 7.539099e-05
    ## AE017221 Thermus thermophilus             0.000000e+00 3.356448e-06
    ## NC_011663 Shewanella baltica              0.000000e+00 0.000000e+00
    ##                                             SRR9217435
    ## CP000139 Bacteroides vulgatus             5.792533e-05
    ## AE015928 Bacteroides thetaiotaomicron     0.000000e+00
    ## NZ_DS996397 Desulfovibrio piger           0.000000e+00
    ## CP001071 Akkermansia muciniphila          0.000000e+00
    ## NZ_KE136524 Enterococcus faecalis         0.000000e+00
    ## CP001472 Acidobacterium capsulatum        0.000000e+00
    ## AE009951 Fusobacterium nucleatum          5.649507e-05
    ## CP001013 Leptothrix cholodnii             1.923693e-04
    ## NZ_JGWU01000001 Bordetella bronchiseptica 3.983260e-04
    ## NZ_KQ961402 Zymomonas mobilis             2.932023e-05
    ## AP009380 Porphyromonas gingivalis         8.581530e-06
    ## NC_007951 Burkholderia xenovorans         4.290765e-06
    ## AE017226 Treponema denticola              0.000000e+00
    ## CP000850 Salinispora arenicola            0.000000e+00
    ## AE006470 Chlorobium tepidum               0.000000e+00
    ## AE017221 Thermus thermophilus             0.000000e+00
    ## NC_011663 Shewanella baltica              1.623339e-04

``` r
TAX
```

    ## Taxonomy Table:     [64 taxa by 8 taxonomic ranks]:
    ##                                                     superkingdom
    ## AE000782 Archaeoglobus fulgidus                     "Archaea"   
    ## NC_000909 Methanocaldococcus jannaschii             "Archaea"   
    ## NC_003272 Nostoc sp. PCC 7120                       "Bacteria"  
    ## AE009441 Pyrobaculum aerophilum                     "Archaea"   
    ## AE009950 Pyrococcus furiosus                        "Archaea"   
    ## AE009951 Fusobacterium nucleatum                    "Bacteria"  
    ## AE010299 Methanosarcina acetivorans                 "Archaea"   
    ## AE009439 Methanopyrus kandleri                      "Archaea"   
    ## NC_003911 Ruegeria pomeroyi                         "Bacteria"  
    ## AE006470 Chlorobaculum tepidum                      "Bacteria"  
    ## AE015928 Bacteroides thetaiotaomicron               "Bacteria"  
    ## AL954747 Nitrosomonas europaea                      "Bacteria"  
    ## BX119912 Rhodopirellula baltica                     "Bacteria"  
    ## BX571656 Wolinella succinogenes                     "Bacteria"  
    ## AE017180 Geobacter sulfurreducens                   "Bacteria"  
    ## AE017226 Treponema denticola                        "Bacteria"  
    ## BX950229 Methanococcus maripaludis                  "Archaea"   
    ## AE017221 Thermus thermophilus                       "Bacteria"  
    ## BA000001 Pyrococcus horikoshii                      "Archaea"   
    ## BA000023 Sulfolobus tokodaii                        "Archaea"   
    ## NC_007951 Paraburkholderia xenovorans               "Bacteria"  
    ## CP000492 Chlorobium phaeobacteroides                "Bacteria"  
    ## NC_008751 Desulfovibrio vulgaris                    "Bacteria"  
    ## CP000568 Ruminiclostridium thermocellum             "Bacteria"  
    ## CP000561 Pyrobaculum calidifontis                   "Archaea"   
    ## CP000609 Methanococcus maripaludis                  "Archaea"   
    ## CP000607 Chlorobium phaeovibrioides                 "Bacteria"  
    ## CP000660 Pyrobaculum arsenaticum                    "Archaea"   
    ## CP000667 Salinispora tropica                        "Bacteria"  
    ## CP000679 Caldicellulosiruptor saccharolyticus       "Bacteria"  
    ## CP000702 Thermotoga petrophila                      "Bacteria"  
    ## CP000139 Bacteroides vulgatus                       "Bacteria"  
    ## NC_009665 Shewanella baltica                        "Bacteria"  
    ## CP000816 Ignicoccus hospitalis                      "Archaea"   
    ## CP000850 Salinispora arenicola                      "Bacteria"  
    ## CP000909 Chloroflexus aurantiacus                   "Bacteria"  
    ## CP000924 Thermoanaerobacter pseudethanolicus        "Bacteria"  
    ## CP000969 Thermotoga sp. RQ2                         "Bacteria"  
    ## CP001013 Leptothrix cholodnii                       "Bacteria"  
    ## CP001071 Akkermansia muciniphila                    "Bacteria"  
    ## AP009380 Porphyromonas gingivalis                   "Bacteria"  
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1          "Bacteria"  
    ## CP001097 Chlorobium limicola                        "Bacteria"  
    ## CP001110 Pelodictyon phaeoclathratiforme            "Bacteria"  
    ## CP001130 Hydrogenobaculum sp. Y04AAS1               "Bacteria"  
    ## NZ_CH959311 Sulfitobacter sp. EE-36                 "Bacteria"  
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1              "Bacteria"  
    ## CP001251 Dictyoglomus turgidum                      "Bacteria"  
    ## NC_011663 Shewanella baltica                        "Bacteria"  
    ## CP000916 Thermotoga neapolitana                     "Bacteria"  
    ## NZ_DS996397 Desulfovibrio piger                     "Bacteria"  
    ## CP001230 Persephonella marina                       "Bacteria"  
    ## CP001472 Acidobacterium capsulatum                  "Bacteria"  
    ## AP009153 Gemmatimonas aurantiaca                    "Bacteria"  
    ## CP001941 Aciduliprofundum boonei                    "Archaea"   
    ## NC_013968 Haloferax volcanii                        "Archaea"   
    ## NZ_KE136524 Enterococcus faecalis                   "Bacteria"  
    ## NZ_KQ961402 Zymomonas mobilis                       "Bacteria"  
    ## NZ_CP015081 Deinococcus radiodurans                 "Bacteria"  
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense "Bacteria"  
    ## NZ_JGWU01000001 Bordetella bronchiseptica           "Bacteria"  
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii         "Bacteria"  
    ## NC_009972 Herpetosiphon aurantiacus                 "Bacteria"  
    ## NC_005213 Nanoarchaeum equitans                     "Archaea"   
    ##                                                     phylum               
    ## AE000782 Archaeoglobus fulgidus                     "Euryarchaeota"      
    ## NC_000909 Methanocaldococcus jannaschii             "Euryarchaeota"      
    ## NC_003272 Nostoc sp. PCC 7120                       "Cyanobacteria"      
    ## AE009441 Pyrobaculum aerophilum                     "Crenarchaeota"      
    ## AE009950 Pyrococcus furiosus                        "Euryarchaeota"      
    ## AE009951 Fusobacterium nucleatum                    "Fusobacteria"       
    ## AE010299 Methanosarcina acetivorans                 "Euryarchaeota"      
    ## AE009439 Methanopyrus kandleri                      "Euryarchaeota"      
    ## NC_003911 Ruegeria pomeroyi                         "Proteobacteria"     
    ## AE006470 Chlorobaculum tepidum                      "Chlorobi"           
    ## AE015928 Bacteroides thetaiotaomicron               "Bacteroidetes"      
    ## AL954747 Nitrosomonas europaea                      "Proteobacteria"     
    ## BX119912 Rhodopirellula baltica                     "Planctomycetes"     
    ## BX571656 Wolinella succinogenes                     "Proteobacteria"     
    ## AE017180 Geobacter sulfurreducens                   "Proteobacteria"     
    ## AE017226 Treponema denticola                        "Spirochaetes"       
    ## BX950229 Methanococcus maripaludis                  "Euryarchaeota"      
    ## AE017221 Thermus thermophilus                       "Deinococcus-Thermus"
    ## BA000001 Pyrococcus horikoshii                      "Euryarchaeota"      
    ## BA000023 Sulfolobus tokodaii                        "Crenarchaeota"      
    ## NC_007951 Paraburkholderia xenovorans               "Proteobacteria"     
    ## CP000492 Chlorobium phaeobacteroides                "Chlorobi"           
    ## NC_008751 Desulfovibrio vulgaris                    "Proteobacteria"     
    ## CP000568 Ruminiclostridium thermocellum             "Firmicutes"         
    ## CP000561 Pyrobaculum calidifontis                   "Crenarchaeota"      
    ## CP000609 Methanococcus maripaludis                  "Euryarchaeota"      
    ## CP000607 Chlorobium phaeovibrioides                 "Chlorobi"           
    ## CP000660 Pyrobaculum arsenaticum                    "Crenarchaeota"      
    ## CP000667 Salinispora tropica                        "Actinobacteria"     
    ## CP000679 Caldicellulosiruptor saccharolyticus       "Firmicutes"         
    ## CP000702 Thermotoga petrophila                      "Thermotogae"        
    ## CP000139 Bacteroides vulgatus                       "Bacteroidetes"      
    ## NC_009665 Shewanella baltica                        "Proteobacteria"     
    ## CP000816 Ignicoccus hospitalis                      "Crenarchaeota"      
    ## CP000850 Salinispora arenicola                      "Actinobacteria"     
    ## CP000909 Chloroflexus aurantiacus                   "Chloroflexi"        
    ## CP000924 Thermoanaerobacter pseudethanolicus        "Firmicutes"         
    ## CP000969 Thermotoga sp. RQ2                         "Thermotogae"        
    ## CP001013 Leptothrix cholodnii                       "Proteobacteria"     
    ## CP001071 Akkermansia muciniphila                    "Verrucomicrobia"    
    ## AP009380 Porphyromonas gingivalis                   "Bacteroidetes"      
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1          "Aquificae"          
    ## CP001097 Chlorobium limicola                        "Chlorobi"           
    ## CP001110 Pelodictyon phaeoclathratiforme            "Chlorobi"           
    ## CP001130 Hydrogenobaculum sp. Y04AAS1               "Aquificae"          
    ## NZ_CH959311 Sulfitobacter sp. EE-36                 "Proteobacteria"     
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1              "Proteobacteria"     
    ## CP001251 Dictyoglomus turgidum                      "Dictyoglomi"        
    ## NC_011663 Shewanella baltica                        "Proteobacteria"     
    ## CP000916 Thermotoga neapolitana                     "Thermotogae"        
    ## NZ_DS996397 Desulfovibrio piger                     "Proteobacteria"     
    ## CP001230 Persephonella marina                       "Aquificae"          
    ## CP001472 Acidobacterium capsulatum                  "Acidobacteria"      
    ## AP009153 Gemmatimonas aurantiaca                    "Gemmatimonadetes"   
    ## CP001941 Aciduliprofundum boonei                    "Euryarchaeota"      
    ## NC_013968 Haloferax volcanii                        "Euryarchaeota"      
    ## NZ_KE136524 Enterococcus faecalis                   "Firmicutes"         
    ## NZ_KQ961402 Zymomonas mobilis                       "Proteobacteria"     
    ## NZ_CP015081 Deinococcus radiodurans                 "Deinococcus-Thermus"
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense "Aquificae"          
    ## NZ_JGWU01000001 Bordetella bronchiseptica           "Proteobacteria"     
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii         "Firmicutes"         
    ## NC_009972 Herpetosiphon aurantiacus                 "Chloroflexi"        
    ## NC_005213 Nanoarchaeum equitans                     "Nanoarchaeota"      
    ##                                                     class                  
    ## AE000782 Archaeoglobus fulgidus                     "Archaeoglobi"         
    ## NC_000909 Methanocaldococcus jannaschii             "Methanococci"         
    ## NC_003272 Nostoc sp. PCC 7120                       ""                     
    ## AE009441 Pyrobaculum aerophilum                     "Thermoprotei"         
    ## AE009950 Pyrococcus furiosus                        "Thermococci"          
    ## AE009951 Fusobacterium nucleatum                    "Fusobacteriia"        
    ## AE010299 Methanosarcina acetivorans                 "Methanomicrobia"      
    ## AE009439 Methanopyrus kandleri                      "Methanopyri"          
    ## NC_003911 Ruegeria pomeroyi                         "Alphaproteobacteria"  
    ## AE006470 Chlorobaculum tepidum                      "Chlorobia"            
    ## AE015928 Bacteroides thetaiotaomicron               "Bacteroidia"          
    ## AL954747 Nitrosomonas europaea                      "Betaproteobacteria"   
    ## BX119912 Rhodopirellula baltica                     "Planctomycetia"       
    ## BX571656 Wolinella succinogenes                     "Epsilonproteobacteria"
    ## AE017180 Geobacter sulfurreducens                   "Deltaproteobacteria"  
    ## AE017226 Treponema denticola                        "Spirochaetia"         
    ## BX950229 Methanococcus maripaludis                  "Methanococci"         
    ## AE017221 Thermus thermophilus                       "Deinococci"           
    ## BA000001 Pyrococcus horikoshii                      "Thermococci"          
    ## BA000023 Sulfolobus tokodaii                        "Thermoprotei"         
    ## NC_007951 Paraburkholderia xenovorans               "Betaproteobacteria"   
    ## CP000492 Chlorobium phaeobacteroides                "Chlorobia"            
    ## NC_008751 Desulfovibrio vulgaris                    "Deltaproteobacteria"  
    ## CP000568 Ruminiclostridium thermocellum             "Clostridia"           
    ## CP000561 Pyrobaculum calidifontis                   "Thermoprotei"         
    ## CP000609 Methanococcus maripaludis                  "Methanococci"         
    ## CP000607 Chlorobium phaeovibrioides                 "Chlorobia"            
    ## CP000660 Pyrobaculum arsenaticum                    "Thermoprotei"         
    ## CP000667 Salinispora tropica                        "Actinobacteria"       
    ## CP000679 Caldicellulosiruptor saccharolyticus       "Clostridia"           
    ## CP000702 Thermotoga petrophila                      "Thermotogae"          
    ## CP000139 Bacteroides vulgatus                       "Bacteroidia"          
    ## NC_009665 Shewanella baltica                        "Gammaproteobacteria"  
    ## CP000816 Ignicoccus hospitalis                      "Thermoprotei"         
    ## CP000850 Salinispora arenicola                      "Actinobacteria"       
    ## CP000909 Chloroflexus aurantiacus                   "Chloroflexia"         
    ## CP000924 Thermoanaerobacter pseudethanolicus        "Clostridia"           
    ## CP000969 Thermotoga sp. RQ2                         "Thermotogae"          
    ## CP001013 Leptothrix cholodnii                       "Betaproteobacteria"   
    ## CP001071 Akkermansia muciniphila                    "Verrucomicrobiae"     
    ## AP009380 Porphyromonas gingivalis                   "Bacteroidia"          
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1          "Aquificae"            
    ## CP001097 Chlorobium limicola                        "Chlorobia"            
    ## CP001110 Pelodictyon phaeoclathratiforme            "Chlorobia"            
    ## CP001130 Hydrogenobaculum sp. Y04AAS1               "Aquificae"            
    ## NZ_CH959311 Sulfitobacter sp. EE-36                 "Alphaproteobacteria"  
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1              "Alphaproteobacteria"  
    ## CP001251 Dictyoglomus turgidum                      "Dictyoglomia"         
    ## NC_011663 Shewanella baltica                        "Gammaproteobacteria"  
    ## CP000916 Thermotoga neapolitana                     "Thermotogae"          
    ## NZ_DS996397 Desulfovibrio piger                     "Deltaproteobacteria"  
    ## CP001230 Persephonella marina                       "Aquificae"            
    ## CP001472 Acidobacterium capsulatum                  "Acidobacteriia"       
    ## AP009153 Gemmatimonas aurantiaca                    "Gemmatimonadetes"     
    ## CP001941 Aciduliprofundum boonei                    ""                     
    ## NC_013968 Haloferax volcanii                        "Halobacteria"         
    ## NZ_KE136524 Enterococcus faecalis                   "Bacilli"              
    ## NZ_KQ961402 Zymomonas mobilis                       "Alphaproteobacteria"  
    ## NZ_CP015081 Deinococcus radiodurans                 "Deinococci"           
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense "Aquificae"            
    ## NZ_JGWU01000001 Bordetella bronchiseptica           "Betaproteobacteria"   
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii         "Clostridia"           
    ## NC_009972 Herpetosiphon aurantiacus                 "Chloroflexia"         
    ## NC_005213 Nanoarchaeum equitans                     ""                     
    ##                                                     order                   
    ## AE000782 Archaeoglobus fulgidus                     "Archaeoglobales"       
    ## NC_000909 Methanocaldococcus jannaschii             "Methanococcales"       
    ## NC_003272 Nostoc sp. PCC 7120                       "Nostocales"            
    ## AE009441 Pyrobaculum aerophilum                     "Thermoproteales"       
    ## AE009950 Pyrococcus furiosus                        "Thermococcales"        
    ## AE009951 Fusobacterium nucleatum                    "Fusobacteriales"       
    ## AE010299 Methanosarcina acetivorans                 "Methanosarcinales"     
    ## AE009439 Methanopyrus kandleri                      "Methanopyrales"        
    ## NC_003911 Ruegeria pomeroyi                         "Rhodobacterales"       
    ## AE006470 Chlorobaculum tepidum                      "Chlorobiales"          
    ## AE015928 Bacteroides thetaiotaomicron               "Bacteroidales"         
    ## AL954747 Nitrosomonas europaea                      "Nitrosomonadales"      
    ## BX119912 Rhodopirellula baltica                     "Planctomycetales"      
    ## BX571656 Wolinella succinogenes                     "Campylobacterales"     
    ## AE017180 Geobacter sulfurreducens                   "Desulfuromonadales"    
    ## AE017226 Treponema denticola                        "Spirochaetales"        
    ## BX950229 Methanococcus maripaludis                  "Methanococcales"       
    ## AE017221 Thermus thermophilus                       "Thermales"             
    ## BA000001 Pyrococcus horikoshii                      "Thermococcales"        
    ## BA000023 Sulfolobus tokodaii                        "Sulfolobales"          
    ## NC_007951 Paraburkholderia xenovorans               "Burkholderiales"       
    ## CP000492 Chlorobium phaeobacteroides                "Chlorobiales"          
    ## NC_008751 Desulfovibrio vulgaris                    "Desulfovibrionales"    
    ## CP000568 Ruminiclostridium thermocellum             "Clostridiales"         
    ## CP000561 Pyrobaculum calidifontis                   "Thermoproteales"       
    ## CP000609 Methanococcus maripaludis                  "Methanococcales"       
    ## CP000607 Chlorobium phaeovibrioides                 "Chlorobiales"          
    ## CP000660 Pyrobaculum arsenaticum                    "Thermoproteales"       
    ## CP000667 Salinispora tropica                        "Micromonosporales"     
    ## CP000679 Caldicellulosiruptor saccharolyticus       "Thermoanaerobacterales"
    ## CP000702 Thermotoga petrophila                      "Thermotogales"         
    ## CP000139 Bacteroides vulgatus                       "Bacteroidales"         
    ## NC_009665 Shewanella baltica                        "Alteromonadales"       
    ## CP000816 Ignicoccus hospitalis                      "Desulfurococcales"     
    ## CP000850 Salinispora arenicola                      "Micromonosporales"     
    ## CP000909 Chloroflexus aurantiacus                   "Chloroflexales"        
    ## CP000924 Thermoanaerobacter pseudethanolicus        "Thermoanaerobacterales"
    ## CP000969 Thermotoga sp. RQ2                         "Thermotogales"         
    ## CP001013 Leptothrix cholodnii                       "Burkholderiales"       
    ## CP001071 Akkermansia muciniphila                    "Verrucomicrobiales"    
    ## AP009380 Porphyromonas gingivalis                   "Bacteroidales"         
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1          "Aquificales"           
    ## CP001097 Chlorobium limicola                        "Chlorobiales"          
    ## CP001110 Pelodictyon phaeoclathratiforme            "Chlorobiales"          
    ## CP001130 Hydrogenobaculum sp. Y04AAS1               "Aquificales"           
    ## NZ_CH959311 Sulfitobacter sp. EE-36                 "Rhodobacterales"       
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1              "Rhodobacterales"       
    ## CP001251 Dictyoglomus turgidum                      "Dictyoglomales"        
    ## NC_011663 Shewanella baltica                        "Alteromonadales"       
    ## CP000916 Thermotoga neapolitana                     "Thermotogales"         
    ## NZ_DS996397 Desulfovibrio piger                     "Desulfovibrionales"    
    ## CP001230 Persephonella marina                       "Aquificales"           
    ## CP001472 Acidobacterium capsulatum                  "Acidobacteriales"      
    ## AP009153 Gemmatimonas aurantiaca                    "Gemmatimonadales"      
    ## CP001941 Aciduliprofundum boonei                    ""                      
    ## NC_013968 Haloferax volcanii                        "Haloferacales"         
    ## NZ_KE136524 Enterococcus faecalis                   "Lactobacillales"       
    ## NZ_KQ961402 Zymomonas mobilis                       "Sphingomonadales"      
    ## NZ_CP015081 Deinococcus radiodurans                 "Deinococcales"         
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense "Aquificales"           
    ## NZ_JGWU01000001 Bordetella bronchiseptica           "Burkholderiales"       
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii         "Thermoanaerobacterales"
    ## NC_009972 Herpetosiphon aurantiacus                 "Herpetosiphonales"     
    ## NC_005213 Nanoarchaeum equitans                     "Nanoarchaeales"        
    ##                                                     family                                             
    ## AE000782 Archaeoglobus fulgidus                     "Archaeoglobaceae"                                 
    ## NC_000909 Methanocaldococcus jannaschii             "Methanocaldococcaceae"                            
    ## NC_003272 Nostoc sp. PCC 7120                       "Nostocaceae"                                      
    ## AE009441 Pyrobaculum aerophilum                     "Thermoproteaceae"                                 
    ## AE009950 Pyrococcus furiosus                        "Thermococcaceae"                                  
    ## AE009951 Fusobacterium nucleatum                    "Fusobacteriaceae"                                 
    ## AE010299 Methanosarcina acetivorans                 "Methanosarcinaceae"                               
    ## AE009439 Methanopyrus kandleri                      "Methanopyraceae"                                  
    ## NC_003911 Ruegeria pomeroyi                         "Rhodobacteraceae"                                 
    ## AE006470 Chlorobaculum tepidum                      "Chlorobiaceae"                                    
    ## AE015928 Bacteroides thetaiotaomicron               "Bacteroidaceae"                                   
    ## AL954747 Nitrosomonas europaea                      "Nitrosomonadaceae"                                
    ## BX119912 Rhodopirellula baltica                     "Planctomycetaceae"                                
    ## BX571656 Wolinella succinogenes                     "Helicobacteraceae"                                
    ## AE017180 Geobacter sulfurreducens                   "Geobacteraceae"                                   
    ## AE017226 Treponema denticola                        "Spirochaetaceae"                                  
    ## BX950229 Methanococcus maripaludis                  "Methanococcaceae"                                 
    ## AE017221 Thermus thermophilus                       "Thermaceae"                                       
    ## BA000001 Pyrococcus horikoshii                      "Thermococcaceae"                                  
    ## BA000023 Sulfolobus tokodaii                        "Sulfolobaceae"                                    
    ## NC_007951 Paraburkholderia xenovorans               "Burkholderiaceae"                                 
    ## CP000492 Chlorobium phaeobacteroides                "Chlorobiaceae"                                    
    ## NC_008751 Desulfovibrio vulgaris                    "Desulfovibrionaceae"                              
    ## CP000568 Ruminiclostridium thermocellum             "Ruminococcaceae"                                  
    ## CP000561 Pyrobaculum calidifontis                   "Thermoproteaceae"                                 
    ## CP000609 Methanococcus maripaludis                  "Methanococcaceae"                                 
    ## CP000607 Chlorobium phaeovibrioides                 "Chlorobiaceae"                                    
    ## CP000660 Pyrobaculum arsenaticum                    "Thermoproteaceae"                                 
    ## CP000667 Salinispora tropica                        "Micromonosporaceae"                               
    ## CP000679 Caldicellulosiruptor saccharolyticus       "Thermoanaerobacterales Family III. Incertae Sedis"
    ## CP000702 Thermotoga petrophila                      "Thermotogaceae"                                   
    ## CP000139 Bacteroides vulgatus                       "Bacteroidaceae"                                   
    ## NC_009665 Shewanella baltica                        "Shewanellaceae"                                   
    ## CP000816 Ignicoccus hospitalis                      "Desulfurococcaceae"                               
    ## CP000850 Salinispora arenicola                      "Micromonosporaceae"                               
    ## CP000909 Chloroflexus aurantiacus                   "Chloroflexaceae"                                  
    ## CP000924 Thermoanaerobacter pseudethanolicus        "Thermoanaerobacteraceae"                          
    ## CP000969 Thermotoga sp. RQ2                         "Thermotogaceae"                                   
    ## CP001013 Leptothrix cholodnii                       ""                                                 
    ## CP001071 Akkermansia muciniphila                    "Akkermansiaceae"                                  
    ## AP009380 Porphyromonas gingivalis                   "Porphyromonadaceae"                               
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1          "Hydrogenothermaceae"                              
    ## CP001097 Chlorobium limicola                        "Chlorobiaceae"                                    
    ## CP001110 Pelodictyon phaeoclathratiforme            "Chlorobiaceae"                                    
    ## CP001130 Hydrogenobaculum sp. Y04AAS1               "Aquificaceae"                                     
    ## NZ_CH959311 Sulfitobacter sp. EE-36                 "Rhodobacteraceae"                                 
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1              "Rhodobacteraceae"                                 
    ## CP001251 Dictyoglomus turgidum                      "Dictyoglomaceae"                                  
    ## NC_011663 Shewanella baltica                        "Shewanellaceae"                                   
    ## CP000916 Thermotoga neapolitana                     "Thermotogaceae"                                   
    ## NZ_DS996397 Desulfovibrio piger                     "Desulfovibrionaceae"                              
    ## CP001230 Persephonella marina                       "Hydrogenothermaceae"                              
    ## CP001472 Acidobacterium capsulatum                  "Acidobacteriaceae"                                
    ## AP009153 Gemmatimonas aurantiaca                    "Gemmatimonadaceae"                                
    ## CP001941 Aciduliprofundum boonei                    ""                                                 
    ## NC_013968 Haloferax volcanii                        "Haloferacaceae"                                   
    ## NZ_KE136524 Enterococcus faecalis                   "Enterococcaceae"                                  
    ## NZ_KQ961402 Zymomonas mobilis                       "Sphingomonadaceae"                                
    ## NZ_CP015081 Deinococcus radiodurans                 "Deinococcaceae"                                   
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense "Hydrogenothermaceae"                              
    ## NZ_JGWU01000001 Bordetella bronchiseptica           "Alcaligenaceae"                                   
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii         "Thermoanaerobacterales Family III. Incertae Sedis"
    ## NC_009972 Herpetosiphon aurantiacus                 "Herpetosiphonaceae"                               
    ## NC_005213 Nanoarchaeum equitans                     "Nanoarchaeaceae"                                  
    ##                                                     genus                 
    ## AE000782 Archaeoglobus fulgidus                     "Archaeoglobus"       
    ## NC_000909 Methanocaldococcus jannaschii             "Methanocaldococcus"  
    ## NC_003272 Nostoc sp. PCC 7120                       "Nostoc"              
    ## AE009441 Pyrobaculum aerophilum                     "Pyrobaculum"         
    ## AE009950 Pyrococcus furiosus                        "Pyrococcus"          
    ## AE009951 Fusobacterium nucleatum                    "Fusobacterium"       
    ## AE010299 Methanosarcina acetivorans                 "Methanosarcina"      
    ## AE009439 Methanopyrus kandleri                      "Methanopyrus"        
    ## NC_003911 Ruegeria pomeroyi                         "Ruegeria"            
    ## AE006470 Chlorobaculum tepidum                      "Chlorobaculum"       
    ## AE015928 Bacteroides thetaiotaomicron               "Bacteroides"         
    ## AL954747 Nitrosomonas europaea                      "Nitrosomonas"        
    ## BX119912 Rhodopirellula baltica                     "Rhodopirellula"      
    ## BX571656 Wolinella succinogenes                     "Wolinella"           
    ## AE017180 Geobacter sulfurreducens                   "Geobacter"           
    ## AE017226 Treponema denticola                        "Treponema"           
    ## BX950229 Methanococcus maripaludis                  "Methanococcus"       
    ## AE017221 Thermus thermophilus                       "Thermus"             
    ## BA000001 Pyrococcus horikoshii                      "Pyrococcus"          
    ## BA000023 Sulfolobus tokodaii                        "Sulfolobus"          
    ## NC_007951 Paraburkholderia xenovorans               "Paraburkholderia"    
    ## CP000492 Chlorobium phaeobacteroides                "Chlorobium"          
    ## NC_008751 Desulfovibrio vulgaris                    "Desulfovibrio"       
    ## CP000568 Ruminiclostridium thermocellum             "Ruminiclostridium"   
    ## CP000561 Pyrobaculum calidifontis                   "Pyrobaculum"         
    ## CP000609 Methanococcus maripaludis                  "Methanococcus"       
    ## CP000607 Chlorobium phaeovibrioides                 "Chlorobium"          
    ## CP000660 Pyrobaculum arsenaticum                    "Pyrobaculum"         
    ## CP000667 Salinispora tropica                        "Salinispora"         
    ## CP000679 Caldicellulosiruptor saccharolyticus       "Caldicellulosiruptor"
    ## CP000702 Thermotoga petrophila                      "Thermotoga"          
    ## CP000139 Bacteroides vulgatus                       "Bacteroides"         
    ## NC_009665 Shewanella baltica                        "Shewanella"          
    ## CP000816 Ignicoccus hospitalis                      "Ignicoccus"          
    ## CP000850 Salinispora arenicola                      "Salinispora"         
    ## CP000909 Chloroflexus aurantiacus                   "Chloroflexus"        
    ## CP000924 Thermoanaerobacter pseudethanolicus        "Thermoanaerobacter"  
    ## CP000969 Thermotoga sp. RQ2                         "Thermotoga"          
    ## CP001013 Leptothrix cholodnii                       "Leptothrix"          
    ## CP001071 Akkermansia muciniphila                    "Akkermansia"         
    ## AP009380 Porphyromonas gingivalis                   "Porphyromonas"       
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1          "Sulfurihydrogenibium"
    ## CP001097 Chlorobium limicola                        "Chlorobium"          
    ## CP001110 Pelodictyon phaeoclathratiforme            "Pelodictyon"         
    ## CP001130 Hydrogenobaculum sp. Y04AAS1               "Hydrogenobaculum"    
    ## NZ_CH959311 Sulfitobacter sp. EE-36                 "Sulfitobacter"       
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1              "Sulfitobacter"       
    ## CP001251 Dictyoglomus turgidum                      "Dictyoglomus"        
    ## NC_011663 Shewanella baltica                        "Shewanella"          
    ## CP000916 Thermotoga neapolitana                     "Thermotoga"          
    ## NZ_DS996397 Desulfovibrio piger                     "Desulfovibrio"       
    ## CP001230 Persephonella marina                       "Persephonella"       
    ## CP001472 Acidobacterium capsulatum                  "Acidobacterium"      
    ## AP009153 Gemmatimonas aurantiaca                    "Gemmatimonas"        
    ## CP001941 Aciduliprofundum boonei                    "Aciduliprofundum"    
    ## NC_013968 Haloferax volcanii                        "Haloferax"           
    ## NZ_KE136524 Enterococcus faecalis                   "Enterococcus"        
    ## NZ_KQ961402 Zymomonas mobilis                       "Zymomonas"           
    ## NZ_CP015081 Deinococcus radiodurans                 "Deinococcus"         
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense "Sulfurihydrogenibium"
    ## NZ_JGWU01000001 Bordetella bronchiseptica           "Bordetella"          
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii         "Caldicellulosiruptor"
    ## NC_009972 Herpetosiphon aurantiacus                 "Herpetosiphon"       
    ## NC_005213 Nanoarchaeum equitans                     "Nanoarchaeum"        
    ##                                                     species                               
    ## AE000782 Archaeoglobus fulgidus                     "Archaeoglobus fulgidus"              
    ## NC_000909 Methanocaldococcus jannaschii             "Methanocaldococcus jannaschii"       
    ## NC_003272 Nostoc sp. PCC 7120                       "Nostoc sp. PCC 7120"                 
    ## AE009441 Pyrobaculum aerophilum                     "Pyrobaculum aerophilum"              
    ## AE009950 Pyrococcus furiosus                        "Pyrococcus furiosus"                 
    ## AE009951 Fusobacterium nucleatum                    "Fusobacterium nucleatum"             
    ## AE010299 Methanosarcina acetivorans                 "Methanosarcina acetivorans"          
    ## AE009439 Methanopyrus kandleri                      "Methanopyrus kandleri"               
    ## NC_003911 Ruegeria pomeroyi                         "Ruegeria pomeroyi"                   
    ## AE006470 Chlorobaculum tepidum                      "Chlorobaculum tepidum"               
    ## AE015928 Bacteroides thetaiotaomicron               "Bacteroides thetaiotaomicron"        
    ## AL954747 Nitrosomonas europaea                      "Nitrosomonas europaea"               
    ## BX119912 Rhodopirellula baltica                     "Rhodopirellula baltica"              
    ## BX571656 Wolinella succinogenes                     "Wolinella succinogenes"              
    ## AE017180 Geobacter sulfurreducens                   "Geobacter sulfurreducens"            
    ## AE017226 Treponema denticola                        "Treponema denticola"                 
    ## BX950229 Methanococcus maripaludis                  "Methanococcus maripaludis"           
    ## AE017221 Thermus thermophilus                       "Thermus thermophilus"                
    ## BA000001 Pyrococcus horikoshii                      "Pyrococcus horikoshii"               
    ## BA000023 Sulfolobus tokodaii                        "Sulfolobus tokodaii"                 
    ## NC_007951 Paraburkholderia xenovorans               "Paraburkholderia xenovorans"         
    ## CP000492 Chlorobium phaeobacteroides                "Chlorobium phaeobacteroides"         
    ## NC_008751 Desulfovibrio vulgaris                    "Desulfovibrio vulgaris"              
    ## CP000568 Ruminiclostridium thermocellum             "Ruminiclostridium thermocellum"      
    ## CP000561 Pyrobaculum calidifontis                   "Pyrobaculum calidifontis"            
    ## CP000609 Methanococcus maripaludis                  "Methanococcus maripaludis"           
    ## CP000607 Chlorobium phaeovibrioides                 "Chlorobium phaeovibrioides"          
    ## CP000660 Pyrobaculum arsenaticum                    "Pyrobaculum arsenaticum"             
    ## CP000667 Salinispora tropica                        "Salinispora tropica"                 
    ## CP000679 Caldicellulosiruptor saccharolyticus       "Caldicellulosiruptor saccharolyticus"
    ## CP000702 Thermotoga petrophila                      "Thermotoga petrophila"               
    ## CP000139 Bacteroides vulgatus                       "Bacteroides vulgatus"                
    ## NC_009665 Shewanella baltica                        "Shewanella baltica"                  
    ## CP000816 Ignicoccus hospitalis                      "Ignicoccus hospitalis"               
    ## CP000850 Salinispora arenicola                      "Salinispora arenicola"               
    ## CP000909 Chloroflexus aurantiacus                   "Chloroflexus aurantiacus"            
    ## CP000924 Thermoanaerobacter pseudethanolicus        "Thermoanaerobacter pseudethanolicus" 
    ## CP000969 Thermotoga sp. RQ2                         "Thermotoga sp. RQ2"                  
    ## CP001013 Leptothrix cholodnii                       "Leptothrix cholodnii"                
    ## CP001071 Akkermansia muciniphila                    "Akkermansia muciniphila"             
    ## AP009380 Porphyromonas gingivalis                   "Porphyromonas gingivalis"            
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1          "Sulfurihydrogenibium sp. YO3AOP1"    
    ## CP001097 Chlorobium limicola                        "Chlorobium limicola"                 
    ## CP001110 Pelodictyon phaeoclathratiforme            "Pelodictyon phaeoclathratiforme"     
    ## CP001130 Hydrogenobaculum sp. Y04AAS1               "Hydrogenobaculum sp. Y04AAS1"        
    ## NZ_CH959311 Sulfitobacter sp. EE-36                 "Sulfitobacter sp. EE-36"             
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1              "Sulfitobacter sp. NAS-14.1"          
    ## CP001251 Dictyoglomus turgidum                      "Dictyoglomus turgidum"               
    ## NC_011663 Shewanella baltica                        "Shewanella baltica"                  
    ## CP000916 Thermotoga neapolitana                     "Thermotoga neapolitana"              
    ## NZ_DS996397 Desulfovibrio piger                     "Desulfovibrio piger"                 
    ## CP001230 Persephonella marina                       "Persephonella marina"                
    ## CP001472 Acidobacterium capsulatum                  "Acidobacterium capsulatum"           
    ## AP009153 Gemmatimonas aurantiaca                    "Gemmatimonas aurantiaca"             
    ## CP001941 Aciduliprofundum boonei                    "Aciduliprofundum boonei"             
    ## NC_013968 Haloferax volcanii                        "Haloferax volcanii"                  
    ## NZ_KE136524 Enterococcus faecalis                   "Enterococcus faecalis"               
    ## NZ_KQ961402 Zymomonas mobilis                       "Zymomonas mobilis"                   
    ## NZ_CP015081 Deinococcus radiodurans                 "Deinococcus radiodurans"             
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense "Sulfurihydrogenibium yellowstonense" 
    ## NZ_JGWU01000001 Bordetella bronchiseptica           "Bordetella bronchiseptica"           
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii         "Caldicellulosiruptor bescii"         
    ## NC_009972 Herpetosiphon aurantiacus                 "Herpetosiphon aurantiacus"           
    ## NC_005213 Nanoarchaeum equitans                     "Nanoarchaeum equitans"               
    ##                                                     strain                                          
    ## AE000782 Archaeoglobus fulgidus                     "Archaeoglobus fulgidus DSM 4304"               
    ## NC_000909 Methanocaldococcus jannaschii             "Methanocaldococcus jannaschii DSM 2661"        
    ## NC_003272 Nostoc sp. PCC 7120                       ""                                              
    ## AE009441 Pyrobaculum aerophilum                     "Pyrobaculum aerophilum str. IM2"               
    ## AE009950 Pyrococcus furiosus                        "Pyrococcus furiosus DSM 3638"                  
    ## AE009951 Fusobacterium nucleatum                    ""                                              
    ## AE010299 Methanosarcina acetivorans                 "Methanosarcina acetivorans C2A"                
    ## AE009439 Methanopyrus kandleri                      "Methanopyrus kandleri AV19"                    
    ## NC_003911 Ruegeria pomeroyi                         "Ruegeria pomeroyi DSS-3"                       
    ## AE006470 Chlorobaculum tepidum                      "Chlorobaculum tepidum TLS"                     
    ## AE015928 Bacteroides thetaiotaomicron               "Bacteroides thetaiotaomicron VPI-5482"         
    ## AL954747 Nitrosomonas europaea                      "Nitrosomonas europaea ATCC 19718"              
    ## BX119912 Rhodopirellula baltica                     "Rhodopirellula baltica SH 1"                   
    ## BX571656 Wolinella succinogenes                     "Wolinella succinogenes DSM 1740"               
    ## AE017180 Geobacter sulfurreducens                   "Geobacter sulfurreducens PCA"                  
    ## AE017226 Treponema denticola                        "Treponema denticola ATCC 35405"                
    ## BX950229 Methanococcus maripaludis                  "Methanococcus maripaludis S2"                  
    ## AE017221 Thermus thermophilus                       "Thermus thermophilus HB27"                     
    ## BA000001 Pyrococcus horikoshii                      "Pyrococcus horikoshii OT3"                     
    ## BA000023 Sulfolobus tokodaii                        "Sulfolobus tokodaii str. 7"                    
    ## NC_007951 Paraburkholderia xenovorans               "Paraburkholderia xenovorans LB400"             
    ## CP000492 Chlorobium phaeobacteroides                "Chlorobium phaeobacteroides DSM 266"           
    ## NC_008751 Desulfovibrio vulgaris                    "Desulfovibrio vulgaris DP4"                    
    ## CP000568 Ruminiclostridium thermocellum             "Ruminiclostridium thermocellum ATCC 27405"     
    ## CP000561 Pyrobaculum calidifontis                   "Pyrobaculum calidifontis JCM 11548"            
    ## CP000609 Methanococcus maripaludis                  "Methanococcus maripaludis C5"                  
    ## CP000607 Chlorobium phaeovibrioides                 "Chlorobium phaeovibrioides DSM 265"            
    ## CP000660 Pyrobaculum arsenaticum                    "Pyrobaculum arsenaticum DSM 13514"             
    ## CP000667 Salinispora tropica                        "Salinispora tropica CNB-440"                   
    ## CP000679 Caldicellulosiruptor saccharolyticus       "Caldicellulosiruptor saccharolyticus DSM 8903" 
    ## CP000702 Thermotoga petrophila                      "Thermotoga petrophila RKU-1"                   
    ## CP000139 Bacteroides vulgatus                       "Bacteroides vulgatus ATCC 8482"                
    ## NC_009665 Shewanella baltica                        "Shewanella baltica OS185"                      
    ## CP000816 Ignicoccus hospitalis                      "Ignicoccus hospitalis KIN4/I"                  
    ## CP000850 Salinispora arenicola                      "Salinispora arenicola CNS-205"                 
    ## CP000909 Chloroflexus aurantiacus                   "Chloroflexus aurantiacus J-10-fl"              
    ## CP000924 Thermoanaerobacter pseudethanolicus        "Thermoanaerobacter pseudethanolicus ATCC 33223"
    ## CP000969 Thermotoga sp. RQ2                         ""                                              
    ## CP001013 Leptothrix cholodnii                       "Leptothrix cholodnii SP-6"                     
    ## CP001071 Akkermansia muciniphila                    "Akkermansia muciniphila ATCC BAA-835"          
    ## AP009380 Porphyromonas gingivalis                   "Porphyromonas gingivalis ATCC 33277"           
    ## NC_010730 Sulfurihydrogenibium sp. YO3AOP1          ""                                              
    ## CP001097 Chlorobium limicola                        "Chlorobium limicola DSM 245"                   
    ## CP001110 Pelodictyon phaeoclathratiforme            "Pelodictyon phaeoclathratiforme BU-1"          
    ## CP001130 Hydrogenobaculum sp. Y04AAS1               ""                                              
    ## NZ_CH959311 Sulfitobacter sp. EE-36                 ""                                              
    ## NZ_CH959317 Sulfitobacter sp. NAS-14.1              ""                                              
    ## CP001251 Dictyoglomus turgidum                      "Dictyoglomus turgidum DSM 6724"                
    ## NC_011663 Shewanella baltica                        "Shewanella baltica OS223"                      
    ## CP000916 Thermotoga neapolitana                     "Thermotoga neapolitana DSM 4359"               
    ## NZ_DS996397 Desulfovibrio piger                     "Desulfovibrio piger ATCC 29098"                
    ## CP001230 Persephonella marina                       "Persephonella marina EX-H1"                    
    ## CP001472 Acidobacterium capsulatum                  "Acidobacterium capsulatum ATCC 51196"          
    ## AP009153 Gemmatimonas aurantiaca                    "Gemmatimonas aurantiaca T-27"                  
    ## CP001941 Aciduliprofundum boonei                    "Aciduliprofundum boonei T469"                  
    ## NC_013968 Haloferax volcanii                        "Haloferax volcanii DS2"                        
    ## NZ_KE136524 Enterococcus faecalis                   "Enterococcus faecalis V583"                    
    ## NZ_KQ961402 Zymomonas mobilis                       ""                                              
    ## NZ_CP015081 Deinococcus radiodurans                 "Deinococcus radiodurans R1"                    
    ## NZ_ABZS01000228 Sulfurihydrogenibium yellowstonense "Sulfurihydrogenibium yellowstonense SS-5"      
    ## NZ_JGWU01000001 Bordetella bronchiseptica           "Bordetella bronchiseptica D989"                
    ## NZ_FWDH01000003 Caldicellulosiruptor bescii         ""                                              
    ## NC_009972 Herpetosiphon aurantiacus                 "Herpetosiphon aurantiacus DSM 785"             
    ## NC_005213 Nanoarchaeum equitans                     "Nanoarchaeum equitans Kin4-M"

``` r
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
metadata
```

    ##                   study_name      sample_id         subject_id  body_site
    ## ERR1136663       ZeeviD_2015   PNP_Main_553       PNP_Main_553      stool
    ## ERR4562131 MetaCardis_2020_a   M0x20MCx1812       M0x20MCx1812      stool
    ## ERR4087400 MetaCardis_2020_a   M0x20MCx1360       M0x20MCx1360      stool
    ## SRR9217435      GhensiP_2019 SP_101SPI_T016            sub_101 oralcavity
    ## SRR059888           HMP_2012      SRS013723 HMP_2012_159268001 oralcavity
    ##            antibiotics_current_use study_condition disease age age_category
    ## ERR1136663                      no         control healthy  51        adult
    ## ERR4562131                      no             T2D     T2D  63        adult
    ## ERR4087400                     yes             IGT  IGT;MS  51        adult
    ## SRR9217435                    <NA>         control healthy  NA        adult
    ## SRR059888                     <NA>         control healthy  27        adult
    ##            gender country location non_westernized sequencing_platform
    ## ERR1136663 female     ISR     <NA>              no       IlluminaHiSeq
    ## ERR4562131   male     DEU  Leipzig              no           IonProton
    ## ERR4087400 female     DEU  Leipzig              no           IonProton
    ## SRR9217435   <NA>     ITA     <NA>              no       IlluminaHiSeq
    ## SRR059888    male     USA     <NA>              no       IlluminaHiSeq
    ##            dna_extraction_kit     pmid number_reads number_bases
    ## ERR1136663               <NA> 26590418     28756422   2515847316
    ## ERR4562131               <NA> 34880489     15374062   2326962418
    ## ERR4087400               <NA> 34880489     20888458   3200734644
    ## SRR9217435               <NA> 33127901     20061461   1987501438
    ## SRR059888              Qiagen 22699609     94436101   8712278548
    ##            minimum_read_length median_read_length ncbi_accession hscrp
    ## ERR1136663                  20                 80     ERR1136663    NA
    ## ERR4562131                  NA                 NA     ERR4562131    NA
    ## ERR4087400                  NA                 NA     ERR4087400  6.78
    ## SRR9217435                  75                100     SRR9217435    NA
    ## SRR059888                   60                101      SRR059888    NA
    ##                           treatment                   curator pregnant      bmi
    ## ERR1136663                       no Jacob_Wirbel;Paolo_Manghi       no 23.70000
    ## ERR4562131 antidiab;metformin;dppiv              Paolo_Manghi     <NA> 27.67959
    ## ERR4087400          antihta;ca2_cbl              Paolo_Manghi     <NA> 50.14402
    ## SRR9217435                     <NA>              Paolo_Manghi     <NA>       NA
    ## SRR059888                      <NA>              Paolo_Manghi     <NA> 23.00000
    ##            smoker         body_subsite
    ## ERR1136663   <NA>          rectal_swab
    ## ERR4562131    yes                 <NA>
    ## ERR4087400     no                 <NA>
    ## SRR9217435     no   subgingival_plaque
    ## SRR059888    <NA> supragingival_plaque
    ##                                                                                                                    uncurated_metadata
    ## ERR1136663 no_immuno_suppressive;no_T2D;no_T1D;no_related_treatments;no_psychiatric_diseases;no_gastro_intestinal_disorder;non_celiac
    ## ERR4562131                                                                                                                       <NA>
    ## ERR4087400                                                                                                                       <NA>
    ## SRR9217435                                                                                                                       <NA>
    ## SRR059888                                                                                                                        <NA>
    ##                 ldl triglycerides hba1c ever_smoker dental_sample_type
    ## ERR1136663       NA            NA    NA        <NA>               <NA>
    ## ERR4562131       NA            NA    NA        <NA>               <NA>
    ## ERR4087400 210.3648      240.0247  5.73        <NA>               <NA>
    ## SRR9217435       NA            NA    NA         yes            implant
    ## SRR059888        NA            NA    NA        <NA>               <NA>
    ##            history_of_periodontitis ppd_m ppd_b ppd_d ppd_l bristol_score
    ## ERR1136663                     <NA>    NA    NA    NA    NA            NA
    ## ERR4562131                     <NA>    NA    NA    NA    NA            NA
    ## ERR4087400                     <NA>    NA    NA    NA    NA             4
    ## SRR9217435                      yes     4     3     4     3            NA
    ## SRR059888                      <NA>    NA    NA    NA    NA            NA
    ##            ncbi_bioproject
    ## ERR1136663      PRJEB11532
    ## ERR4562131      PRJEB38742
    ## ERR4087400      PRJEB37249
    ## SRR9217435     PRJNA547717
    ## SRR059888       PRJNA48479

Finally, merge the metadata with the sourmash phyloseq object

``` r
pseq_final <- merge_phyloseq(physeq, metadata)
pseq_final
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 15 taxa and 5 samples ]
    ## sample_data() Sample Data:       [ 5 samples by 41 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 15 taxa by 8 taxonomic ranks ]
