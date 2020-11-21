DADA2
================
Madison Wright
11/11/2020

This script processes trimmed (w/o primers) sequences through the [DADA2
pipline (v 1.16)](https://benjjneb.github.io/dada2/tutorial.html), which
can be installed following these
[steps](https://benjjneb.github.io/dada2/dada-installation.html)

# Install and Load DADA2 and ShortRead from Bioconductor

``` r
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.11")
# BiocManager::install("ShortRead")
```

``` r
# BiocManager::install(version = '3.11')
```

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.4     ✓ dplyr   1.0.2
    ## ✓ tidyr   1.1.2     ✓ stringr 1.4.0
    ## ✓ readr   1.4.0     ✓ forcats 0.5.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(dada2)
```

    ## Loading required package: Rcpp

``` r
library(ShortRead)
```

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Loading required package: BiocParallel

    ## Loading required package: Biostrings

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'XVector'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: Rsamtools

    ## Loading required package: GenomeInfoDb

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomicAlignments

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

    ## 
    ## Attaching package: 'GenomicAlignments'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     last

    ## 
    ## Attaching package: 'ShortRead'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     id

    ## The following object is masked from 'package:purrr':
    ## 
    ##     compose

    ## The following object is masked from 'package:tibble':
    ## 
    ##     view

# Import file names

``` r
path <- "~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq/"

#store the names of the forward and rev files as lists
fnFs <- list.files(path, pattern = "_R1_001.fastq", full.names = TRUE)
fnRs <- list.files(path, pattern = "_R2_001.fastq", full.names = TRUE)
```

# Retrieve orientation of primers

The primers targeted the V4 region and are known 514F-Y and 806RB
primers (see Apprill et al.,
2015)\[<http://www.int-res.com/articles/ame_oa/a075p129.pdf>\]

``` r
#store the  forward and reverse primers
FWD = "GTGYCAGCMGCCGCGGTAA"
REV = "GGACTACNVGGGTWTCTAAT"

#now store all the orientations of your forward and reverse  primers
allOrients <- function(primer) {
  # The Biostrings works w/ DNAString objects rather than character vectors
  require(Biostrings)
  dna <- DNAString(primer) 
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  # Convert back to character vector
  return(sapply(orients, toString))  
}
#store the fwd and reverse oreintations separately
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

#view the orientations of the primers
FWD.orients
```

    ##               Forward            Complement               Reverse 
    ## "GTGYCAGCMGCCGCGGTAA" "CACRGTCGKCGGCGCCATT" "AATGGCGCCGMCGACYGTG" 
    ##               RevComp 
    ## "TTACCGCGGCKGCTGRCAC"

``` r
REV.orients
```

    ##                Forward             Complement                Reverse 
    ## "GGACTACNVGGGTWTCTAAT" "CCTGATGNBCCCAWAGATTA" "TAATCTWTGGGVNCATCAGG" 
    ##                RevComp 
    ## "ATTAGAWACCCBNGTAGTCC"

# search for Primers

``` r
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
```

    ##                  Forward Complement Reverse RevComp
    ## FWD.ForwardReads       0          0       0       0
    ## FWD.ReverseReads       0          0       0     283
    ## REV.ForwardReads       0          0       0    1195
    ## REV.ReverseReads       0          0       0       0

# Inspect read quality profiles

## forward reads

``` r
plotQualityProfile(fnFs[1:12])
```

![](Dada2_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Reverse reads

``` r
plotQualityProfile(fnRs[1:12])
```

![](Dada2_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

# Filtering and Trimming

``` r
sample.names <- sapply(strsplit(basename(fnFs),"_L"), `[`,1)
sample.names
```

    ##  [1] "144_A0_S6"  "144_A4_S7"  "144_A8_S8"  "144_B0_S9"  "144_B4_S10"
    ##  [6] "144_B8_S11" "144_C0_S12" "144_C4_S13" "144_C8_S14" "144_D0_S15"
    ## [11] "144_D4_S16" "144_D8_S17" "144_E0_S18" "144_E4_S19" "144_E8_S20"
    ## [16] "144_F0_S21" "144_F4_S22" "144_F8_S23" "144_G0_S24" "144_G4_S25"
    ## [21] "144_G8_S26" "144_H0_S27" "144_H4_S28" "144_H8_S29"

``` r
#create a "filtered" folder in the working directory as a place to put all the new filtered fastQ files
filt_path <- file.path(path,"filtered")
#add the appropriate designation string to any new files made that will be put into the "filtered" folder
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
```

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(200,150),  maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE) 
out
```

    ##                              reads.in reads.out
    ## 144_A0_S6_L001_R1_001.fastq     79521     70868
    ## 144_A4_S7_L001_R1_001.fastq     36911     33172
    ## 144_A8_S8_L001_R1_001.fastq     50351     45283
    ## 144_B0_S9_L001_R1_001.fastq     14331      3286
    ## 144_B4_S10_L001_R1_001.fastq    50889     45432
    ## 144_B8_S11_L001_R1_001.fastq    62383     56216
    ## 144_C0_S12_L001_R1_001.fastq    33962     25214
    ## 144_C4_S13_L001_R1_001.fastq    69734     62942
    ## 144_C8_S14_L001_R1_001.fastq    60793     54965
    ## 144_D0_S15_L001_R1_001.fastq    62933     57502
    ## 144_D4_S16_L001_R1_001.fastq    49383     44769
    ## 144_D8_S17_L001_R1_001.fastq    61144     55119
    ## 144_E0_S18_L001_R1_001.fastq    53714     49019
    ## 144_E4_S19_L001_R1_001.fastq    41686     37467
    ## 144_E8_S20_L001_R1_001.fastq    34947     30413
    ## 144_F0_S21_L001_R1_001.fastq    54554     49975
    ## 144_F4_S22_L001_R1_001.fastq    32800     29889
    ## 144_F8_S23_L001_R1_001.fastq    33312     30792
    ## 144_G0_S24_L001_R1_001.fastq    40935     37648
    ## 144_G4_S25_L001_R1_001.fastq    40109     36506
    ## 144_G8_S26_L001_R1_001.fastq    35610     32774
    ## 144_H0_S27_L001_R1_001.fastq    63711     59242
    ## 144_H4_S28_L001_R1_001.fastq    27892     25062
    ## 144_H8_S29_L001_R1_001.fastq    36860     33766

# Learn the error rates

``` r
errF <- learnErrors(filtFs, multithread = TRUE)
```

    ## 110953600 total bases in 554768 reads from 12 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread = TRUE)
```

    ## 100750050 total bases in 671667 reads from 15 samples will be used for learning the error rates.

<img src="Dada2_files/figure-gfm/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

<img src="Dada2_files/figure-gfm/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

# Dereplication

``` r
derepFs <- derepFastq(filtFs, verbose = TRUE)
```

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_A0_S6_F_filt.fastq

    ## Encountered 17611 unique sequences from 70868 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_A4_S7_F_filt.fastq

    ## Encountered 8523 unique sequences from 33172 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_A8_S8_F_filt.fastq

    ## Encountered 11776 unique sequences from 45283 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_B0_S9_F_filt.fastq

    ## Encountered 1394 unique sequences from 3286 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_B4_S10_F_filt.fastq

    ## Encountered 10927 unique sequences from 45432 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_B8_S11_F_filt.fastq

    ## Encountered 13739 unique sequences from 56216 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_C0_S12_F_filt.fastq

    ## Encountered 7426 unique sequences from 25214 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_C4_S13_F_filt.fastq

    ## Encountered 12758 unique sequences from 62942 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_C8_S14_F_filt.fastq

    ## Encountered 12894 unique sequences from 54965 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_D0_S15_F_filt.fastq

    ## Encountered 13967 unique sequences from 57502 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_D4_S16_F_filt.fastq

    ## Encountered 8744 unique sequences from 44769 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_D8_S17_F_filt.fastq

    ## Encountered 13616 unique sequences from 55119 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_E0_S18_F_filt.fastq

    ## Encountered 12580 unique sequences from 49019 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_E4_S19_F_filt.fastq

    ## Encountered 8678 unique sequences from 37467 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_E8_S20_F_filt.fastq

    ## Encountered 7830 unique sequences from 30413 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_F0_S21_F_filt.fastq

    ## Encountered 12981 unique sequences from 49975 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_F4_S22_F_filt.fastq

    ## Encountered 7337 unique sequences from 29889 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_F8_S23_F_filt.fastq

    ## Encountered 7303 unique sequences from 30792 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_G0_S24_F_filt.fastq

    ## Encountered 9886 unique sequences from 37648 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_G4_S25_F_filt.fastq

    ## Encountered 7586 unique sequences from 36506 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_G8_S26_F_filt.fastq

    ## Encountered 7592 unique sequences from 32774 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_H0_S27_F_filt.fastq

    ## Encountered 13584 unique sequences from 59242 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_H4_S28_F_filt.fastq

    ## Encountered 5261 unique sequences from 25062 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_H8_S29_F_filt.fastq

    ## Encountered 8537 unique sequences from 33766 total sequences read.

``` r
derepRs <- derepFastq(filtRs, verbose = TRUE)
```

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_A0_S6_R_filt.fastq

    ## Encountered 24175 unique sequences from 70868 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_A4_S7_R_filt.fastq

    ## Encountered 11316 unique sequences from 33172 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_A8_S8_R_filt.fastq

    ## Encountered 16625 unique sequences from 45283 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_B0_S9_R_filt.fastq

    ## Encountered 1770 unique sequences from 3286 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_B4_S10_R_filt.fastq

    ## Encountered 15567 unique sequences from 45432 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_B8_S11_R_filt.fastq

    ## Encountered 20652 unique sequences from 56216 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_C0_S12_R_filt.fastq

    ## Encountered 12162 unique sequences from 25214 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_C4_S13_R_filt.fastq

    ## Encountered 19438 unique sequences from 62942 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_C8_S14_R_filt.fastq

    ## Encountered 18981 unique sequences from 54965 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_D0_S15_R_filt.fastq

    ## Encountered 20947 unique sequences from 57502 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_D4_S16_R_filt.fastq

    ## Encountered 14888 unique sequences from 44769 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_D8_S17_R_filt.fastq

    ## Encountered 21171 unique sequences from 55119 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_E0_S18_R_filt.fastq

    ## Encountered 18785 unique sequences from 49019 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_E4_S19_R_filt.fastq

    ## Encountered 14665 unique sequences from 37467 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_E8_S20_R_filt.fastq

    ## Encountered 15767 unique sequences from 30413 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_F0_S21_R_filt.fastq

    ## Encountered 19823 unique sequences from 49975 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_F4_S22_R_filt.fastq

    ## Encountered 10390 unique sequences from 29889 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_F8_S23_R_filt.fastq

    ## Encountered 10140 unique sequences from 30792 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_G0_S24_R_filt.fastq

    ## Encountered 14019 unique sequences from 37648 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_G4_S25_R_filt.fastq

    ## Encountered 12399 unique sequences from 36506 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_G8_S26_R_filt.fastq

    ## Encountered 10835 unique sequences from 32774 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_H0_S27_R_filt.fastq

    ## Encountered 20646 unique sequences from 59242 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_H4_S28_R_filt.fastq

    ## Encountered 10356 unique sequences from 25062 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/Github/144l_students/Input_Data/week5/EEMB144L_2018_fastq//filtered/144_H8_S29_R_filt.fastq

    ## Encountered 12193 unique sequences from 33766 total sequences read.

``` r
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

# Infer the sequence variants

``` r
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
```

    ## Sample 1 - 70868 reads in 17611 unique sequences.
    ## Sample 2 - 33172 reads in 8523 unique sequences.
    ## Sample 3 - 45283 reads in 11776 unique sequences.
    ## Sample 4 - 3286 reads in 1394 unique sequences.
    ## Sample 5 - 45432 reads in 10927 unique sequences.
    ## Sample 6 - 56216 reads in 13739 unique sequences.
    ## Sample 7 - 25214 reads in 7426 unique sequences.
    ## Sample 8 - 62942 reads in 12758 unique sequences.
    ## Sample 9 - 54965 reads in 12894 unique sequences.
    ## Sample 10 - 57502 reads in 13967 unique sequences.
    ## Sample 11 - 44769 reads in 8744 unique sequences.
    ## Sample 12 - 55119 reads in 13616 unique sequences.
    ## Sample 13 - 49019 reads in 12580 unique sequences.
    ## Sample 14 - 37467 reads in 8678 unique sequences.
    ## Sample 15 - 30413 reads in 7830 unique sequences.
    ## Sample 16 - 49975 reads in 12981 unique sequences.
    ## Sample 17 - 29889 reads in 7337 unique sequences.
    ## Sample 18 - 30792 reads in 7303 unique sequences.
    ## Sample 19 - 37648 reads in 9886 unique sequences.
    ## Sample 20 - 36506 reads in 7586 unique sequences.
    ## Sample 21 - 32774 reads in 7592 unique sequences.
    ## Sample 22 - 59242 reads in 13584 unique sequences.
    ## Sample 23 - 25062 reads in 5261 unique sequences.
    ## Sample 24 - 33766 reads in 8537 unique sequences.

``` r
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
```

    ## Sample 1 - 70868 reads in 24175 unique sequences.
    ## Sample 2 - 33172 reads in 11316 unique sequences.
    ## Sample 3 - 45283 reads in 16625 unique sequences.
    ## Sample 4 - 3286 reads in 1770 unique sequences.
    ## Sample 5 - 45432 reads in 15567 unique sequences.
    ## Sample 6 - 56216 reads in 20652 unique sequences.
    ## Sample 7 - 25214 reads in 12162 unique sequences.
    ## Sample 8 - 62942 reads in 19438 unique sequences.
    ## Sample 9 - 54965 reads in 18981 unique sequences.
    ## Sample 10 - 57502 reads in 20947 unique sequences.
    ## Sample 11 - 44769 reads in 14888 unique sequences.
    ## Sample 12 - 55119 reads in 21171 unique sequences.
    ## Sample 13 - 49019 reads in 18785 unique sequences.
    ## Sample 14 - 37467 reads in 14665 unique sequences.
    ## Sample 15 - 30413 reads in 15767 unique sequences.
    ## Sample 16 - 49975 reads in 19823 unique sequences.
    ## Sample 17 - 29889 reads in 10390 unique sequences.
    ## Sample 18 - 30792 reads in 10140 unique sequences.
    ## Sample 19 - 37648 reads in 14019 unique sequences.
    ## Sample 20 - 36506 reads in 12399 unique sequences.
    ## Sample 21 - 32774 reads in 10835 unique sequences.
    ## Sample 22 - 59242 reads in 20646 unique sequences.
    ## Sample 23 - 25062 reads in 10356 unique sequences.
    ## Sample 24 - 33766 reads in 12193 unique sequences.

``` r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, trimOverhang = T)
```

    ## 67488 paired-reads (in 337 unique pairings) successfully merged out of 70016 (in 919 pairings) input.

    ## 31521 paired-reads (in 166 unique pairings) successfully merged out of 32548 (in 495 pairings) input.

    ## 43338 paired-reads (in 229 unique pairings) successfully merged out of 44672 (in 637 pairings) input.

    ## 2834 paired-reads (in 70 unique pairings) successfully merged out of 3053 (in 180 pairings) input.

    ## 43547 paired-reads (in 189 unique pairings) successfully merged out of 44815 (in 563 pairings) input.

    ## 53945 paired-reads (in 249 unique pairings) successfully merged out of 55532 (in 682 pairings) input.

    ## 23399 paired-reads (in 207 unique pairings) successfully merged out of 24632 (in 533 pairings) input.

    ## 61134 paired-reads (in 173 unique pairings) successfully merged out of 62173 (in 545 pairings) input.

    ## 53096 paired-reads (in 229 unique pairings) successfully merged out of 54429 (in 574 pairings) input.

    ## 54656 paired-reads (in 298 unique pairings) successfully merged out of 56615 (in 722 pairings) input.

    ## 43396 paired-reads (in 123 unique pairings) successfully merged out of 44304 (in 344 pairings) input.

    ## 53165 paired-reads (in 220 unique pairings) successfully merged out of 54536 (in 555 pairings) input.

    ## 46486 paired-reads (in 265 unique pairings) successfully merged out of 48192 (in 614 pairings) input.

    ## 35900 paired-reads (in 163 unique pairings) successfully merged out of 37001 (in 436 pairings) input.

    ## 28472 paired-reads (in 161 unique pairings) successfully merged out of 29793 (in 525 pairings) input.

    ## 46737 paired-reads (in 255 unique pairings) successfully merged out of 48889 (in 655 pairings) input.

    ## 28793 paired-reads (in 138 unique pairings) successfully merged out of 29584 (in 299 pairings) input.

    ## 29694 paired-reads (in 190 unique pairings) successfully merged out of 30478 (in 355 pairings) input.

    ## 35711 paired-reads (in 226 unique pairings) successfully merged out of 37011 (in 451 pairings) input.

    ## 35287 paired-reads (in 138 unique pairings) successfully merged out of 36111 (in 377 pairings) input.

    ## 31371 paired-reads (in 158 unique pairings) successfully merged out of 32263 (in 325 pairings) input.

    ## 54517 paired-reads (in 302 unique pairings) successfully merged out of 58549 (in 630 pairings) input.

    ## 23994 paired-reads (in 118 unique pairings) successfully merged out of 24764 (in 319 pairings) input.

    ## 32510 paired-reads (in 162 unique pairings) successfully merged out of 33297 (in 320 pairings) input.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                        sequence
    ## 1 TACGGAGGGTGCAAGCGTTACTCGGAATCACTGGGCGTAAAGAGCGTGTAGGCGGATAGTTAAGTTTGAAGTGAAATCCTATGGCTCAACCATAGAACTGCTTTGAAAACTGATTATCTAGAATATGGGAGAGGTAGATGGAATTTCTGGTGTAGGGGTAAAATCCGTAGAGATCAGAAGGAATACCGATTGCGAAGGCGATCTACTGGAACATTATTGACGCTGAGACGCGAAAGCGTGGGGAGCAAACAGG
    ## 2 TACGGAGGGTCCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTTGTCAGTCAGATGTGAAAGCCCAGGGCTCAACCTTGGAACTGCACCTGATACTGCAAGACTAGAGTACAATAGAGGGGAGTGGAATTTCCGGTGTAGCGGTGAAATGCGTAGAGATCGGAAGGAACACCAGTGGCGAAGGCGACTCCCTGGATTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACAGG
    ## 3 TACGGAAGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCTGGGCTCAACCTAGGAACTGCATCCAAAACTAACTCACTAGAGTACGATAGAGGGAGGTAGAATTCATAGTGTAGCGGTGGAATGCGTAGATATTATGAAGAATACCAGTGGCGAAGGCGGCCTCCTGGATCTGTACTGACACTGAGGTGCGAAAGCGTGGGTAGCGAACAGG
    ## 4 TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTAGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATTAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGG
    ## 5 TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTAGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATTAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAGCGAAAGCATGGGTAGCGAAGAGG
    ## 6 TACGGGAGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGTCCGCAGGCGGCCCTTCAAGTCTGCTGTTAAAAAGTGGAGCTTAACTCCATCATGGCAGTGGAAACTGAGGGGCTTGAGTGTGGTAGGGGCAGAGGGAATTCCCGGTGTAGCGGTGAAATGCGTAGATATCGGGAAGAACACCAGTGGCGAAGGCGCTCTGCTGGGCCATCACTGACGCTCATGGACGAAAGCCAGGGGAGCGAAAGGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1      3773       2       1     97         0      0      1   TRUE
    ## 2      3731       3       3     97         0      0      1   TRUE
    ## 3      3139       4       4     97         0      0      1   TRUE
    ## 4      3038       1       2     97         0      0      1   TRUE
    ## 5      2715       1       5     97         0      0      1   TRUE
    ## 6      2390       5       6     97         0      0      1   TRUE

``` r
saveRDS(mergers, "~/Desktop/Github/144l_students/Output_Data/week 5/dada_merged.rds")
saveRDS(mergers, "~/Desktop/Github/144l_students/Input_Data/week 6/dada_merged.rds")
```

# Sequence Table

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # samples by unique sequence
```

    ## [1]  24 970

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 252 253 254 255 256 257 258 265 266 270 
    ##  25 873  52   4   5   2   4   1   1   3

# Remove the Chimeras

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose = TRUE)
```

    ## Identified 81 bimeras out of 970 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  24 889

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.9929833

# Assign taxonomy using a reference database

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/Github/144l_students/Input_Data/week5/Reference_Database/silva_nr_v138_train_set.fa", multithread = TRUE)
```

``` r
saveRDS(t(seqtab.nochim), "~/Desktop/Github/144l_students/Output_Data/week 5/seqtab-nochimtaxa.rds")
saveRDS(taxa,"~/Desktop/Github/144l_students/Output_Data/week 5/taxa.rds")

saveRDS(t(seqtab.nochim), "~/Desktop/Github/144l_students/Input_Data/week 6/seqtab-nochimtaxa.rds")
saveRDS(taxa,"~/Desktop/Github/144l_students/Input_Data/week 6/taxa.rds")
```
