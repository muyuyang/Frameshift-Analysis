# Frameshift-Analysis
Frameshift-Analysis is a python-based tool for analyzing frameshift mutations between genomes. 


## Requirements & Installation
* Python >= 3.3
* LAST ([Install](http://last.cbrc.jp/doc/last.html))
* Biopython ([Install](https://biopython.org/wiki/Download))
* ete3 ([Install](http://etetoolkit.org/download/))



## Features
* Pairwise comparison

  Given two genomes, assuming one is ancestral and the other is derived, we can find:
  * the statistics of mutations(e.g. Number of frameshifts)
  * all frameshifts on protein-coding regions
  * information of proteins affected by frameshifts(e.g. function, annotation of domains)


* Timing of frameshift mutation

  Given genomes and their phylogenetic relationship, we can find where each frameshift happened.


## License

MIT license

