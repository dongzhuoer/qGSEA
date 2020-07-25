# qGSEA
[![Build Status](https://travis-ci.com/dongzhuoer/qGSEA.svg?branch=master)](https://travis-ci.com/dongzhuoer/qGSEA)

## Overview

[GSEA](http://software.broadinstitute.org/gsea/index.jsp) is a powerful tool, but preparing its input files, escpecially from GSE raw data, is a painful process. I know that from very personal experience, that's why I develop thistool. You just need a GSE accession and qGSEA would take care of everything else for you.


## Installation

```r
if (!('devtools' %in% .packages(T))) install.packages('devtools')
remotes::install_github('dongzhuoer/qGSEA')
```


## Usage

refer to `vignette('qGSEA')`.


## Tips

For using GSEA in non-human organisms, you can find gene sets [here](http://ge-lab.org/gskb/) (`.gmt` files ready for downloading).

Especially, for mouse, you shouldn't miss the Bioconductor package [gskb](https://mirrors.tuna.tsinghua.edu.cn/bioconductor/packages/release/data/experiment/html/gskb.html), it collect gene sets from [various sources](http://ge-lab.org/gskb/Table%201-sources.pdf).

