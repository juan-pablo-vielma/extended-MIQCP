#!/bin/bash
source activate py3k
perprof --tikz --tex CPLEXLP-prof CPLEXCP-prof CPLEXSepLazy-prof CPLEXSepLP-prof GurobiLP-prof GurobiCP-prof GurobiSepLazy-prof GurobiSepLP-prof LiftedLP-prof -o all  --semilog  -f --maxtime 3600
perprof --tikz --tex CPLEXSepLP-prof CPLEXTowerLP-prof CPLEXTowerSepLP-prof GurobiSepLP-prof GurobiTowerLP-prof GurobiTowerSepLP-prof -o formulations  --semilog  -f --maxtime 3600

sed -e 's/mark=none/mark\ repeat=600/g' -e 's/\\end{center}//g' -e 's/\\begin{center}//g' -e 's/plot,/plot,cycle\ list\ name=mylist,/g' all.tex > delme.tex
cat header.tex delme.tex footer.tex > all.tex
lualatex all.tex

sed -e 's/mark=none/mark\ repeat=600/g' -e 's/\\end{center}//g' -e 's/\\begin{center}//g' -e 's/plot,/plot,cycle\ list\ name=mylist,/g' formulations.tex > delme.tex
cat header.tex delme.tex footer.tex > formulations.tex
lualatex formulations.tex
