#!/bin/bash

## Este script agrega el simbolo ^ al principo del cada nombre de muestra
## y concatena todos los archivos phy (un archivo phy por locus)

for i in `ls ../data/phy`; do ## enlista todos los archivos phy

  ## Agrega el simbolo ^ a cada nombre de muestra de cada locus
  cat ../data/phy/$i | sed -e '2,$ s/^/^/'  >> ../out/mamm_opt_70_bpp.txt

done
