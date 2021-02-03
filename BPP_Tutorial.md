## Tutorial para delimitación de especies con BBP (Bayesian Phylogenetics and Phylogeography)

___


Por José Rubén Montes y Cristian Cervantes


-

[BPP](https://github.com/bpp/bpp) es un programa Bayesiano publicado por [Yang, 2015](http://abacus.gene.ucl.ac.uk/ziheng/pdf/2015YangCZv61p854.pdf) que emplea Cadenas Markoviana (MCMC) y analiza datos multi-locus bajo el modelo coalescente de especies (MSC) que toma en cuenta el polimorfismo ancestral compartido y el sorteo incompleto de limajes (ILS). BBP construye árboles de especies y permite poner a prueba hipótesis de delimitación de especies bajo un árbol guía o sin *priors*. BPP permite correr 4 tipos de análisis (A00, A01, A10, A11). Cada tipo de análisis permite estimar diferentes aspectos. Por ejemplo, parámetros poblacionales (A00), árboles de especies (A01, A10) y delimatación de especies (A10, A11). 

**NOTA:** El **0** y **1** en los archivos **A** son comandos que permiten activar o desactivar un parámetro en el archivo control o ctl. En los tipos análisis BPP veremos cómo se utilizan el 0 y el 1.

 

## **1. Requerimentos primera parte**
|Archivos | inputs |
| ------   | ------ | 
|**A.** Secuencias de genes en formato PHYLIP | cem_opt _20 _bpp.txt |
|**B.** Hipótesis de delimitación (asignación de individuos a especies)| cem_opt _20 _bppmap.txt  |
|**C.** Árbol de especies en formato NEWICK| species_tree.nwk |

Ejemplo del archivo **species_bpp.txt.**

```

20 118
^ariKF30      GCTAAAATGCAAATATATTATGAAATAAATGTTACAAATAGTTAATTGCTTAACTGTTTGTGTGCTGTACGACATAATGATATTATCATCTGTTCTCTAGGAA
^balfOM11     GCTAAAATGCAAATATATTATGAAAGACATGTTACAAATAGTTAATTGCTTAACTG-TTGTGTGCTGTAATACATAA--ATGATATCATCTTTTCTCTAGGAA
^bun03s3A     ---------------TATTATGAAATACATGTTACAAATGGTTAATTGCTTAACTG-TTGTGTGGTGTAATACATAA--ATGATATCATCTGTTCTCTAGGAA
^calD1537     GCTAAAATGCAAATATATTATGAAATAAATGTTACAAATAGTTAATTGCTTAACTGTTTGTGTGCTGTACTACATAATGATATTATCATCTGTTCTCTAGGAA
...

20 119
^ariKF30      CGTAAAATGCAAATATATTATGAAATAAATGTTACAAATAGTTAATTGCTTAACTGTTTGTGTGCTGTCCGACATAATGATATTATCATCTGTTCTCTAGGAA
^balfOM11     CGTAAAATGCAAATATATTATGAAAGACATGTTACAAATAGTTAATTGCTTAACTG-CTGTGTGCTGTAATACATAA--ATGATATCATCTTTTCTCTAGGAA
^bun03s3A     ATTGCTTAACTGCCCCATTATGAAATACATGTTACAAATGGTTAATTGCTTAACTG-TTGTGTGGTGTAATACATAA--CGGATATCATCTGTTCTCTAGGAA
^calD1537     CGTAAAATGCAAATATATTATGAAATAAATGTTACAAATAGTTAATTGCTTAACTGTTTGTGTGCTGTCTTACATAATGATATTATCATCTGTTCTCTAGGAA
...



```

**NOTA**: El archivo debe contener todos los genes que vas a analizar y al inicio de cada individuo en todos los genes **debes** agregar el símbolo "^"



Ejemplo del archivo **individuals_Imap.txt.**

```
ariKF30   aristata
balfOM11  balfouriana
bun03s3A  bungeana
kremp371  krempfii
lamb1195  lambertiana
calD1537  californiarum
```


Ejemplo del árbol de especies **species_tree.nwk.**

```
(((californiarum,(aristata, balfouriana)),((bungeana),(krempfii,lambertiana)))); 
    
```
**NOTA**: El árbol debe editarse para eliminar los valores de probabilidad posterior o bootstrap y longitudes de rama.

Si necesitas estimar tu propio árbol de especies te recomiendo visitar el tutorial de [ASTRAL](https://github.com/JR-Montes/Tutorial_ASTRAL/blob/master/Tutorial_ASTRAL_Multi_and_Single_locus.md) o [SVDquartets](https://github.com/JR-Montes/Tutorial-SVDquartets)


## **2.** Elegir tipo de análisis en BPP.

#### Opción A00: Calcula los parámetros poblacionales bajo el modelo coalescente 

Esta opcíon permite calcular los valores de *theta* y *tau* (tamaño poblacional y grado de divergencia) a partir del tipo de datos o secuencias que el porgrama recibe. Esta opción funciona cuando no desconoces sobre parámetros genéticos y poblacionales de tus especies. 


```
speciesdelimitation = 0 * fixed species tree
speciestree = 0 * speciestree pSlider ExpandRatio ShrinkRatio
```
  

#### Opción A01: Estima el árbol de especies

Esta opción solo permite estimar un árbol de especies y se asume la asignación de individuos a especies y la delimitación de especies.
 
```
speciesdelimitation = 0 * fixed species tree
speciestree = 1 * speciestree pSlider ExpandRatio ShrinkRatio
```


#### Opción A10: Estima la delimitación de especies a partir de un árbol de especies guía.

Esta opción permite estimar la delimitación de especies utilizando las Cadenas Markovianas con salto reversible (rjMCMC) para muestrear o visitar ditintos esquemas de delimitación. Las rjMCMC se mueven entre varios modelos de delimitació a partir de un árbol guía. El árbol guía se considera un *prior*.


```
speciesdelimitation = 1 * fixed species tree
speciestree = 0 * speciestree pSlider ExpandRatio ShrinkRatio
```


#### Opción A11: Estima la delimitación de especies sin árbol de especies guía.

Esta opción no necesita un árbol guía pero tampoco se genera un árbol guía a partir de este análisis. 



```
speciesdelimitation = 1 * fixed species tree
speciestree = 1 * speciestree pSlider ExpandRatio ShrinkRatio
```


## **3.** Generar archivo control (archivo.ctl) 

El archvio control se puede generar de manera fácil una vez que cuentas con el archivo de secuencias.txt. Para generar el archivo vamos a utilizar el programa en línea [minimalist BPP](https://brannala.github.io/bpps/#/) desarrollado por B. Rannala. Minimalist BPP solo necesita el archivo de las secuencias para crear el individual_Imap.txt y el arcivo.ctl.

Interface gráfica de Minimalist BPP

![Minimalist BPP](https://github.com/JR-Montes/Tutorial_BPP_Delimitation/blob/main/Captura%20de%20pantalla%202021-01-31%20a%20la(s)%2022.52.02.png)



#### **Paso 1** Cargar archivo de secuencias 

![Step 1](https://github.com/JR-Montes/Tutorial_BPP_Delimitation/blob/main/Step_1.png)

#### **Paso 2** Asignar individuos a especies

![Step 2](https://github.com/JR-Montes/Tutorial_BPP_Delimitation/blob/main/Step_2.png)


**Nota**: Aquí debes asignar el epíteto de las especies a los individuos. Ejemplo, cembroides a cem. Eso significa que el programa va a asignar el nombre cembroides a todos aquellos individuos que inicien con las lestras cem en tu lista de individuos





Ejemplo del archivo control generado por Minimalist BPP.

```
	seed = -1 
    
    seqfile = cem_opt_20_bpp.txt * secuencias
    Imapfile = cem_opt_20_bppmap.txt * asignación de individuos a especies
    outfile = out.txt 
    mcmcfile = mcmc.txt 

    speciesdelimitation = 1 1 2 1 * delimitación de especies rjMCMC algoritmo_1 finetune (a m)
    speciestree = 0 * árbol de especies estimado
    
    speciesmodelprior = 1  * 0:  Verosimilitud uniforme; 1: árboles enraizados balanceados; 2: uniforme SLH; 3: uniformeSEnraizados
    
    species&tree = 6  aristata  balfouriana  bungeana  californiarum  krempfii lambertiana  * especies
                      5  6  4  2  2  1    * individuos por especie
                    (((californiarum,(aristata, balfouriana)),((bungeana),(krempfii,lambertiana)))); *árbol estimado por BPP
    
    usedata = 1  * 0: no usar datos (prior); 1:seq like
    
    cleandata = 0 * remover sitios con ambiguedades (1:si, 0:no)?

 
    nloci = 20
    thetaprior = 3  0.0075 e * invgamma(a, b) para theta
    tauprior = 3  1.5.  * invgamma(a, b) para tau & Dirichlet (a) for other tau's
    finetune = 1: 0.02 0.02 0.02 0.02 0.02 0.02 0.02
    
    print = 1 0 0 0 * MCMC, locusrate, heredityscalars, Genetrees
    burnin = 5000 * porcentaje de datos eliminados
    sampfreq = 2  * muestreo de frecuencias
    nsample = 500000 * número de generaciones
    threads = 14 1 1 * numero de núcleos

```
 

#### I. Se recomienda leer el [manual](https://hal.archives-ouvertes.fr/hal-02536475/document) para entender esta opción:

```
speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted
```

#### II. Se recomienda leer el [manual](https://hal.archives-ouvertes.fr/hal-02536475/document) para entender esta opción:

```
finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

```



## **4.** Correr el análisis

Una vez que se ha generado el archivo control y los has descargado, ahora puedes ejecutar el programa. Ejecuta la siguiente línea de comando en el directorio donde vive el archivo control, las secuencas y el archivo de asignación:

```
~ ./bpp --cfile cem_opt_20_bpp.ctl & 

```
  
  
**Nota**: Se recomiendo correr una réplica por cada análisis para corroborar la estabilidad y convergencia de las MCMC (Ver [manual](https://hal.archives-ouvertes.fr/hal-02536475/document)). También se recomiendo cambiar los valores de *thera* y *tau* para valorar la sensibilidad de BPP con tus datos. Si las hipótesis no cambian entonces puede confiar en los resultados de BPP. 


## **5.** Archivos de salida

Al final BPP genera dos archivos de salida. El archivo de las cadenas y el archivo con los resultados de los parámetros elejidos.


**NOTA**: El archivo MCMC solo se puede leer en el programa Tracer cuando elegiste las opción A00 pero para el resto de opciones no se puede leer ese archivo. 


Ejemplo del archivo MCMC:

```
Gen     np      tree                    theta_1         theta_2          theta_3        theta_6         theta_7         theta_8         theta_10
2       67      1111111111111111111111  0.001477        0.001907        0.001768        0.002835        0.003550        0.001290        0.00$
4       67      1111111111111111111111  0.001459        0.001884        0.001747        0.002801        0.003507        0.001274        0.00$
6       67      1111111111111111111111  0.001485        0.001918        0.001778        0.002852        0.003570        0.001297        0.00$
8       67      1111111111111111111111  0.001485        0.001918        0.001778        0.002852        0.003570        0.001297        0.00$
10      67      1111111111111111111111  0.001485        0.001918        0.001778        0.002852        0.003570        0.001297        0.00$
12      67      1111111111111111111111  0.001485        0.001918        0.001778        0.002852        0.003570        0.001297        0.00$
14      67      1111111111111111111111  0.001485        0.001952        0.001778        0.002852        0.003570        0.001297        0.00$
16      67      1111111111111111111111  0.001485        0.001610        0.001510        0.002852        0.003570        0.001297        0.00$
18      67      1111111111111111111111  0.001485        0.001610        0.001510        0.002852        0.003570        0.001297        0.00$
20      67      1111111111111111111111  0.001485        0.001610        0.001510        0.002852        0.003570        0.001297        0.00$
22      67      1111111111111111111111  0.001476        0.001600        0.001535        0.002833        0.003548        0.001246        0.00$
24      67      1111111111111111111111  0.001476        0.001600        0.001535        0.002833        0.003548        0.001246        0.00$
26      67      1111111111111111111111  0.001476        0.001618        0.001535        0.002833        0.003548        0.001246        0.00$
28      67      1111111111111111111111  0.001462        0.001603        0.001521        0.002807        0.003514        0.001234        0.00$
30      67      1111111111111111111111  0.001486        0.001629        0.001546        0.002852        0.003571        0.001254        0.00$

```

Ejemplo del archivo OUT:

```
  maximartineziipinceanarzedowskii
  maximartineziipinceana
  longaevaaristatabalfouriana
  longaevaaristata
  bungeanagerardianakrempfiilambertiana
  bungeanagerardiana
  krempfiilambertiana

Guide tree with posterior probability for presence of nodes:
(((((((((johannis, discolor)#1.000000, culminicola)#1.000000, remota)#1.000000, ((cembroides, orizabensis)#1.000000, lagunae)#1.000000)#1.000000, ((((monophylla, quadrifolia)#1.000000, californiarum)#1.000000, fallax)#1.000000, edulis)#1.000000)#1.000000, ((maximartinezii, pinceana)#1.000000, rzedowskii)#1.000000)#1.000000, ((longaeva, aristata)#1.000000, balfouriana)#1.000000)#1.000000, nelsonii)#1.000000, ((bungeana, gerardiana)#1.000000, (krempfii, lambertiana)#1.000000)#1.000000)#1.000000;;

```

## **6.** Análisis de los resultados


Para analizar los resultados recomendamos leer el [manual](https://hal.archives-ouvertes.fr/hal-02536475/document) de BPP
  


