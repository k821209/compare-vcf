# compare-vcf
```python
#!/usr/bin/python
from __future__ import print_function
import sys
sys.path.append('./')
import kang
import pandas as pd
import numpy as np
from tqdm import tqdm
import textwrap
```

This code is for retrieve SNP in gene regions between two vcf files 

- genic SNP
- CDS SNP
- non-syn/syn SNP
- Alignment text visualization


Result 1
```
# CA05g07370 + Amidophosphoribosyltransferase
# pep_seq
YCM334 TVSVHNKVLNLATGIGLFSDVFNQSNLDKLPGDMAIGHVRIPVEGCDVVIEVPDSGVVATLGYAAKSGVPYQQGLIRSHYVGRTFIEPSQKIRNFGVKLK
TAEAHN MVSVHNKVLNLATGIGLVSDVFNQSNLDKLPGDMAIGHVRIPVEGCDVVIEVPDSGVVATLGYAAKSGVPYQQGLIRSHYVGRTFIEPSQKIRNFGVKLK
       *................*..................................................................................
YCM334 LSPVKAVLEGKTVVVVDDSIVRGTTSSKIVRLLKEAGAKEVHMRIPSPPLITSCYYGIDTPSTEELISNSMNVEEIGKFIGADSLAFLPTDSFCKLLSSD
TAEAHN LSPVKAVLEGKTVVVVDDSIVRGTTSSKIVRLLKEAGAKEVHMRIPSPPLITSCYYGIDTPSTEELISNSMNVEEIGKFIGADSLAFLPTDSFCKLLSSD
       ....................................................................................................
YCM334 YTNFCYACFSGRYPVYPTTAMEECIDDKMMSAVPAN*
TAEAHN YTNFCYACFSGRYPVYPTTAMEECIDDKMMSAVPAN*
       .....................................
# cds_seq
YCM334 ACGGTTTCTGTCCACAATAAAGTTCTCAATTTGGCTACAGGTATTGGGTTGTTCTCTGATGTCTTCAACCAGTCGAATCTTGATAAACTTCCAGGTGACA
TAEAHN ATGGTTTCTGTCCACAATAAAGTTCTCAATTTGGCTACAGGTATTGGGTTGGTCTCTGATGTCTTCAACCAGTCGAATCTTGATAAACTTCCAGGTGACA
       .*.................................................*................................................
YCM334 TGGCTATTGGTCATGTCAGAATCCCCGTGGAGGGCTGTGATGTTGTGATAGAAGTACCAGATTCAGGTGTTGTAGCAACACTTGGCTATGCAGCCAAATC
TAEAHN TGGCTATTGGTCATGTCAGAATCCCCGTGGAGGGCTGTGATGTTGTGATAGAAGTACCAGATTCAGGTGTTGTAGCAACACTTGGCTATGCAGCCAAATC
       ....................................................................................................
YCM334 AGGAGTACCATATCAACAAGGCTTAATAAGGTCACACTATGTTGGCAGGACATTTATCGAGCCATCCCAGAAGATTAGGAATTTTGGGGTGAAGCTTAAG
TAEAHN AGGAGTACCATATCAACAAGGCTTAATAAGGTCACACTATGTTGGCAGGACATTTATCGAGCCATCCCAGAAGATTAGGAATTTTGGGGTGAAGCTTAAG
       ....................................................................................................
YCM334 CTCTCACCGGTTAAAGCAGTATTGGAAGGGAAAACAGTGGTAGTTGTGGATGATTCAATCGTCAGAGGGACGACATCGTCAAAGATTGTGAGGCTGTTAA
TAEAHN CTCTCACCGGTTAAAGCAGTATTGGAAGGGAAAACAGTGGTAGTTGTGGATGATTCAATCGTCAGAGGGACGACATCGTCAAAGATTGTGAGGCTGTTAA
       ....................................................................................................
YCM334 AGGAGGCAGGTGCTAAAGAGGTCCACATGAGGATTCCAAGCCCTCCACTTATAACTTCTTGTTACTATGGGATAGATACTCCTAGTACAGAGGAATTGAT
TAEAHN AGGAGGCAGGTGCTAAAGAGGTCCACATGAGGATTCCAAGCCCTCCACTTATAACTTCTTGTTACTATGGGATAGATACTCCTAGTACAGAGGAATTGAT
       ....................................................................................................
YCM334 ATCTAACAGCATGAACGTGGAGGAGATTGGGAAGTTTATTGGGGCGGATTCTCTAGCCTTTCTACCAACTGATAGCTTTTGTAAGCTATTAAGCAGTGAT
TAEAHN ATCTAACAGCATGAACGTGGAGGAGATTGGGAAGTTTATTGGGGCGGATTCTCTAGCCTTTCTACCAACTGATAGCTTTTGTAAGCTATTAAGCAGTGAT
       ....................................................................................................
YCM334 TATACAAACTTTTGTTATGCTTGCTTTTCAGGTAGGTATCCAGTATACCCAACCACAGCTATGGAGGAGTGTATAGATGATAAAATGATGTCAGCTGTAC
TAEAHN TATACAAACTTTTGTTATGCTTGCTTTTCAGGTAGGTATCCAGTATACCCAACCACAGCTATGGAGGAGTGTATAGATGATAAAATGATGTCAGCTGTAC
       ....................................................................................................
YCM334 CTGCTAATTAA
TAEAHN CTGCTAATTAA
...
```

Result 2

```
Genomic SNP ref Pepper1.55ch02  170238955       CA02g30750      -       A       TAAGGCCTCCCTCCATAACC(G/G)CCTAGCCACCTTACAATGCA
Genomic SNP among       Pepper1.55ch03  244580286       CA03g29430      +       C       GATATCATCATCAATCTACC(C/T)ATTTTCGTGGATTATAGGTC
Genomic SNP among       Pepper1.55ch03  244582159       CA03g29430      +       A       AGAAGAATGGGCAAAGTGCT(A/G)CACATATTTCATGATTAAAG
Genomic SNP among       Pepper1.55ch03  244583070       CA03g29430      +       T       TTACTTGAACAAAAAAGTGT(T/G)TTCAGTCATAGTAACTGTAT
Genomic SNP among       Pepper1.55ch03  244583192       CA03g29430      +       A       AAAAGAGGAGAAGAAATAAG(A/G)ATACTTCATTAAAGTTTATT
Genomic SNP among       Pepper1.55ch03  244583408       CA03g29430      +       C       ACGAAAAGGTGGTTTATTTG(C/T)TTAGACTTTAGATGGGTGTC
Genomic SNP among       Pepper1.55ch05  36994031        CA05g07370      +       T       CTGTCAGGAAGGTGTTGGCA(C/T)GGTTTCTGTCCACAATAAAG
Genomic SNP among       Pepper1.55ch05  36994081        CA05g07370      +       G       TGGCTACAGGTATTGGGTTG(T/G)TCTCTGATGTCTTCAACCAG
Genomic SNP ref Pepper1.55ch05  36994606        CA05g07370      +       A       TCCCTGAGCCTAAATCTTGC(G/G)TCTTTGAGCACATTTACTTT
Genomic SNP ref Pepper1.55ch05  36995093        CA05g07370      +       A       GAGGAATTGATATCTAACAG(C/C)ATGAACGTGGAGGAGATTGG
Genomic SNP ref Pepper1.55ch05  36995128        CA05g07370      +       T       GATTGGGAAGTTTATTGGGG(C/C)GGATTCTCTAGCCTTTCTAC
CDS SNP among   Pepper1.55ch05  36994031        CA05g07370      +       T       CTGTCAGGAAGGTGTTGGCA(C/T)GGTTTCTGTCCACAATAAAG
CDS SNP among   Pepper1.55ch05  36994081        CA05g07370      +       G       TGGCTACAGGTATTGGGTTG(T/G)TCTCTGATGTCTTCAACCAG
CDS SNP ref     Pepper1.55ch05  36995093        CA05g07370      +       A       GAGGAATTGATATCTAACAG(C/C)ATGAACGTGGAGGAGATTGG
CDS SNP ref     Pepper1.55ch05  36995128        CA05g07370      +       T       GATTGGGAAGTTTATTGGGG(C/C)GGATTCTCTAGCCTTTCTAC
nonsyn SNP among        Pepper1.55ch05  36994031        CA05g07370      +       T       CTGTCAGGAAGGTGTTGGCA(C/T)GGTTTCTGTCCACAATAAAG
nonsyn SNP among        Pepper1.55ch05  36994081        CA05g07370      +       G       TGGCTACAGGTATTGGGTTG(T/G)TCTCTGATGTCTTCAACCAG
```
