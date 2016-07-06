#!/usr/bin/python
# Genomic SNP ref Pepper1.55ch02  170238955       CA02g30750      -       A       TAAGGCCTCCCTCCATAACC(G/G)CCTAGCCACCTTACAATGCA
# Genomic SNP among       Pepper1.55ch03  244580286       CA03g29430      +       C       GATATCATCATCAATCTACC(C/T)ATTTTCGTGGATTATAGGTC
import pandas as pd
from tqdm import tqdm  
file_in = 'pepper.allgenes.snp.pos.txt'


dic = {'genename':[],
       'num_genic_snp':[],
       'num_cds_snp':[],
       'num_nonsyn_snp':[] }

for line in tqdm(open(file_in)):
    if line[0] == '#':
        continue
    cell = line.strip().split('\t')
    SNPtype  = cell[0]
    genename = cell[3]
    try:
        index = dic['genename'].index(genename)
        if SNPtype == 'Genomic SNP among':
            dic['num_genic_snp'][index]  += 1
        elif SNPtype == 'CDS SNP among':
            dic['num_cds_snp'][index]    += 1
        elif SNPtype == 'nonsyn SNP among':
            dic['num_nonsyn_snp'][index] += 1
    except ValueError:
        if SNPtype == 'Genomic SNP among':
            dic['genename'].append(genename)
            dic['num_genic_snp'].append(1)
            dic['num_cds_snp'].append(0)
            dic['num_nonsyn_snp'].append(0)
        elif SNPtype == 'CDS SNP among':
            dic['genename'].append(genename)
            dic['num_genic_snp'].append(0)
            dic['num_cds_snp'].append(1)
            dic['num_nonsyn_snp'].append(0)
        elif SNPtype == 'nonsyn SNP among':
            dic['genename'].append(genename)
            dic['num_genic_snp'].append(0)
            dic['num_cds_snp'].append(0)
            dic['num_nonsyn_snp'].append(1)
df = pd.DataFrame(dic)
df.to_csv(file_in+'.gene2snp.txt',sep='\t')
