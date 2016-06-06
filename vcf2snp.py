#!/usr/bin/python

from __future__ import print_function 
import sys
sys.path.append('./')
import kang
import pandas as pd
import numpy as np
from tqdm import tqdm 
import textwrap
# Reference index
print('indexing references..')


file_ref_fa     = '/ref/analysis/juntaehwan/ref/Cannuum/Pepper.v.1.55.total.chr.fa.seeders.fa'
file_gff       	= '/ref/analysis/juntaehwan/ref/Cannuum/Pepper_1.55.gene_models.gff3'
file_ref_annotation = '/ref/analysis/juntaehwan/ref/Pepper.v.1.55.proteins.annotated.fasta'
list_vcf       = ['/ref/analysis/juntaehwan/data/YCM334_samtools.raw.vcf','/ref/analysis/juntaehwan/data/Taeahn_samtools.raw.vcf']
                  #'/ref/analysis/juntaehwan/data/Perennial.SRR2751913.cv.vcf','/ref/analysis/juntaehwan/data/Dempsey.SRR2751914.cv.sambamba.vcf']
list_vcf_label = ['YCM334','TAEAHN']#,'PERENN','DEMPSE'] 


dic_annotation_fa   = kang.Fasta2dic_all(file_ref_annotation)
dic_ref_fa      = kang.Fasta2dic(file_ref_fa)
dic_annot = {}
for line in dic_annotation_fa.keys():
    cell = line.split()
    try:
        dic_annot[cell[0]] = ' '.join(cell[1:])
    except IndexError:
        dic_annot[cell[0]] = 'None'
print('gff parsing')
df_gff         	= pd.read_csv(file_gff,sep='\t',header=None)
mask        	= (df_gff[2] == 'CDS')
df_gff_cds  	= df_gff[mask]
mask        	= (df_gff[2] == 'gene')
df_gff_gene 	= df_gff[mask]
df_gff_cds['ID'] 	= df_gff_cds[8].apply(lambda x : x.split(';')[0].split(':')[1])
df_gff_gene['ID'] 	= df_gff_gene[8].apply(lambda x : x.split(';')[0].replace('ID=',''))
df_gff_gene_ix 		= df_gff_gene.set_index('ID')
df_gff_cds_ix 		= df_gff_cds.set_index('ID')
# Reference index done


# definition vcf2snpdic  
dic_base     = {'A':1,'T':2,'G':3,'C':4}
dic_base_rev = {1:'A',2:'T',3:'G',4:'C'}
def get_snp_dic(file_vcf):
    dic = {}
    for line in tqdm(open(file_vcf)):
        #print line
        if line[0] == '#':
            continue
        cell       = line.strip().split('\t')
        if cell[7][0:5] == 'INDEL':
            continue
        if cell[9].split(':')[0] == '1/1': # Only homovar catch! 
            pass
        else : continue 
        chromosome = cell[0]
        intloc     = int(cell[1]) - 1
        var        = cell[4].split(',')[0]
        info       = cell[7]
        fqual      = float(cell[5])
        dicinfo    = dict(zip([x.split('=')[0] for x in info.split(';')],[x.split('=')[1] for x in info.split(';')]))
        try:
            if dic_ref_fa[chromosome]:
                pass
        except KeyError:
            chromosome = 'Pepper1.55'+chromosome.split('.')[-1].replace('r','') # in the case of reference_fa's headers are not consistent with gff file. 
        try:
            dic[chromosome][intloc] = dic_base[var]
        except KeyError:
            print(chromosome)
            dic[chromosome]         = np.zeros(len(dic_ref_fa[chromosome]),dtype=np.int8)
            dic[chromosome][intloc] = dic_base[var]
    return dic

def seq_comp(list_label,list_seq,Outfile,linecut=100):
    list_array_seq = [np.array(list(x)) for x in list_seq]
    mask = np.zeros(len(list_seq[0]),dtype=bool)
    for n,array_seq in enumerate(list_array_seq):
        if n == len(list_array_seq)-1:
            continue
        array_seq2 = list_array_seq[n+1]
        mask += (array_seq != array_seq2)
    list_wrap = [textwrap.wrap(x,linecut) for x in list_seq]
    c  = textwrap.wrap(''.join(['.' if x == False else '*' for x in mask]),linecut)
    #print(list_label,list_wrap)
    #print(mask.nonzero()[0])
    for n,line in enumerate(list_wrap[0]):
        for m,seq in enumerate(list_wrap):
            print (list_label[m],seq[n],file=Outfile)
        print (' '*len(list_label[0]),c[n],file=Outfile)
    return(mask.nonzero()[0])




# VCF index
print('vcf to snp dictionary..')
list_chromosome = dic_ref_fa.keys()
list_dic_snp = [get_snp_dic(vcf) for vcf in list_vcf]

# CDS array indexing
print('cds indexing')
dic_genename_cds = {}
for genename in tqdm(df_gff_cds_ix.index):
    df = df_gff_cds_ix.loc[genename]
    if str(type(df)) == "<class 'pandas.core.frame.DataFrame'>":
        df.reset_index(inplace=True)
        chromosome = df[0][0]
        strand     = df[6][0]
        array = np.zeros(len(dic_ref_fa[chromosome]))
        for i in df.index:
            left       = df.loc[i][3]
            right      = df.loc[i][4]
            array[left-1:right] = 1
    else:
        chromosome = df[0]
        strand     = df[6]
        array = np.zeros(len(dic_ref_fa[chromosome]))
        left       = df[3]
        right      = df[4]
        array[left-1:right] = 1
    dic_genename_cds[genename] = (chromosome,strand, array)

# SNPs gathering
print('start SNP gathering')
Outfile_nonsynalign = open('pepper.allgenes.non-syn.snp.align.txt','w')
Outfile_snp_pos     = open('pepper.allgenes.snp.pos.txt','w')
for genename in tqdm(dic_genename_cds):
    g_left  = df_gff_gene_ix.loc[genename][3]
    g_right = df_gff_gene_ix.loc[genename][4]
    
    chromosome,strand, cds_array = dic_genename_cds[genename]
    
    matrix_snp      = np.matrix([dic_snp[chromosome][g_left-1:g_right] for dic_snp in list_dic_snp[0:2]])
    snp_mask = (np.sum(matrix_snp,axis=0) > 0)
    
    # SNP genomic postions
    for pos in snp_mask.nonzero()[1]:
        refpos  = pos+g_left-1 # zerobase
        refbase = dic_ref_fa[chromosome][refpos]
        bases   = [refbase if x == 0 else dic_base_rev[x] for x in np.array(matrix_snp[:,pos].T)[0]]
        context = dic_ref_fa[chromosome][refpos-20:refpos] +'(%s)'%('/'.join(bases)) + dic_ref_fa[chromosome][refpos+1:refpos+21]
        if len(set(bases)) == 1:
            print ('Genomic SNP ref',chromosome,pos+g_left,genename,strand,refbase,context,sep='\t',file=Outfile_snp_pos)
        else:
            print ('Genomic SNP among',chromosome,pos+g_left,genename,strand,refbase,context,sep='\t',file=Outfile_snp_pos)
    # SNP cds positions
    cds_mask = (cds_array[g_left-1:g_right]>0)
    cds_snp_mask = (snp_mask & cds_mask)
    for pos in cds_snp_mask.nonzero()[1]:
        refpos  = pos+g_left-1 # zerobase
        refbase = dic_ref_fa[chromosome][refpos]
        bases = [refbase if x == 0 else dic_base_rev[x] for x in np.array(matrix_snp[:,pos].T)[0]]
        context = dic_ref_fa[chromosome][refpos-20:refpos] +'(%s)'%('/'.join(bases)) + dic_ref_fa[chromosome][refpos+1:refpos+21]
        if len(set(bases)) == 1:
            print ('CDS SNP ref',chromosome,pos+g_left,genename,strand,refbase,context,sep='\t',file=Outfile_snp_pos)
        else:
            print ('CDS SNP among',chromosome,pos+g_left,genename,strand,refbase,context,sep='\t',file=Outfile_snp_pos)
    
    mask = cds_snp_mask
    
    if len(mask.nonzero()[0]) > 0:
        chr_seq     = dic_ref_fa[chromosome]
        r,c = matrix_snp.shape
        list_cds_seq_snp = [] # SNP. ... . vcf.. .. .... .. 
        for n,snp_array in enumerate(matrix_snp):
            chr_seq_snp  = chr_seq[g_left-1:g_right]
            cds_snp_mask_acc = cds_mask & (snp_array > 0)
            for each in cds_snp_mask_acc.nonzero()[1]:
                #print(each)
                chr_seq_snp = chr_seq_snp[:each] + dic_base_rev[matrix_snp[n,each]] + chr_seq_snp[each+1:]
            cds_seq_snp = ''
            for base in cds_mask.nonzero()[0]:
                cds_seq_snp += chr_seq_snp[base]
            list_cds_seq_snp.append(cds_seq_snp)
            
        
        
        if strand == '+':
            list_cds_seq_snp_pep = []
            for cds_seq_snp in list_cds_seq_snp:
                cds_seq_snp_pep = kang.translation(cds_seq_snp)
                list_cds_seq_snp_pep.append(cds_seq_snp_pep)
            
            if len(set(list_cds_seq_snp_pep)) > 1: # ... .. .. .... ... 1. ... 
                pass
            else: continue
        else:
            list_cds_seq_snp_pep = []
            for cds_seq_snp in list_cds_seq_snp:
                cds_seq_snp_pep = kang.translation(kang.rev_comp(cds_seq_snp))
                list_cds_seq_snp_pep.append(cds_seq_snp_pep)
            
            if len(set(list_cds_seq_snp_pep)) > 1: # ... .. .. .... ... 1. ... 
                pass
            else: continue
            
        print ('#',genename,strand,dic_annot[genename],file=Outfile_nonsynalign)
        print ('# pep_seq',file=Outfile_nonsynalign)
        pos_pepvar_array = seq_comp(list_vcf_label,list_cds_seq_snp_pep,Outfile_nonsynalign)
        print ('# cds_seq',file=Outfile_nonsynalign)
        if strand == '+':
            pos_snp_array = seq_comp(list_vcf_label,list_cds_seq_snp,Outfile_nonsynalign)
        else:
            pos_snp_array = seq_comp(list_vcf_label,[kang.rev_comp(x) for x in list_cds_seq_snp],Outfile_nonsynalign)
        for n,pos_snp in enumerate(pos_snp_array):
            if pos_snp/3 in pos_pepvar_array:
                if strand == '+':
                    pos = cds_snp_mask.nonzero()[1][n]
                    realpos = g_left + pos
                    refpos  = realpos-1 # zerobase
                    refbase = dic_ref_fa[chromosome][refpos]
                    bases = [refbase if x == 0 else dic_base_rev[x] for x in np.array(matrix_snp[:,pos].T)[0]]
                    context = dic_ref_fa[chromosome][refpos-20:refpos] +'(%s)'%('/'.join(bases)) + dic_ref_fa[chromosome][refpos+1:refpos+21]
                    if len(set(bases)) == 1:
                        print ('nonsyn SNP ref',chromosome, realpos,genename,strand,refbase,context,sep='\t',file=Outfile_snp_pos)
                    else:
                        print ('nonsyn SNP among',chromosome, realpos,genename,strand,refbase,context,sep='\t',file=Outfile_snp_pos)
                else:
                    pos = cds_snp_mask.nonzero()[1][-(n+1)]
                    realpos = g_left + pos
                    refpos  = realpos-1 # zerobase
                    refbase = dic_ref_fa[chromosome][refpos]
                    bases = [refbase if x == 0 else dic_base_rev[x] for x in np.array(matrix_snp[:,pos].T)[0]]
                    context = dic_ref_fa[chromosome][refpos-20:refpos] +'(%s)'%('/'.join(bases)) + dic_ref_fa[chromosome][refpos+1:refpos+21]
                    if len(set(bases)) == 1:
                        print ('nonsyn SNP ref',chromosome, realpos,genename,strand,refbase,context,sep='\t',file=Outfile_snp_pos)
                    else:
                        print ('nonsyn SNP among',chromosome, realpos,genename,strand,refbase,context,sep='\t',file=Outfile_snp_pos)
print('#Done',file=Outfile_nonsynalign)
print('#Done',file=Outfile_snp_pos)
Outfile_nonsynalign.close()
Outfile_snp_pos.close()
