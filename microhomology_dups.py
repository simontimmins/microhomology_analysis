#microhomology INSERTION DUPLICATION quantification in python:
import re
import pandas as pd
from bgreference import hg38
import sys
#input: VCF of insertions
#take breakpoint coordinates
#also extract insertion identity
#
#caluclate deletion length, therefore can work out endpoint of deletion (i.e. len(REF)-len(ALT))
#expand e.g. 10bp either side from startpoint and endpoint of deletion
#use bgreference to obtain sequence for these on either side
#thus obtain a datafrmae with:
    #chr, start, end. 10bp up flank, 10bp down flank, deleted sequence, deleted sequence length. 

#on this dataframe, define a function which cpunts the number of identical/complementary nucelotides between the deleted sequence (REF) and the flanks. 
#summing this across genomic positions to identify microhomology patters. 
vcf_path = sys.argv[1] #needs to be minimised vcf
name = sys.argv[2]

#determine if duplication


def duplication_up_down(vcf_path):
    vcf = pd.read_table(vcf_path)
    vcf.columns = ['chr', 'bkpt', 'ID', 'ref', 'alt', '0', 'filt']
    vcf.drop(['ID', '0', 'filt',], axis=1, inplace=True)
    vcf.dropna(inplace=True)
    vcf['length'] = vcf['alt'].astype(str).map(len)-vcf['ref'].astype(str).map(len)
    expanded = vcf
    expanded['length'] = vcf['length'].astype(int)
    expanded['up_coor_start'] = vcf['bkpt'].astype(int) - vcf['length'].astype(int)
    expanded['down_coor_start'] = vcf['bkpt'].astype(int)
    expanded['up_coor_start'] = expanded['up_coor_start']-1
    expanded['down_coor_start'] = expanded['down_coor_start'].astype(int)
    expanded['up_seq'] =  expanded.apply(lambda x: hg38(x['chr'], x['up_coor_start']-1,x['length']+1), axis=1)#needs +1 to have correct lengths of sequence
    expanded['down_seq'] = expanded.apply(lambda x: hg38(x['chr'], x['down_coor_start'],x['length']+1), axis=1)
    expanded['dup'] = expanded['chr']
    expanded['compare_seq'] = expanded['dup']
    for i in range(len(expanded['down_seq'] )):
        print(['up', expanded['up_seq'].iloc[i]])
        print(['down', expanded['down_seq'].iloc[i]])
        print(['insert', expanded['alt'].iloc[i]])
	#remove N regions
        if re.search('N',expanded['down_seq'].iloc[i])is not None or re.search('N',expanded['up_seq'].iloc[i]) is not None:
            expanded.loc[i, 'dup']= 'N region'
            expanded.loc[i, 'compare_seq'] = 'na'
        elif expanded['down_seq'].iloc[i] == expanded['alt'].iloc[i]: #ie if down flank is duplicated
            #add to output dataframe
            print("duplication")
            expanded.loc[i, 'dup'] = 'TD'
            expanded.loc[i, 'compare_seq'] = hg38(expanded.loc[i, 'chr'], expanded.loc[i, 'down_coor_start'] + expanded.loc[i, 'length'], expanded.loc[i, 'length'])
        elif str(''.join(reversed(expanded['up_seq'].iloc[i]))) == expanded['alt'].iloc[i]: #ie if up flank is duplicated
            expanded.loc[i, 'dup'] = 'TD'
            expanded.loc[i, 'compare_seq'] = ''.join(reversed(hg38(expanded.loc[i, 'chr'], expanded.loc[i, 'up_coor_start'] - expanded.loc[i, 'length'], expanded.loc[i, 'length'])))
        else:
            expanded.loc[i, 'dup']= 'non-TD'
            expanded.loc[i, 'compare_seq'] = 'na' #comapre seq is the sequence which we then use to comapre for microhomology. 
    print(expanded)
    dups_df = expanded.loc[expanded['dup'] == 'TD']
    return(dups_df)


duplications = duplication_up_down(vcf_path)


#combine horizontally to make dataframe, 4 by 20

#theh need to find a way in which to apply this to every single line
hom=[0]*20
pos_index=[0]*20
for j in range( len(duplications['bkpt'])):
    compare_seq = duplications['compare_seq'].iloc[j]
    length = duplications['length'].iloc[j]
    alt = duplications['alt'].iloc[j]
    for i in range(1, length):
        pos_index[i] +=1
        if alt[i] == compare_seq[i]:
            hom[i]+=1

df = pd.DataFrame(zip(hom, pos_index), columns=['hom', 'count'])
df['rate_hom'] = df['hom']/df['count']
df.drop(0, axis=0, inplace=True)
df.to_csv(name)
