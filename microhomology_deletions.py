#python microhomology deletions
import pandas as pd
from bgreference import hg38
import sys
#input: VCF of deletions
#take breakpoint coordinates
#also extract deletion identity, ie ref
vcf_path = sys.argv[1] #needs to be minimised vcf, ie no INFO
name = sys.argv[2]
def deletion_up_down(vcf_path):
    vcf = pd.read_table(vcf_path)
    vcf.columns = ['chr', 'bkpt', 'ID', 'ref', 'alt', '0', 'filt']
    vcf.drop(['ID', '0', 'filt',], axis=1, inplace=True)
    vcf.dropna(inplace=True)
    vcf['length'] = vcf['ref'].astype(str).map(len)-vcf['alt'].astype(str).map(len)
    expanded = vcf
    expanded['length'] = vcf['length'].astype(int)
    expanded['up_coor_start'] = vcf['bkpt'].astype(int) - vcf['length'].astype(int)
    expanded['down_coor_start'] = vcf['bkpt'].astype(int)
    expanded['up_coor_start'] = expanded['up_coor_start']-1
    expanded['down_coor_start'] = expanded['down_coor_start'].astype(int)
    expanded['up_seq'] =  expanded.apply(lambda x: hg38(x['chr'], x['up_coor_start']-1,x['length']+1), axis=1)#needs +1 to have correct lengths of sequence
    expanded['down_seq'] = expanded.apply(lambda x: hg38(x['chr'], x['down_coor_start'],x['length']+1), axis=1)
    expanded['del'] = expanded['chr']
    expanded['compare_seq'] = expanded['del']
    for i in range(len(expanded['down_seq'] )):
        print(['up', expanded['up_seq'].iloc[i]])
        print(['down', expanded['down_seq'].iloc[i]])
        print(['insert', expanded['ref'].iloc[i]])
        if re.search('N',expanded['down_seq'].iloc[i])is not None or re.search('N',expanded['up_seq'].iloc[i]) is not None:
            expanded.loc[i, 'del']= 'N region'
            expanded.loc[i, 'compare_seq'] = 'na'
        elif expanded['down_seq'].iloc[i] == expanded['ref'].iloc[i]: #ie if down flank is duplicated
            #add to output dataframe
            print("duplication")
            expanded.loc[i, 'del'] = 'del'
            expanded.loc[i, 'compare_seq'] = hg38(expanded.loc[i, 'chr'], expanded.loc[i, 'down_coor_start'] + expanded.loc[i, 'length'], expanded.loc[i, 'length'])
        elif str(''.join(reversed(expanded['up_seq'].iloc[i]))) == expanded['ref'].iloc[i]: #ie if up or down of the breakpoint in the reference is the deleted region
            expanded.loc[i, 'del'] = 'del'  ##NOTE- if upstream, then the sequence for comparison is the REVERSE of the sequence
            expanded.loc[i, 'compare_seq'] = ''.join(reversed(hg38(expanded.loc[i, 'chr'], expanded.loc[i, 'up_coor_start'] - expanded.loc[i, 'length'], expanded.loc[i, 'length'])))
        else:
            expanded.loc[i, 'del']= 'non'
            expanded.loc[i, 'compare_seq'] = 'na' #comapre seq is the sequence which we then use to comapre for microhomology. nb not sure whether i need to look at both for this?
    print(expanded)
    dups_df = expanded.loc[expanded['del'] == 'del']
    #returns dataframe with compare seq that is the sequence flanking the side of the deleted region away from the breakpoint
    return(dups_df)

deletions = deletion_up_down(vcf_path)

hom=[0]*20
pos_index=[0]*20
for j in range( len(deletions['bkpt'])):
    compare_seq = deletions['compare_seq'].iloc[j]
    length = deletions['length'].iloc[j]
    ref = deletions['ref'].iloc[j]
    for i in range(1, length):
        pos_index[i] +=1
        if ref[i] == compare_seq[i]:
            hom[i]+=1


deletion_hom_df = pd.DataFrame(zip(hom, pos_index), columns=['hom', 'count'])
deletion_hom_df['rate_hom'] = deletion_hom_df['hom']/deletion_hom_df['count']
deletion_hom_df.drop(0, axis=0, inplace=True)
deletion_hom_df.to_csv(str(name + '_5-20bp_deletions_microhom.csv'))
