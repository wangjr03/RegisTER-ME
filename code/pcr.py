##2. tag store (forward and reverse strand respectively) 
##3. do the filter (according to the max tolerate
##score comparation: the reads map to the same 

import argparse
import re
from scipy.stats import poisson
import subprocess
import os




parser=argparse.ArgumentParser()
parser.add_argument('--pvalue',help='the p_value, before corrected',type=float,default=0.001)
parser.add_argument('--genome',help='the whole genome size. default the human whole genome size',type=int,default=3088286401)
parser.add_argument('--pre',help='the output prefix',required=True)
args=parser.parse_args()


chrom_size=args.genome
prefix=args.pre
pvalue=args.pvalue

input_file=prefix+'_first_filter'
output1=prefix+'_all_reads'
output2=prefix+'_uniq_reads'


###sort the input file first
os.system('sort -k 1 '+input_file+' -o' +input_file)



###calculate the number of unique reads
f=open(input_file)
label={}
uniq_tag=0
line=next(f)
line=line.split('\t')
index=line[0]
mark='uniq'

###print out the ambiguous tags
for line in f:
   line=line.split('\t')
   if ( line[0]==index):
      label[index]='ambiguous'
      mark='ambiguous'
   else:
      label[index]=mark
      if (mark=='ambiguous'):
         mark='uniq'
      else:
         uniq_tag+=1
   index=line[0]


label[index]=mark
if mark=='uniq':
   uniq_tag+=1

f.close()



###calculate the max num of tags allowed
mu=float(uniq_tag)/chrom_size
pois_score=poisson.ppf(1-pvalue/chrom_size,mu)  ###the max number of tags mapped to a same position






##2. score[chrom][start]
forward={}
reverse={}
score_forward={}
score_reverse={}
chromosome=['chr'+str(i) for i in list(range(1,23))+['X','Y']]
for i in chromosome:
   forward[i]={}
   reverse[i]={}
   score_forward[i]={}
   score_reverse[i]={}

f=open(input_file)
for line in f:
   line=line.split('\t')
   index=line[0]
   chrom=line[2]
   start=int(line[3])
   strand=int(line[1])
   score=int(line[11].split(':')[-1])
   if (label[index]=='ambiguous'):
      continue
   else:
      if ( (strand/16)%2 ):
         try:
            if ( len(score_reverse[chrom][start])< pois_score ):
               (reverse[chrom][start]).append(index)
               score_reverse[chrom][start].append(score)
            else:
               min_score=min(score_reverse[chrom][start])
               if ( score> min_score ):
                  position=(score_reverse[chrom][start]).index(min_score)
                  score_reverse[chrom][start][position]=score
                  reverse[chrom][start][position]=index
         except:
            reverse[chrom][start]=[index]
            score_reverse[chrom][start]=[score]
      else:
         try:
            if ( len(score_forward[chrom][start])<pois_score ):
               score_forward[chrom][start].append(score)
               (forward[chrom][start]).append(index)
            else:
               min_score=min(score_forward[chrom][start])
               if ( score> min_score ):
                  position=(score_forward[chrom][start]).index(min_score)
                  score_forward[chrom][start][position]=score
                  forward[chrom][start][position]=index
         except:
            forward[chrom][start]=[index]
            score_forward[chrom][start]=[score]
f.close()


###free the memory for 2 dictionary
del score_forward
del score_reverse


######################################
###print out all the uniq reads
f=open(input_file)
w1=open(output1,'w')
w2=open(output2,'w')
for line in f:
   tem=line.split('\t')
   index=tem[0]
   chrom=tem[2]
   start=int(tem[3])
   if (label[index]=='ambiguous'):
      w1.writelines(index+'\t'+chrom+'\t'+str(start)+'\n')
   else:
      try:
         if (index in forward[chrom][start]):  ##in forward chain
            w1.writelines(index+'\t'+chrom+'\t'+str(start)+'\n')
            w2.writelines(line)
      except:
         try:
            if (index in reverse[chrom][start]): ##in reverse chain
               w1.writelines(index+'\t'+chrom+'\t'+str(start)+'\n')
               w2.writelines(line)
         except:      ##filtered out
            pass
         




f.close()
w1.close()
w2.close()



#########################################

print('finish: pcr.py\n')
