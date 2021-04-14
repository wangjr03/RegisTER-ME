import argparse
import os
import subprocess
import sys


parser=argparse.ArgumentParser('This code is used to write out 2 files, one containing the distance matrix for enriched kmers, one for depleted kmers')
parser.add_argument('--pre',help='the output prefix',required=True)
parser.add_argument('--codepath',help='the path for the code',required=True)
args=parser.parse_args()

################################################################
###define a function for both enriched and depleted group
def cat_kmer_group(kmer_file,kmer_list,output_file,distance,local_path):
   if len(kmer_list)<20:
      os.system('echo [1] 0 >>'+output_file)
      os.system('echo group1 >>'+output_file)
      sen=kmer_list[0]
      for i in range(1,len(kmer_list)):
         sen=sen+'\t'+kmer_list[i]
      #sen=sen+'\n'
      os.system('echo '+sen+'>>'+output_file)
   else: 
      w=open(kmer_file,'w')
      w.writelines(kmer_list[0])
      for i in range(1,len(kmer_list)):
         w.writelines('\t'+kmer_list[i])
      w.writelines('\n')
      for kmer1 in kmer_list:
         kmer2=kmer_list[0]
         w.writelines(str(distance[(kmer1,kmer2)]))
         for i in range(1,len(kmer_list)):
            kmer2=kmer_list[i]
            w.writelines('\t'+str(distance[(kmer1,kmer2)]))
         w.writelines('\n')
      w.close()
      group_number=int(len(kmer_list)/20)+1
      os.system('Rscript '+local_path+'kmeans.R '+kmer_file+' '+str(group_number)+'  >>'+output_file)
#################################################################3


local_path=args.codepath
prefix=args.pre

discriminative_kmer_file=prefix+'_discriminative_kmer'
enriched_output=prefix+'_enriched_kmer_matrix'
depleted_output=prefix+'_depleted_kmer_matrix'
kmer_group_file=prefix+'_kmer_group'



###record the 6mer 
kmer=[]
f=open(local_path+'6_mer.bed')
for line in f:
   kmer.append(line.strip())
f.close()


###record the distance between every 6mer
distance={}
pointer1=0
f=open(local_path+'6_mer_bitscore.matrix')
for line in f:
   line=line.strip().split('\t')
   kmer1=kmer[pointer1]
   for pointer2 in range(len(line)):
      kmer2=kmer[pointer2]
      distance[(kmer1,kmer2)]=1/float(line[pointer2])
   pointer1+=1

f.close()


###get the enriched kmer list and depleted kmer list
enriched=[]
depleted=[]
f=open(discriminative_kmer_file)
for line in f:
   line=line.strip().split('\t')
   fraction=float(line[-1])
   if (fraction>1):
      enriched.append(line[0])
   else:
      depleted.append(line[0])
f.close()



#######print out the kmer groups
if len(depleted)==0 and len(enriched)==0:
   print('No discriminative kmers\n')
   sys.exit()
else:
   if len(enriched)!=0:
      cat_kmer_group(enriched_output,enriched,kmer_group_file,distance,local_path)
   if len(depleted)!=0:
      cat_kmer_group(depleted_output,depleted,kmer_group_file,distance,local_path)
      

print('finish: kmeans.py\n')
