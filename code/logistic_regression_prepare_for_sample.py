from pyfasta import Fasta
import re
import argparse
import numpy as np



parser=argparse.ArgumentParser()
parser.add_argument('--flanking',help='the sequence length up/downstream of the central pointer of a site, which is used to calculate kmer frequency. Totally 2f bp sequence would be used then',type=int,default=200)
parser.add_argument('--weight',help='the weight for ambiguous reads',type=float,default=0)
parser.add_argument('--pre',help='the output path',required=True)
parser.add_argument('--ref_path',help='the path to the reference genome',required=True)
args=parser.parse_args()


prefix=args.pre
flanking_length=args.flanking
ambiguous_weight=args.weight
ref_path=args.ref_path


k=6
discriminative_file=prefix[:-1]+'_discriminative_kmer'
sample_file=prefix+'_sample'
output_file=prefix+'_sample_kmer_feature'
kmer_group_file=prefix[:-1]+'_kmer_group'
kmer_frequency_file=prefix[:-1]+'_kmer_frequency'


############### function to delete overlapping kmers
def calculate_occurance(array,k):
   if (len(array)==0):
      return 0
   else:
      tem=[array[0]]
      for i in array:
         if i>=(tem[-1]+k):
            tem.append(i)
      return len(tem)

#######################################################


####store kmer group information######
f=open(kmer_group_file)
group_foldchange={}
kmer_belonging={}
pointer=0
for line in f:
   if (re.search('[ACTG]',line)):
      pointer+=1
      group_foldchange[pointer]=[0,0]
      line=line.strip().split()
      for i in line:
         kmer_belonging[i]=pointer
f.close()
group_keys=group_foldchange.keys()

####store the frequency for each kmer group
f=open(kmer_frequency_file)
next(f)
within_count=0
outside_count=0
for line in f:
   line=line.split('\t')
   a_kmer=line[0]
   within_count+=int(line[1])
   outside_count+=int(line[3])
   if (a_kmer in kmer_belonging.keys()):
      group_foldchange[kmer_belonging[line[0]]][0]+=int(line[1])
      group_foldchange[kmer_belonging[line[0]]][1]+=int(line[3])

f.close()
foldchange=float(outside_count)/within_count
for i in group_keys:
   group_foldchange[i]=group_foldchange[i][0]*foldchange/group_foldchange[i][1]


#########################################################
   
######################################################################
sequence=Fasta(ref_path)
#sequence=Fasta('/mnt/research/compbio/wanglab/huangbi4/k_mer/all.fasta')
f=open(sample_file)
w=open(output_file,'w')
group_occurance={}
for line in f:
   tem=line.split('\t')
   chrom=tem[0]
   start=int(tem[1])
   end=int(tem[2])
   middle=start+int((end-start)/2)
   uniq_count=int(tem[4])
   weighted_count=float(tem[3])
   for i in group_keys:
      group_occurance[i]=[]
###go through the sequence
   seq=sequence[chrom][middle-flanking_length-1:middle+flanking_length]
   for s in range(len(seq)-k+1):
      tem=seq[s:s+k]
      try:
         group=kmer_belonging[tem]
         group_occurance[group].append(s)
      except:
         pass
   w.writelines(line.strip())
   for i in group_keys:
      group_occurance[i]=np.log10(group_foldchange[i]**(calculate_occurance(group_occurance[i],k)))
      w.writelines('\t'+str(group_occurance[i]))
   w.writelines('\n')
f.close()
w.close()


print('finish logistic_regression_prepare_for_sample.py\n')
