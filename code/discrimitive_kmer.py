import argparse
import numpy as np
import scipy.stats


parser=argparse.ArgumentParser('')
parser.add_argument('--pre',help='the output prefix',required=True)
args=parser.parse_args()

###get the input file name and output file name
prefix=args.pre

input_file=prefix+'_kmer_frequency'
output_file=prefix+'_discriminative_kmer'


final={}
within={}
within_frequency={}
outside={}
outside_frequency={}
foldchange_threshold=1.5
p_value=0.001
k=6
kmer_type_count=4**k
binomial_p=p_value/kmer_type_count

f=open(input_file)
next(f)
for  line in f:
   line=line.strip().split('\t')
   kmer=line[0]
   within_frequency[kmer]=float(line[2])
   within[kmer]=int(line[1])
   outside_frequency[kmer]=float(line[4])
   outside[kmer]=int(line[3])

f.close()


within_median=np.median(list(within_frequency.values()))
within_total=sum(list(within.values()))
outside_median=np.median(list(outside_frequency.values()))
outside_total=sum(list(outside.values()))

for kmer in within.keys():
   foldchange=within_frequency[kmer]/outside_frequency[kmer]
   if foldchange > foldchange_threshold:
      p_value=scipy.stats.binom_test(within[kmer],within_total,outside_frequency[kmer],alternative='greater')
      if (within_frequency[kmer]> within_median and p_value<binomial_p):
         final[kmer]=foldchange
   elif (foldchange < (1/foldchange_threshold) ):
      p_value=scipy.stats.binom_test(within[kmer],within_total,outside_frequency[kmer],alternative='less')
      if (outside_frequency[kmer]>outside_median and p_value<binomial_p):
         final[kmer]=foldchange



w=open(output_file,'w')
for i in final:
   w.writelines(i+'\t'+str(final[i])+'\n')

w.close()


print('finish:discrimitive_kmer.py\n')

