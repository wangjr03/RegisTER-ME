from pyfasta import Fasta
import re
import argparse

parser=argparse.ArgumentParser('')
parser.add_argument('--flanking',help='the sequence length up/downstream of the central pointer of a site, which is used to calculate kmer frequency. Totally 2f bp sequence would be used then',type=int,default=200) 
parser.add_argument('--pre',help='the output prefix',required=True)
parser.add_argument('--ref_path',help='the path to the reference genome',required=True)
args=parser.parse_args()


###The input file name
prefix=args.pre
flanking_length=args.flanking
k=6
ref_path=args.ref_path

sample_file=prefix+'_sample'
output_file=prefix+'_kmer_frequency'

###count the kmer frequency within and outside peaks
sequence=Fasta(ref_path)
#sequence=Fasta('/mnt/research/compbio/wanglab/huangbi4/k_mer/all.fasta')
within={}
outside={}



f=open(sample_file)
for line in f:
   line=line.strip().split('\t')
   chrom=line[0]
   start=int(line[1])
   end=int(line[2])
   middle=start+int((end-start)/2)
   seq=sequence[chrom][middle-flanking_length-1:middle+flanking_length]
   label=line[-1]
   if (label=='1'):
      for s in range(len(seq)-k+1):
         tem=seq[s:s+k]
         if ( re.search('N',tem)):
            continue
         else:
            try:
               within[tem]+=1
            except:
               within[tem]=1
   else:
      for s in range(len(seq)-k+1):
         tem=seq[s:s+k]
         if ( re.search('N',tem)):
            continue
         else:
            try:
               outside[tem]+=1
            except:
               outside[tem]=1



total_within=float(sum(within.values()))
total_outside=float(sum(outside.values()))
w=open(output_file,'w')
w.writelines('kmer\twithin_counts\twithin_freq\toutside_counts\toutside_freq\n')
for key in within.keys():
   if key in outside.keys():
      w.writelines(key+'\t'+str(within[key])+'\t'+str(within[key]/total_within)+'\t'+str(outside[key])+'\t'+str(outside[key]/total_outside)+'\n')
   else:
      pass


w.close()

print('finish:kmer_group.py\n')
