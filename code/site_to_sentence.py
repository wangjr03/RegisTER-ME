import argparse
import time
import re

current_time=time.time()
parser=argparse.ArgumentParser('this code can turn the site file into the original bowtie result file')
parser.add_argument('-i',help='the input file',required=True)
parser.add_argument('--pre',help='the output prefix',required=True)
args=parser.parse_args()




prefix=args.pre
source_file=args.i

input_file=prefix+'_reads'
outputfile=input_file

##store the site information
site={}
chromosome=['chr'+str(i) for i in list(range(1,23))+['X','Y']]
for chrom in chromosome:
   site[chrom]={}
f=open(input_file)
for line in  f:                           ###line[0]-->chrom, line[1]-->start, line[2]--> index
   line=line.strip().split('\t')
   site[line[0]][(line[2],line[1])]=1     ###site[chrom][(index,start)]=1
f.close()


##write out based on the source file 
f=open(source_file)
w=open(outputfile,'w')
for line in f:
   if ( re.search('^@',line)):
      continue
   else:
      tem=line.split('\t')
      index=tem[0]
      chrom=tem[2]
      start=tem[3]
      if (chrom in chromosome):
         if ( (index,start) in site[chrom]):
            w.writelines(line)

w.close()
f.close()


print('finish site_to_sentence.py\n')
