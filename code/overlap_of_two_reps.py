import argparse
import os
import re

parser=argparse.ArgumentParser('This script is used to find out the overlap between two input files. Only lines from file1 would be kept')
parser.add_argument('--pre',help='the output prefix',required=True)
args=parser.parse_args()




####get the output file name (output file would be stored in current folder
prefix=args.pre
input_file1=prefix+'1_peaks.bed'
input_file2=prefix+'2_peaks.bed'
output_file=prefix+'_overlap_peaks'



###initialization
chromosome=['chr'+str(i) for i in list(range(1,23))+['X','Y']]

file1,file2={},{}
for chrom in chromosome:
   file1[chrom],file2[chrom]=[],[]


###store the peaks from both files
f=open(input_file1)
for line in f:
   line=line.split('\t')
   chrom=line[0]
   start=int(line[1])
   end=int(line[2])
   file1[chrom].append((start,end))


f.close()


f=open(input_file2)
for line in f:
   line=line.split('\t')
   chrom=line[0]
   start=int(line[1])
   end=int(line[2])
   file2[chrom].append((start,end))


f.close()



####do the overlapping
overlap={}
for chrom in chromosome:
   if len(file1[chrom])>0 and len(file2[chrom])>0:
      overlap[chrom]={}
      peak=sorted(file2[chrom])
      pointer=0
      total_peak=len(peak)
      for pos in sorted(file1[chrom]):
         start=pos[0]
         end=pos[1]
         if (end<peak[pointer][0]):
            overlap[chrom][(start,end)]=0
         elif (start>peak[pointer][1]):
            if (pointer==(total_peak-1)):
               overlap[chrom][(start,end)]=0
            else:
               while (start>peak[pointer][1] and pointer< (total_peak-1)):
                  pointer+=1
               if (end<peak[pointer][0]):
                  overlap[chrom][(start,end)]=0
               elif (start >peak[pointer][1] ):
                  overlap[chrom][(start,end)]=0
               else:
                  overlap[chrom][(start,end)]=1
         else: 
            overlap[chrom][(start,end)]=1
   

###output the result
f=open(input_file1)
w=open(output_file,'w')
for line in f:
   tem=line.split('\t')
   chrom=tem[0]
   start=int(tem[1])
   end=int(tem[2])
   if chrom in overlap.keys():
      if overlap[chrom][(start,end)]:
         w.writelines(line)


f.close()
w.close()



print('finish: overlap_of_two_reps.py\n')
