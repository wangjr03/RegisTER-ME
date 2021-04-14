import argparse
import re


parser=argparse.ArgumentParser()
parser.add_argument('-i',help='the input file',required=True)
parser.add_argument('-a',help='the total mismatches and insertion/deletion allowed',default=2,type=int)
parser.add_argument('--pre',help='the output path',required=True)
args=parser.parse_args()


prefix=args.pre
output=prefix+'_first_filter'




### Initialize the chromosome list
chromosome=['chr'+str(i) for i in list(range(1,23))+['X','Y']]


### begin filtering
f=open(args.i)
w=open(output,'w')

#store_line=[]
for line in f:
   if re.search('^@',line):      ### delete the header part
      continue
   else:
      tem=line.strip().split('\t')
      if tem[3]=='0':           ###delete the non-mapping part
         continue
      elif tem[2] not in chromosome:
         continue 
      else:
         total=0
         for i in tem[12:]:
            if re.search('XM',i):
               n_mismatch=int(i.split(':')[-1])
               break
         label_cigar=tem[5]
         label_number=[int(x) for x in re.split(r'[MDI]',label_cigar) if x]
         label=[x for x in re.split(r'[0-9]',label_cigar) if x]
         insert_n_delete=0
         for i in range(len(label)):
            if ( label[i]=='I' or label[i]=='D'):
               insert_n_delete+=label_number[i]
         if ( (n_mismatch+insert_n_delete)> args.a ):
            continue
         else:
            #store_line.append(line)
            w.writelines(line)
      

f.close()
w.close()


print('finish: filter.py\n')
