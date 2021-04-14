import argparse
import re
from scipy.stats import poisson
import subprocess
import os
import generate_site
import peaks_mapping
import random
import scipy.stats
import numpy as np
from pyfasta import Fasta
import time


parser=argparse.ArgumentParser('')
parser.add_argument('--t1',help='The first input sam file')
parser.add_argument('--t2',help='The second input sam file')
parser.add_argument('--pre',help='the output prefix',required=True)
parser.add_argument('-a',help='the total mismatches and insertion/deletion allowed',default='2')
parser.add_argument('--genome',help='the whole genome size. default the human whole genome size',default='3088286401')
parser.add_argument('--pvalue',help='the p_value, before corrected',default='0.001')
parser.add_argument('--gap',help='the gap size allowed between two consecutive tag to form a site',default='1')
parser.add_argument('--weight',help='the weight for ambiguous reads')
parser.add_argument('--ref_path',help='the path to the reference genome')
parser.add_argument('--codepath',help='the path for the code')
parser.add_argument('--flanking',help='the sequence length up/downstream of the central pointer of a site, which is used to calculate kmer frequency. Totally 2f bp sequence would be used then',default='200')
parser.add_argument('--threshold',help='the threshold',default='0.01')
parser.add_argument('--fraction',help='the fraction of change between every two iters',default='0.1')
parser.add_argument('--max_iter',help='the maximum iteration number',default='10')
parser.add_argument('--num_bin',help='the number of bins to grid the prior',type=int,default=100)
parser.add_argument('-p', help = 'if or not output the probability an ambiguous read to its possible mapping positions', action = 'store_true')
args=parser.parse_args()



####get the parameters
rep1_input=args.t1
rep2_input=args.t2
prefix=args.pre
mismatch_allowed=args.a
genome=args.genome
pvalue=args.pvalue
gap=args.gap
code_path=args.codepath
flanking=args.flanking
weight=args.weight
ref_path=args.ref_path
threshold=args.threshold
fraction=args.fraction
max_iter=args.max_iter
num_bin = args.num_bin



output_path='/'.join(prefix.split('/')[:-1])
output_path='./' if output_path=='' else output_path

########The code path (make sure its format is /a/b/c/###
code_path=code_path.split('/')
code_path=[x for x in code_path if x]
code_path='/'+'/'.join(code_path)+'/' if code_path[0]!='.' else './'



####delete the existing kmer group file
ls_content=subprocess.check_output(['ls',output_path])
ls_content=ls_content.decode('utf-8')
kmer_file=prefix+'_kmer_group'
kmer_file_name=kmer_file.split('/')[-1]
if re.search(kmer_file_name,ls_content):
   os.system('rm '+prefix+'_kmer_group')

####the first filter on the two replicates
os.system('python '+code_path+'filter.py -i '+rep1_input+' -a '+mismatch_allowed+' --pre '+prefix+'1')
os.system('python '+code_path+'filter.py -i '+rep2_input+' -a '+mismatch_allowed+' --pre '+prefix+'2')

####delete pcr duplications
os.system('python '+code_path+'pcr.py --pre '+prefix+'1 --genome '+genome+' --pvalue '+pvalue)
os.system('python '+code_path+'pcr.py --pre '+prefix+'2 --genome '+genome+' --pvalue '+pvalue)

####calculate the tag size
uniq_read_file=prefix+'1_uniq_reads'
f=open(uniq_read_file)
line=next(f)
line=line.strip().split('\t')
label_cigar=line[5]
label_number=[int(x) for x in re.split(r'[MDI]',label_cigar) if x]
label_character=[x for x in re.split('[0-9]',label_cigar) if x]
tag_length=0
for i in range(len(label_number)):
   if label_character[i]!='D':
      tag_length+=label_number[i]


tag_length=str(tag_length)

####run macs
os.system('python '+code_path+'macs.py --pre '+prefix+'1 --tagsize '+tag_length+' --genome '+genome)
os.system('python '+code_path+'macs.py --pre '+prefix+'2 --tagsize '+tag_length+' --genome '+genome)

####get overlap peaks of the two replicates
os.system('python '+code_path+'overlap_of_two_reps.py --pre '+prefix)

####calculate median(peak length), use 0.5median as site size
peak_file=prefix+'_overlap_peaks'
peak_size=[]
f=open(peak_file)
for line in f:
   line=line.strip().split('\t')
   start=int(line[1])
   end=int(line[2])
   peak_size.append(end-start)
f.close()

site_size=int(np.median(peak_size)/2)


####get the flanking region size
flanking_size=str(int(float(flanking.split('site')[0])*site_size)) if re.search('site',flanking) else str(int(flanking))


####print out some intermediate files
os.system('python '+code_path+'print_out_median_file.py --tagsize '+tag_length+' --pre '+prefix+'1 --weight '+weight+' --gap '+gap+' --maxsite '+str(site_size))
os.system('python '+code_path+'print_out_median_file.py --tagsize '+tag_length+' --pre '+prefix+'2 --weight '+weight+' --gap '+gap+' --maxsite '+str(site_size))


####gather two sample files
os.system('cat '+prefix+'1_sample >'+prefix+'_sample')
os.system('cat '+prefix+'2_sample >>'+prefix+'_sample')


####kmer frequency within and outside peak sites
os.system('python '+code_path+'kmer_group.py --flanking '+flanking_size+' --pre '+prefix+' --ref_path '+ref_path)

####get the discriminative kmers
os.system('python '+code_path+'discrimitive_kmer.py --pre '+prefix)


####use kmeans to get kmer groups
os.system('python '+code_path+'kmeans.py --pre '+prefix+' --codepath '+code_path)


####feature values for sample file
os.system('python '+code_path+'logistic_regression_prepare_for_sample.py --flanking '+flanking_size+' --pre '+prefix+'1 --ref_path '+ref_path+' --weight '+weight)

os.system('python '+code_path+'logistic_regression_prepare_for_sample.py --flanking '+flanking_size+' --pre '+prefix+'2 --ref_path '+ref_path+' --weight '+weight)

####call logistic regression
os.system('Rscript '+code_path+'logist.R '+prefix+'1')
os.system('Rscript '+code_path+'logist.R '+prefix+'2')


####calculate the prior values for all the sites
os.system('python '+code_path+'logistic_regression_for_all_site.py  --flanking '+flanking_size+' --pre '+prefix+'1 --ref_path '+ref_path+' --weight '+weight)

os.system('python '+code_path+'logistic_regression_for_all_site.py  --flanking '+flanking_size+' --pre '+prefix+'2 --ref_path '+ref_path+' --weight '+weight)


####run Gibbs sampling
os.system('python '+code_path+'gibbs.py --pre '+prefix+'1 --threshold '+str(threshold)+' --fraction '+str(fraction)+' --max_iter '+str(max_iter)+' --num_bin '+str(num_bin)+' -p ')

os.system('python '+code_path+'gibbs.py --pre '+prefix+'2 --threshold '+str(threshold)+' --fraction '+str(fraction)+' --max_iter '+str(max_iter)+' --num_bin '+str(num_bin)+' -p ')



####turn bed to sam format
os.system('python '+code_path+'site_to_sentence.py -i '+rep1_input+' --pre '+prefix+'1')
os.system('python '+code_path+'site_to_sentence.py -i '+rep2_input+' --pre '+prefix+'2')




###cat the unique file to the reads file
os.system('cat '+prefix+'1_uniq_reads >> '+prefix+'1_reads')
os.system('cat '+prefix+'2_uniq_reads >> '+prefix+'2_reads')
