import argparse
import subprocess
import os
import re

parser=argparse.ArgumentParser('')
parser.add_argument('--tagsize',help='tag size, the batch file should calculate this at the begining and use it in later scripts',required=True,type=int)
parser.add_argument('--genome',help='the whole genome size. default the human whole genome size',type=int,default=3088286401)
parser.add_argument('--pre',help='the output prefix',required=True)
args=parser.parse_args()


##the input file always refers to the original bowtie output
prefix=args.pre

uniq_tag_file=prefix+'_uniq_reads'
peak_file=prefix+'_peaks.bed'
peak_file_name=peak_file.split('/')[-1]


genome_size=args.genome
tag_size=args.tagsize

output_path='/'.join(prefix.split('/')[:-1])
output_path='.' if output_path=='' else output_path

os.system('macs -t '+uniq_tag_file+' --tsize '+str(tag_size)+' --format SAM --gsize '+str(genome_size)+' --name '+prefix)

ls_content=subprocess.check_output(['ls',output_path])
ls_content=ls_content.decode('utf-8')
#if not (re.search(peak_file_name,ls_content)):
#   os.system('macs -t '+uniq_tag_file+' --tsize '+str(tag_size)+' --format SAM --gsize '+str(genome_size)+' --name '+prefix+' --mfold 10')
   
if not (re.search(peak_file_name,ls_content)):
   os.system('macs -t '+uniq_tag_file+' --tsize '+str(tag_size)+' --format SAM --gsize '+str(genome_size)+' --name '+prefix+' --mfold 5,30')


print('print: macs.py\n')
