import random
import numpy as np
import generate_site
import argparse
import re
import time
from collections import Counter


parser=argparse.ArgumentParser('the gibbs sampling code for chromosome2')
parser.add_argument('--threshold',help='the threshold',type=float,default=0.01)
parser.add_argument('--fraction',help='the fraction of change between every two iters',type=float,default=0.1)
parser.add_argument('--max_iter',help='the maximum iteration number',default=10,type=int)
parser.add_argument('--pre',help='the prefix',required=True)
parser.add_argument('--num_bin',help='the number of bins to grid the prior',type=int,default=100)
parser.add_argument('-p', help = 'if or not output the probability an ambiguous read to its possible mapping positions', action = 'store_true')
args=parser.parse_args()




#####file name
prefix=args.pre
max_iteration=args.max_iter
threshold_we_set=args.threshold
fraction_change_we_set=args.fraction
bin_num=args.num_bin


final_read_file=prefix+'_reads'
ps_pn_file=prefix+'_ps_pn'
uniq_ps_pn_file=prefix+'_uniq_ps_pn'
all_sites_file=prefix+'_all_sites'
ambiguous_file=prefix+'_ambiguous'
prior_file=prefix+'_prior'
bin_region=1.0/bin_num

##############################
small_value=1e-15
#############################



######functions
def sampling_method(possibility):
   add_up=sum(possibility)
   cumulative_possibility=list(np.cumsum(possibility))
   value=random.uniform(0,add_up)
   index=cumulative_possibility.index(next(i for i in cumulative_possibility if i>value))
   return index



def count_probability(round_decial,max_density,length,counts,ps,pn,prior_value):
   density=round(float(counts+1)/length*100,round_decial)
   if (density>max_density):
      density=max_density
   probability=ps[density]*prior_value/(1-prior_value+small_value)
   density=round(float(counts)/length*100,round_decial)
   if (density>max_density):
      density=max_density
   probability/=pn[density]
   return probability



def modify_possibility(possibility):
   if sum(possibility)==0:
      possibility=[x+1 for x in possibility]
   return possibility



def binned_prior(prior_value,bin_region):
   pos=int(prior_value/bin_region)
   pos=(pos+0.5)*bin_region
   if pos>1:
      pos-=bin_region
   return pos

###########################################################3


###ps, pn part
ps,pn={},{}
density_list=[]
with open(ps_pn_file) as f:
   for line in f:
      line=line.strip().split('\t')
      index=float(line[0])
      density_list.append(index)
      ps[index]=float(line[1])
      pn[index]=float(line[2])
      max_density=float(line[0])
   



density_gap=density_list[1]-density_list[0]
round_decimal=int(-np.log10(density_gap))
ps[density_list[0]]=ps[density_list[1]]
pn[density_list[0]]=pn[density_list[1]]


####the ps, pn for unique reads
uniq_ps,uniq_pn={},{}
with open(uniq_ps_pn_file) as f:
   for line in f:
      line=line.strip().split('\t')
      index=float(line[0])
      uniq_ps[index]=float(line[1])
      uniq_pn[index]=float(line[2])





###prior, sites, kmer_values
chromosome=['chr'+str(i) for i in list(range(1,23))+['X','Y']]
time1=time.time()
prior,site={},{}
for chrom in chromosome:
   prior[chrom],site[chrom]={},{}

with open(all_sites_file) as f:
   for line in f:
      line=line.strip().split('\t')
      chrom=line[0]
      start=int(line[1])
      end=int(line[2])
      prior[chrom][(start,end)]=binned_prior(float(line[-1]),bin_region)
      site[chrom][(start,end)]=int(line[4])



######for ambiguous parts
ambiguous={}
ambiguous_link={}
with open(ambiguous_file) as f:
   for line in f:
      line=line.strip().split('\t')
      index=line[0]
      chrom=line[1]
      site_start=int(line[2])
      site_end=int(line[3])
      tag_start=int(line[4])
      tag_end=int(line[5])
      try:
         ambiguous[index].append((chrom,site_start,site_end))
      except:
         ambiguous[index]=[(chrom,site_start,site_end)]
      try:
         ambiguous_link[index][(chrom,site_start,site_end)]=(chrom,tag_start,tag_end)
      except:
         ambiguous_link[index]={(chrom,site_start,site_end):(chrom,tag_start,tag_end)}



time2=time.time()
print('time for input file and generate data structure:'+str(time2-time1)+'\n')


total_ambiguous=len(ambiguous.keys())       
final_mapping={}


time3=time.time()
threshold=1
iter=0
fraction_change=1
track_mapping={}
while (iter<=max_iteration    and   (threshold>threshold_we_set   or  fraction_change>fraction_change_we_set) ):
   iter+=1
   if (iter==1):
      for every_tag in ambiguous:
         track_mapping[every_tag]=[]
         possibility=[]
         for every_site in ambiguous[every_tag]:
            chrom=every_site[0]
            start=every_site[1]
            end=every_site[2]
            length_of_site=end-start+1
            possibility.append(count_probability(round_decimal,max_density,length_of_site,site[chrom][(start,end)],uniq_ps,uniq_pn,prior[chrom][(start,end)]))
         possibility=modify_possibility(possibility)  ###new added, in case 
         index=sampling_method(possibility)
         chrom=ambiguous[every_tag][index][0]
         start=ambiguous[every_tag][index][1]
         end=ambiguous[every_tag][index][2]
         val=[track_mapping[every_tag].append(ambiguous[every_tag][index]) for i in range(10)]
      time4=time.time()
      ####begin to map all the reads after the initialization step
      for every_tag in ambiguous:
         mapped_pos=track_mapping[every_tag][-1]
         chrom=mapped_pos[0]
         start=mapped_pos[1]
         end=mapped_pos[2]
         site[chrom][(start,end)]+=1 
   if (iter!=1):
      threshold_count=0
      fraction_count=0
      fraction_all=0
      for every_tag in ambiguous:
         possibility=[]
         val=track_mapping[every_tag].pop(0)                       ###pop out the one enters the first
         mapped_tag=track_mapping[every_tag][-1]              ###select the last mapped position
         site[mapped_tag[0]][(mapped_tag[1],mapped_tag[2])]-=1
         for every_site in ambiguous[every_tag]:
            chrom=every_site[0]
            start=every_site[1]
            end=every_site[2]
            length_of_site=end-start+1
            possibility.append(count_probability(round_decimal,max_density,length_of_site,site[chrom][(start,end)],ps,pn,prior[chrom][(start,end)]))
         possibility=modify_possibility(possibility)
         index=sampling_method(possibility)
         chrom=ambiguous[every_tag][index][0]
         start=ambiguous[every_tag][index][1]
         end=ambiguous[every_tag][index][2]
         site[chrom][(start,end)]+=1       ###this site++
         track_mapping[every_tag].append(ambiguous[every_tag][index])
         if (ambiguous[every_tag][index] != mapped_tag):
            threshold_count+=1
         if (mapped_tag!=track_mapping[every_tag][-3]):
            fraction_all+=1
            if (ambiguous[every_tag][index] != mapped_tag):
               fraction_count+=1   
      threshold=float(threshold_count)/total_ambiguous
      if fraction_all==0:
         fraction_change=0
      else:
         fraction_change=float(fraction_count)/fraction_all
      print('Iteration:'+str(iter)+'\n')
      print('Delta:'+str(threshold)+'\n')
      print('Delta(Delta):'+str(fraction_change)+'\n\n')




#if args.p: w = open(prefix+'_probability','w')
for index in track_mapping:
   site_info=max(set(track_mapping[index]),key=track_mapping[index].count)
   final_mapping[index]=site_info

if args.p:  
   with open(prefix+'_probability','w') as w:
      for index in track_mapping:
         pos_count = Counter(track_mapping[index])
         for a_site in pos_count:
            chrom, start, end = ambiguous_link[index][a_site]
            w.writelines('\t'.join([index,chrom,str(start),str(pos_count[a_site]*10)])+'\n')
   



time5=time.time()
print('time for all iteration:'+str(time5-time3)+'\n')



w=open(final_read_file,'w')

for index in final_mapping:
   site_info=final_mapping[index]
   site_info=ambiguous_link[index][site_info]
   w.writelines(site_info[0]+'\t'+str(site_info[1])+'\t'+index+'\n')


w.close()
