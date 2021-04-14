import argparse
import generate_site
import peaks_mapping
import random
import scipy.stats
import numpy as np
import time


parser=argparse.ArgumentParser('')
parser.add_argument('--tagsize',help='tag size',required=True,type=int)
parser.add_argument('--gap',help='the gap size allowed between two consecutive tag to form a site',default=1,type=int)
parser.add_argument('--maxsite',help='maximum length of a site',default=300,type=int)
#parser.add_argument('--maxdensity',help='in gamma distribution, when the tag density exceeds this value, it would be change to this maximun density value',default=30,type=int)
parser.add_argument('--pre',help='the output prefix',required=True)
parser.add_argument('--weight',help='the weight for ambiguous reads',required=True,type=float)
args=parser.parse_args()


random.seed(int(time.time()))


prefix=args.pre
tag_size=args.tagsize
gap=args.gap
max_site_length=args.maxsite
#max_density=args.maxdensity
weight=args.weight

peak_file=prefix[:-1]+'_overlap_peaks'
all_read_file=prefix+'_all_reads'
within_peak_sites=prefix+'_within_peak_sites'
outside_peak_sites=prefix+'_outside_peak_sites'
ambiguous_file=prefix+'_ambiguous'
sample_file=prefix+'_sample'
ps_pn_file=prefix+'_ps_pn'
uniq_read_ps_pn_file=prefix+'_uniq_ps_pn'  ##one new output





chromosome=['chr'+str(i) for i in list(range(1,23))+['X','Y']]

a=generate_site.generate_site(all_read_file)
tag,ambiguous_link=a._store_file()
site,all_reads,ambiguous_link=a._generating_site(tag,ambiguous_link,tag_size,gap,max(105,max_site_length),weight)


a=peaks_mapping.file_store(peak_file)
peaks=a._record_file()

a=peaks_mapping.map_site_to_peak()
label0,label1=a._map_site_to_peak(peaks,site)

########smooth-->normalization-->new mean, variance--> new gamma
def fit_gamma(index,real):
   smoothed_real=[real[0]]
   for i in range(1,len(real)-1):
      tem=[real[i-1],real[i],real[i+1]]
      if real[i-1]==0 and real[i+1]==0:
         smoothed_real.append(real[i])
      else:
         smoothed_real.append(np.median(tem))
   smoothed_real.append(real[-1])
   smoothed_real=[i/sum(smoothed_real) for i in smoothed_real]
###calculate mean
   mean_value=0
   for i in range(len(smoothed_real)):
      mean_value=mean_value+(index[i]+0.5)*smoothed_real[i]
###calculate variance
   variance=0
   for i in range(len(smoothed_real)):
      variance=variance+smoothed_real[i]*((index[i]+0.5)**2)
   variance=variance-mean_value**2
###refit the gamma distribution
   lam=mean_value/variance
   alp=(mean_value**2)/variance
   fitted_gamma={}
   #for i in range(len(index)):
   for i in index:
      fitted_gamma[i]=scipy.stats.gamma.pdf(i,a=alp,scale=1.0/lam)
      #fitted_gamma.append(scipy.stats.gamma.pdf(i,a=alp,scale=1.0/lam))
   return fitted_gamma

#####



def cal_gap(site_size):
   gap=1.0
   density=1.0/site_size/2*100
   round_decimal=0
   while density<1:
      gap/=10
      density*=10
      round_decimal+=1
   return gap,round_decimal






###seperate the data: input1 (only contains 0), input2 (others)
def seperate_data(input_list):
   input1=[x for x in input_list if x==0]
   input2=[x for x in input_list if x!=0]
   return input1,input2




####get the density gap and round decimal
density_gap,round_decimal=cal_gap(max_site_length)

if max_site_length<100:
   density_gap=1
   round_decimal=0


###get the sample files (to calculate discriminative kmers)
len_label0=len(label0)
len_label1=len(label1)

###select all of the within_peak sites, and the same number of outside_peak sites
sampling_num=int(len_label1)
sample0_index=random.sample(range(1,len_label0+1),sampling_num)
'''
for i in range(sampling_num):
   number0=random.randint(1,len_label0)
   while (number0 in sample0_index):
      number0=random.randint(1,len_label0)
   sample0_index.append(number0)
'''


sample0_index=sorted(sample0_index)


ws=open(sample_file,'w')
w=open(within_peak_sites,'w')
within_density=[]
uniq_within_density=[]
for every_site in label1:
   chrom=every_site[0]
   start=every_site[1]
   end=every_site[2]
   position=(start,end)
   all_read_count=all_reads[(chrom,start,end)]
   uniq_read_count=site[chrom][position]
   density=round(float(all_read_count)/(end-start)*100,round_decimal)
   uniq_density=round(float(uniq_read_count)/(end-start)*100,round_decimal)
   within_density.append(density)
   uniq_within_density.append(uniq_density)
   w.writelines(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+str(all_read_count)+'\t'+str(uniq_read_count)+'\t1\n')
   ws.writelines(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+str(all_read_count)+'\t'+str(uniq_read_count)+'\t1\n')
   
w.close()

   
w=open(outside_peak_sites,'w')
pointer=0
line_index=1
outside_density=[]
uniq_outside_density=[]
for every_site in label0:
   chrom=every_site[0]
   start=every_site[1]
   end=every_site[2]
   position=(start,end)
   all_read_count=all_reads[(chrom,start,end)]
   uniq_read_count=site[chrom][position]
   density=round(float(all_read_count)/(end-start)*100,round_decimal)
   uniq_density=round(float(uniq_read_count)/(end-start)*100,round_decimal)
   outside_density.append(density)
   uniq_outside_density.append(uniq_density)
   w.writelines(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+str(all_read_count)+'\t'+str(uniq_read_count)+'\t0\n')
   if (pointer<sampling_num):
      if (line_index==sample0_index[pointer]):
         ws.writelines(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+str(all_read_count)+'\t'+str(uniq_read_count)+'\t0\n')
         pointer+=1
   line_index+=1 

w.close()
ws.close()

###output format: chrom site_start site_end all_reads uniq_reads site_label


####print out the ambiguous information
w=open(ambiguous_file,'w')
for index in ambiguous_link:
   for every_site in ambiguous_link[index]:
      for every_possible_read in ambiguous_link[index][every_site]:
         site_chrom=every_site[0]
         site_start=every_site[1]
         site_end=every_site[2]
         tag_start=every_possible_read[1]
         tag_end=every_possible_read[2]
         w.writelines(index+'\t'+site_chrom+'\t'+str(site_start)+'\t'+str(site_end)+'\t'+str(tag_start)+'\t'+str(tag_end)+'\n')

w.close()



##############################
####first get the max density
max_within_density=np.percentile(within_density,95)
max_uniq_within_density=np.percentile(uniq_within_density,95)
max_outside_density=np.percentile(outside_density,95)
max_uniq_outside_density=np.percentile(uniq_outside_density,95)


max_density=round(max([max_within_density,max_uniq_within_density,max_outside_density,max_uniq_outside_density]),round_decimal)
#max_density=int(max([max_within_density,max_uniq_within_density,max_outside_density,max_uniq_outside_density]))

#####modify the density list based on the max_density
within_density=[min(x,max_density) for x in within_density]
uniq_within_density=[min(x,max_density) for x in uniq_within_density]
outside_density=[min(x,max_density) for x in outside_density]
uniq_outside_density=[min(x,max_density) for x in uniq_outside_density]




###for ps
ps_prob=[]
index=[round(x*density_gap,round_decimal) for x in range(int(max_density/density_gap)+1)]
total_ps=float(len(within_density))
for i in index:
   ps_prob.append(within_density.count(i)/total_ps)

fitted_ps=fit_gamma(index,ps_prob)

###for pn
pn_prob=[]
total_pn=float(len(outside_density))
for i in index:
   pn_prob.append(outside_density.count(i)/total_pn)

fitted_pn=fit_gamma(index,pn_prob)
   

####output ps,pn
w=open(ps_pn_file,'w')
for i in index:
   w.writelines(str(i)+'\t'+str(fitted_ps[i])+'\t'+str(fitted_pn[i])+'\n')

w.close() 


###for unique ps
###first seperate the data
unique_ps0,unique_ps1=seperate_data(uniq_within_density)
###calculate the data for density=0
ps0_density=len(unique_ps0)/total_ps
###calculate the distribution for density!=0 part
total_ps1=float(len(unique_ps1))
unique_ps1_prob=[]
for i in index:
   unique_ps1_prob.append(unique_ps1.count(i)/total_ps1)

fitted_ps1=fit_gamma(index,unique_ps1_prob)

###combine the two parts together
fitted_uniq_ps={}
for i in index:
   fitted_uniq_ps[i]=fitted_ps1[i]*(1-ps0_density)


fitted_uniq_ps[index[0]]=ps0_density
#fitted_uniq_ps=[(1-ps0_density)*i for i in fitted_ps1]
#fitted_uniq_ps[0]+=ps0_density


###for unique pn
###first seperate the data
unique_pn0,unique_pn1=seperate_data(uniq_outside_density)
###calculate the data for density=0
pn0_density=len(unique_pn0)/total_pn
###calculate the distribution for density!=0 part
total_pn1=float(len(unique_pn1))
unique_pn1_prob=[]
for i in index:
   unique_pn1_prob.append(unique_pn1.count(i)/total_pn1)

fitted_pn1=fit_gamma(index,unique_pn1_prob)

###combine the two parts together
fitted_uniq_pn={}
for i in index:
   fitted_uniq_pn[i]=fitted_pn1[i]*(1-pn0_density)

fitted_uniq_pn[index[0]]=pn0_density
   

#fitted_uniq_pn=[(1-pn0_density)*i for i in fitted_pn1]
#fitted_uniq_pn[0]+=pn0_density


w=open(uniq_read_ps_pn_file,'w')
for i in index:
   w.writelines(str(i)+'\t'+str(fitted_uniq_ps[i])+'\t'+str(fitted_uniq_pn[i])+'\n')

w.close()




'''
####print out the ps_pn value
w=open(ps_pn_file,'w')
within_mean=np.mean(within_density)
within_variance=np.std(within_density)**2
within_lam=within_mean/within_variance
within_alp=within_mean**2/within_variance
outside_mean=np.mean(outside_density)
outside_variance=np.std(outside_density)**2
outside_lam=outside_mean/outside_variance
outside_alp=outside_mean**2/outside_variance
###print out the ps, pn value  (output format: density ps pn)
for i in range(max_density+1):
   within_likelihood=scipy.stats.gamma.pdf(i,a=within_alp,scale=1./within_lam)
   outside_likelihood=scipy.stats.gamma.pdf(i,a=outside_alp,scale=1./outside_lam)
   w.writelines(str(i)+'\t'+str(within_likelihood)+'\t'+str(outside_likelihood)+'\n')

w.close() 


########print out the ps_pn value for unique reads
w=open(uniq_read_ps_pn_file,'w')
within_mean=np.mean(uniq_within_density)
within_variance=np.std(uniq_within_density)**2
within_lam=within_mean/within_variance
within_alp=within_mean**2/within_variance
outside_mean=np.mean(uniq_outside_density)
outside_variance=np.std(uniq_outside_density)**2
outside_lam=outside_mean/outside_variance
outside_alp=outside_mean**2/outside_variance
###print out the ps, pn value  (output format: density ps pn)
for i in range(max_density+1):
   within_likelihood=scipy.stats.gamma.pdf(i,a=within_alp,scale=1./within_lam)
   outside_likelihood=scipy.stats.gamma.pdf(i,a=outside_alp,scale=1./outside_lam)
   w.writelines(str(i)+'\t'+str(within_likelihood)+'\t'+str(outside_likelihood)+'\n')

w.close()
'''


print('finish: print_out_median_file.py\n')
