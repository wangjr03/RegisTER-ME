class file_store():
   def __init__(self,file_name):
      self.file_name=file_name
      self.chromosome=['chr'+str(i) for i in list(range(1,23))+['X','Y']]

   def _record_file(self):
      peaks={}
      for i in self.chromosome:
         peaks[i]=[]
      f=open(self.file_name)
      for line in f:
         line=line.split('\t')
         chrom=line[0]
         start=int(line[1])
         end=int(line[2])
         peaks[chrom].append([start,end])
      f.close()
      return peaks




class map_site_to_peak():
   def __init__(self):
      self.chromosome=['chr'+str(i) for i in list(range(1,23))+['X','Y']]

   def _map_site_to_peak(self,peaks,sites):
      label1=[]  ##within peak
      label0=[]  ##outside peak
      for chrom in sites:
         if peaks[chrom]==[] :
            for i in sites[chrom]: label0.append([chrom,i[0],i[1]])
         else:
            peak=sorted(peaks[chrom])
            site=sorted(sites[chrom])
            pointer=0
            len_peak=len(peak)
            for i in site:
               if ( i[1]<peak[pointer][0]):
                  label0.append([chrom,i[0],i[1]])
               elif (i[0] >peak[pointer][1]):
                  if (pointer==(len_peak-1)):
                     label0.append([chrom,i[0],i[1]])
                  else:
                     while (i[0]>peak[pointer][1] and pointer<(len_peak-1)):
                        pointer+=1
                     if (i[1]<peak[pointer][0]):
                        label0.append([chrom,i[0],i[1]])
                     elif (i[0] >peak[pointer][1] ):
                        label0.append([chrom,i[0],i[1]])
                     else:
                        label1.append([chrom,i[0],i[1]])
               else:
                  label1.append([chrom,i[0],i[1]])
      return label0,label1
      
   def _sampling_sites(self,input_array,size):
      import random
      import time
      random.seed(int(time.time()))
      end=len(input_array)-1
      sample=[] 
      position_record=[]
      for i in range(size):
         random_value=random.randint(0,end)
         while ( random_value in position_record):
            random_value=random.randint(0,end)
         position_record.append(random_value)
         sample.append(input_array[random_value])
         i+=1
      return sample
