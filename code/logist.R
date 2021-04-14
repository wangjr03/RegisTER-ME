###input the file name
args=commandArgs(T)
prefix=args[1]
sample_file=paste(prefix,'_sample_kmer_feature',sep='')
output_file=paste(prefix,'_prior',sep='')


a=read.table(file=sample_file,sep='\t',header=F)
n_col=ncol(a)
colnames(a)=c('chrom','start','stop','weighted_reads','uniq_reads','label',paste0('kmer',seq(1,n_col-6,1)))


a$length=a$stop-a$start+1
#a$density=a$weighted_reads/a$length*100
#b='label~density'
b='label~'
for (i in seq(1,n_col-6,1)) { b=paste0(b,'+kmer',i) } 
c=glm(formula=b,family=binomial,data=a)
#write.table(exp(c$coefficients),file=output)
write.table(c[1],file=output_file)
