args=commandArgs(T)
file=args[1]
k=as.numeric(args[2])


a=read.table(file,header=T,sep='\t')
length=min(nrow(a),ncol(a))
b=a[1:length,1:length]
col_name=colnames(b)


group=kmeans(b,k)
labels=group$cluster
percentage=group$betweenss/group$totss
print(percentage)

for (i in 1:k) {
   cat(paste('group',i))
   cat('\n')
   for (j in 1:length) {
      if (labels[j]==i) {
         cat(col_name[j])
         cat('\t')  }  }
   cat('\n\n\n')  }

