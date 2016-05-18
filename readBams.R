library(dnar)
library(ampcountr)
library(parallel)

bamFiles<-list.files('data','bam$',full.names=TRUE)

cover<-list()
for(ii in bamFiles){
	outFile<-sub('data/','work/',sub('bam$','cover',ii))
	sortPrefix<-sub('data/','work/',sub('\\.bam','_sort',ii))
	sortFile<-sprintf('%s.bam',sortPrefix)
	if(!file.exists(sortFile)){
		cmd<-sprintf('samtools sort -@ 6 -m 8G %s %s',ii,sortPrefix)
		message(cmd)
		system(cmd)
		cmd<-sprintf('samtools index %s',sortFile)
		message(cmd)
		system(cmd)
	}
	if(!file.exists(outFile)){
		cmd<-sprintf('~/installs/bedCount/bam2depth %s>%s',sortFile,outFile)
		message(cmd)
		system(cmd)
	}
	cover[[ii]]<-read.table(outFile)
	colnames(cover[[ii]])<-c('chr','pos','count')
}

ref<-read.fa('data/Pf3D7_v3.fasta')
ref$seq<-toupper(ref$seq)
primerData<-read.csv('data/ChosenSets.csv',stringsAsFactors=FALSE,row.names=1)
primerSets<-strsplit(primerData$primers,',')
names(primerSets)<-primerData$X_id
predictCover<-mclapply(primerSets,function(primers){
	fs<-gregexpr(paste(primers,collapse='|'),ref$seq)
	rs<-gregexpr(paste(revComp(primers),collapse='|'),ref$seq)
	rs<-lapply(rs,function(r)r+attr(r,'match.length')-1)
	out<-mcmapply(predictAmplifications,fs,rs,SIMPLIFY=FALSE,MoreArgs=list(maxLength=30000),mc.cores=length(primers))
	names(out)<-ref$name
	return(list('pred'=out,'for'=fs,'rev'=rs))
})

for(ii in names(predictCover)){
	message(ii)
	outFile<-file.path('out',sprintf('%s.png',ii))
	thisCovers<-cover[grep(sprintf('data/%s',ii),names(cover))]
	thisPredicts<-predictCover[[ii]][['pred']]
	png(outFile,width=3000,height=3000,res=250)
	par(mfrow=c(4,4))
	ylim<-c(0,max(sapply(thisCovers,function(xx)max(xx$count))))
	maxPredict<-max(sapply(thisPredicts,function(xx)max(xx$amplifications)))
	for(jj in names(thisPredicts)){
		thisCover<-lapply(thisCovers,function(x)x[x$chr==jj,])
		thisPredict<-thisPredicts[[jj]]
		plot(1,1,xlim=c(1,max(thisPredict$end)),ylim=ylim,main=jj,xlab='Position',ylab='Cover')
		sapply(thisCover,function(xx)points(xx$pos,xx$count))
		segments(thisPredict$start,thisPredict$amplifications/maxPredict*ylim[2],thisPredict$end,thisPredict$amplifications/maxPredict*ylim[2],col='#FF000099')
		segments(thisPredict$end[-nrow(thisPredict)],thisPredict$amplifications[-nrow(thisPredict)]/maxPredict*ylim[2],thisPredict$start[-1],thisPredict$amplifications[-1]/maxPredict*ylim[2],col='#FF000099')
	}
	dev.off()
}
