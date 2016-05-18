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
	return(out)
})

