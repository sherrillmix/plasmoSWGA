library(dnar)

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
primers<-strsplit(primerData$primers,'.')

