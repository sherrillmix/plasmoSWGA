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
predictCover<-lapply(primerSets,function(primers){
	fs<-gregexpr(paste(primers,collapse='|'),ref$seq)
	rs<-gregexpr(paste(revComp(primers),collapse='|'),ref$seq)
	rs<-lapply(rs,function(r)r+attr(r,'match.length')-1)
	names(rs)<-names(fs)<-ref$name
	out<-mapply(predictAmplifications,fs,rs,SIMPLIFY=FALSE,MoreArgs=list(maxLength=30000))
	names(out)<-ref$name
	return(list('pred'=out,'for'=fs,'rev'=rs))
})

mclapply(names(predictCover),function(ii){
	message(ii)
	outFile<-file.path('out',sprintf('%s.png',ii))
	thisCovers<-cover[grep(sprintf('data/%s',ii),names(cover))]
	thisPredicts<-predictCover[[ii]][['pred']]
	png(outFile,width=9000,height=3000,res=250)
	par(mfrow=c(4,4))
	ylim<-c(0,max(sapply(thisCovers,function(xx)max(xx$count))))
	maxPredict<-max(sapply(thisPredicts,function(xx)max(xx$amplifications)))
	for(jj in names(thisPredicts)){
		thisCover<-lapply(thisCovers,function(x)x[x$chr==jj,])
		thisPredict<-thisPredicts[[jj]]
		thisF<-predictCover[[ii]][['for']][[jj]]
		thisR<-predictCover[[ii]][['rev']][[jj]]
		plot(1,1,xlim=c(1,max(thisPredict$end)),ylim=ylim,main=jj,xlab='Position',ylab='Cover')
		abline(v=thisF,col='#FF000055',lty=2)
		#abline(v=thisR,col='#0000FF55',lty=2)
		sapply(thisCover,function(xx)points(xx$pos,xx$count))
		segments(thisPredict$start,thisPredict$amplifications/maxPredict*ylim[2],thisPredict$end,thisPredict$amplifications/maxPredict*ylim[2],col='#00FF0099')
		segments(thisPredict$end[-nrow(thisPredict)],thisPredict$amplifications[-nrow(thisPredict)]/maxPredict*ylim[2],thisPredict$start[-1],thisPredict$amplifications[-1]/maxPredict*ylim[2],col='#00FF0099')
	}
	dev.off()
},mc.cores=10)


mclapply(names(predictCover),function(ii){
	message(ii)
	outFile2<-file.path('out',sprintf('repro_%s.png',ii))
	png(outFile2,height=3000,width=3000,res=300)
		par(mfrow=c(2,2))
		targets<-c('A','B','C')
		for(jj in 1:2){
			for(kk in (jj+1):3){
				cover1<-cover[[grep(sprintf('%s%s_',ii,targets[jj]),names(cover))]]
				cover2<-cover[[grep(sprintf('%s%s_',ii,targets[kk]),names(cover))]]
				count1Col<-sprintf('count%s',targets[jj])
				count2Col<-sprintf('count%s',targets[kk])
				colnames(cover1)[colnames(cover1)=='count']<-count1Col
				colnames(cover2)[colnames(cover2)=='count']<-count2Col
				compare<-merge(cover1,cover2,all=TRUE)
				compare[is.na(compare)]<-0
				plot(compare[,count1Col],compare[,count2Col],xlab=count1Col,ylab=count2Col,col='#00000033',cex=.5,main=ii)
			}
		}
	dev.off()
},mc.cores=10)

png('out/control.png',width=9000,height=3000,res=250)
	par(mfrow=c(4,4))
	thisCovers<-cover[grep('data/Control',names(cover))]
	ylim<-c(0,max(sapply(thisCovers,function(xx)max(xx$count))))
	for(jj in names(predictCover[[1]][['pred']])){
		thisCover<-lapply(thisCovers,function(x)x[x$chr==jj,])
		plot(1,1,xlim=c(1,max(thisCover[[1]]$pos)),ylim=ylim,main=jj,xlab='Position',ylab='Cover')
		sapply(thisCover,function(xx)points(xx$pos,xx$count))
	}
dev.off()




pdf('out/regionOfInterest.pdf',width=15)
	targetPrimer<-'36562'
	targetChr<-'Pf3D7_05_v3'
	targetPos<-141000
	window<-30000
	covers<-lapply(cover[grep(sprintf('data/%s',targetPrimer),names(cover))],function(x)x[x$chr==targetChr&abs(x$pos-targetPos)<window,])
	thisF<-predictCover[[targetPrimer]][['for']][[targetChr]]
	thisR<-predictCover[[targetPrimer]][['rev']][[targetChr]]
	thisPredict<-predictCover[[targetPrimer]][['pred']][[targetChr]]
	#maxPredict<-max(thisPredict$amplifications)
	maxPredict<-1e3
	yMax<-max(sapply(covers,function(xx)max(xx$count)))
	plot(1,1,type='n',xlab=sprintf('%s position',targetChr),main=targetPrimer,ylim=c(0,yMax/1000),xlim=targetPos+window*c(-1,1),ylab='Coverage (x1000)',las=1)
	sapply(cover,function(xx)lines(xx$pos,xx$count/1000))
	abline(v=thisF,col='#FF000055',lty=2)
	abline(v=thisR,col='#0000FF55',lty=2)
	topQuarter<-par('usr')[4]-c(0.05,.25)*diff(par('usr')[3:4])
	secondQuarter<-par('usr')[4]-c(.25,.45)*diff(par('usr')[3:4])
	fStack<-stackRegions(thisF,thisF+30000)
	rStack<-stackRegions(thisR-30000,thisR)
	fPos<-seq(topQuarter[1],topQuarter[2],length.out=max(fStack))[fStack]
	rPos<-seq(secondQuarter[1],secondQuarter[2],length.out=max(rStack))[rStack]
	arrows(thisF,fPos,thisF+30000,fPos,col='#FF000099',length=.1)
	arrows(thisR,rPos,thisR-30000,rPos,col='#0000FF99',length=.1)
	segments(thisPredict$start,thisPredict$amplifications/maxPredict*yMax/1000,thisPredict$end,thisPredict$amplifications/maxPredict*yMax/1000,col='#00FF0099')
	segments(thisPredict$end[-nrow(thisPredict)],thisPredict$amplifications[-nrow(thisPredict)]/maxPredict*yMax/1000,thisPredict$start[-1],thisPredict$amplifications[-1]/maxPredict*yMax/1000,col='#00FF0099')
dev.off()

