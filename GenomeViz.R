#!/usr/bin/env Rscript

args=commandArgs(trailingOnly = TRUE)
if (length(args)<2){stop("2 parameters needed: 1. Sample ID which is part of the VCF file name; 2. Absolute path to analysis directory containing a folder 'VCF' which in turn contains the VCF file.")}

library(data.table)
library(ggplot2)
options(warn=-1)

#define hilbert curve function
Hilbert <- function(level=5, x=0, y=0, xi=1, xj=0, yi=0, yj=1) {
  if (level <= 0) {
    return(c(x + (xi + yi)/2, y + (xj + yj)/2))
  } else {
    return(rbind(
      Hilbert(level-1, x,           y,           yi/2, yj/2,  xi/2,  xj/2),
      Hilbert(level-1, x+xi/2,      y+xj/2 ,     xi/2, xj/2,  yi/2,  yj/2),
      Hilbert(level-1, x+xi/2+yi/2, y+xj/2+yj/2, xi/2, xj/2,  yi/2,  yj/2),
      Hilbert(level-1, x+xi/2+yi,   y+xj/2+yj,  -yi/2,-yj/2, -xi/2, -xj/2)
    ))
  }
}

sample=args[1]  #e.g. sample="PGA_0002"
dir=args[2]    #e.g. VCF=fread(paste0("/fhgfs/groups/lab_bock/jklughammer/projects/otherProjects/PGA/)

message("Now reading the VCF file")
file=system(paste0("ls VCF/*",sample,"_parsed.vcf"),intern=TRUE)
system(paste0("mkdir -p ",sample))
VCF=fread(file)
setnames(VCF,names(VCF),c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","EFFECT","FUNCTIONAL_CLASS","IMPACT"))

message("Now filtering for high quality variants")
VCF_pass=VCF[FILTER=="PASS"]
VCF_pass[,loc:=c(1:nrow(VCF_pass)),]

message("Now calculating or loading hilbert curve")
if (file.exists("h11.tsv")){h11=read.table("h11.tsv",sep="\t",header=TRUE)}else{
  h11=Hilbert(level=11)
  write.table(h11,"h11.tsv",quote=FALSE,sep="\t",row.names=FALSE)
}

VCF_pass[,hilb_x:=h11[,1],]
VCF_pass[,hilb_y:=h11[,2],]

message("Now producing color code and legend")
Effects_col=VCF_pass[,.N,by=EFFECT][order(N,decreasing=TRUE)]
Effects_col[,col:=c("#B0C0E1","#6285DA","gray70","#29C68F","#D8C443","#506F9C","#A28AFB","#DD8643","#93DB3F","#C890C0","#C9CB8C","#C28879","#7FD9B1","#474323","#618532","#7F312B","#D150CC","#C08035","#65D170","#723F7E","coral2","#C8D11A","#CF477B"),]
Effects_col[,pos_x:=1,]
Effects_col[,pos_y:=c(nrow(Effects_col):1),]

png(paste0(sample,"/",sample,"_legend.png"),width=600,height=800)
ggplot(Effects_col,aes(x=pos_x,y=pos_y,label=paste0(EFFECT,": ",N)))+geom_text(col=Effects_col$col,hjust=0)+scale_x_continuous(limits=c(0.5,10))+theme(panel.grid=element_blank(),panel.background=element_rect(fill="white"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank())
dev.off()

#select mutation to be used for zoom
#PGA_0002: rs2228570 = VDR
#PGA_0001: rs11549407 = HBB

print(VCF_pass[EFFECT%in%Effects_col[N<100]$EFFECT,list(CHROM,POS,ID,EFFECT)][order(EFFECT)],nrows=Inf)
message("Please select ID of focus mutation from the displayed list and press enter: ")
mut=readLines(file("stdin"),1)

while (grepl("rs",mut)!=TRUE){
  message("Incorrect mutation ID. Mutation ID must start with rs.")
  message("Please select ID of focus mutation from the displayed list and press enter: ")
  mut=readLines(file("stdin"),1)
}

#assign coloring to the variants
setkey(Effects_col,"EFFECT")
setkey(VCF_pass,"EFFECT")
VCF_col=merge(VCF_pass,Effects_col,"EFFECT")
VCF_col=VCF_col[order(loc)]

#actually plot the hilbert curve zooming
message("Start plotting...")
mut_loc=VCF_col[ID==mut]$loc
mut_x=VCF_col[ID==mut]$hilb_x
mut_y=VCF_col[ID==mut]$hilb_y
for (i in 0:11){
  message(paste0("Plotting hilbert ",i))
  n=4^i
  anchor=floor(mut_loc/n)*n
  if(i==0){sub=VCF_col[(anchor):(anchor)]}else{sub=VCF_col[(anchor+1):(anchor+n)] }
  x_min=min(sub$hilb_x)
  x_max=max(sub$hilb_x)
  y_min=min(sub$hilb_y)
  y_max=max(sub$hilb_y)
  
  fac_x=2.0/((i+0.00000001)*2.3)^1.8*(x_max-x_min)
  fac_y=2.0/((i+0.00000001)*2.3)^1.8*(y_max-y_min)

#plot without zooming square  
  png(paste0(sample,"/",sample,"_hilb_",11-i,"_1.png"),width=800,height=800)
  plot(sub$hilb_x,sub$hilb_y,col="black",type="l",xlab="",ylab="",axes=F,lwd=20/(i*1.7),xlim=c(x_min-fac_x,x_max+fac_x),ylim=c(y_min-fac_y,y_max+fac_y))
  
  if (i==0){par(cex=2)}
  points(sub$hilb_x,sub$hilb_y,col=sub$col,type="p",pch=".",cex=(200+50*i)/(2^(i)))
  if(i==0){
    text(x=sub$hilb_x,y=sub$hilb_y,lab=paste0(sub$ID,"\n","Mut:",sub$REF,">",sub$ALT,"\n","chr",sub$CHROM,":",sub$POS,"\n",sub$EFFECT))
    dev.off()
    next
  }
  dev.off()
  
#plot with zooming square
  png(paste0(sample,"/",sample,"_hilb_",11-i,"_2.png"),width=800,height=800)  
  plot(sub$hilb_x,sub$hilb_y,col="black",type="l",xlab="",ylab="",axes=F,lwd=20/(i*1.7),xlim=c(x_min-fac_x,x_max+fac_x),ylim=c(y_min-fac_y,y_max+fac_y))
  points(sub$hilb_x,sub$hilb_y,col=sub$col,type="p",pch=".",cex=(200+50*i)/(2^(i)))
  
  xmin <- x_min-fac_x
  xmax <- x_max+fac_x
  ymin <- y_min-fac_y
  ymax <- y_max+fac_x
  
  if (mut_x>xmin+(xmax-xmin)/2&mut_y>ymin+(ymax-ymin)/2) {
    rect(xmin+(xmax-xmin)/2,ymin+(ymax-ymin)/2,xmax,ymax,col=NA,border="red",lwd=2/((i+2)*0.1))
  }
  if (mut_x<xmin+(xmax-xmin)/2&mut_y>ymin+(ymax-ymin)/2) {
    rect(xmin,ymin+(ymax-ymin)/2,xmin+(xmax-xmin)/2,ymax,col=NA,border="red",lwd=2/((i+2)*0.1))
  }
  if (mut_x>xmin+(xmax-xmin)/2&mut_y<ymin+(ymax-ymin)/2) {
    rect(xmin+(xmax-xmin)/2,ymin,xmax,ymin+(ymax-ymin)/2,col=NA,border="red",lwd=2/((i+2)*0.1))
  }
  if (mut_x<xmin+(xmax-xmin)/2&mut_y<ymin+(ymax-ymin)/2) {
    rect(xmin,ymin,xmin+(xmax-xmin)/2,ymin+(ymax-ymin)/2,col=NA,border="red",lwd=2/((i+2)*0.1))
  }
  dev.off()
}
message("Done.")