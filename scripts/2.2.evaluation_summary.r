panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use="complete.obs"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}


#########################
######### summary #######
#########################
library(ggplot2)
load("s2.assign_eval.polycat.Rdata")
load("s2.assign_eval.rsem.Rdata")
load("s2.assign_eval.hylite.Rdata")
load("s2.assign_eval.salmon.Rdata")
load("s2.assign_eval.kallisto.Rdata")



## check sample order
polycatSummary$info$sample == hyliteSummary$info$sample
rsemSummary$info$sample == hyliteSummary$info$sample
salmonSummary$info$sample == hyliteSummary$info$sample
kallistoSummary$info$sample == hyliteSummary$info$sample

## double check libSize
size<-data.frame(polycat=polycatSummary$info$ADsLibSize, hylite=hyliteSummary$info$ADsLibSize, rsem=rsemSummary$info$ADsLibSize, salmon=salmonSummary$info$ADsLibSize, kallisto=kallistoSummary$info$ADsLibSize )
rownames(size)=polycatSummary$info$sample
pairs( size, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, main ="Library size")
pairs( size[1:22,], lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, main ="Library size: PE")
pairs( size[23:33,], lower.panel=panel.smooth, upper.panel=panel.cor, pch=20,  main ="Library size: SE")
# bar plots not very useful
# library(reshape)
# df<-size
# df$sample<-rownames(df)
# df<-melt(df,id="sample")
# ggplot(df, aes(x=sample,y=value, fill=variable))+ geom_bar(stat="identity", position=position_dodge())

## summary
sumTbl=list()
metrics = c( "Efficiency" , "Discrepancy", "Accuracy",    "Precision", "Fmeasure","MCC"  )
methods = c("polycat","hylite","rsem", "salmon", "kallisto")
for(i in metrics){
    if(exists("temp")){rm(temp)}
    for(m in methods)
    {
        x = get(paste0(m,"Summary"))[[i]]["overall",]
        if(!exists("temp")) {temp=x}else {temp=rbind(temp,x)}
        rownames(temp)[nrow(temp)] =m
    }
    sumTbl[[i]] = temp
}
print(sumTbl)

## gene-wise metrics comparison, noting that sample order is different in hylite
message("Compare metrics between methods: ")
print("Compare metrics between methods: ")
pw<-combn(length(methods),2)
for(i in 1:ncol(pw)){
    print(paste0(methods[pw[1,i]]," vs ",methods[pw[2,i]]))
    res1<- get(methods[pw[1,i]])
    res2<- get(methods[pw[2,i]])
    for(m in c( "efficiency" , "discrepancy", "accuracy"))
    {
        select <- is.finite(res1[[m]]) & is.finite(res2[[m]])
        
        print(paste0("Testing ",m))
        print( t.test(res1[[m]][select],res2[[m]][select], paired=TRUE) )
        print( wilcox.test(res1[[m]][select],res2[[m]][select], paired=TRUE) )
    }
    
}

# compare At vs Dt
message("Compare metrics between homoeologs within each method: ")
print("Compare metrics between homoeologs within each method: ")
for(i in paste0(methods,"H")){
    resH<-get(i)
    for(m in c("efficiency","discrepancy","accuracy","precision","Fstat"))
    {
        # m<-"discrepancy"
        # res<-polycatH
        col<-grep(m,names(resH))
        select <- is.finite(resH[,col[1]]) & is.finite(resH[,col[2]])
        print(paste0(i,": ",m))
        print( t.test(resH[select,col[1]],resH[select, col[2]], paired=TRUE) )
        print( wilcox.test(resH[select,col[1]],resH[select, col[2]], paired=TRUE) )
    }
}

# multiple measurements as potential explainational variables to performance metrics

message("Linear regression analysis of method")
for(i in methods){
    print(paste0("Linear regression analysis of method: ",i))
    res<- get(i)
    
    for(m in c("accuracy","efficiency","discrepancy"))
    {
        message(paste0("...",m))
        res<-na.omit(res)
        
        model.null = glm(get(m) ~ 1,
        data=res)
        
        model.full = glm(get(m) ~ geneLenM + percentageEffectM + expression.log2 +sample,
        data=res)
        
        step(model.null,
        scope = list(upper=model.full),
        direction="both",
        test="Chisq",
        data=Data)
        
    }
}

# bar plots not very useful
library(reshape)
plotMetric=function(x, main=""){
    x$method=rownames(x)
    df<-melt(x,ids="method")
    df$method<- factor(df$method,levels = methods) # my order
    p<-ggplot(df, aes(x=method,y=value, fill=method))+ geom_bar(stat="identity", position=position_dodge()) +  geom_text( aes(label = round(value,4), y = value - 0.05), position = position_dodge(0.9), vjust = 0)+ facet_grid(variable~.) + ggtitle(main)
    print(p)
}
pdf("s2.assign_eval.summary.pdf")
for(i in names(sumTbl)){
    plotMetric(sumTbl[[i]],i)
}
dev.off()
