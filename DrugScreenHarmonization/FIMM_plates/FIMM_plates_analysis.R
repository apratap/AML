#boxplot : control variation across plates
df <- filter(k562_FO3D, assay == "CTG")
p1 <- ggplot(df, aes(x=content, y=value, color=plateNum)) + geom_boxplot() 
p1 <- p1 + theme_bw() + xlab('well types') + ylab('CTG raw values') + ggtitle("CTG: plate variation")
df <- filter(k562_FO3D, assay == "CTxG")
p2 <- ggplot(df, aes(x=content, y=value, color=plateNum)) + geom_boxplot() 
p2 <- p2 + theme_bw() + xlab('well types') + ylab('CTxG raw values') + ggtitle("CTxG: plate variation")
df <- filter(k562_FO3D, assay == "MTS")
p3 <- ggplot(df, aes(x=content, y=value, color=plateNum)) + geom_boxplot() 
p3 <- p3 + theme_bw() + xlab('well types') + ylab('MTS norm vals') + ggtitle("MTS: plate variation")
multiplot(p1,p2,p3)



# x <- dcast(filter(k562_FO3D, ProductName != ''), ProductName + center + cellLine + assay + set_number ~ dose  , 
#            value.var = "normViab" )
# 
# m <- x[,c(6:10)]
# head(m)
# m <- t(scale(t(m)))
# library("NMF")
# aheatmap(m,
#          scale='none',
#          distfun="euclidean",
#          Colv=NA,
#          Rowv =  T,
#          annRow = x[,c(4)],
#          info=TRUE,
#          cexRow=0)


# filter(k562_FO3D_fit, is.na(IC50) == T & goodNess_of_fit > .8 ) %>% 
#   select(center, cellLine, set_number,assay, ProductName)


#correlating assays and replicates
k562_FO3D_fit_cor <- dcast(k562_FO3D_fit, ProductName ~ center + assay + set_number , 
                           value.var = "Simpson" )
k562_FO3D_fit_cor$ProductName <- NULL
k562_FO3D_fit_cor <- apply(k562_FO3D_fit_cor,2,as.numeric)

corPlot(k562_FO3D_fit_cor, title="CellLine K562 comparison (based on AUC)")

#scatter : IC50 v/s goodness
min_conc <- min(k562_FO3D$conc, na.rm=T)
max_conc <- max(k562_FO3D$conc, na.rm=T)
p1 <- ggplot(data=k562_FO3D_fit, aes(x=log10(IC50), y=goodNess_of_fit, color=assay)) + geom_point(alpha=0.8)
p1 <- p1 + geom_vline(xintercept=c(log10(min_conc), log10(max_conc)), color="red")
p1





df <- filter(k562_FO3D, ProductName == "Anagrelide" & assay =='CTG')
View(k562_FO3D_fit)
get_drugResponse_stats()
View(k562_FO3D)
x <- filter(molm_FO3D, ProductName == 'Dabrafenib', assay =='CTG')
View(x)
ggplot(data=x, aes(x=dose, y=normViab,group=test)) + geom_point() + geom_line(group=test)






#correlating assays and replicates
molm_FO3D_fit_cor <- dcast(molm_FO3D_fit, ProductName ~  center + assay + set_number , 
                           value.var = "Simpson" )
molm_FO3D_fit_cor$ProductName <- NULL
corPlot(molm_FO3D_fit_cor, title="CellLine MOLM14 comparison (AUC based)")

min_conc <- min(molm_FO3D$conc__num_)
max_conc <- max(molm_FO3D$conc__num_)


#scatter : IC50 v/s AUC
p2 <- ggplot(data=k562_FO3D_fit, aes(y=Simpson, x=log10(IC50))) + geom_point(alpha=0.4) + ylab('AUC')
p2 <- p2 + geom_vline(xintercept=c(log10(min_conc), log10(max_conc)), color="red")
p2


#boxplot : control variation across plates
df <- filter(molm_FO3D, assay == "CTG")
p1 <- ggplot(df, aes(x=content, y=value, color=plateNum)) + geom_boxplot() 
p1 <- p1 + theme_bw() + xlab('well types') + ylab('CTG raw values') 
p1 <- p1 + ggtitle("MOLM14 plates variation (CTG assay)")
p1



molm_FO3D_fit_flt <- filter(molm_FO3D_fit, goodNess_of_fit > 0.6)


#AUC variation across assays
ggplot(data=molm_FO3D_fit,
       aes(x=paste(assay,set_number), y=Simpson)) + geom_boxplot()

#IC50 variation across assays
ggplot(data=molm_FO3D_fit,
       aes(x=paste(assay,set_number), y=log10(IC50))) + geom_boxplot()

#goodNess of fit v/s IC50 
ggplot(data=molm_FO3D_fit,
       aes(x=log10(IC50), y=goodNess_of_fit)) + geom_point(alpha=0.3)

# AUC v/s IC50 
ggplot(data=molm_FO3D_fit,
       aes(x=Simpson, y=log10(IC50))) + geom_point(alpha=0.3)

#comparing raw CTG values 
df <- molm_FO3D %>%
  filter(assay == "CTG")
ggplot(df, aes(x=content, y=value, color=plateNum)) + geom_boxplot()



selected_cols <- c('center', 'assay', 'ProductName', 'set_number', 'cellLine',
                   'dose', 'normViab')
combined_data <- rbind(k562_FO3D[,selected_cols], molm_FO3D[,selected_cols])
combined_data <- filter(combined_data, ProductName != "")
x <- dcast(combined_data, ProductName + center + cellLine + assay + set_number ~ dose  , 
           value.var = "normViab" )


mean(x$D5 - x$D1)

m <- x[,c(6:10)]
#m <- m[sample(1:nrow(m),1000),]
library("NMF")
aheatmap(m,
         scale='none',
         distfun="euclidean",
         Colv=NA,
         Rowv =  T,
         annRow = F,
         info=TRUE,
         cexRow=0)



combined_fit_data <- rbind(molm_FO3D_fit,k562_FO3D_fit)
View(combined_fit_data)


#summary 
View(ddply(.data = combined_fit_data,
           .variables=c('center','cellLine','assay','set_number'),
           .fun = function(x) nrow(x)))







