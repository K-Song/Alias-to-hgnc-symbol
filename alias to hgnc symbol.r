
alias_to_hgnc_symbol<-function(df,metric){
  library(limma)
  require(reshape2)
  df$gene<-rownames(df)
  list<-lapply(df$gene, FUN = function (x) alias2Symbol(x, species = "Hs"))
  names(list)<-df$gene
  indices_alias_not_found<-which(sapply(list, length) == 0)
  list[indices_alias_not_found]<-names(list)[indices_alias_not_found]
  df_symbol_gene<-setNames(melt(list), c("hgnc_symbol", "gene"))
  df_value_symbol_gene<-merge(df, df_symbol_gene, by.x = "gene", by.y = "gene") #have muliple rows for one gene
  columns_numeric<-lapply(df_value_symbol_gene, class) == "numeric" #pick numeric columns
  df_only_numeric<-df_value_symbol_gene[,columns_numeric] #pick numeric columns
  if(metric == "variance"){
    df_value_symbol_gene$variance<-apply(df_only_numeric, 1, var)
    df_value_symbol_gene$hgnc_symbol<-as.character(df_value_symbol_gene$hgnc_symbol)
    df_value_symbol_gene_ordered<-df_value_symbol_gene[order(df_value_symbol_gene$hgnc_symbol, -df_value_symbol_gene$variance), ]
  }#rank varaince from highest to lowest
  else {
    df_value_symbol_gene$mean<-rowSums(df_only_numeric)/apply(df_only_numeric, 1, function(x) sum(!is.na(x)))
    df_value_symbol_gene$hgnc_symbol<-as.character(df_value_symbol_gene$hgnc_symbol)
    df_value_symbol_gene_ordered<-df_value_symbol_gene[order(df_value_symbol_gene$hgnc_symbol, -df_value_symbol_gene$mean), ]
  }#rank mean from highest to lowest
  df_value_symbol_gene_duplicated_deleted<-df_value_symbol_gene_ordered[!duplicated(df_value_symbol_gene_ordered$hgnc_symbol),]
  df_value_symbol_gene_duplicated_deleted$variance<-NULL
  df_value_symbol_gene_duplicated_deleted$gene<-NULL
  return(df_value_symbol_gene_duplicated_deleted)
}