#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(pheatmap))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(MESS))

#######################################

## parameter grid 
mist.train.getParamGrid = function(steps=0.1) {
  param_weight_detail = seq(from=0,to=1, by=steps)
  param_grid = expand.grid(param_weight_detail, param_weight_detail, param_weight_detail)
  param_grid = param_grid[apply(param_grid[,1:3],1,sum) == 1,]
  colnames(param_grid) = c("Reproducibility","Abundance","Specificity")
  param_grid = data.frame(param_grid, ID=paste(param_grid$Reproducibility, param_grid$Abundance, param_grid$Specificity,sep="|"))
  param_grid
}

mist.train.getPredictionRates = function(prediction, thresholds = seq(from=0,to=0.99, by=.01)){
  pos = nrow(prediction[prediction$label==1,])
  neg = nrow(prediction[prediction$label==0,])
  pos_false = c()
  pos_true = c()
  neg_true = c()
  neg_false = c()
  ID=unique(prediction$ID)
  for(t in thresholds){
    pos_predicted = prediction[prediction$predicted_scores >= t ,]
    neg_predicted = prediction[prediction$predicted_scores < t ,]
    pf = nrow(pos_predicted[pos_predicted$label==0,])
    pos_false = c(pos_false, pf)
    pos_true = c(pos_true, nrow(pos_predicted)-pf)
    nf = nrow(neg_predicted[neg_predicted$label==1,])
    neg_false = c(neg_false, nf)
    neg_true = c(neg_true, nrow(neg_predicted)-nf)
  }
  res = data.frame(R_A_S=ID, threshold=thresholds, tpr=pos_true/pos, fpr=pos_false/neg, specificity=neg_true/neg, precision=pos_true/(pos_true+pos_false), fdr=pos_false/(pos_false+pos_true), acc=(pos_true+neg_true)/(pos+neg), f1=(2*pos_true)/((2*pos_true)+pos_false+neg_false), total_pos=pos_true+pos_false)
  res
}

# 
# plotROC = function(prediction_rates){
#   dim_predictions = round(sqrt(length(unique(prediction_rates$ID))))
#   prediction_rates = prediction_rates[with(prediction_rates,order(ID,tpr,fpr)),]
#   prediction_auc_summ = unique(prediction_rates[,c("ID","auc")])
#   p = ggplot(data=prediction_rates, aes(x=fpr,y=tpr))
#   p + geom_line() + 
#     facet_wrap(~ID, nrow=dim_predictions, ncol=dim_predictions) + 
#     geom_text(aes(x = 0.6, y = 0.1, label = round(auc,digits=3)), parse=T) +
#     theme(axis.text.x=element_text(angle=-90))
# }

## generic grid search 
mist.train.gridSearch = function(metrics, param_grid, predictionFun=mist.train.getPredictionRates){
  prediction_rates = c()
  for(i in 1:nrow(param_grid)){
    param_set = param_grid[i,]
    print(sprintf('%s/%s : R %s A %s S %s',i,nrow(param_grid),param_set$Reproducibility,param_set$Abundance ,param_set$Specificity))
    predicted_scores = (metrics$Reproducibility * param_set$Reproducibility) + (metrics$Abundance * param_set$Abundance) + (metrics$Specificity * param_set$Specificity)
    prediction = data.frame(metrics, param_set, predicted_scores=predicted_scores)
    prediction_rates = rbind(prediction_rates, predictionFun(prediction))
  }
  prediction_rates
}

mist.train.label = function(metrics, training_set){
  training_set = data.frame(training_set, label=1)
  metrics = merge(metrics, training_set, by=c('Bait','Prey'), all.x=T)
  metrics[is.na(metrics$label),]$label=0
  metrics
}

mist.train.main = function(metrics, training_set, output_file, training_steps=0.1){
  metrics = mist.train.label(metrics, training_set)  
  param_grid = mist.train.getParamGrid(training_steps)
  prediction_rates = mist.train.gridSearch(metrics, param_grid)
  #plotROC(prediction_rates)
  write.table(prediction_rates, file=output_file, sep="\t", quote=F, eol="\n", row.names=F)
  prediction_rates = prediction_rates[order(prediction_rates$f1, decreasing=T),]
  
  prediction_weights = as.numeric(unlist(strsplit(as.character(prediction_rates$R_A_S[1]),"\\|")))
  prediction_weights = data.frame(R=prediction_weights[1],A=prediction_weights[2],S=prediction_weights[3], stringsAsFactors=F)
  return(prediction_weights)
}

# metrics = read.delim('tests/entero/processed/preprocessed_MAT_MIST_METRICS.txt',stringsAsFactors=F)
# training_set = read.delim('tests/entero/input/PV_goldset.txt',stringsAsFactors=F, header=F)
# colnames(training_set) = c('Bait','Prey')
# mist.train.main(metrics, training_set)
# 
# 
