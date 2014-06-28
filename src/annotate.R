suppressMessages(library(biomaRt))

annotate.queryBMUniprot = function(results, attributes=c("accession","name","gene_name","protein_name","organism")){
  unique_uniprot_ac = unique(results$Prey)
  uniProt <- useMart("unimart", dataset="uniprot")
  annotations = getBM(attributes=attributes,filter="accession",values=uniprot_ac,mart=uniProt)
  results_annotated = merge(results, annotations, all.x=T, by.x='Prey',by.y='accession')
  results_annotated
}

annotate.queryFile = function(results, species='HUMAN', uniprot_dir="~/Projects/HPCKrogan/Scripts/MSPipeline/files/"){
  species_split = unlist(strsplit(species, "-"))
  
  Uniprot = NULL
  for(org in species_split){
    cat(sprintf("\tLOADING %s\n",org))
    tmp = read.delim2(sprintf("%s/uniprot_protein_descriptions_%s.txt",uniprot_dir,org), stringsAsFactors=F, quote="")    
    if(is.null(Uniprot)){
      Uniprot = as.data.frame(tmp)
    }else{
      Uniprot = rbind(Uniprot, tmp)  
    }
  }
  results_annotated = merge(results, Uniprot, all.x=T, by.x='Prey',by.y='Entry')
  results_annotated
}
