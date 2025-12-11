
suppressPackageStartupMessages({library(Maaslin2);library(tidyverse);library(readr)})
args<-commandArgs(trailingOnly=TRUE)
metadata_path<-if(length(args)>=1) args[1] else "metadata.tsv"
virome_dir<-if(length(args)>=2) args[2] else "out"
analysis_out<-if(length(args)>=3) args[3] else "maaslin_results"
group_var<-if(length(args)>=4) args[4] else "Group"
covariates<-if(length(args)>=5) strsplit(args[5],",")[[1]] else c()
levels_to_run<-if(length(args)>=6) strsplit(args[6],",")[[1]] else c("order","family","species")
normalization<-if(length(args)>=7) args[7] else "TSS"
transform<-if(length(args)>=8) args[8] else "LOG"
dir.create(analysis_out,showWarnings=FALSE,recursive=TRUE)
md<-read_tsv(metadata_path,show_col_types=FALSE)
stopifnot("SampleID"%in%colnames(md))
md[[group_var]]<-as.factor(md[[group_var]])
if(length(covariates)>0){for(v in covariates){if(is.numeric(md[[v]])) md[[v]]<-scale(md[[v]])%>%as.numeric()}}
build_feature_matrix<-function(level_name,virome_dir,sample_ids=md$SampleID){
  files<-file.path(virome_dir,paste0(sample_ids,".virome_",level_name,".tsv"))
  names(files)<-sample_ids
  tbls<-lapply(names(files),function(sid){
    f<-files[[sid]]
    x<-read_tsv(f,col_names=c("Taxon","Count"),show_col_types=FALSE)
    x<-x%>%group_by(Taxon)%>%summarize(Count=sum(Count),.groups="drop")
    x%>%pivot_wider(names_from=Taxon,values_from=Count,values_fill=0)%>%mutate(.sample=sid)
  })
  mat<-bind_rows(tbls)%>%relocate(.sample)
  features<-mat%>%select(-.sample)
  rownames(features)<-mat$.sample
  t(as.matrix(features))
}
run_maaslin_level<-function(level){
  feat_mat<-build_feature_matrix(level,virome_dir)
  common_samples<-intersect(colnames(feat_mat),md$SampleID)
  feat_mat<-feat_mat[,common_samples,drop=FALSE]
  md_sub<-md%>%filter(SampleID%in%common_samples)%>%arrange(match(SampleID,common_samples))
  stopifnot(identical(colnames(feat_mat),md_sub$SampleID))
  feat_mat[feat_mat<0]<-0
  res_dir<-file.path(analysis_out,paste0("maaslin_",level))
  dir.create(res_dir,showWarnings=FALSE,recursive=TRUE)
  Maaslin2(input_data=feat_mat,input_metadata=md_sub%>%column_to_rownames("SampleID"),output=res_dir,normalization=normalization,transform=transform,fixed_effects=c(group_var,covariates),random_effects=c(),min_prevalence=0.1,min_abundance=0.0,standardize=TRUE,plot_heatmap=TRUE,plot_scatter=TRUE,cores=max(1,parallel::detectCores()-1),max_significance=0.05)
}
fits<-lapply(levels_to_run,run_maaslin_level)
summaries<-lapply(levels_to_run,function(level){read_tsv(file.path(analysis_out,paste0("maaslin_",level),"all_results.tsv"),show_col_types=FALSE)%>%filter(qval<=0.05,metadata==group_var)%>%mutate(level=level)})
sig_all<-bind_rows(summaries)
write_tsv(sig_all,file.path(analysis_out,"significant_associations.tsv"))
for(lv in levels_to_run){
  df<-read_tsv(file.path(analysis_out,paste0("maaslin_",lv),"all_results.tsv"),show_col_types=FALSE)%>%filter(metadata==group_var)
  p<-df%>%ggplot(aes(x=coef,y=-log10(qval),color=qval<=0.05))+geom_point(alpha=0.8)+scale_color_manual(values=c("FALSE"="grey70","TRUE"="#D55E00"))+labs(title=paste("MaAsLin2 associations:",lv),x=paste("Effect size (",group_var,")",sep=""),y=expression(-log[10](q)))+theme_bw()
  ggsave(filename=file.path(analysis_out,paste0("volcano_",lv,".png")),plot=p,width=7,height=5,dpi=300)
  