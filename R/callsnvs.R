#' callsnvs
#'
#' @param sample.info.file sample info file listing cases and controls. tab-delimeted file
#' @param targetbp folder with RData for each annotated positions
#' @param pbem_dir folder with pbem data. default: file.path(outdir, BaseErrorModel)
#' @param outdir output folder for this step analysis
#' @param outdir.calls.name folder will be created in outdir. default: "Results"
#' @export
callsnvs <- function(sample.info.file,
                     targetbp,
                     pbem_dir = file.path(outdir,"BaseErrorModel"),
                     PBEsim = "~/datasetoy/RData/PBEsim.RData",
                     AFbycov = TRUE,
                     minaf_cov_corrected,
                     mincov = 10,
                     mincovgerm = 10,
                     minalt = 1,
                     maxafgerm = 0.2,
                     outdir,
                     outdir.calls.name="Results"){

  cat(paste("[",Sys.time(),"]\tReading the sample.info.file","\n"))
  sif <- import_sif(main_sif = sample.info.file)

  cat(paste("[",Sys.time(),"]\tReading chromosomes from bpcovered.tsv","\n"))
  chromosomes = read.delim(file = file.path(targetbp,"bpcovered.tsv"),as.is=T)
  chromosomes = sort(paste0("chr",chromosomes[-nrow(chromosomes),1]))

  cat(paste("[",Sys.time(),"]\tDetection of somatic SNVs in case samples","\n"))
  if(!file.exists(file.path(outdir, outdir.calls.name))){
    dir.create(file.path(outdir, outdir.calls.name), showWarnings = T)
  }

  # Apply filters to define a set of putative SNVs
  cat(paste("[",Sys.time(),"]\tApply basic filters"),"\n")
  fpam = data.frame(AFbycov = as.character(AFbycov),
                    spec = as.character(minaf_cov_corrected[1,1]),
                    mincov = as.character(mincov),
                    minalt = as.character(minalt),
                    mincovgerm = as.character(mincovgerm),
                    maxafgerm = as.character(maxafgerm),
                    filtering_date = Sys.time(),
                    stringsAsFactors = F)
  write.table(fpam,file = file.path(outdir, outdir.calls.name, "filtering_criteria.txt"),col.names = T,row.names = F,sep="\t",quote = F)

  # Import background pbem
  tab_bg_pbem = read.delim(file = file.path(pbem_dir, "pbem_background.tsv"),as.is=T,header=T,stringsAsFactors = F)

  # Import matrix for coverages and pbem
  load(PBEsim)
  tab_cov_pbem = tab.list[[6]]

  chrom.in.parallel = 1 # to be corrected
  mc.cores = 1 # to be corrected

  TableSif <- sif[[2]]
  for(id in 1:nrow(TableSif)){
    this = TableSif[id,]
    name.patient = TableSif$patient[id]
    name.plasma = gsub(basename(this$plasma.bam),pattern = ".bam",replacement = "")
    name.germline = gsub(basename(this$germline.bam),pattern = ".bam",replacement = "")
    cat(paste("[",Sys.time(),"]\tPatient:",name.patient,"\tCase:",name.plasma,"\tControl:",name.germline),"\n")
    out1 = paste0("pmtab_F1_",name.plasma,".tsv")
    out2 = paste0("pmtab_F2_",name.plasma,".tsv")
    out3 = paste0("pmtab_F3_",name.plasma,".tsv")
    # create patient sub-folder
    patient_folder = file.path(outdir, outdir.calls.name, name.patient)
    dir.create(patient_folder)
    germline.folder = list.files(pacbamfolder, pattern = paste0(name.germline,"$"),full.names = T)
    if(length(germline.folder)==0){
      message("[ ERROR ] Cannot find folder:\t",file.path(pacbamfolder, name.germline))
      quit()
    }
    plasma.folder = list.files(pacbamfolder, pattern = paste0(name.plasma,"$"),full.names = T)
    if(length(plasma.folder)==0){
      message("[ ERROR ] Cannot find folder:\t",file.path(pacbamfolder, name.plasma))
      quit()
    }

    # run in parallel on chromosomes
    mclapply(seq(1,length(chromosomes),1),
             filter,
             chromosomes=chromosomes,
             patient_folder=patient_folder,
             plasma.folder=plasma.folder,
             germline.folder=germline.folder,
             out1=out1,
             out2=out2,
             out3=out3,
             mc.cores = chrom.in.parallel)

    # collapse all chromosome outs into a single table
    setwd(patient_folder)
    tabs_list = list.files(patient_folder,full.names = T,recursive = T,pattern = 'chrpm_f1.tsv')
    if(length(tabs_list)>0){
      cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out1)
      system(cmd)
    } else {
      cat(file = "f1_table_WARNING.txt","No calls found in chrpm_f1.tsv")
    }
    tabs_list = list.files(patient_folder,full.names = T,recursive = T,pattern = 'chrpm_f2.tsv')
    if(length(tabs_list)>0){
      cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out2)
      system(cmd)
    } else {
      cat(file = "f2_table_WARNING.txt","No calls found in chrpm_f2.tsv")
    }
    tabs_list = list.files(patient_folder,full.names = T,recursive = T,pattern = 'chrpm_f3.tsv')
    if(length(tabs_list)>0){
      cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out3)
      system(cmd)
    } else {
      cat(file = "f3_table_WARNING.txt","No calls found in chrpm_f3.tsv")
    }
  }




}
