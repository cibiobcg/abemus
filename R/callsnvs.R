#' callsnvs
#'
#' @param sample.info.file sample info file listing cases and controls. tab-delimeted file
#' @param targetbp folder with RData for each annotated positions
#' @param pbem_dir folder with pbem data. default: file.path(outdir, BaseErrorModel)
#' @param minaf_cov_corrected the output of adjust_af_threshold_table()
#' @param PBEsim default: "~/datasetoy/RData/PBEsim.RData"
#' @param mincov Minimum locus coverage in case sample. default: 10
#' @param minalt Minimum number of reads supporting the alternative allele in case sample. default: 1
#' @param mincovgerm Minimum locus coverage in germline sample. default: 10
#' @param maxafgerm Maximum allelic fraction observed in matched germline locus. default: 0.2
#' @param outdir output folder for this step analysis
#' @param outdir.calls.name folder will be created in outdir. default: "Results"
#' @param AFbycov Apply afth coverage based. default: AFbycov=TRUE
#' @param pacbamfolder folder with pileups
#' @param coverage_binning Bins of coverage into which divide allelic fractions. default: 50
#' @param chrom.in.parallel number of chromosomes to run in parallel. default: 1
#' @export
callsnvs <- function(sample.info.file,
                     targetbp,
                     pacbamfolder,
                     pbem_dir = file.path(outdir,"BaseErrorModel"),
                     PBEsim = "~/datasetoy/RData/PBEsim.RData",
                     AFbycov = TRUE,
                     minaf_cov_corrected,
                     mincov = 10,
                     mincovgerm = 10,
                     minalt = 1,
                     maxafgerm = 0.2,
                     coverage_binning = 50,
                     outdir,
                     outdir.calls.name="Results",
                     chrom.in.parallel=1){

  cat(paste("[",Sys.time(),"]\tReading the sample.info.file","\n"))
  sif <- import_sif(main_sif = sample.info.file)

  cat(paste("[",Sys.time(),"]\tReading chromosomes from bpcovered.tsv","\n"))
  chromosomes = read.delim(file = file.path(targetbp,"bpcovered.tsv"),as.is=T)
  chromosomes = sort(paste0("chr",chromosomes[-nrow(chromosomes),1]))

  cat(paste("[",Sys.time(),"]\tDetection of somatic SNVs in case samples","\n"))
  if(!file.exists(file.path(outdir, outdir.calls.name))){
    dir.create(file.path(outdir, outdir.calls.name), showWarnings = T)
  }

  # summary of thresholds applied
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
  xbg <- as.numeric(tab_bg_pbem$background_pbem)

  # Import matrix for coverages and pbem
  load(PBEsim)
  tab_cov_pbem = tab.list[[6]]

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
    # create sub-folder to save outputs per each case
    caseout_folder = file.path(outdir, outdir.calls.name, name.plasma)
    dir.create(caseout_folder,showWarnings = T)

    germline.folder = list.files(pacbamfolder, pattern = paste0(name.germline,"$"),full.names = T)
    if(length(germline.folder)==0){
      message("[ ERROR ] Cannot find folder:\t",file.path(pacbamfolder, name.germline))
      stop()
    }
    plasma.folder = list.files(pacbamfolder, pattern = paste0(name.plasma,"$"),full.names = T)
    if(length(plasma.folder)==0){
      message("[ ERROR ] Cannot find folder:\t",file.path(pacbamfolder, name.plasma))
      stop()
    }

    # run in parallel on chromosomes
    mclapply(seq(1,length(chromosomes),1),
             filter,
             name.patient=name.patient,
             name.plasma=name.plasma,
             name.germline=name.germline,
             chromosomes=chromosomes,
             caseout_folder=caseout_folder,
             plasma.folder=plasma.folder,
             germline.folder=germline.folder,
             pbem_dir=pbem_dir,
             out1=out1,
             out2=out2,
             out3=out3,
             mincov=mincov,
             minalt=minalt,
             mincovgerm=mincovgerm,
             maxafgerm=maxafgerm,
             AFbycov=AFbycov,
             minaf_cov_corrected=minaf_cov_corrected,
             coverage_binning=coverage_binning,
             xbg=xbg,
             tab_cov_pbem=tab_cov_pbem,
             afs=afs,
             covs=covs,
             mc.cores = chrom.in.parallel)

    # collapse all chromosome outs into a single table
    setwd(caseout_folder)
    tabs_list = list.files(caseout_folder,full.names = T,recursive = T,pattern = 'chrpm_f1.tsv')
    if(length(tabs_list)>0){
      cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out1)
      system(cmd)
    } else {
      cat(file = "f1_table_WARNING.txt","No calls found in chrpm_f1.tsv")
    }
    tabs_list = list.files(caseout_folder,full.names = T,recursive = T,pattern = 'chrpm_f2.tsv')
    if(length(tabs_list)>0){
      cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out2)
      system(cmd)
    } else {
      cat(file = "f2_table_WARNING.txt","No calls found in chrpm_f2.tsv")
    }
    tabs_list = list.files(caseout_folder,full.names = T,recursive = T,pattern = 'chrpm_f3.tsv')
    if(length(tabs_list)>0){
      cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out3)
      system(cmd)
    } else {
      cat(file = "f3_table_WARNING.txt","No calls found in chrpm_f3.tsv")
    }
  }

  # summary table index for calls
  tabsnvs_index <- get_tabcalls_path(TableSif = TableSif,
                                     outdir = outdir,
                                     outdir.calls.name = outdir.calls.name)

  cat(paste("[",Sys.time(),"]\talright.","\n"))
  return(tabsnvs_index)
}
