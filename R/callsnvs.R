#' Call somatic SNVs by using both global and local error models
#' @param sample.info.file sample info file listing cases and controls. tab-delimeted file
#' @param targetbed targeted regions in BED format.
#' @param pacbamfolder_bychrom folder with pileups
#' @param pbem_dir folder with pbem data. default: file.path(outdir, "BaseErrorModel")
#' @param controls_dir folder with afth data. default: file.path(outdir,"Controls")
#' @param detection.specificity The quantile of the GSE distribution(s) to use to compute the afth. default: 0.995
#' @param replicas Relpica sampling to define stability of afth in bins of coverage. default: 1000
#' @param replicas.in.parallel default: 1
#' @param coeffvar.threshold Consider a bin as stable if the Coeff of variations after n replicas is lower than coeffvar.threshold. default: 0.01
#' @param AFbycov Apply afth coverage based. default: AFbycov=TRUE
#' @param mincov Minimum locus coverage in case sample. default: 10
#' @param mincovgerm Minimum locus coverage in germline sample. default: 10
#' @param minalt Minimum number of reads supporting the alternative allele in case sample. default: 1
#' @param maxafgerm Maximum allelic fraction observed in matched germline locus. default: 0.2
#' @param coverage_binning Bins of coverage into which divide allelic fractions. default: 50
#' @param outdir output folder for this step analysis
#' @param outdir.calls.name folder will be created in outdir. default: "Results"
#' @param chrom.in.parallel number of chromosomes to run in parallel. default: 1
#' @examples
#' sample.info.file <- system.file("extdata", "test_sif_toy.tsv", package = "abemus")
#' outdir <- tempdir()
#' targetbed <- system.file("extdata", "regions_toy.bed", package = "abemus")
#' pacbamfolder_bychrom <- system.file("extdata", "pacbam_data_bychrom", package = "abemus")
#' pbem_dir <- system.file("extdata", "BaseErrorModel", package = "abemus")
#' controls_dir <- system.file("extdata", "Controls", package = "abemus")
#' calls <- callsnvs(sample.info.file = sample.info.file,outdir=outdir,targetbed=targetbed,pbem_dir=pbem_dir,controls_dir=controls_dir,pacbamfolder_bychrom=pacbamfolder_bychrom,replicas = 2)
#' @return list of data frames with path to files with SNV calls.
#' @export
callsnvs <- function(sample.info.file,
                     targetbed,
                     pacbamfolder_bychrom,
                     pbem_dir = file.path(outdir,"BaseErrorModel"),
                     controls_dir = file.path(outdir,"Controls"),
                     detection.specificity = 0.995,
                     replicas = 1000,
                     replicas.in.parallel = 1,
                     coeffvar.threshold = 0.01,
                     AFbycov = TRUE,
                     mincov = 10,
                     mincovgerm = 10,
                     minalt = 1,
                     maxafgerm = 0.2,
                     coverage_binning = 50,
                     outdir,
                     outdir.calls.name = "Results",
                     chrom.in.parallel = 1){

  old_wd <- getwd()
  on.exit(setwd(old_wd))

  message(paste("[",Sys.time(),"]\tReading the sample.info.file"))
  sif <- import_sif(main_sif = sample.info.file)

  message(paste("[",Sys.time(),"]\tReading chromosomes from 'targetbed'"))
  chromosomes <- bed2positions(targetbed = targetbed,get_only_chromosomes = TRUE)[[1]]

  message(paste("[",Sys.time(),"]\tDetection of somatic SNVs in case samples"))
  if(!file.exists(file.path(outdir, outdir.calls.name))){
    dir.create(file.path(outdir, outdir.calls.name), showWarnings = TRUE)
  }

  # compute minaf_cov_corrected
  afthcor <- adjust_af_threshold_table(controls_dir = controls_dir,
                                       pbem_dir = pbem_dir,
                                       detection.specificity = detection.specificity,
                                       replicas = replicas,
                                       coverage_binning = coverage_binning,
                                       replicas.in.parallel = replicas.in.parallel,
                                       coeffvar.threshold = coeffvar.threshold)
  minaf_cov_corrected <- afthcor$minaf_cov_corrected

  # summary of thresholds applied
  fpam = data.frame(AFbycov = as.character(AFbycov),
                    spec = as.character(minaf_cov_corrected[1,1]),
                    mincov = as.character(mincov),
                    minalt = as.character(minalt),
                    mincovgerm = as.character(mincovgerm),
                    maxafgerm = as.character(maxafgerm),
                    filtering_date = Sys.time(),
                    stringsAsFactors = FALSE)
  write.table(fpam,file = file.path(outdir, outdir.calls.name, "filtering_criteria.txt"),col.names = TRUE,row.names = FALSE,sep="\t",quote = FALSE)

  # Import background pbem
  load(file.path(pbem_dir, "pbem_background.RData"))
  xbg <- as.numeric(bgpbem)

  TableSif <- sif$df_cases
  for(id in 1:nrow(TableSif)){
    this = TableSif[id,]
    name.patient = TableSif$patient[id]
    name.plasma = gsub(basename(this$plasma.bam),pattern = ".bam",replacement = "")
    name.germline = gsub(basename(this$germline.bam),pattern = ".bam",replacement = "")
    message(paste("[",Sys.time(),"]\tPatient:",name.patient,"\tCase:",name.plasma,"\tControl:",name.germline))
    out1 = paste0("pmtab_F1_",name.plasma,".tsv")
    out2 = paste0("pmtab_F2_",name.plasma,".tsv")
    out3 = paste0("pmtab_F3_",name.plasma,".tsv")
    # create sub-folder to save outputs per each case
    caseout_folder = file.path(outdir, outdir.calls.name, name.plasma)
    dir.create(caseout_folder,showWarnings = TRUE)

    if(is.na(name.germline)){

      germline.folder <- NA

      plasma.folder = list.files(pacbamfolder_bychrom, pattern = paste0(name.plasma,"$"),full.names = TRUE)
      if(length(plasma.folder)==0){
        message("[ ERROR ] Cannot find folder:\t",file.path(pacbamfolder_bychrom, name.plasma))
        stop()
      }

    } else {

      germline.folder = list.files(pacbamfolder_bychrom, pattern = paste0(name.germline,"$"),full.names = TRUE)
      if(length(germline.folder)==0){
        message("[ ERROR ] Cannot find folder:\t",file.path(pacbamfolder_bychrom, name.germline))
        stop()
      }
      plasma.folder = list.files(pacbamfolder_bychrom, pattern = paste0(name.plasma,"$"),full.names = TRUE)
      if(length(plasma.folder)==0){
        message("[ ERROR ] Cannot find folder:\t",file.path(pacbamfolder_bychrom, name.plasma))
        stop()
      }

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
             tab_cov_pbem=bombanel_tab_cov_pbem,
             afs=bombanel_afs,
             covs=bombanel_covs,
             mc.cores = chrom.in.parallel)

    # collapse all chromosome outs into a single table
    setwd(caseout_folder)
    tabs_list = list.files(caseout_folder,full.names = TRUE,recursive = TRUE,pattern = 'chrpm_f1.tsv')
    if(length(tabs_list)>0){
      cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out1)
      system(cmd)
    } else {
      cat(file = "f1_table_WARNING.txt","No calls found in chrpm_f1.tsv")
    }
    tabs_list = list.files(caseout_folder,full.names = TRUE,recursive = TRUE,pattern = 'chrpm_f2.tsv')
    if(length(tabs_list)>0){
      cmd = paste('cat',paste(tabs_list,collapse = ' '),'>>',out2)
      system(cmd)
    } else {
      cat(file = "f2_table_WARNING.txt","No calls found in chrpm_f2.tsv")
    }
    tabs_list = list.files(caseout_folder,full.names = TRUE,recursive = TRUE,pattern = 'chrpm_f3.tsv')
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

  message(paste("[",Sys.time(),"]\talright."))
  return(list(tabsnvs_index=tabsnvs_index))
}
