#' add_names
#'
add_names <- function(pm,
                      name.patient,
                      name.plasma,
                      name.germline){
  pm$PatientID = name.patient
  pm$CaseID = name.plasma
  pm$GermlineID = name.germline
  pm = pm[with(pm,order(chr,pos)),]
  return(pm)
}
