#' Get mitochondrial genes (names) for species
#' @export
GetMitoGenes <- function(species) {
  species <- tolower(x = species)
  if (species == "salamanders") {
    mt.genes <- c(
      "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "MT-TA", "MT-TC", "MT-TD", "MT-TE", "MT-TF", "MT-TG", "MT-TH", "MT-TI", "MT-TK", "MT-TL1", "MT-TL2", "MT-TM", "MT-TN", "MT-TP", "MT-TQ", "MT-TR", "MT-TS1", "MT-TS2", "MT-TT", "MT-TV", "MT-TW", "MT-TY", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "RNR1", "RNR2"
    )
  } else if (species == "human") {
    mt.genes <- c("MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "MT-RNR1", "MT-RNR2", "MT-TA", "MT-TC", "MT-TD", "MT-TE", "MT-TF", "MT-TG", "MT-TH", "MT-TI", "MT-TK", "MT-TL1", "MT-TL2", "MT-TM", "MT-TN", "MT-TP", "MT-TQ", "MT-TR", "MT-TS1", "MT-TS2", "MT-TT", "MT-TV", "MT-TW", "MT-TY")
  } else if (species == "mouse") {
    mt.genes <- c("mt-Atp6", "mt-Atp8", "mt-Co1", "mt-Co2", "mt-Co3", "mt-Cytb", "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd4l", "mt-Nd5", "mt-Nd6", "mt-Rnr1", "mt-Rnr2", "mt-Ta", "mt-Tc", "mt-Td", "mt-Te", "mt-Tf", "mt-Tg", "mt-Th", "mt-Ti", "mt-Tk", "mt-Tl1", "mt-Tl2", "mt-Tm", "mt-Tn", "mt-Tp", "mt-Tq", "mt-Tr", "mt-Ts1", "mt-Ts2", "mt-Tt", "mt-Tv", "mt-Tw", "mt-Ty")
  } else if (species == "ferret") {
    mt.genes <- c("MT-ATP8", "MT-CO1", "MT-CYTB", "MT-ND2", "MT-ND5", "MT-RNR1", "MT-RNR2", "MT-TM", "MT-TS1", "MTATP6P8", "MTCO1P12", "MTCO3P43", "MTCYBP40", "MTND1P15", "MTND4P14", "MTND5P20", "MTRNR2L10")
  } else if (species == "pig") {
    mt.genes <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6")
  } else if (species == "chicken") {
    mt.genes <- c("ATP6", "ATP8", "COII", "COX3", "CYTB", "MT-CO1", "MT-ND2", "ND1", "ND3", "ND4", "ND4L", "ND5", "ND6")
  } else {
    stop(paste0("Unsupported species", species))
  }

  return(mt.genes)
}
