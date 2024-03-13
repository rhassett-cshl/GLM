#### load packages ####
suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.UCSC.hg38)
})

### function to generate .bw of predicted zeta based on epigenomic model
### plus with + values, minus with - values
pred_ep_zeta = function(shared_gb, rc, ep, ep_kappa){
  # prepare gb: gene body coordinates + RC + EP
  # DON'T FORGET to remove 'width' if there is any 
  gb = shared_gb %>%
    dplyr::bind_cols(rc, ep) %>% 
    dplyr::select(-dplyr::any_of(c('width')))

  # Yji contains gene_id, xji and features
  Yji <- gb %>%
    dplyr::select(5:last_col())

  # the order of ep column
  ep_order = Yji %>%
    dplyr::select(3:last_col()) %>%
    colnames()

  # get kappa based on cell line and make the order follow Yji col order
  k = ep_kappa %>%
    dplyr::filter(stringr::str_detect(type, cell)) %>%
    dplyr::mutate(feature = dplyr::case_when(
      feature == '3\' spl' ~ 'sj3',
      feature == '5\' spl' ~ 'sj5',
      feature == 'CTCF' ~'ctcf',
      feature == 'DNAm' ~'wgbs',
      startsWith(feature, 'H') == T ~ tolower(feature),
      feature == 'low complx' ~ 'rpts'
    )) %>%
    dplyr::arrange(match(feature, ep_order)) %>% # feature order should match between Yji and kappa
    dplyr::pull(kappa)


  ## get pred zeta
  power <- Yji %>%
    dplyr::select(3:last_col()) %>%
    as.matrix(.) %*% k %>%
    as.vector()

  ## get zeta
  pred_zeta <- exp(power)


  ########### transfer to grng and prepare it for saving to bigwig file ########
  ## prepare gb+zeta
  gb_zeta = shared_gb %>%
    dplyr::select(-ensembl_gene_id) %>%
    dplyr::mutate(score = pred_zeta)

  ## change sequence type as 'chr1' 
  gb_zeta = plyranges::as_granges(gb_zeta)
  GenomeInfoDb::seqlevelsStyle(gb_zeta) = 'UCSC'

  ## map seqlevels of annotation hg38 and give seqinfo to zeta grng
  hs = BSgenome.Hsapiens.UCSC.hg38
  GenomeInfoDb::seqlevels(gb_zeta) = GenomeInfoDb::seqlevels(hs)
  GenomeInfoDb::seqinfo(gb_zeta) = GenomeInfoDb::seqinfo(hs)
  
  ## split strands: '+' and '-'
  gb_zeta_plus = gb_zeta[strand(gb_zeta) == '+']
  gb_zeta_minus = gb_zeta[strand(gb_zeta) == '-']

  ## assign zeta on '-' strand with negative values
  gb_zeta_minus$score = gb_zeta_minus$score * -1
  
  ## remove the positions that appears in plus strand from minus strand
  dedup_zeta_minus = gb_zeta_minus[! gb_zeta_minus %over% gb_zeta_plus]
  
  ## merge plus and minus
  dedup_gb_zeta = c(gb_zeta_plus, dedup_zeta_minus)

  ## return "deduplicated" zeta prediction with grng 
  return(dedup_gb_zeta)
  
}




### function to generate .bw of predicted zeta based on combined model
### plus with + values, minus with - values
pred_epAllmer_zeta = function(shared_gb, rc, ep, allmer, epAllmer_kappa){
  # prepare gb: gene body coordinates + RC
  # DON'T FORGET to remove 'width' if there is any 
  gb = shared_gb %>%
    dplyr::bind_cols(rc) %>% 
    dplyr::select(-dplyr::any_of(c('width')))
  
  # get kappa based on cell line and make the order follow Yji col order
  k = epAllmer_kappa %>%
    dplyr::filter(stringr::str_detect(type, cell)) %>%
    dplyr::pull(kappa)
  
  ## start the combination calculation of two sets of features
  # y1 as the allmer features
  y1  = allmer
  
  # y2 as pre-normalized epi features, needs to be transformed into matrix
  y2 = ep %>% as.matrix()
  
  # y1 linear transformation
  norm_item = calculate_norm_item(y1)
  c1 = norm_item$c1
  c2 = norm_item$c2
  
  # separate k into k1 and k2 for allmer and ep, respectively
  k1 = k[1:ncol(y1)]
  k2 = k[(ncol(y1)+1):length(k)]
  
  # y1 based calculation
  expNdot1 <- calculate_expNdot_norm(k1, y1, c1, c2)
  
  # y2 based calculation
  expNdot2 <- calculate_expNdot(k2, y2)
  
  # combined expNdot
  expNdot = expNdot1 * expNdot2
  
  ## get pred zeta, which should be 1/expNdot_norm
  pred_zeta <- 1/expNdot %>%
    as.vector()
  
  ## prepare gb+zeta
  gb_zeta = shared_gb %>%
    dplyr::select(-ensembl_gene_id) %>%
    dplyr::mutate(score = pred_zeta)
  
  ## change sequence type as 'chr1' 
  gb_zeta = plyranges::as_granges(gb_zeta)
  GenomeInfoDb::seqlevelsStyle(gb_zeta) = 'UCSC'
  
  ## map seqlevels of annotation hg38 and give seqinfo to zeta grng
  hs = BSgenome.Hsapiens.UCSC.hg38
  GenomeInfoDb::seqlevels(gb_zeta) = GenomeInfoDb::seqlevels(hs)
  GenomeInfoDb::seqinfo(gb_zeta) = GenomeInfoDb::seqinfo(hs)
  
  ## split strands: '+' and '-'
  gb_zeta_plus = gb_zeta[strand(gb_zeta) == '+']
  gb_zeta_minus = gb_zeta[strand(gb_zeta) == '-']
  
  ## assign zeta on '-' strand with negative values
  gb_zeta_minus$score = gb_zeta_minus$score * -1
  
  ## remove the positions that appears in plus strand from minus strand
  dedup_zeta_minus = gb_zeta_minus[! gb_zeta_minus %over% gb_zeta_plus]
  
  ## merge plus and minus
  dedup_gb_zeta = c(gb_zeta_plus, dedup_zeta_minus)
  
  ## return "deduplicated" zeta prediction with grng 
  return(dedup_gb_zeta)
}
