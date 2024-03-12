### function to generate plus.bw/minus.bw based on epigenomic model
### plus with + values, minus with - values
pred_ep_zeta = function(shared_gbsg_rc){
  ## prepare gb: gene body coordinates + RC + EP
  gb = shared_gbsg %>%
    dplyr::bind_cols(rc, ep)

  # Yji contains gene_id, xji and features
  Yji <- gb %>%
    dplyr::select(5:last_col())

  # the order of ep column
  ep_order = Yji %>%
    dplyr::select(3:last_col()) %>%
    colnames()

  # get kappa based on cell line and make the order follow Yji col order
  k = all_oriK_df %>%
    dplyr::filter(stringr::str_detect(type, cell)) %>%
    dplyr::mutate(feature = dplyr::case_when(
      feature == '3\' spl' ~ 'sj3',
      feature == '5\' spl' ~ 'sj5',
      feature == 'CTCF' ~'ctcf',
      feature == 'DNAm' ~'wgbs',
      startsWith(feature, 'H') == T ~ tolower(feature),
      feature == 'low complx' ~ 'rpts'
    )) %>%
    dplyr::arrange(match(feature, ep_order)) %>%
    dplyr::pull(mean_coef)


  # ep_order = c("ctcf", "h4k20me1", "h3k79me2", "h3k4me1", "h3k9me3", "h3k36me3",
  #              "sj5", "sj3", "rpts", "wgbs")


  ## get pred zeta
  power <- Yji %>%
    dplyr::select(3:last_col()) %>%
    as.matrix(.) %*% k %>%
    as.vector()

  ## get zeta
  pred_zeta <- exp(power)

  ## print see
  print("finish predicting")


  ## prepare gb+zeta
  gb_zeta = shared_gbsg %>%
    dplyr::select(-ensembl_gene_id) %>%
    dplyr::mutate(score = pred_zeta)

  ## change sequence type
  gb_zeta$seqnames = paste0('chr', gb_zeta$seqnames)
  gb_zeta = plyranges::as_granges(gb_zeta)

  ## see print
  print(head(gb_zeta))

  ## give seqinfo
  library(BSgenome.Hsapiens.UCSC.hg38)
  hs = BSgenome.Hsapiens.UCSC.hg38
  seqinfo(gb_zeta) = seqinfo(hs)

  ## output bigwig path
  bw_plus_out = paste0(comp_dir, '/pred_zeta/', cell, '_predZeta_plus.bw')
  bw_minus_out = paste0(comp_dir, '/pred_zeta/', cell, '_predZeta_minus.bw')
  # bw_out = paste0(comp_dir, '/pred_zeta/', cell, '_predZeta.bw')

  gb_zeta_plus = gb_zeta[strand(gb_zeta) == '+']
  gb_zeta_minus = gb_zeta[strand(gb_zeta) == '-']

  gb_zeta_minus$score = gb_zeta_minus$score * -1


  ## export
  rtracklayer::export.bw(gb_zeta_plus, bw_plus_out)
  rtracklayer::export.bw(gb_zeta_minus, bw_minus_out)

  ## export
  # gb_zeta = c(gb_zeta_plus, gb_zeta_minus)
  # rtracklayer::export.bw(gb_zeta, bw_out)
}
