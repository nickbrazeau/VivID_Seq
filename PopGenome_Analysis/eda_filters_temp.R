smpls <- readxl::read_excel("~/Documents/GitHub/VivID_Seq/scrape_pubseqs/vivid_seq_public_NGS.xlsx") %>%
  dplyr::select(c("iid", "acc", "country", "vividregion"))

drc <- data.frame(iid = c("D9U3K", "O6Y4K", "Q8J6O"),
                  acc = c("D9U3K", "O6Y4K", "Q8J6O"),
                  country = "CD",
                  vividregion = "AF")

smpls <- rbind.data.frame(smpls, drc) %>%
  dplyr::mutate(country = ifelse(vividregion == "NHA", "NHA", country)) %>%
  dplyr::filter(vividregion != "Lab")

vcf <- vcfR::read.vcfR(file = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/vcfs_hardfilt_variants/passed.joint.vcf.gz")
vcf.snp <- vcfR::extract.indels(vcf, return.indels = F)
vcf.snp

#........................................................
# Find samples with a high level of heterozygosity
# which mimicks COI > 1 but is really heteroplasmy
# or sequencing error
#........................................................
smplhet <- apply(vcfR::is.het(vcfR::extract.gt(vcf.snp)), 2, mean)
smplhet <- tibble::tibble(smpls = names(smplhet),
                          hetprop = smplhet)
smplhetpass <- smplhet$smpls[smplhet$hetprop < 0.2] # less than 20%

vcf.snp <- vcfRmanip::select_samples(vcfRobject = vcf.snp,
                                     smplvctr = smplhetpass)

saveRDS(smplhetpass, file = "data/derived_data/final_smpl_list.RDS")


#........................................................
# Remove loci that have deleted allele alternative
# which is encoded in vcf as a *
#........................................................
alts <- vcfR::getALT(vcf.snp)
asterickfilt <- grepl(pattern = "\\*", x = alts)
asterickfilt <- data.frame(seqname = vcfR::getCHROM(vcf.snp),
                           start = vcfR::getPOS(vcf.snp),
                           end = vcfR::getPOS(vcf.snp),
                           geneid = 1:length(vcfR::getPOS(vcf.snp))) %>%
  dplyr::filter(asterickfilt)

vcf.snp <- vcfRmanip::vcffilter_ChromPos(vcfRobject = vcf.snp,
                                     chromposbed = asterickfilt)


#........................................................
# Clean het calls
#........................................................
vcf.snp.hom <- vcfR_force_hets_to_homo(vcf.snp)
vcf.snp.hom.ss <- vcfR_filter_unique_sites(vcf.snp.hom)

#........................................................
# Remove Private Alleles within a country
#........................................................
vcf.snp.hom.ss.long <- vcf.snp.hom.ss %>%
  vcfR::extract_gt_tidy(., format_fields = "GT") %>%
  dplyr::group_by(Indiv) %>%
  dplyr::mutate(POS = vcfR::getPOS(vcf.snp.hom.ss)) %>%
  dplyr::rename(acc = Indiv) %>%
  dplyr::left_join(., smpls, by = "acc") %>%
  dplyr::ungroup()

countryn <- vcf.snp.hom.ss.long %>%
  dplyr::select(c("acc", "country")) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(countryn = n())

privateloci.withincntry <- vcf.snp.hom.ss.long %>%
  dplyr::left_join(., countryn, by = "country") %>%
  dplyr::mutate(I = 1,
                limit = 0.1 * (countryn)/(countryn - 1), # note correction for small sample size
                limit2 = 1/countryn
                ) %>%
  dplyr::group_by(country, POS, gt_GT) %>%
  dplyr::summarise(
    locicount = sum(I)/mean(countryn), # countryn all same
    limit = mean(limit), # limit all same
    limit2 = mean(limit2) # limit all same
  ) %>%
  dplyr::mutate(rm = ifelse(locicount <= limit & !is.na(gt_GT), T, F),
                rm = ifelse(locicount == limit2 & !is.na(gt_GT) & limit2 != 1, T, rm)
                ) %>%
  dplyr::select(c("country", "POS", "gt_GT", "rm"))


vcf.snp.hom.ss.long.nopriv <- vcf.snp.hom.ss.long %>%
  dplyr::left_join(., y = privateloci.withincntry, by = c("country", "POS", "gt_GT")) %>%
  dplyr::filter(rm != T) %>%
  dplyr::select(-c("rm"))


#........................................................
# Final Object
#........................................................
vcf.snp.hom.ss.long.nopriv


jpeg("~/Desktop/mtdna_explore/tempfig.jpg", height = 20, width = 20, units = "in", res = 300)

vcf.snp.hom.ss.long.nopriv %>%
  dplyr::mutate(POS = factor(POS, levels = c(1:5989), ordered = T)) %>%
  ggplot() +
  geom_tile(aes(x = iid, y = POS, fill = gt_GT)) +
  facet_wrap(~country, scales = "free", nrow = 1) +
  theme(axis.text.x = element_text(angle = 90))

graphics.off()

jpeg("~/Desktop/mtdna_explore/cambodia.jpg", height = 20, width = 50, units = "in", res = 300)
vcf.snp.hom.ss.long.nopriv %>%
  dplyr::mutate(POS = factor(POS, levels = c(1:5989), ordered = T)) %>%
  dplyr::filter(country == "KH") %>%
  ggplot() +
  geom_tile(aes(x = iid, y = POS, fill = gt_GT)) +
  facet_wrap(~country, scales = "free", nrow = 1) +
  theme(axis.text.x = element_text(angle = 90))

graphics.off()





jpeg("~/Desktop/mtdna_explore/peru.jpg", height = 20, width = 50, units = "in", res = 300)
vcf.snp.hom.ss.long.nopriv %>%
  dplyr::mutate(POS = factor(POS, levels = c(1:5989), ordered = T)) %>%
  dplyr::filter(country == "PE") %>%
  ggplot() +
  geom_tile(aes(x = iid, y = POS, fill = gt_GT)) +
  facet_wrap(~country, scales = "free", nrow = 1) +
  theme(axis.text.x = element_text(angle = 90))

graphics.off()




jpeg("~/Desktop/mtdna_explore/thailand.jpg", height = 20, width = 50, units = "in", res = 300)
vcf.snp.hom.ss.long.nopriv %>%
  dplyr::mutate(POS = factor(POS, levels = c(1:5989), ordered = T)) %>%
  dplyr::filter(country == "TH") %>%
  ggplot() +
  geom_tile(aes(x = iid, y = POS, fill = gt_GT)) +
  facet_wrap(~country, scales = "free", nrow = 1) +
  theme(axis.text.x = element_text(angle = 90))

graphics.off()









vcf.noape.ss <- vcfR2segsites_gt(vcf.noape[is.biallelic(vcf.noape),])

#........................................................
# Viz and Critique SNPs no ape
#........................................................
vcf.noape.ss.gt <- vcfR::extract_gt_tidy(vcf.noape.ss)

vcf.noape.ss.gt <- vcf.noape.ss.gt %>%
  dplyr::group_by(Indiv) %>%
  dplyr::mutate(POS = vcfR::getPOS(vcf.noape.ss),
                POS = factor(POS, levels = c(1:5898), ordered = T)) %>%
  dplyr::rename(acc = Indiv) %>%
  dplyr::left_join(., smpls, by = "acc")

vcf.noape.ss.gt %>%
  ggplot() +
  geom_tile(aes(x = acc, y = POS, fill = gt_GT)) +
  facet_wrap(~country, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))




#........................................................
# Viz and Critique SNPs no bad hets
#........................................................

smplhet <- apply(vcfR::is.het(vcfR::extract.gt(vcf.noape)), 2, mean)
smplhet <- tibble::tibble(smpls = names(smplhet),
                          hetprop = smplhet)
smplhetpass <- smplhet$smpls[smplhet$hetprop < 0.2] # less than 20%

vcf.noape.nohet <- vcfRmanip::select_samples(vcfRobject = vcf.noape,
                                             smplvctr = smplhetpass)

vcf.noape.nohet.ss <- vcfRmanip::vcfR2segsites_gt(vcf.noape.nohet[is.biallelic(vcf.noape.nohet), ])
vcf.noape.nohet.ss.gt <- vcfR::extract_gt_tidy(vcf.noape.nohet.ss)

vcf.noape.nohet.ss.gt <- vcf.noape.nohet.ss.gt %>%
  dplyr::group_by(Indiv) %>%
  dplyr::mutate(POS = vcfR::getPOS(vcf.noape.nohet.ss),
                POS = factor(POS, levels = c(1:5898), ordered = T)) %>%
  dplyr::rename(acc = Indiv) %>%
  dplyr::left_join(., smpls, by = "acc")

vcf.noape.nohet.ss.gt %>%
  ggplot() +
  geom_tile(aes(x = acc, y = POS, fill = gt_GT)) +
  facet_wrap(~country, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))



# cambodia??
vcf.noape.nohet.ss.gt %>%
  dplyr::filter(country == "KH") %>%
  ggplot() +
  geom_tile(aes(x = acc, y = POS, fill = gt_GT)) +
  facet_wrap(~country, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))


# what about reading in polarized only sites
anc <- Biostrings::readDNAStringSet("analyses/00-find_ancestral_alleles/ancestral.fa")


maskpos <- which( unlist( stringr::str_split(anc, pattern = "") ) == "N")
maskdf <- data.frame(seqname = "PvP01_MIT_v1",
                     start = maskpos,
                     end = maskpos,
                     geneid = 1:length(maskpos))


vcf.noape.nohet.ss.mask <- vcfRmanip::vcffilter_ChromPos(vcfRobject = vcf.noape.nohet.ss,
                                                          chromposbed = maskdf)


vcf.noape.nohet.ss.mask %>%
  vcfR::extract_gt_tidy(.) %>%
  dplyr::group_by(Indiv) %>%
  dplyr::mutate(POS = vcfR::getPOS(vcf.noape.nohet.ss.mask),
                POS = factor(POS, levels = c(1:5898), ordered = T)) %>%
  dplyr::rename(acc = Indiv) %>%
  dplyr::left_join(., smpls, by = "acc") %>%
  ggplot() +
  geom_tile(aes(x = acc, y = POS, fill = gt_GT)) +
  facet_wrap(~country, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))

