import magenpy as mgp
import viprs as vp


# Load genotype and GWAS summary statistics data (chromosome 22):
gdl = mgp.GWADataLoader(bed_files=mgp.tgp_eur_data_path(),
                        sumstats_files=mgp.ukb_height_fastGWA_path(),
                        sumstats_format="fastGWA")
