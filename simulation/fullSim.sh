# Make the json by altering the prefix field
PREFIX=/home/christian/Research/Integrative/prsPPMx/temp/simulation

# in test.json replace the "prefix" field with the above prefix
sed -i "s|\"prefix\":.*|\"prefix\": \"$PREFIX\"|g" test.json


for i in {0..4}
do
    cd /home/christian/Research/Stat_gen/tools/MASH
    # Simulate

    python -m Simulate \
     --argfile /home/christian/Research/Integrative/prsPPMx/test.json
    cd /home/christian/Research/Integrative/prsPPMx

    # GWAS
    # Mkae it space delimited 
    awk '{$7=""; print $0}' ${PREFIX}.fam > ${PREFIX}.fam2
    mv ${PREFIX}.fam2 ${PREFIX}.fam

    # ldmatrix 
    plink --bfile ${PREFIX} --r2 --out ${PREFIX} 

    # GWAS
    plink2 --bfile ${PREFIX} --glm --out ${PREFIX} --covar ${PREFIX}.covar --pheno ${PREFIX}.pheno

    # PCs
    plink2 --bfile ${PREFIX} --pca 10 --out ${PREFIX}

    for i in {0..4}
    do
      awk -F'\t' 'NR==1 || $7 == "ADD" {print}' ${PREFIX}.Y${i}.glm.linear > ${PREFIX}.Y${i}.glm.linear.filtered
      PRSice_linux --base ${PREFIX}.Y${i}.glm.linear.filtered --score avg --target ${PREFIX} --out ${PREFIX}.Y${i} \
        --A1 A1 --A2 ALT --stat BETA --snp ID --bp POS --pvalue P --keep-ambig --pheno ${PREFIX}.pheno --pheno-col Y${i}
    done

    awk -F'\t' 'NR==1 || $7 == "ADD" {print}' ${PREFIX}.Z.glm.linear > ${PREFIX}.Z.glm.linear.filtered
    PRSice_linux --base ${PREFIX}.Z.glm.linear.filtered --score avg --target ${PREFIX} --out ${PREFIX}.Z \
      --A1 A1 --A2 ALT --stat BETA --snp ID --bp POS --pvalue P --keep-ambig --pheno ${PREFIX}.pheno --pheno-col Z

    # Run estimation
    Rscript simulation/ppmx.R --filepath ${PREFIX}
done

