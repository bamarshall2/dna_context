#!/bin/bash


blast_genomes () {
    for species in $(<bacteroides_manifest.txt)
    do hmmsearch --tblout ./results/${species}_${protein}_results.tsv ./proteins/$protein ../orfs/extracted_orfs/${species}
        echo -n $species","$protein"," >> output_bacteroides_nusG.csv
        cat ./results/${species}_${protein}_results.tsv | awk 'NR==4{print$5}' >> output_bacteroides_nusG.csv
done
}


for protein in $(<GOI_manifest.txt)
    do blast_genomes
done
