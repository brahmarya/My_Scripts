
# filter the pdb names which contains anomalies (e.g. MODRES, Flip, Kink)

grep -i "MODRES \| kink \| flip " RCSB/fltr_pdb/* | cut -c15-22 | sort | uniq  > list_of_DNA_containing_anomalies.txt

file=`cat list_of_DNA_containing_anomalies.txt`
for f in $file; do rm /storage/ICMR/RCSB/fltr_pdb/$f; done
for f in $file; do name="${f%%.pdb}"; echo $name; rm /storage/ICMR/RCSB/fltr_fasta/${name^^}.fasta.txt; done


# cluster the filtered fasta sequnces

cat /storage/ICMR/RCSB/fltr_fasta/* > fltr_fasta.fasta
dnaclust fltr_fasta.fasta  > fltr_cluster.txt


# Shell script to calculate ASA of DNA using NACCESS tool.

rm -rf /storage/ICMR/RCSB/ASA
file=`ls /storage/ICMR/RCSB/fltr_pdb`
for f in $file; do name="${f%%.pdb}"; echo $name; `mkdir -p /storage/ICMR/RCSB/ASA/$name`; cd /storage/ICMR/RCSB/ASA/$name; `pwd` ;  /raid/ajay/Tools/Naccess/naccess  /storage/ICMR/RCSB/fltr_pdb/$f ; done

file=`cat /storage/ICMR/protein_list.txt`
for f in $file; do name="${f%%.pdb}"; echo $name; `mkdir -p /storage/ICMR/data/ASA/$name`; cd /storage/ICMR/data/ASA/$name; `pwd` ;  /raid/ajay/Tools/Naccess/naccess  /storage/ICMR/RCSB/PDB/$f ; done

# trim the *.rsa file and combine all DNA sequences

cd /storage/ICMR/data/ASA
grep "^RES" */*.rsa | cut --complement -c5-17 > ../all_dna.rsa # This command should be run within the folder.

cd /storage/ICMR/RCSB/ASA
grep "^RES" */*.rsa | cut --complement -c5-17 > ../all_dna.rsa # This command should be run within the folder.


