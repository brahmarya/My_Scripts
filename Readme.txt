


Step 1: Download the data from NDB and RCSB-PDB database.
	NDB: only DNA sequnces has been downloaded
	RCSB: Advance Search options.
		Search Parameter:
		"Chain Type: there is a DNA chain but not any Protein or RNA or Hybrid and Experimental Method is X-RAY and Resolution is between 0.0 and 2.99"
		On 20 Feb, 2020

	Note: A DNA list from a previous work is also selected.( Supplymentry Data ) 
		"Andrabi, M., Hutchins, A.P., Miranda-Saavedra, D. et al. Predicting conformational ensembles and genome-wide transcription factor binding sites from DNA sequences. Sci Rep 7, 4071 (2017). https://doi.org/10.1038/s41598-017-03199-6"
		
############################ Further steps are for RCSB data #####################################
		
Step 2: Filtering the downloaded data.
	- Only those DNA were selected which had only two chains.
	- Any DNA sequences that were less than 10 in length and did not have equal length of both chains were removed.
	- Any DNA sequence/chain containing following characters "x","n","X","N","I","i","U","u" were filtered.
	
	-- "/storage/obj2-29March2019/src/filter_dna.R" # R code to filter DNA sequences.
	
Step 3: Filtered PDB files still contains anomalies like: MODRES, Flip, Kink in DNA sequences.
	- These anomalies were greped and removed the pdb files contains.
	-- grep -i "MISSING\|MODRES\|kink\|flip" /storage/obj2-29March2019/RCSB/fltr_pdb/* | cut -c41-48 | sort | uniq  > /storage/obj2-29March2019/RCSB/list_of_DNA_containing_anomalies.txt
	-- file=`cat /storage/obj2-29March2019/RCSB/list_of_DNA_containing_anomalies.txt`
	-- for f in $file; do rm /storage/obj2-29March2019/RCSB/fltr_pdb/$f; done
	-- for f in $file; do name="${f%%.pdb}"; echo $name; rm /storage/obj2-29March2019/RCSB/fltr_fasta/${name^^}.fasta.txt; done
	
Step 4: Calculate ASA of DNA using NACCESS tool
	-- rm -rf /storage/obj2-29March2019/RCSB/ASA
	-- file=`ls /storage/obj2-29March2019/RCSB/fltr_pdb`
	-- for f in $file; do name="${f%%.pdb}"; echo $name; `mkdir -p /storage/obj2-29March2019/RCSB/ASA/$name`; cd /storage/obj2-29March2019/RCSB/ASA/$name; `pwd` ;  /raid/ajay/Tools/Naccess/naccess  /storage/obj2-29March2019/RCSB/fltr_pdb/$f ; done

Step 5: Trim the *.rsa file and combine all DNA sequences
	-- cd /storage/obj2-29March2019/RCSB/ASA
	-- grep "^RES" */*.rsa | cut --complement -c5-17 > ../all_dna.rsa

Step 6: Manually insert space between the chain and residue positions in all_dna.rsa file.

Step 7: Cluster the filtered fasta sequences
	-- cat /storage/obj2-29March2019/RCSB/fltr_fasta/* > /storage/obj2-29March2019/RCSB/fltr_fasta.fasta
	-- dnaclust /storage/obj2-29March2019/RCSB/fltr_fasta.fasta  > /storage/obj2-29March2019/RCSB/fltr_cluster.txt
	
Step 8: Calculate the Dynamic and Static sequence features using dynaseq tool
	- Note: To add DNA names as rownames use "/storage/obj2-29March2019/src/mydynaseqR.R" function.
	-- /storage/obj2-29March2019/src/cr_dynaseq_features.R   # To create dynaseq ( dynamic and static ) features for further use
	-- /storage/obj2-29March2019/src/cr_pattern.R 			# To create pattern files for all window sizes. ( One-hot encoding + dynaseq features )
	-- /storage/obj2-29March2019/src/cr_pattern_pair.R   	# To create pattern files for base pairs. (No window size is taken)
	-- /storage/obj2-29March2019/src/cr_pattern_step.R		# To create pattern files for base steps. (No window size is taken)
	
	Note: "/storage/obj2-29March2019/src/cr_features.R"      # This script to calculate dynseq features for only window size 5 DNA, directly from the static dyna-seq feature matrix 
	
Step 9: Write diffrent DNA sequence patterns into single file.
	-- awk 'FNR==1 && NR!=1 { while (/^ptrn/) getline; } 1 {print}' /storage/obj2-29March2019/RCSB/pat/win_1/*:A.tsv > /storage/obj2-29March2019/RCSB/pat/win_1_A.tsv  # similarly for all window size
	-- awk 'FNR==1 && NR!=1 { while (/^ptrn/) getline; } 1 {print}' /storage/obj2-29March2019/RCSB/pat/win_1/*:B.tsv > /storage/obj2-29March2019/RCSB/pat/win_1_B.tsv	# similarly for all window size
	-- awk 'FNR==1 && NR!=1 { while (/^"name/) getline; } 1 {print}' /storage/obj2-29March2019/RCSB/pat_pair/* > /storage/obj2-29March2019/RCSB/pattern_pair.tsv
	-- awk 'FNR==1 && NR!=1 { while (/^"base_step/) getline; } 1 {print}' /storage/obj2-29March2019/RCSB/pat_step/* > /storage/obj2-29March2019/RCSB/pattern_step.tsv
	
Step 10: 
	
	
	



























bedops --everything file1.bed file2.bed ... fileN.bed  | bedmap --echo-map - | awk '(split($0, a, ";") == 3)' - | sed 's/\;/\n/g' - | sort-bed - | uniq - > answer.bed


bedops -u file1.bed file2.bed ... fileN.bed | bedmap --echo --echo-map-id-uniq --fraction-both 0.5 - | awk -F"|" '(split($2, a, ";") > 1)' > answer.bed
















	
	
