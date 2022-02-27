
# convert Fastq to VCF and Fasta files.  ----------------------------------

path <- dir("/storage/Collab/Mohit/Clinical_Data", full.names = TRUE)
file <- sapply(path,function(x) strsplit(x,"[/-]")[[1]][6])
fullPath <- paste(path, file,sep = "/")


lapply(1:28, function(n){
   fq <- list.files(fullPath[n], pattern = "gz")
   show(fq)
   Path <- paste(fullPath[n],"/", sep = "")

# Fastq file alignment using bwa ------------------------------------------
   
   bwa_cmd <- paste("bwa mem /storage/Collab/Mohit/Reference_genome/Mycobacterium_tuberculosis_H37Rv.fasta ",Path,fq[1], " ", Path,fq[2], " > ",Path,"Aligned_PE.sam", sep = "")
   system(bwa_cmd)
   system("PID=$!; wait $PID")
   
# SAM to BAM conversion ---------------------------------------------------
   sam_cmd <- paste("samtools view -S -b ",Path,"Aligned_PE.sam > ",Path,"Aligned_PE.bam", sep = "/")
   system(sam_cmd)
   system("PID=$!; wait $PID")

# Sorting the BAM file ----------------------------------------------------
   samt_cmd <- paste("samtools sort ",Path,"Aligned_PE.bam -o ",Path,"Aligned_PE_sorted.bam", sep = "")
   system(samt_cmd)
   system("PID=$!; wait $PID")

# Mark and remove the PCR duplicates using PICARD -------------------------
   pica_cmd <- paste("picard MarkDuplicates I=",Path,"Aligned_PE_sorted.bam O=",Path,"Aligned_PE_sorted_marked.bam M=",Path,"picard_info.txt REMOVE_DUPLICATES=true AS=true", sep = "")
   system(pica_cmd)
   system("PID=$!; wait $PID")

# Indexing the sorted marked file --------------------------------------------
   sami_cmd <- paste("samtools index ",Path,"Aligned_PE_sorted_marked.bam > ", Path,"Aligned_PE_sorted_marked.bam.bai" , sep = "")
   system(sami_cmd)
   system("PID=$!; wait $PID")

# Correction of alignments in marked files -----------------------------------
   picam_cmd <- paste("picard CollectMultipleMetrics I=",Path,"Aligned_PE_sorted_marked.bam O=",Path,"alignmentsummary R=/storage/Collab/Mohit/Reference_genome/Mycobacterium_tuberculosis_H37Rv.fasta", sep = "")
   system(picam_cmd)
   system("PID=$!; wait $PID")

# Generating vcf files using freebayes ------------------------------------
   fb_cmd <- paste("freebayes -f /storage/Collab/Mohit/Reference_genome/Mycobacterium_tuberculosis_H37Rv.fasta ",Path,"Aligned_PE_sorted_marked.bam > ",Path,"Rv_freebayes.vcf", sep = "")
   system(fb_cmd)
   system("PID=$!; wait $PID")

# VCF to fasta conversion -------------------------------------------------
   bg_cmd <- paste("bgzip -f ",Path,"Rv_freebayes.vcf > ",Path,"Rv_freebayes.vcf.gz", sep = "")
   system(bg_cmd)
   system("PID=$!; wait $PID")
   
   tb_cmd <- paste("tabix -f ",Path,"Rv_freebayes.vcf.gz > ",Path,"Rv_freebayes.vcf.gz.tbi", sep = "")
   system(tb_cmd)
   system("PID=$!; wait $PID")
   
   fst_cmd <- paste("bcftools consensus -f /storage/Collab/Mohit/Reference_genome/Mycobacterium_tuberculosis_H37Rv.fasta ",Path,"Rv_freebayes.vcf.gz > /storage/Collab/Mohit/Fasta/", file[n], ".fasta", sep = "")
   system(fst_cmd)
})
