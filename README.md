# Using bioinformatics to identify and understand microbial communities. BIOL6220 Project by Kelsey Coates

## About the project: Plants and bacteria form mutualistic symbiotic relationships to process nutrients. There are bacteria capable of fixing nitrogen without plant hosts (free-living nitrogen fixers (FLNF)), but these bacteria are understudied outside of agricultural research. Human nutrient addition to wetland soils may alter how bacteria process nitrogen by disrupting mutualistic relationships and impacting the prevalence of different functional groups in soil communities (mutualists versus FLNF). Considering our lab's interest in microbial communities, their functions, and human impacts on these communities we want to know if there are FLNF in our local wetland soils? We can test this at our West Research Campus (WRC) long-term fertilization and disturbance experimental site. Our lab's workflow involves collecting soil samples from WRC where vegetated plots are fertilized, disturbed, both, or neither (control), extracting DNA from those soil samples, amplifying the 16S region rRNA region, and using bioinformatics to generate sequences. By combining this workflow with published data on FLNF genes, we can determine if bacteria present in WRC soils also have the ability to fix nitrogen independently of plant hosts.  

'''
### The first step was using Windowd Command Terminal to navigate onto the virtual computer/ server and navigating into the ECU bio department directory -> Peralta lab -> Kelsey. 
 cd kelsey
 ls
### We used wget to get reference bacterial sequence database from online (Mothur, software for processing bioinformatic information for microbial communities). Then unzipping the reference database file. We needed to get arb software onto the server to work with the .arb file.
gunzip SILVA_132_SSURef_NR99_13_12_17_opt.arb.gz
ls
less -S SILVA_132_SSURef_NR99_13_12_17_opt.arb
 arb SILVA_132_SSURef_NR99_13_12_17_opt.arb
 mamba install arb
 mamba install arb -c bioconda -c conda-forge
 micromamba install arb -c bioconda -c conda-forge
 arb SILVA_132_SSURef_NR99_13_12_17_opt.arb
 micromamba activate /home/peraltalab/miniconda3/
 eval "$(micromamba shell hook --shell bash)"
 micromamba activate
 arb
 micromamba activate /home/peraltalab/miniconda3/
 arb
 su
 conda deactivate
 conda activate base
  ### Arb and conda were having a hard time working together, so we decided to try a different reference database from Mothur
  wget https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_LSURef_tax_silva.fasta.gz
  gunzip SILVA_132_LSURef_tax_silva.fasta.gz
  ### We used SILVA database 138 to create a BLAST file containing only the 16S sequences from the NFix database. (In the near future, I will try this command again with SILVA 132, and see if that fixes downstream compatibility issues)
  blastall -p blastn -m 8 -d ../SILVA_dbs/SILVA_138.2_LSURef_tax_silva.BLAST.db -i 16S_only_seqs.fasta -o 16S_only_seqs.fasta_BLASTm8.SILVA.out
  ### We copied the SILVA reference database into the Peralta lab directory for future use 
  cp SILVA_132_LSURef_tax_silva.fasta ../SILVA_dbs/
  ### We made a BLAST file from the reference database and continue to put the file in a format that will be compatible with Mothur's software and our established workflow. 
  makeblastdb -in SILVA_132_LSURef_tax_silva.fasta -out SILVA_132_LSURef_tax_silva.fasta_BLASTdb -dbtype nucl
  history | grep blastall
  blastall -p blastn -m 8 -d SILVA_132_LSURef_tax_silva.fasta_BLASTdb -i 16S_only_seqs.fasta -o 16S_only_seqs.fasta_BLASTm8.SILVA132.out
  ls
 ### We used the sed command to remove the quotes that interrupted the sequences then piped that into a csv. In R, we used tidyverse to filter the csv for sequences that were >95% similar to the sequences in our SILVA 132 reference database. Back in windows command shell, we take the sequences with >95% similarity and pipe them into a text file. 
sed "s/\"//g" 16S_only_seqs.fasta_BLASTm8.SILVA132.out.95percent.csv
 sed "s/\"//g" 16S_only_seqs.fasta_BLASTm8.SILVA132.out.95percent.csv > 16S_only_seqs.fasta_BLASTm8.SILVA132.out.95percent_noQuotes.csv
  history | grep "noQuotes"
  cut -f2 -d"," 16S_only_seqs.fasta_BLASTm8.SILVA132.out.95percent_noQuotes.csv > SILVA132_Nfixers95.txt
  ### We wanted to confirm the text file looked how we wanted. 
 less -S SILVA132_Nfixers95.txt
 tail -1 SILVA132_Nfixers95.txt
 head 2+ SILVA132_Nfixers95.txt
tail +2 SILVA132_Nfixers95.txt
 tail +2 SILVA132_Nfixers95.txt > SILVA132_ Nfixers95.txt
 ### We selected the second column of the text file containing the sequence IDs
 cut -f2 -d"," 16S_only_seqs.fasta_BLASTm8.SILVA132.out.95percent_noQuotes.csv | tail +2 > SILVA132_Nfixers95.txt
 less -S SILVA132_Nfixers95.txt
 history | grep "awk"
 ### The text file was interwoven with text and line breaks we didn't want so we used the awk command to remove the junk, and have continuous sequences running with no line breaks
  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < SILVA_132_LSURef_tax_silva.fasta > SILVA_132_LSURef_tax_silva.fasta_singleline.fasta
  less -S SILVA_132_LSURef_tax_silva.fasta_singleline.fasta
  history | grep "noQuotes"
  ### We created the final fasta file, preparing for comparing our sequences generated from WRC soils to sequences in the NFix database, that have been BLASTed against the SILVA database. 
  while read line; do grep -A1 "$line" SILVA_138.2_LSURef_tax_silva.fasta_singleline.fasta ; done < SILVA132_Nfixers95.txt
  while read line; do grep -A1 "$line" SILVA_138.2_LSURef_tax_silva.fasta_singleline.fasta ; done < SILVA132_Nfixers95.txt > SILVA132_Nfixers95.fasta
  while read line; do grep -A1 "$line" SILVA_132_LSURef_tax_silva.fasta_singleline.fasta ; done < SILVA132_Nfixers95.txt > SILVA132_Nfixers95.fasta
  grep -c "^>" SILVA132_Nfixers95.fasta
### Next Steps: 1) Use fasta files in Mothur software to determine if our lab’s soil microbes sequences from WRC share FLNF genes from the NFix database. 2) Identify FLNF taxon in our microbial community data set via sequences and/or names. 3) Graph the prevalence of these FLNF taxon in fertilized versus unfertilized plots to understand their role in nitrogen cycling using R. 4) Graph the prevalence of these FLNF taxon over time to understand how fertilization and disturbance (mowing) treatments might impact their presence and function. 

