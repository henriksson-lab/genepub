getgenome:
	cd input; wget ftp://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz


getuniprot:
	wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000589_10090.fasta.gz
	#scp UP000000589_10090.fasta.gz  beagle.henlab.org:/data/henlab/ref_genome/mouse/uniprot.fa.gz   to /homology/uniprot.fa

put_red_uniprot:
	#the renamed uniprot file for BLASTing
	#scp uniprot.red.fa  beagle.henlab.org:/data/henlab/ref_genome/mouse/uniprot.fa

get_blastp:
	#scp beagle.henlab.org:/data/henlab/ref_genome/mouse/results_prot.out.gz ./mouse_all2lall_blast_prot.csv.gz

make_blastdb:
	makeblastdb -in uniprot.fa -title "mouse_protein" -dbtype prot

runblast:
	blastp -db uniprot.fa -query uniprot.fa -out results_prot.out -outfmt 6



mount:
	sshfs beagle.henlab.org:/data/henlab/project/bias/new/greta ./x

uploadfeature:
	scp greta/feature/* beagle.henlab.org:/data/henlab/project/bias/new/greta/feature/
	#scp totfeature.csv beagle.henlab.org:/data/henlab/project/bias/new/
	#scp features/* beagle.henlab.org:/data/henlab/project/bias/new/features/

downloadfeature:
	scp beagle.henlab.org:/data/henlab/project/bias/new/features/* feature/
	scp beagle.henlab.org:/data/henlab/project/bias/new/totfeature.csv ./


uploadwebsite:
	scp website/data/*sqlite beagle.henlab.org:/data/henlab/project/bias/new/website/data/

getgreta:
	scp -r beagle.henlab.org:/data/henlab/project/bias/new/greta/feature greta/
	scp -r beagle.henlab.org:/data/henlab/project/bias/new/greta/out greta/
