# Roary 
roary gff3-input*.gff -f output/roary -p 4 -v 

# Ptolemy 
java -jar ptolemy.jar extract -g ptolemy-input-list.txt -o output/ptolemy
java -jar ptolemy.jar syntenic-anchors --db output -o output/ptolemy
java -jar ptolemy.jar canonical-quiver -s output/ptolemy/syntenic_anchors.txt --db output/ptolemy -o output/ptolemy --dump

# PPanGGoLin
ppanggolin workflow --anno ppanggolin-input-list.txt -output output/ppanggolin
ppanggolin write -p output/ppanggolin/pangenome.h5 --families_tsv -o output/ppanggolin
ppanggolin write -p output/ppanggolin/pangenome.h5 --all_gene_families -o output/ppanggolin -f
ppanggolin write -p output/ppanggolin/pangenome.h5 --all_genes -o output/ppanggolin -f

# PIRATE
./PIRATE -i gff3-input -s "50,70,95" -o output/pirate
perl pirate-scripts/subsample_outputs.pl -i output/pirate/PIRATE.gene_families.tsv -g
output/pirate/modified_gffs/ -o output/pirate/gene_families.prev_locus.tsv --field "prev_locus"
perl pirate-scripts/gene_cluster_to_binary_fasta.pl -i output/pirate/PIRATE.gene_families.tsv output/pirate/binary_presence_absence.fasta
FastTree -fastest -nocat -nome -noml -nosupport -nt output/pirate/binary_presence_absence.fasta > output/pirate/binary_presence_absence.nwk 2>/dev/null

# Panaroo
python convert_refseq_to_prokka_gff.py -g input.gff -f input.fasta -o output.gff # repeat for all inputs
panaroo -i gff3-input/*.gff -o output/panaroo -t 4 -verbose # strict mode
panaroo -i gff3-input/*.gff -o output/panaroo --mode relaxed -t 4 -verbose # relaxed mode

# Combining Panaroo and Ptolemy
panaroo -i gff3-input/*.gff -o output/panaroo --mode relaxed -t 4 -verbose # obtain gene clusters 
java -jar ptolemy.jar extract -g ptolemy-input-list.txt -o output # index graph 
python createSA.py output/panaroo output # create syntenic anchor file input to Ptolemy 
java -jar ptolemy.jar canonical-quiver -s output/syntenic_anchors.txt --db output -o output -dump

# Perform GO enrichment study
python runGOE.py go.obo orf2go.tsv studyinput.txt goestudy.out
