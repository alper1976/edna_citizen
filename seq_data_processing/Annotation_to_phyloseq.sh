
module load BLAST+/2.13.0-gompi-2022a

RAWDIR=$RAWDIR/Figs
DATABASE_DIR=/cluster/projects/nn9745k/03_databases/fish/ScandiFish_12s_v1.4/
Script_path
cd $RAWDIR

RESDIR=$RAWDIR


for f in *rep-seqs.fna
   do
      FILE=${f#$DIRS}
      if [[ $FILE =~ seqs ]]
      OUTPUT_FILE=${FILE//rep-seqs.fna/output_blast_results}
      then
      	blastn -max_target_seqs 100 -evalue 0.01 -query $FILE -out $OUTPUT_FILE -db $DATABASE_DIR/ScandiFish_12s_v1.4.fasta_db -outfmt 6 -num_threads 2 # BLASTN compares the sequences and keeps up to 100 reference sequences. # Input is asv_table in fasta format, DADA2; uniquesToFasta(seqtab.nochim) function
      fi
   done

echo " "
echo " "
echo "Finished with blast annotation"
echo " "
echo " "
module purge

module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2



Rscript $Script_path/Blast_annotation.R -i $RAWDIR -o $RESDIR -e $Error_correction_model # starts r script to create sequence table in dada2 and phyloseq compatible format.

#works!




wget -nc https://raw.githubusercontent.com/naturalis/galaxy-tool-lca/master/lca.py -P $Script_path #downloads lca.py if it do not exist

rm $RESDIR/*lca_03_98_out.tabular # Removes table if existing LCA analysis exist.

## LCA analysis
## Best hit goes to species, taxonomic sorting is conducted if top hit <99%
Input_FILE=$RAWDIR/*.lca.tabular
cd $RAWDIR

for f in *.lca.tabular
   do
      FILE=${f#$DIRS}
      echo $FILE
      if [[ $FILE =~ $Error_correction_model ]]
		OUTPUT_FILE=${FILE//.lca.tabular/.lca.03_98.tabular}
		rm $OUTPUT_FILE
		echo $OUTPUT_FILE
      then
      	python $Script_path/lca.py -i $FILE -o $OUTPUT_FILE -b 8 -id 98 -cov 95 -t best_hit -tid 99 -tcov 100 -flh unknown
		fi
	done

echo " "
echo " "
echo "Finished with LCA annotation"
echo " "
echo " "

#RAWDIR=$Input_path
cd $RAWDIR
#RESDIR=$Input_path/Figs

Rscript $Script_path/Create_phyloseq_object.R -i $RESDIR -e $Error_correction_model
echo " "
echo " "
echo "Finished to make Phyloseq object"
echo " "
echo " "