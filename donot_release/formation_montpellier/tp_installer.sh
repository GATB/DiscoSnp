
# GET THE DATA
wget http://igm.univ-mlv.fr/~peterlon/humch1_00096_reads.fasta.gz
wget http://igm.univ-mlv.fr/~peterlon/humch1_00100_reads.fasta.gz
wget http://igm.univ-mlv.fr/~peterlon/humch1_first_10M.fasta
wget http://igm.univ-mlv.fr/~peterlon/ref_human

# GET THE TOOLS
mkdir tp_tools
cd tp_tools
wget http://igm.univ-mlv.fr/~peterlon/tools.zip

# INSTALL THE TOOLS
unzip tools.zip 
cd tools/gassst
pwd
make #make osx=1 with mac
cd ../../..
ln -s tp_tools/tools/* .

 