# Get and decompress the data
rm -f data_test_disco.zip
wget http://bioinformatique.rennes.inria.fr/disco_tests/data_test_disco.zip
# wget http://gatb-discosnp.gforge.inria.fr/data_test_disco.zip 

# create: sh generate_refs.sh
unzip data_test_disco.zip

 
sh local_large_test.sh


rm -f created ref discoRes*
rm -f data_test_disco.zip
# rm -f humch1_*
# rm -f fof.txt
# rm -f ref*



