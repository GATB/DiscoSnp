
# GET ALLL
wget http://bioinformatique.rennes.inria.fr/data_and_tools_tp_disco.zip

# INSTALL THE TOOLS
unzip data_and_tools_tp_disco.zip 
mv human_5M/* .
cd gassst
pwd
make #make osx=1 with mac
cd ..
chmod a+x validator.sh
 
