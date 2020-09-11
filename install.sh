#bin/bash
#
# Script to install third-party softwares
#------------------------
# Pre-requisites needed
#------------------------
# cmake :
# zlib :
# gcc version:
#
#
#

export DROGOM_HOME=`pwd`

echo $DROGOM_HOME
echo -e "\nInstalling SGA..."
cd $DROGOM_HOME/lib
tar xvzf sga.tar.gz

echo -e "\n\nInstalling SPAdes..."
cd $DROGOM_HOME/lib
if [ -d SPAdes.tar.gz ]; then
        rm -rf SPAdes.tar.gz
fi
tar xvzf SPAdes.tar.gz
[ $? -ne 0 ] && exit $?

echo -e "\n\nInstalling BWA..."
cd $DROGOM_HOME/lib
tar xvzf bwa.tar.gz

echo -e "\n\nInstalling CDHIT..."
cd $DROGOM_HOME/lib
tar xvzf cd-hit-est.tar.gz

echo -e "\n\nInstalling CMPress..."
cd $DROGOM_HOME/lib
tar xvzf cmpress.tar.gz

echo -e "\n\nInstalling CMSearch..."
cd $DROGOM_HOME/lib
tar xvzf cmsearch.tar.gz

echo -e "\n\nInstalling Samtools..."
cd $DROGOM_HOME/lib
tar xvzf samtools.tar.gz

echo -e "\nCompleted installing all third-party tools."
echo -e "\n----------------------------------------------\n"
## IMPP installation
echo -e "\nInstalling iMPP ..."
if [ ! -e $DROGOM_HOME/bin ]; then
        mkdir $DROGOM_HOME/bin
fi
cd $DROGOM_HOME/src
make clean
make

cd $DROGOM_HOME
echo -e "\n-----------------------------------------------\n"
echo -e "\nSuccessfully installed DRAGoM.\n"
echo -e "\n-----------------------------------------------\n"
