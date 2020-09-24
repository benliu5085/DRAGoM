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
rm -rf "env.config"
echo "home " $DROGOM_HOME >> env.config

echo -e "\nInstalling SGA..."
cd $DROGOM_HOME/lib
tar xvzf sga.tar.gz
echo "sga " $DROGOM_HOME"/lib/sga" >> $DROGOM_HOME/env.config

echo -e "\n\nInstalling SPAdes..."
cd $DROGOM_HOME/lib
if [ -d SPAdes.tar.gz ]; then
        rm -rf SPAdes.tar.gz
fi
tar xvzf SPAdes.tar.gz
[ $? -ne 0 ] && exit $?
echo "spades " $DROGOM_HOME"/lib/SPAdes/bin/spades.py"  >> $DROGOM_HOME/env.config

echo -e "\n\nInstalling BWA..."
cd $DROGOM_HOME/lib
tar xvzf bwa.tar.gz
echo "bwa " $DROGOM_HOME"/lib/bwa"  >> $DROGOM_HOME/env.config

echo -e "\n\nInstalling CDHIT..."
cd $DROGOM_HOME/lib
tar xvzf cd-hit-est.tar.gz
echo "cdhit " $DROGOM_HOME"/lib/cd-hit-est"  >> $DROGOM_HOME/env.config

echo -e "\n\nInstalling CMPress..."
cd $DROGOM_HOME/lib
tar xvzf cmpress.tar.gz
echo "cmpress " $DROGOM_HOME"/lib/cmpress"  >> $DROGOM_HOME/env.config

echo -e "\n\nInstalling CMSearch..."
cd $DROGOM_HOME/lib
tar xvzf cmsearch.tar.gz
echo "cmsearch " $DROGOM_HOME"/lib/cmsearch"  >> $DROGOM_HOME/env.config

echo -e "\n\nInstalling Samtools..."
cd $DROGOM_HOME/lib
if test `cat /etc/os-release | grep -c 'Ubuntu'` -gt 0 ; then
  tar xvzf samtools_ubuntu.gz.tar
elif test `cat /etc/os-release | grep -c 'Red Hat'` -gt 0; then
  tar xvzf samtools_redhat.tar.gz
fi
echo "samtools " $DROGOM_HOME"/lib/samtools"  >> $DROGOM_HOME/env.config

echo -e "\nCompleted installing all third-party tools."
echo -e "\n----------------------------------------------\n"
## DRAGoM installation
echo -e "\nInstalling DRAGoM ..."
if [ ! -e $DROGOM_HOME/bin ]; then
        mkdir $DROGOM_HOME/bin
fi
cd $DROGOM_HOME/src
make clean
make
echo "dragom " $DROGOM_HOME"/bin/dragom.exe"  >> $DROGOM_HOME/env.config

cd $DROGOM_HOME
echo -e "\n-----------------------------------------------\n"
echo -e "\nSuccessfully installed DRAGoM.\n"
echo -e "\n-----------------------------------------------\n"
