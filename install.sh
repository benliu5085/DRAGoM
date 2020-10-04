#bin/bash
#
# Script to install third-party softwares
#------------------------
# Pre-requisites needed
#------------------------
# cmake : 3.12
# make : 3.82
# gcc version: 4.8.5
#
#
#

export DROGOM_HOME=`pwd`
rm -rf "env.config"
echo "home " $DROGOM_HOME >> env.config

echo -e "Installing SGA..."
cd $DROGOM_HOME/lib
tar xvzf sga.tar.gz
echo "sga " $DROGOM_HOME"/lib/sga" >> $DROGOM_HOME/env.config

echo -e "\nInstalling SPAdes..."
cd $DROGOM_HOME/lib
tar xvzf SPAdes.tar.gz
echo "spades " $DROGOM_HOME"/lib/SPAdes/bin/spades.py"  >> $DROGOM_HOME/env.config

echo -e "\nInstalling BWA..."
cd $DROGOM_HOME/lib
tar xvzf bwa.tar.gz
echo "bwa " $DROGOM_HOME"/lib/bwa"  >> $DROGOM_HOME/env.config

echo -e "\nInstalling CDHIT..."
cd $DROGOM_HOME/lib
tar xvzf cd-hit-est.tar.gz
echo "cdhit " $DROGOM_HOME"/lib/cd-hit-est"  >> $DROGOM_HOME/env.config

echo -e "\nInstalling CMPress..."
cd $DROGOM_HOME/lib
tar xvzf cmpress.tar.gz
echo "cmpress " $DROGOM_HOME"/lib/cmpress"  >> $DROGOM_HOME/env.config

echo -e "\nInstalling CMSearch..."
cd $DROGOM_HOME/lib
tar xvzf cmsearch.tar.gz
echo "cmsearch " $DROGOM_HOME"/lib/cmsearch"  >> $DROGOM_HOME/env.config

echo -e "\nInstalling Samtools..."
cd $DROGOM_HOME/lib
if test `cat /etc/os-release | grep -c 'Ubuntu'` -gt 0 ; then
  tar xvzf samtools_ubuntu.gz.tar
elif test `cat /etc/os-release | grep -c 'Red Hat'` -gt 0; then
  tar xvzf samtools_redhat.tar.gz
fi
echo "samtools " $DROGOM_HOME"/lib/samtools"  >> $DROGOM_HOME/env.config

echo -e "\nCompleted installing all third-party tools."
echo -e "----------------------------------------------"
## DRAGoM installation
echo -e "Installing DRAGoM ..."
mkdir -p $DROGOM_HOME/bin
cd $DROGOM_HOME/src
cmake .

if test -f "Makefile"; then
	make
	mv dragom.exe ../bin
	echo "dragom " $DROGOM_HOME"/bin/dragom.exe"  >> $DROGOM_HOME/env.config

	cd $DROGOM_HOME
	echo -e "\n-----------------------------------------------\n"
	echo -e "Successfully installed DRAGoM."
	echo -e "\n-----------------------------------------------\n"
else
	echo -e "\n-----------------------------------------------\n"
	echo -e "Failed to installed DRAGoM."
	echo -e "Check the error massage above"
	echo -e "\n-----------------------------------------------\n"
fi
