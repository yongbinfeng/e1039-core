#!/bin/bash

#sudo mount -t cvmfs seaquest.opensciencegrid.org /cvmfs/seaquest.opensciencegrid.org
#sudo yum install sqlite-devel
#source  /cvmfs/seaquest.opensciencegrid.org/seaquest/software/current/setup.sh

#git clone https://github.com/cmantill/e1039-core.git
git clone https://github.com/wpmccormack/e1039-core.git
#git clone git@github.com:wpmccormack/e1039-core.git
cd e1039-core/
git checkout patrick_new_tracking_withEmbedding
cd script/

cp /home/submit/pcharris/DQ/setup-install.sh
cd ..
./script/setup-install.sh auto
source ../core-inst/this-e1039.sh
./build.sh

#git clone https://github.com/DarkQuest-FNAL/DarkQuest.git
git clone https://github.com/wpmccormack/DarkQuest.git
#git clone git@github.com:wpmccormack/DarkQuest.git
cd DarkQuest
#git checkout add_ana_flag_withEmbedding
cd e1039-analysis/SimHits/
sed "s@XXX@$PWD@g" /home/submit/pcharris/DQ/setup_mye1039_tmp.sh > setup_mye1039.sh
#sed "s@XXX@$PWD@g" /data/t3home000/pharris/DQ/DarkQuest/e1039-analysis/SimHits/setup_mye1039_tmp.sh > setup_mye1039.sh

source setup_mye1039.sh
mkdir work
mkdir install
cd work
cmake ../src/ -DCMAKE_INSTALL_PREFIX=../install
make clean
make
make install
cd ../../
pwd

cd HitEmbedding
sed "s@XXX@$PWD@g" /work/submit/wmccorma/DQ_Embed/setup_my1039_tmpHitEmbed.sh > setup_my1039.sh
source setup_my1039.sh
cmake-this
make-this
cd ..
