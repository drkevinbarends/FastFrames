asetup StatAnalysis,0.5.3
mkdir -p build install
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install/ ../FastFrames/
make -j4
make install
source setup.sh
cd ../tWZClass
mkdir -p build install
cd build
cmake -DCMAKE_PREFIX_PATH=/eos/user/k/kebarend/tWZ/FastFrames/install -DCMAKE_INSTALL_PREFIX=../install ../
make -j4
make install
source setup.sh
cd ../../