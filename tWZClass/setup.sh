cd build
cmake -DCMAKE_PREFIX_PATH=/eos/user/k/kebarend/tWZ/FastFrames/install -DCMAKE_INSTALL_PREFIX=../install ../
make -j4
make install
source setup.sh
cd ../../