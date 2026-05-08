#!/usr/bin/env bash

install_openfhe() {
  export CC=/usr/bin/clang;
  export CXX=/usr/bin/clang++;
  mkdir -p ./tmp/openfhe/;
  git clone https://github.com/openfheorg/openfhe-development ./tmp/openfhe/
  mkdir -p ./tmp/openfhe/build
  cmake ./tmp/openfhe -B./tmp/openfhe/build -DCMAKE_INSTALL_PREFIX="$1" -DWITH_NATIVEOPT=ON -DWITH_OPENMP=OFF
  cmake --build ./tmp/openfhe/build -j 8
  cmake --install ./tmp/openfhe/build
}

install_hexl() {
  export CC=/usr/bin/clang;
  export CXX=/usr/bin/clang++;
  mkdir -p ./tmp/hexl;
  git clone https://github.com/intel/hexl ./tmp/hexl;
  mkdir -p ./tmp/hexl/build;
  cmake ./tmp/hexl -B./tmp/hexl/build -DCMAKE_INSTALL_PREFIX="$1" -DHEXL_TESTING=OFF;
  cmake --build ./tmp/hexl/build -j 8;
  cmake --install ./tmp/hexl/build;
}

install_openfhe "$PWD"/dependencies/openfhe
install_hexl "$PWD"/dependencies/hexl
