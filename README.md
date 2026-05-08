# Pirouette
Implementation of our lightweight PIR method in ia.cr/2025/680 by Jiayi Kang and Leonard Schild, to appear in PoPETs 2026, Issue 2.

# Build instructions
- Run ```install_dependencies.sh``` to install openfhe and intel-hexl
- ```mkdir build```
- ```cd build```
- cmake ..
- make -j 16
- run the gtest suite as needed
- Everything is nicer if one uses ```CLion``` 