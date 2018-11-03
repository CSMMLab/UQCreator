# Dependencies
- Compiler with C++14 support e.g. g++ >= v4.9
- cmake >= v3.5
- LAPACK
- git
- vtk 
 
# Build instructions:
Note that an active internet connection is required for the first build in order to download the latest versions of the libraries/tools blaze, cpptoml, gnuplot-iostream and matplotlib-cpp!
If an if an active internet connection is present, run:

     git submodule init
     git submodule update

     cd code/build/release
     cmake -DCMAKE_BUILD_TYPE=Release ../../src
     make
 
The binary will afterwards be placed in the code/bin folder.
 
# Run instructions
Execute the compiled binary and hand over a valid TOML styled config file by using the parameter '-c'.
Example:

     ./bin/UQCreator -c input/example.toml

