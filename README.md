# Dependecies
- Compiler with C++14 support e.g. g++ >= v4.9
- BLAS
- LAPACK
- boost
- cmake >= v3.5
- git
- gnuplot
 
# Build instructions:
Note that an active internet connection is required for the first build in order to download the latest versions of the header-only libs blaze, cpptoml and gnuplot-iostream!
If an active connection is present, enter the code/build folder and run:

     cmake ../src
     make
 
The binary will afterwards be placed in the code/bin folder.
 
# Run instructions
Execute the compiled binary and hand over a valid TOML styled config file by using the '-c' keyword.
Example:

     ./bin/UQCreator -c input/example.toml

