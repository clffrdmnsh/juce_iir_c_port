# juce_iir_c_port
A simple port of all the useful function from the IIR Filter class in juce to c

1) Thou shall compile and build an executable if you got gcc on you command line. 

gcc -0 name_for_the_Lib main.c IIRFilter.c

if this dosent work. get gcc. 

2) Thou shall create a project file incluing this API using cmake

mkdir build
cd build
cmake .. -G"Unix Makefiles"or -G"XCode" or -G"Visual Studio 16 2022" -A x64

3) Simply add the .h and .c files to a pre existing project and import functions from it. 

bless
