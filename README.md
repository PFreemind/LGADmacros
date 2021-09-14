# LGADmacros
Macros for LGAD measurements 
Partick Freeman
pmfreeman@ucsb.edu


Script to convert from .dat binary files from DRS4, modified version of code at http://er.jinr.ru/git/vratislav.chudoba/NeuRad_tests/commit/7b3acf5f30be57cc5d5c15e1ade81df559a44b8c

build command with:
g++ -o read_binary read_binary.cpp -lm `root-config --cflags --libs`

run with: 
./read_binary ../Sr90KU4_run2 ../test.root
// form: ./read_binary input_file.dat outpur_file.root
//reads .dat file from DRS4 output
//modified verision has a TBranch for board serial number for use in analysis code
//allows for arbitrary number of boards (assumed to have 4 channels read out each)


Waveform analysis macro: analysisDRSN.C
in a ROOT session

root -l 
.L analysisDRSN.C
 analyze("../data/Sr90/processed/run0000.root",1,2,3,4,2)
 //that is, command of the form: analyze("input file", 1,2,3,4,nBoards)
 //assumes same channels on each board
// 1,2,3,4, indicates which channels read out, code needs to be modified to use less than all 4

analyzed root files have vectors of various pulse measurements



