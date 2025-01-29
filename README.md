# Artifact for "Designing an 8-bit General-Purpose (T)FHE Processor Abstraction"
Dear referees,

The submitted artifact is the source code of the 8-bit TFHE processor abstraction proposed in the paper ``Designing an 8-bit General-Purpose (T)FHE Processor Abstraction''. 
It relies on version 1.1 of the TFHElib library (which runs over C++). Our artifact implements each instruction given in the paper (the complete list is in Appendix A of the paper) as well as the algorithms designed as a benchmark for the new approach (see Table 6 in Section 8).
All multivariate instructions are implemented in the instr.cpp file. The corresponding header instr.h defines each function of instr.cpp and give the corresponding section in the article. 
In the same fashion, we have respectively created the instri.cpp (resp. algos.cpp and instr_16.cpp) files to implement the univariate instructions (resp. the algorithms and 16-bit fixed point arithmetic operations) and their corresponding headers instri.h (resp. algos.h and instr_16.h) that link the functions to the correct sections of the paper.
The main is composed of an exhaustive set of tests that allows to execute each instruction and algorithm on inputs that are given by the user.

## Install TFHElib
To run the code, you first have to install the TFHElibrary (https://tfhe.github.io/tfhe/) and apply the patch provided in the file patch_fft.patch. The git clone command should automatically install version 1.1 of TFHElib, wich we used to implement our work.
1/ git clone https://github.com/tfhe/tfhe.git   
2/ cd tfhe
3/ cp path/to/patch_fft.patch ./
4/ git apply patch_fft.patch
5/ make


## Run the code
Once TFHElib version 1.1 and the provided patch are installed, you should:
1/ Go to the newly created directory tfhe_instructions_set.
2/ Update the path to the TFHElib directory in the provided sources/CMakeList.txt.
3/ In the tfhe_instructon_set directory, run the "cmake -S . -B ./build" command. 
4/ Run "cd build ; make" to finalize compilation.
5/ Command "../bin/tches n1 n2 n3 n4 n5" will then run an exhaustive test of all our instructions with timing measurements on inputs n1 and n2. 
Then, it will run an exhaustive test of all our algorithm implementations on the array composed of n1, n2, n3, n4 and n5. (The n_i must be byte values).
Finally, the programm will execute some 16-bit fixed-point arithmetic operations, including the execution of a neuron with sigmoid.
Note that you can enter as many byte as you want, as long as you give at least two of them.

## Instructions on how to interpret the outputs
The execution timings as well as the instruction/algorithm being executed and the corresponding results are printed on the terminal in real time. 
The user will easily recognize the instructions being tested as they are denoted with the notations given in the paper. For instance the addition of two 8-bit payload ciphertexts encrypting 45 and 103 will output:
ADD(45 + 103) = 148
timings : 0.49377s

This allows the user to verify if the homomorphic results are correct and to check the timings. 
For comparison, the execution timings obtained on our machine with the same code are given in Appendix A and Section 8 of the article.

Many thanks in advance for the time you will spend evaluating our work.
