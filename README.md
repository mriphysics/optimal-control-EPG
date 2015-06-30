##Optimal Control EPG

This Matlab code may be used to reproduce results from an upcoming paper (link to be posted soon) on optimal control design of MRI pulse sequences. The work builds on prior work (see [this](http://dx.doi.org/10.1002/mrm.25192) paper, and [this](http://dx.doi.org/10.1002/mrm.25192) one) but treats the design problem using optimal control methods, rather than simple numerical optimization. The result is highly accelerated and much stabler calculation.

**Releases:**

The [releases](https://github.com/mriphysics/optimal-control-EPG/releases) page contains downloadable archives of the code along with test data (B1 field maps in .mat format). The Windows version also contains .exe compiled binaries.

###Contents
  1. [Introduction](#intro)
  2. [Compile C++ code (optional)](#compile)
  3. [Choose the input parameters](#input)
  4. [Run the main script](#run)

<a name="intro"></a>

### 1. Introduction

This code runs on Matlab with the optimization toolbox (necessary for the function fmincon.m). The EPG simulator is implemented in C++ (Linux and Mac users) and in Matlab (all operative systems). The user can choose which version. The C++ implementation is faster but it requires the open-source linear algebra library ARMADILLO version 5.100 or higher. For more info and download see: http://arma.sourceforge.net
For information about the compilation of the c++ code see section 2.

<a name="compile"></a>

###2. COMPILE C++ code (optional)

The EPG simulator is implemented in both Matlab as well as in c++. The user can choose which version by indicating it in the info file (see below). To improve performance use the c++ version: obj_EPG13_cpp.
To compile the C++ file, download ARMADILLO and save in a directory of your choice (indicated here as ARMADILLO). Then type in matlab:

```system('g++ obj_EPG13.cpp -o obj_EPG13 -O2 -I DIRECTORY/armadillo-5.100.1/include -lblas -llapack');```

If you have troubles compiling or running obj_EPG13_cpp.m you can choose for the analogous matlab version obj_EPG13.m (slower). Just replace "obj_EPG13_cpp.m" by "obj_EPG13.m" in the "OPTIMIZATION" cell of the MAIN.m file

NOTE: with some extra effort it is possible to compile obj_EPG13.cpp on Windows and modify the line

`system('./obj_EPG13');`

in "obj_EPG13_cpp.m" to

`system('obj_EPG13.exe');`

The C++ EPG simulator should then run also on Windows machines.

<a name="input"></a>

###3. Choose the input parameters

The optimal control EPG design requires input parameters. Their are given in the file named info_NAME.m where they are listed. The B1 maps have to be presented as .mat file.
Two examples are already given. You can open them to have more detailed instructions.

<a name="run"></a>

###4. The main function is MAIN.m .

After commenting the relevant sections, run the script all the way through.
