# README #

**Quick compile:**

       mkdir bin
       cd bin
       cmake ..
       make
       cp -r ../data data
Then you will find and run the executable file named “\*D.ex”.

**Structure:**

* main.cpp: main file, which gives the initial and boundary conditions 
* include/: ufl and main header file
* src/: source file
* data/: 
*       mesh/: mesh files
*       parameter_fsi.xml: FSI parameter file
* unittest/: test codes for different functions
