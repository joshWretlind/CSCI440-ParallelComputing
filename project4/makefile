#Compiler
#--------------------

mpiCC = mpiCC

#Compiler Switches
#--------------------

options = -std=c++0x

#Files
#--------------------
main = main.cpp



#Commands
#--------------------

all: main

main: $(main)
	$(mpiCC) $(main) -o main_exe $(options)

clean: rm main_exe 

