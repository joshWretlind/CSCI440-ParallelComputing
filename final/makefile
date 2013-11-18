#Compiler
#--------------------

mpiCC = mpiCC

#Files
#--------------------
main = main.cpp
#Compiler Flags
#--------------------
flags = -std=c++0x

#Commands
#--------------------

all: main_exe

main_exe: $(main)
	$(mpiCC) $(flags) $(main) -o $@
		  
clean: rm main_exe
