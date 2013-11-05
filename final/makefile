#Compiler
#--------------------

mpiCC = mpiCC

#Files
#--------------------
main = main.cpp

#Commands
#--------------------

all: main_exe

main_exe: $(main)
	$(mpiCC) $(main) -o $@
		  
clean: rm main_exe
