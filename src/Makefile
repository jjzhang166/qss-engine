#Target OS variables
OS ?= unix
DEBUG ?= False
CORES ?= 1

#The Directories, Source, Includes, Objects, Binary 
ROOT   := ..

all: qss-compiler qss-engine qss-gui qss-tools qss-user

qss-tools: create-dirs
	@echo BUILDING SBML TRANSLATOR
	@cd tools/sbml &&  $(MAKE) OS=$(OS) DEBUG=$(DEBUG) -j $(CORES)
	@echo Done

qss-engine: create-dirs
	@echo BUILDING QSS LIBRARY
	@cd libs && $(MAKE) OS=$(OS) -j $(CORES)
	@echo Done
	@echo BUILDING QSS  DEBUG LIBRARY
	@cd libs && $(MAKE) DEBUG=True OS=$(OS) -j $(CORES)
	@echo Done

qss-user: create-dirs
	@echo BUILDING USER LIBS
	@cd ../usr && $(MAKE) OS=$(OS) DEBUG=$(DEBUG) -j $(CORES)      
	@echo Done

qss-compiler: create-dirs
	@echo BUILDING MICROMODELICA C COMPILER
	@cd mmo && $(MAKE) OS=$(OS) DEBUG=$(DEBUG) -j $(CORES) 
	@echo Done

qss-gui: create-dirs
	@echo BUILDING MMO EDITOR
ifeq ($(OS),unix)
	@cd gui && qmake-qt4
else
	@cd gui && qmake
endif
	@cd gui && $(MAKE)      
	@echo Done

qss-doc:
	@echo BUILDING PROJECT DOCUMENTATION
	@mkdir -p $(ROOT)/doc/engine
	@mkdir -p $(ROOT)/doc/engine/html
	@mkdir -p $(ROOT)/doc/engine/latex
	@cd engine && doxygen QSSSolver.doxyfile 
	@mkdir -p $(ROOT)/doc/mmoc
	@mkdir -p $(ROOT)/doc/mmoc/html
	@mkdir -p $(ROOT)/doc/mmoc/latex
	@cd mmo && doxygen MMOCompiler.doxyfile 
	@mkdir -p $(ROOT)/doc/gui
	@mkdir -p $(ROOT)/doc/gui/html
	@mkdir -p $(ROOT)/doc/gui/latex
	@cd gui && doxygen QSSSolverGUI.doxyfile 
	@mkdir -p $(ROOT)/doc/sbml
	@mkdir -p $(ROOT)/doc/sbml/html
	@mkdir -p $(ROOT)/doc/sbml/latex
	@cd tools/sbml && doxygen sbml.doxyfile
	
create-dirs:
	@mkdir -p $(ROOT)/bin
	@mkdir -p $(ROOT)/build
	@mkdir -p $(ROOT)/output
	@mkdir -p $(ROOT)/usr/libs	
	
.PHONY: clean all 

clean: 
	@cd mmo && $(MAKE) OS=$(OS) clean 
	@cd libs && $(MAKE) clean 
	@cd libs && $(MAKE) DEBUG=True clean 
	@cd ../usr && $(MAKE) clean
	@cd gui && $(MAKE) clean    
	@cd tools/sbml && $(MAKE) clean     

help:
	@echo "make DEBUG=<True|False> OS=<unix|win32|mac> CORES=<n>"
	@echo "Default values:"
	@echo ""
	@echo "DEBUG=False"
	@echo "OS=unix"
	@echo "CORES=1"
