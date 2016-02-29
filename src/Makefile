OS ?= unix
DEBUG ?= False
CORES ?= 1

all: qss-compiler qss-engine qss-gui qss-tools qss-user

qss-tools:
	@echo BUILDING SBML TRANSLATOR
	@cd tools/sbml &&  $(MAKE) OS=$(OS) DEBUG=$(DEBUG) -j $(CORES)
	@echo Done

qss-engine:
	@echo BUILDING QSS LIBRARY
	@cd libs && $(MAKE) OS=$(OS) -j $(CORES)
	@echo Done
	@echo BUILDING QSS  DEBUG LIBRARY
	@cd libs && $(MAKE) DEBUG=True OS=$(OS) -j $(CORES)
	@echo Done

qss-user:
	@echo BUILDING USER LIBS
	@cd ../usr && $(MAKE) OS=$(OS) DEBUG=$(DEBUG) -j $(CORES)      
	@echo Done

qss-compiler:
	@echo BUILDING MICROMODELICA C COMPILER
	@cd mmo && $(MAKE) OS=$(OS) DEBUG=$(DEBUG) -j $(CORES) 
	@echo Done

qss-gui:
	@echo BUILDING MMO EDITOR
ifeq ($(OS),unix)
	@cd gui && qmake-qt4
else
	@cd gui && qmake
endif
	@cd gui && $(MAKE)      
	@echo Done

.PHONY: clean all 

clean: 
	cd mmo && $(MAKE) OS=$(OS) clean 
	cd libs && $(MAKE) clean 
	cd libs && $(MAKE) DEBUG=True clean 
	cd ../usr && $(MAKE) clean
	cd gui && $(MAKE) clean    
	cd tools/sbml && $(MAKE) clean     

help:
	echo "make DEBUG=<True|False> OS=<unix|win32|mac> CORES=<n>"
	echo "Default values:"
	echo ""
	echo "DEBUG=False"
	echo "OS=unix"
	echo "CORES=1"
