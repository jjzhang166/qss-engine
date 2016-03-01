#!/bin/bash
#===================================================================================
#
# 				 FILE: doc.sh
#
# 				USAGE: doc.sh 
#
# 	DESCRIPTION: Generate the framework documentation and upload the files to the
#                CIFASIS server. 
#
#    PARAMETERS: ---  
#  REQUIREMENTS: ---
#         NOTES: ---
#        AUTHOR: Joaquin Fernandez, joaquin.f.fernandez@gmail.com
#       PROJECT: QSS Solver
#       VERSION: 3.1
#===================================================================================



cd $MMOC_SRC/engine

doxygen QSSSolver.doxyfile

cd $MMOC_SRC/mmo

doxygen QSSSolver.doxyfile

cd $MMOC_SRC/sbml

doxygen sbml.doxyfile

cd $MMOC_SRC/gui

doxygen QSSSolverGUI.doxyfile

cd $MMOC_DOC/engine/html

SSH_AUTH_SOCK=0 scp -r * jfernandez@www.cifasis-conicet.gov.ar:/home/jfernandez/public_html/engine/

cd $MMOC_DOC/mmoc/html

SSH_AUTH_SOCK=0 scp -r * jfernandez@www.cifasis-conicet.gov.ar:/home/jfernandez/public_html/mmoc/

cd $MMOC_DOC/gui/html

SSH_AUTH_SOCK=0 scp -r * jfernandez@www.cifasis-conicet.gov.ar:/home/jfernandez/public_html/gui/

cd $MMOC_DOC/sbml/html

SSH_AUTH_SOCK=0 scp -r * jfernandez@www.cifasis-conicet.gov.ar:/home/jfernandez/public_html/sbml/
