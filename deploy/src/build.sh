#!/bin/bash
#===================================================================================
#
# 				 FILE: build.sh  
#
# 				USAGE: build.sh 
#
# 	DESCRIPTION: Build the tar.gz file that contains source files and libraries.
#
#    PARAMETERS: ---
#       OPTIONS: ---
#  REQUIREMENTS: ---
#         NOTES: --- 
#        AUTHOR: Joaquin Fernandez, joaquin.f.fernandez@gmail.com
#       PROJECT: QSS Solver
#       VERSION: 3.0
#===================================================================================

rm -rf qss-solver-3.1.tar.gz
cd ../../
mkdir qss-solver-3.1
svn export ./bin qss-solver-3.1/bin
svn export ./build qss-solver-3.1/build
svn export ./deploy qss-solver-3.1/deploy
svn export ./doc qss-solver-3.1/doc
svn export ./models qss-solver-3.1/models
svn export ./output qss-solver-3.1/output
svn export ./packages qss-solver-3.1/packages
svn export ./src qss-solver-3.1/src
svn export ./testsuite qss-solver-3.1/testsuite
svn export ./usr qss-solver-3.1/usr
tar -cvzf qss-solver-3.1.tar.gz qss-solver-3.1
mv qss-solver-3.1.tar.gz ./deploy/src/
rm -rf qss-solver-3.1
