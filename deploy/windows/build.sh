#!/bin/bash
#===================================================================================
#
# 				 FILE: build.sh  
#
# 				USAGE: build.sh 
#
# 	DESCRIPTION: Build the Windows installer including the QSS Solver engine, 
# 							 the QSS Solver GUI, the MicroModelica C Compiler with the 
# 						   corresponding user libraries and the SBML translator tool.
#
#    PARAMETERS: ---
#       OPTIONS: ---
#  REQUIREMENTS: ---
#         NOTES: --- 
#        AUTHOR: Joaquin Fernandez, joaquin.f.fernandez@gmail.com
#       PROJECT: QSS Solver
#       VERSION: 3.1
#===================================================================================

cd ../../
echo Building Windows Package
echo Retrieving latest from SVN
svn update 
svnversion | head -c 4 > rvn
REV=`cat rvn`
./deploy/linux/setRevision.sh ./deploy/windows/qss-solver.ini $REV 
echo Done.
echo Creating directories
rm -rf tmp-win-installer 
echo Done.
echo Exporting files
svn export ./deploy/windows/export ./tmp-win-installer
svn export ./models ./tmp-win-installer/qss-solver/models
svn export ./build ./tmp-win-installer/qss-solver/build
svn export ./src ./tmp-win-installer/qss-solver/src
svn export ./output ./tmp-win-installer/qss-solver/output
svn export ./usr ./tmp-win-installer/qss-solver/usr
svn export ./packages ./tmp-win-installer/qss-solver/packages
svn export ./doc ./tmp-win-installer/qss-solver/doc
cp doc/COPYING ./tmp-win-installer/qss-solver
cp doc/INSTALL ./tmp-win-installer/qss-solver
cp doc/README.txt ./tmp-win-installer/qss-solver
cp deploy/images/integrator.svg ./tmp-win-installer/qss-solver/bin
cp deploy/windows/qss-solver.ini ./tmp-win-installer/qss-solver/bin
cp deploy/windows/qss-solver.gpr ./tmp-win-installer
echo Done.
echo Building binaries
cd src
make OS=win32 
cd ..
echo Done.
echo Copying executables files 
cd bin
cp qss-solver.exe ../tmp-win-installer/qss-solver/bin/
cp mmoc.exe ../tmp-win-installer/qss-solver/bin/
cp translate-sbml.exe ../tmp-win-installer/qss-solver/bin/
echo Done.
echo Copying libraries
cd ../src/libs
cp *.a ../../tmp-win-installer/qss-solver/src/libs/
cp ../../usr/libs/*.a ../../tmp-win-installer/qss-solver/usr/libs/
cd ../../tmp-win-installer
mv ./qss-solver/bin/*.a ./qss-solver/src/libs/
echo Done.
echo Creating installer
gibuild qss-solver.gpr
cp ./Output/qss-solver-install.exe ../
cd ..
rm rvn
rm -rf tmp-win-installer 
mv qss-solver-install.exe deploy/windows/
echo Done.
