#!/bin/bash
#===================================================================================
#
# 				 FILE: build.sh
#
# 				USAGE: build.sh 
#
# 	DESCRIPTION: Build the MAC bundle that contains the QSS Solver engine, 
# 							 the QSS Solver GUI, the MicroModelica C Compiler with the 
# 						   corresponding user libraries and the SBML translator tool.
#
#    PARAMETERS: ---
#       OPTIONS: ---
#  REQUIREMENTS: ---
#         NOTES: --- 
#         NOTES: ---
#        AUTHOR: Joaquin Fernandez, joaquin.f.fernandez@gmail.com
#       PROJECT: QSS Solver
#       VERSION: 3.1
#===================================================================================

cd ../../
echo "Retrieving latest from SVN";
svn update > rev
tail rev | awk '/At/ {print $NF}' | head -c 4 > rvn
REV=`cat rvn`
./deploy/linux/setRevision.sh ./deploy/mac/qss-solver.ini $REV 
echo "Done"
echo "Building QSS Solver MAC bundle.";
echo "Building Binaries";
rm -rf ./bin/qss-solver.app
cd ./src
make OS=mac  
cd ..
svn export deploy/mac/scripts ./bin/qss-solver.app/Contents/MacOS/scripts
chmod +x ./bin/qss-solver.app/Contents/MacOS/scripts/*
cp bin/mmoc  ./bin/qss-solver.app/Contents/MacOS/ 
cp bin/translate-sbml  ./bin/qss-solver.app/Contents/MacOS/ 
cp deploy/mac/qss-solver.ini ./bin/qss-solver.app/Contents/MacOS/qss-solver.ini
cp doc/COPYING  ./bin/qss-solver.app/Contents/MacOS/ 
cp doc/INSTALL  ./bin/qss-solver.app/Contents/MacOS/ 
cp doc/README.txt  ./bin/qss-solver.app/Contents/MacOS/ 
svn export build ./bin/qss-solver.app/Contents/Resources/build
svn export doc ./bin/qss-solver.app/Contents/Resources/doc
svn export models ./bin/qss-solver.app/Contents/Resources/models
svn export output ./bin/qss-solver.app/Contents/Resources/output
svn export usr ./bin/qss-solver.app/Contents/Resources/usr
svn export packages ./bin/qss-solver.app/Contents/Resources/packages
svn export src ./bin/qss-solver.app/Contents/Resources/src
cp src/libs/*.a ./bin/qss-solver.app/Contents/Resources/src/libs
cp usr/libs/*.a ./bin/qss-solver.app/Contents/Resources/usr/libs
cp /usr/local/lib/libsbml.5.10.2.dylib ./bin/qss-solver.app/Contents/Resources/usr/libs/libsbml.5.dylib
cd ./bin
rm -rf qss-solver.dmg
rm -rf qss-solver.app/Contents/Frameworks/*
rm -rf qss-solver.app/Contents/Resources/qt.conf
macdeployqt qss-solver.app -verbose=2 -dmg
cd ..
mv bin/qss-solver.dmg deploy/mac/
rm ./src/rvn ./src/rev
