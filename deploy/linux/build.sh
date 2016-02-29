#!/bin/bash
#===================================================================================
#
# 				 FILE: build.sh  
#
# 				USAGE: build.sh 
#
# 	DESCRIPTION: Build the Linux deb package including the QSS Solver engine, 
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

rm -rf qss-solver-i386.deb
rm -rf qss-solver-amd64.deb

cd ../../
ARCH=`uname -m`
echo "Retrieving latest from SVN";
svn update
head ./doc/version.major -c 4 > ./src/vm
cd ./src
svn update > rev
tail rev | awk '/At/ {print $NF}' | head -c 4 > rvn
cat vm rvn > version
REV=`cat rvn`
../deploy/linux/setRevision.sh ../deploy/linux/qss-solver.ini $REV
VER=`cat version`
echo "Building QSS Solver DEB package for $ARCH version $VER";
echo "Building Binaries";
make clean
make 
cd ..
rm -rf tmp_deb
svn export deploy/linux/deb tmp_deb
chmod 0755 tmp_deb/DEBIAN/post*
mkdir ./tmp_deb/opt/qss-solver
svn export deploy/linux/scripts ./tmp_deb/opt/qss-solver/bin
if [ "$ARCH" == "i686" ]; then
  cat ./tmp_deb/DEBIAN/control.i386 | awk -v VERSION="$VER" '{ if(index($0,"Version: ")>=1) print "Version: " VERSION ; else print $0;}' >  ./tmp_deb/DEBIAN/control
  rm ./tmp_deb/DEBIAN/control.amd64; 
  rm ./tmp_deb/DEBIAN/control.i386; 
fi
if [ "$ARCH" == "x86_64" ]; then
  cat ./tmp_deb/DEBIAN/control.amd64 | awk -v VERSION="$VER" '{ if(index($0,"Version: ")>=1) print "Version: " VERSION ; else print $0;}' >  ./tmp_deb/DEBIAN/control
  rm ./tmp_deb/DEBIAN/control.amd64; 
  rm ./tmp_deb/DEBIAN/control.i386; 
fi
cp bin/mmoc  ./tmp_deb/opt/qss-solver/bin/ 
cp bin/qss-solver ./tmp_deb/opt/qss-solver/bin/ 
cp bin/translate-sbml  ./tmp_deb/opt/qss-solver/bin/
cp src/tools/partitioners/hmetis/khmetis ./tmp_deb/opt/qss-solver/bin/
cp deploy/linux/qss-solver.ini ./tmp_deb/opt/qss-solver/bin/qss-solver.ini
cp deploy/images/integrator.svg ./tmp_deb/opt/qss-solver/bin/
chmod 0755 `find tmp_deb/opt/qss-solver/bin`
cp doc/COPYING ./tmp_deb/opt/qss-solver/
cp doc/INSTALL ./tmp_deb/opt/qss-solver/
cp doc/README.txt ./tmp_deb/opt/qss-solver/
cp version ./tmp_deb/opt/qss-solver/
svn export build ./tmp_deb/opt/qss-solver/build
svn export doc ./tmp_deb/opt/qss-solver/doc
svn export models ./tmp_deb/opt/qss-solver/models
svn export output ./tmp_deb/opt/qss-solver/output
svn export usr ./tmp_deb/opt/qss-solver/usr
svn export packages ./tmp_deb/opt/qss-solver/packages
svn export src ./tmp_deb/opt/qss-solver/src
cp src/libs/*.a ./tmp_deb/opt/qss-solver/src/libs
cp usr/libs/*.a ./tmp_deb/opt/qss-solver/usr/libs
if [ "$ARCH" == "i686" ]; then
	cp /usr/lib/libsbml.so.5.6.0 ./tmp_deb/opt/qss-solver/src/libs/libsbml.so.5
	cp src/tools/partitioners/patoh/Linux-i386/libpatoh.a ./tmp_deb/opt/qss-solver/src/libs/libpatoh.a
fi
if [ "$ARCH" == "x86_64" ]; then
	cp /usr/lib64/libsbml.so.5.10.2 ./tmp_deb/opt/qss-solver/src/libs/libsbml.so.5
	cp src/tools/partitioners/patoh/Linux-x86_64/libpatoh.a ./tmp_deb/opt/qss-solver/src/libs/libpatoh.a
fi
chmod 0644 `find tmp_deb/ -iname *.cpp`
chmod 0644 `find tmp_deb/ -iname *.c`
chmod 0644 `find tmp_deb/ -iname *.hpp`
chmod 0644 `find tmp_deb/ -iname *.h`
chmod 0644 `find tmp_deb/ -iname *.ini`
chmod 0644 `find tmp_deb/ -iname *.png`
chmod 0644 `find tmp_deb/ -iname *.tex`
chmod 0644 `find tmp_deb/opt/qss-solver/build/ -type f`
chmod 0644 `find tmp_deb/opt/qss-solver/doc/ -type f`
chmod 0644 `find tmp_deb/opt/qss-solver/src/ -type f`
chmod 0644 `find tmp_deb/opt/qss-solver/packages/ -type f`
chmod 0644 `find tmp_deb/opt/qss-solver/usr/ -type f`
chmod 0755 `find tmp_deb/ -type d`
fakeroot dpkg -b tmp_deb qss-solver.deb
if [ "$ARCH" == "i686" ]; then
  mv qss-solver.deb ./deploy/linux/qss-solver-i386.deb
fi
if [ "$ARCH" == "x86_64" ]; then
  mv qss-solver.deb ./deploy/linux/qss-solver-amd64.deb
fi
rm ./src/rvn ./src/rev ./src/version ./src/vm
rm -rf tmp_deb
cd deploy/linux
