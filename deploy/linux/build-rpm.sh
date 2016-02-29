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

cd ../../
ARCH=`uname -m`
echo "Retrieving latest from SVN";
#svn up
svnversion > rev
head version.major -c 6 >vm
cat vm rev > version
VER=`cat version`
echo "Building QSS Solver DEB package for $ARCH";
echo "Building Binaries";
cd ./src
#make clean
#make SBML=True
cd ..
RPMDIR=~/rpmbuild/BUILDROOT/qss-solver-3.0-1.el7.centos.x86_64 
rm -rf ${RPMDIR} 
mkdir ${RPMDIR} 
mkdir ${RPMDIR}/opt
mkdir ${RPMDIR}/opt/qss-solver
svn export deploy/linux/scripts ${RPMDIR}/opt/qss-solver/bin
cp bin/mmoc  ${RPMDIR}/opt/qss-solver/bin/ 
cp bin/qss-solver ${RPMDIR}/opt/qss-solver/bin/ 
cp bin/translate-sbml  ${RPMDIR}/opt/qss-solver/bin/
cp deploy/linux/qss-solver.ini ${RPMDIR}/opt/qss-solver/bin/qss-solver.ini
chmod 0755 `find ${RPMDIR}/opt/qss-solver/bin`
cp doc/COPYING ${RPMDIR}/opt/qss-solver/
cp doc/INSTALL ${RPMDIR}/opt/qss-solver/
cp doc/README.txt ${RPMDIR}/opt/qss-solver/
cp version ${RPMDIR}/opt/qss-solver/
svn export build ${RPMDIR}/opt/qss-solver/build
svn export doc ${RPMDIR}/opt/qss-solver/doc
svn export models ${RPMDIR}/opt/qss-solver/models
svn export output ${RPMDIR}/opt/qss-solver/output
svn export usr ${RPMDIR}/opt/qss-solver/usr
svn export packages ${RPMDIR}/opt/qss-solver/packages
svn export src ${RPMDIR}/opt/qss-solver/src
cp src/libs/*.a ${RPMDIR}/opt/qss-solver/src/libs
cp usr/libs/*.a ${RPMDIR}/opt/qss-solver/usr/libs
cp /usr/lib64/libsbml.so.5.11.2 ${RPMDIR}/opt/qss-solver/src/libs/libsbml.so.5
chmod 0644 `find ${RPMDIR}/ -iname *.cpp`
chmod 0644 `find ${RPMDIR}/ -iname *.c`
chmod 0644 `find ${RPMDIR}/ -iname *.hpp`
chmod 0644 `find ${RPMDIR}/ -iname *.h`
chmod 0644 `find ${RPMDIR}/ -iname *.ini`
chmod 0644 `find ${RPMDIR}/ -iname *.png`
chmod 0644 `find ${RPMDIR}/ -iname *.tex`
chmod 0644 `find ${RPMDIR}/opt/qss-solver/build/ -type f`
chmod 0644 `find ${RPMDIR}/opt/qss-solver/doc/ -type f`
chmod 0644 `find ${RPMDIR}/opt/qss-solver/src/ -type f`
chmod 0644 `find ${RPMDIR}/opt/qss-solver/packages/ -type f`
chmod 0644 `find ${RPMDIR}/opt/qss-solver/usr/ -type f`
chmod 0755 `find ${RPMDIR}/ -type d`
rpmbuild -bb ~/rpmbuild/SPECS/qss-solver.spec
mv ~/rpmbuild/RPMS/x86_64/qss-solver-3.0-1.el7.centos.x86_64.rpm ./deploy/linux/qss-solver-3.0.x86_64.rpm
rm rev version vm
rm -rf ${RPMDIR}
cd deploy/linux
