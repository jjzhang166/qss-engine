#!/bin/bash
#===================================================================================
#
# 				 FILE: simulate.sh
#
# 				USAGE: simulate.sh <FILE> <PARALLEL> <FOLDER>
#
# 	DESCRIPTION: Creates the outut directory for the simulation in the default 
# 							 location, copy the configuration files in the output directory,
# 							 run the simulation executable and finally remove the
# 							 configuration files. 
#
#    PARAMETERS: <FILE> simulation executable file.
#    						 <PARALLEL> True for parallel simulations, False otherwise.
#    						 <FOLDER> Location of the weights files for graph partition algorithms.
#       OPTIONS: --- 
#  REQUIREMENTS: MMOC_OUTPUT MMOC_BUILD must point to the corresponding default 
# 							 directories used by the QSS Solver GUI.
#         NOTES: ---
#        AUTHOR: Joaquin Fernandez, joaquin.f.fernandez@gmail.com
#       PROJECT: QSS Solver
#       VERSION: 3.1
#===================================================================================

FILE=$1

PARALLEL=$2

FOLDER=$3

mkdir -p $MMOC_OUTPUT/$FILE

cp $MMOC_BIN/khmetis $MMOC_OUTPUT/$FILE/

cd $MMOC_OUTPUT/$FILE

rm -rf *.log	

cp	$MMOC_BUILD/$FILE/$FILE.ini .

if [[	-e $MMOC_BUILD/$FILE/$FILE.part ]]; 
then
	cp $MMOC_BUILD/$FILE/$FILE.part .
fi

if [ "$PARALLEL" == "true" ]; then
	if [[ -e ${FILE}.ewgts ]];
	then
  	cp ${FILE}.ewgts ${FILE}.eweights
	else
		if [[ -e $MMOC_BUILD/$FILE/${FILE}.eweights ]];
		then
			cp $MMOC_BUILD/$FILE/${FILE}.eweights . 
		fi
	fi
	if [[ -e ${FILE}.hewgts ]];
	then
  	cp ${FILE}.hewgts ${FILE}.heweights
	else
		if [[ -e $MMOC_BUILD/$FILE/${FILE}.heweights ]];
		then
			cp $MMOC_BUILD/$FILE/${FILE}.heweights . 
		fi
	fi
	if [[ -e $MMOC_BUILD/$FILE/${FILE}.graph ]];
	then
		cp $MMOC_BUILD/$FILE/${FILE}.graph . 
	fi
	if [[ -e $MMOC_BUILD/$FILE/${FILE}.hgraph ]];
	then
		cp $MMOC_BUILD/$FILE/${FILE}.hgraph . 
	fi
fi

$MMOC_BUILD/$FILE/$FILE

rm -rf hkmetis
rm -rf *.part
rm -rf *.eweights	
rm -rf *.heweights	
rm -rf *.graph	
rm -rf *.hgraph	
rm -rf $FILE.ini  
