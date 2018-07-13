#!/bin/bash
      
scramdir=$1
outputDir=$2
inputDir=$3
inputFile=$4
xsec=$5
id=$6
runMacro=$7
soFile=$8
pcmFile=$9
vomsProxy=${10}


workDir=`pwd`
echo `hostname`
echo "args:    $*"

cd ${scramdir}
eval `scramv1 runtime -sh`
cd $workDir

export X509_USER_PROXY=${scramdir}/src/${vomsProxy}
echo ${X509_USER_PROXY}

xrdfs dcache-cms-xrootd.desy.de ls /store/group/upgrade/delphes_output/

#cp ${scramdir}/setRootEnv.C . 
cp ${scramdir}/src/rootlogon.C .
cp ${scramdir}/src/$runMacro  .
cp ${scramdir}/src/$soFile  .
cp ${scramdir}/src/$pcmFile  .

echo root -l -b -q ${runMacro}+\(\"${inputDir}/${inputFile}\",${xsec},${id},\"${inputFile}\"\)
root -l -b -q ${runMacro}+\(\"${inputDir}/${inputFile}\",${xsec},${id},\"${inputFile}\"\)

echo cp ${inputFile} ${outputDir}/
cp ${inputFile} ${outputDir}/

status=`echo $?`
echo "Status - $status"

exit $status
