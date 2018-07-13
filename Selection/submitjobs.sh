#!/bin/bash

#master submission script
conf_dir=/afs/cern.ch/work/j/jlawhorn/public/forXLL/delphes-intro/Selection/config/
runMacro=selection.C
outputBase=/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-test/

echo "Checking for new samples and"
echo "remaking config files."

while read line #loop over lines in ${conf_file}                                                                                                         
do
    array=($line)
    if [[ "${array[0]:0:1}" != "#" ]]; then
	echo "hi"
        if [[ ${array[3]} -eq "1" ]]; then
            testvar="$( xrdfs dcache-cms-xrootd.desy.de ls ${array[5]} 2>&1 > /dev/null )"
            if [[ ${testvar} != "" ]]; then
                echo "Input folder does not exist or no VOMS proxy... skipping" ${array[0]}
                continue
            fi
            filelist=(`xrdfs dcache-cms-xrootd.desy.de ls ${array[5]} | grep root`)
        elif [[ ${array[3]} -eq "0" ]]; then
            filelist=(`ls ${array[5]} | grep root`)
        fi
        size=${#filelist[@]}

	if [ ! -e "${conf_dir}/${array[0]}.CONF" ]; then
            echo "Creating filelist for" ${array[0]}
            if [ ${array[3]} -eq "1" ]; then
                echo "root://131.169.191.218:1094/"${array[5]} ${array[1]} ${array[4]} > ${conf_dir}/${array[0]}.CONF
            elif [ ${array[3]} -eq "0" ]; then
                echo ${array[5]} ${array[1]} ${array[4]} > ${conf_dir}/${array[0]}.CONF
            fi
            for ((i=0; i<${#filelist[@]}; i++))
            do
		if [ ${array[3]} -eq "1" ]; then
		    echo ${filelist[i]##*/} >> ${conf_dir}/${array[0]}.CONF
		elif [ ${array[3]} -eq "0" ]; then
                    echo ${filelist[i]} >> ${conf_dir}/${array[0]}.CONF
		fi
            done

	elif [[ -e "${conf_dir}/${array[0]}.CONF" ]] && [[ `wc -l ${conf_dir}/${array[0]}.CONF | cut -d' ' -f1` != $((size+1)) ]]; then
            smallsize=`wc -l ${conf_dir}/${array[0]}.CONF | cut -d' ' -f1`
            echo "Regenerating filelist for" ${array[0]}":"
            echo ${size} "files found in" ${array[5]} "but only" $((smallsize-1)) "listed in file"
            if [[ ${array[3]} -eq "1" ]]; then
                echo ${array[5]} ${array[1]} ${array[4]} > ${conf_dir}/${array[0]}.CONF
            elif [[ ${array[3]} -eq "0" ]]; then
                echo ${array[5]} ${array[1]} ${array[4]} > ${conf_dir}/${array[0]}.CONF
            fi
            for ((i=0; i<${#filelist[@]}; i++))
            do
                echo ${filelist[i]} >> ${conf_dir}/${array[0]}.CONF
            done
        fi
    fi
done < "${conf_dir}/xsecs.dat"

echo " "
echo "Done checking config files."
echo "Submitting jobs."
echo " "

for file in ${conf_dir}/*CONF
do
    info=( $(head -n 1 $file) )
    outname=${file%.*}
    outname=${outname##*/}
    outputDir=${outputBase}/${outname}
    workDir=${CMSSW_BASE}/src/
    soFile=`echo $runMacro | sed 's/\./_/'`.so
    pcmFile=`echo $runMacro | sed 's/\./_/'`*pcm
    script=runjobs.sh
    vomsProxy=`voms-proxy-info | grep -oh "\w*x509up_u\w*"`

    # check a few things                                                                                      
    if [ ! "$CMSSW_BASE" ]; then
        echo "-------> error: define cms environment."
        exit 1
    fi
    if [ macros/$runMacro -nt macros/$soFile ]; then
        echo "-------> error: forgot to recompile run macro."
        exit 1
    fi
    if [[ "${vomsProxy:0:1}" != "x" ]]; then
	echo "-------> error: no voms-proxy found, jobs for files outside CERN will fail"
    else
	cp /tmp/${vomsProxy} $workDir
    fi

    mkdir -p ${outputDir}
    
    cp rootlogon.C             $workDir
    cp $runMacro               $workDir
    cp $soFile                 $workDir
    cp $pcmFile                $workDir

    sed 1d ${file} | while read line
    do
	#if [[ ! -e ${outputDir} ]]; then
	    #continue
	    #echo "Output directory doesn't exist. Not submitting."
	#if [[ `bjobs -w 2> /dev/null | grep ${line}` ]]; then
	    #continue
	    #echo "Job is currently running. Not submitting."
	#if [[ `grep "File broken" ${outputDir}/out.${line%.*}.txt` ]]; then
	    #continue
	    #echo "Input file is broken. Not submitting."
	#if [[ -e ${outputDir}/${line} ]] && [[ `grep "Selection complete" ${outputDir}/out.${line%.*}.txt` ]]; then
	if [[ -e ${outputDir}/${line} ]]; then
	    continue
	    #echo "Output file exists. Not submitting."
	else
	    echo $script $workDir $outputDir ${info[0]} ${line} ${info[1]} ${info[2]} $runMacro $soFile $pcmFile $vomsProxy
	    #bsub -o ${outputDir}/out.${line%.*}.txt -e ${outputDir}/err.${line%.*}.txt -q 8nm ${script} $workDir $outputDir ${info[0]} ${line} ${info[1]} ${info[2]} $runMacro $soFile $pcmFile
	fi
    done 
done

echo " " 
echo "Done submitting jobs."