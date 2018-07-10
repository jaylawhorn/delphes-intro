# Now let's just run selection.C interactively as you need to have .so files to submit the lxbatch jobs (I don't like to compile in lxbatch as it slows things down and before it wouldn't compile in slc6 anyways when you submitted the jobs)
root -l -b -q selection.C+
#Open conf.txt file and select the samples you would like to analyze
#Open submitjobs.sh and modify the hardcoded cond_dir and outputBase
#copy the conf.txt file to the conf_dir as the script looks for it over there and not in your current directory
./submitjobs.sh
#and enjoy the physics (optional)

*    selection.C: ROOT macro that actually performs selection

See file for further documentation.

*    conf.txt: configuration file that specifies sample name, cross section, flag for if is it in eos or not, and the base file path for each sample

For example, the first line is "HHToTTBB_14TeV 2.92 0 /afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV/", which means it's going to look at files matching "/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV/HHToTTBB_14TeV*" using normal "ls" (the EOS/not flag is to toggle between "ls" and "eos ls", as well as toggle the "root://eoscms.cern.ch/" prefix.
