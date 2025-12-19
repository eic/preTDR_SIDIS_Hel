#!/bin/bash
runlist=$1
NAME=18x275minQ1
STEP=1

nrun=0
gr=0
run_first=0
run_last=0
counter=0

cat ${runlist} | while read input
do
    file=$(echo $input| awk '{print($1)}')
    echo processing ${file}
    #SUBFILES+=($input)
    SUBFILES=$input
    ((nrun++))
    if [ `expr $nrun % $STEP` = 0 ]; then
	((gr++))
	((counter++))
        #echo counter $counter
        #if [ $counter -eq 10 ]; then
        #    sleep 100s
        #    counter=0
        #fi
        mkdir condor_scripts/${NAME}_${gr}
	mkdir /gpfs/mnt/gpfs02/eic/cvhulse/epic/analysis/${NAME}_${gr}
	cp run_${NAME}.sh condor_scripts/${NAME}_${gr}
	cp con_${NAME}.sub condor_scripts/${NAME}_${gr}
	echo "Output=/gpfs/mnt/gpfs02/eic/cvhulse/epic/analysis/${NAME}_${gr}/ev${gr}.out" >> condor_scripts/${NAME}_${gr}/con_${NAME}.sub
        echo "Error=/gpfs/mnt/gpfs02/eic/cvhulse/epic/analysis/${NAME}_${gr}/ev${gr}.err" >> condor_scripts/${NAME}_${gr}/con_${NAME}.sub
        echo "Log=/gpfs/mnt/gpfs02/eic/cvhulse/epic/analysis/${NAME}_${gr}/ev${gr}.log" >> condor_scripts/${NAME}_${gr}/con_${NAME}.sub
	echo  'Arguments  = ' ${gr} >> condor_scripts/${NAME}_${gr}/con_${NAME}.sub
	echo  '+JobFlavour = "workday"' >> condor_scripts/${NAME}_${gr}/con_${NAME}.sub
	echo "Queue" >> condor_scripts/${NAME}_${gr}/con_${NAME}.sub
	cd condor_scripts/${NAME}_${gr}
	ln -s ../../event_eval.exe ./
	touch runlist.txt
	for com in $SUBFILES; do
            echo ${com} >> runlist.txt
        done
	SUBFILES=""
	condor_submit ./con_${NAME}.sub
	cd ../../
    fi
done
