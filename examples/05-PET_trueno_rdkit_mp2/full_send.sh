#!/bin/bash

# NJOBS          --> Number of jobs sent to the slurm system
# MAXJOBSINSLURM --> Maximum number of jobs in the slurm system
# JOBSEND        --> Number of jobs finished in the slurm system
# TOTALJOBS      --> Jobs to be sent to the slurm system
# jobs.txt       --> Info of the jobs sent or finished in the slurm system

MAXJOBSINSLURM=50

NJOBS=`squeue -h |wc -ll`

while true; do
    for ifile in `ls *.com`; do
        base="${ifile%.*}"
        echo $idir, $base, ${NJOBS}, ${JOBSEND}, ${TOTALJOBS}

        if [[ ! -e ./jobs.txt ]]; then
            echo -n >./jobs.txt
        fi

        if [[ $NJOBS -lt $MAXJOBSINSLURM  ]]; then
            if [[  -z `egrep $ifile ./jobs.txt` ]]; then
                sbatch ${base}.sh 1>tmp.txt;
                jobid=`awk '{print $NF}' tmp.txt`
                echo "${jobid} ${base} ${base}.log" >>./jobs.txt
                rm tmp.txt
            fi
        fi

        NJOBS=`squeue -h |wc -ll`
        #echo $NJOBS
        JOBSEND=`cat ./jobs.txt | wc -ll`
        TOTALJOBS=`ls -ld ./*_[0-9]*com |wc -ll`
        echo "JOBSEND: ${JOBSEND}, TOTALJOBS: ${TOTALJOBS}"

    done
    if [[ ${JOBSEND} -ge ${TOTALJOBS} ]]; then
        break
    fi
    # Each 15 seconds checks the jobs
    sleep 15
done
