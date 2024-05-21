#!/bin/bash
rm -f "summary.txt"
rm -f "summary_energy.txt"
echo -n >summary.txt
echo -n >summary_energy.txt
input="jobs.txt"
index=0
if [[ -f "$input" ]]; then
    while IFS= read -r line
    do
        pid=`echo "$line" | awk '{print $1}'`
        idir=`echo "$line" | awk '{print $2}'`
        output=`echo "$line" | awk '{print $3}'`
        p=`sacct -j $pid --format=JobID,JobName%60,State| egrep -v "batch|----|JobID" | head -1`
        echo $p >>"summary.txt"
        if [[ ! -z `echo $p | egrep "COMPLETED"` ]]; then
            en=`egrep "E\("  $output | tail -1 |awk '{print $5}'`
            enmp2=`egrep "UMP2" $output | tail -1 | awk '{{print $6}}'|tr 'D' 'E'`
            [[ -z "$enmp2" ]] && en=$en || tmp=$enmp2
            [[ ! -z "$tmp" ]] && en=$(printf "%.14f" $tmp)
            time_s=`egrep 'Elapsed' $output  | awk '{print $3*86400+$5*3600+$7*60+$9}' | tail -1`

            if [ $index -eq 0 ]; then
                 e0=$en
            fi
            index=`echo $index+1| bc -l`
            erel=`printf "%.3f" $(echo "scale=3; (($en)-($e0))*627.5095"|bc) `
            echo "$idir $pid $en $erel $time_s" >>summary_energy.txt
         #echo "$idir $pid $en $erel"
        fi
    done < "$input"

    [[ ! -z `egrep "PENDING|RUNNING" summary.txt` ]]  || echo -n >done
fi
#echo -n >done
