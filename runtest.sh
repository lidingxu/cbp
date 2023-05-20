#!/bin/bash
timelimit=3600
timebound=3800
parallelscplex=0 # 1: use parallelism of cplex; 0: not use 
gnuparalleltest=1 # 1: use GNU parallel to speed up test; 0: not use
algorithms=("bcsocp" "bppl" "bpsocp"  "bphybrid") 
datapath="data"
testpath="test"
settingpath="bp_setting"

# create and clear test files
#python3 ./testdir.py

export PYTHONPATH=/home/lxu/software/docplex/lib/python/:$PYTHONPATH

runInstance() {
    parallelscplex=$1
    timelimit=$2
    benchmark=$3
    instance=$4
    algo=$5
    datapath="data"
    testpath="test"
    settingpath="bp_setting"

    benchmarkpath="$datapath/$benchmark"
    logpath="$testpath/$benchmark/logs"
    resultpath="$testpath/$benchmark/results"

    echo $timelimit $instance $benchmark $benchmarkpath $logpath $resultpath $algo
    
  

    if [ $algo == "bcsocp" ]
    then
        python3.8 ./bc.py $timelimit $parallelscplex $instance  $benchmarkpath $logpath/$algo $resultpath/$algo
        #echo $timelimit $parallelscplex $instance  $benchmarkpath $logpath/$algo $resultpath/$algo
        #echo =
    else
	echo ""
        #echo $logpath $instance $algo $instance $algo
        build/cbp -c "set limits time $timelimit" -c  "set cbp is_parallelscplex $parallelscplex" -c "set load $settingpath/$algo.set" -c  "read $benchmarkpath/$instance" -c "opt write statistics" -c "$logpath/$algo/${instance}_${algo}.log"  -c "quit"
    fi

}
export -f runInstance


benchmarks=$(ls ${datapath})


for benchmark in $benchmarks
do
    instances=$(ls $datapath/$benchmark)
    if [ $gnuparalleltest == 0 ]
    then
        for instance in  $instances
        do
            for algo in ${algorithms[@]}
                do
                    runInstance "$parallelscplex" "$timelimit" "$benchmark"  "$instance" "$algo"
            done
        done
    else
        parallel --will-cite --jobs 85% --timeout $timebound runInstance  "$parallelscplex" "$timelimit" "$benchmark"  ::: "$instances" :::  "${algorithms[@]}"
        #$instances | parallel --will-cite   --dryrun  "printls {}"
        #parallel --will-cite  printls0 para ::: 1
        #parallel --will-cite  printls "$benchmark" para  ::: "$instances" :::  "${algorithms[@]}"
        #break 
        #parallel --will-cite -j 4 --dryrun  printls ::: "${instances[@]}" :::  "${algorithms[@]}"  
        #parallel --will-cite  runInstance ::: "$algorithms"  "${instances[@]}"  "$benchmark" "para"
    fi
    #find   $datapath/$benchmark -name *.cbp
done

