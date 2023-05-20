#!/bin/bash
timelimit=3600
parallelscplex=0 # 1: use parallelism of cplex; 0: not use 
gnuparalleltest=1 # 1: use GNU parallel to speed up test; 0: not use
algorithms=("bppl"  "bpknnk3r25u"  "bpknnk3r2u"  "bpknnk5r15u"  "bpknnk5r25u"  "bpknnk5r2u"  "bpknnk1r15d"    "bpknnk1r25d"    "bpknnk1r2d"    "bpknnk3r15d"    "bpknnk3r25d"    "bpknnk3r2d"    "bpknnk5r15d"    "bpknnk5r25d"    "bpknnk5r2d") 
datapath="train_data"
testpath="train_results"
settingpath="bp_train_setting"

# create and clear test files
#python3 ./testdir.py

runInstance() {
    parallelscplex=$1
    timelimit=$2
    benchmark=$3
    instance=$4
    algo=$5
    settingpath="bp_train_setting"

    benchmarkpath=$benchmark
    logpath="train_results"
    resultpath="$testpath/$benchmark/results"

    echo $timelimit $instance $benchmark $benchmarkpath $logpath $resultpath $algo
    

    python3 ./checkexec.py  $benchmark $algo $instance

    if [ $algo == "bcsocp" ]
    then
        python ./bc.py $timelimit $parallelscplex $instance  $benchmarkpath $logpath/$algo $resultpath/$algo
        #echo $timelimit $parallelscplex $instance  $benchmarkpath $logpath/$algo $resultpath/$algo
        #echo =
    else
        #echo $logpath $instance $algo $instance $algo
        build/cbp -c "set limits time $timelimit" -c  "set cbp is_parallelscplex $parallelscplex" -c "set load $settingpath/$algo.set" -c  "read $benchmarkpath/$instance" -c "opt write statistics" -c "$logpath/$algo/${instance}_${algo}.log"  -c "quit"
    fi

}
export -f runInstance


benchmarks=$(ls ${datapath})
benchmarks=("train_data")

for benchmark in $benchmarks
do
    instances=$(ls $benchmark)
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
        parallel --will-cite --jobs 90% runInstance  "$parallelscplex" "$timelimit" "$benchmark"  ::: "$instances" :::  "${algorithms[@]}"
    fi
    #find   $datapath/$benchmark -name *.cbp
done

