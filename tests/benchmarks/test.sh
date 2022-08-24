#HG38_SIM=hg38/simulation
#echo metrics/${HG38_SIM}/error_rate_trials_sensitivity.txt \

errs=( 001 002 003 004 005 )
lengths=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )
params=( )
for err in ${errs[@]}; do
    for length in ${lengths[@]}; do
        param="${length}_${err}"
        params+=( "$param" )
    done
done

for value in "${params[@]}"
do
     echo $value
done
