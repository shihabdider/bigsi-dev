N=4
for j in {1..5};
do
    ((i=i%N)); ((i++==0)) && wait
    sleep 6 &
done
