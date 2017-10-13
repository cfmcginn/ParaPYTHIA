#!/bin/bash

#runVar=(240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400)
runVar=(0 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400)
#runVar=(1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000 4100 4200 4300 4400 4500 4600 4700 4800 4900 5000 5100 5200 5300 5400 5500 5600 5700 5800 5900 6000 6100 6200 6300 6400 6500 6600 6700 6800 6900 7000 7100 7200 7300 7400 7500)

#runVar=(0 10)
maxjobs=10

for i in `seq 0 $((${#runVar[@]}-2))`
do
    if (( $(($i % $maxjobs)) == 0 )) ; then
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $i wait
    fi

    valStart=${runVar[$i]}
    valEnd=$((${runVar[$(($i+1))]} - 1))

    sh $PWD/runLoopTemplate.sh $valStart $valEnd > loop_$valStart\_$valEnd.log &

    sleep 5
done

wait

echo "ParaLoop complete"