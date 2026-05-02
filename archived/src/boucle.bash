#!/bin/bash

# compilation
make -f makeadaptation.make
make -f makeerror.make

# make dir and initialize mesh
rm -rf ../dir
rm -rf ../output
mkdir -p ../output
mkdir -p ../dir
cp ../data/joconde.lowres.mesh ../dir/maillage.mesh
cp ../data/joconde.lowres.sol ../dir/maillage.niveaugris.sol

# Initialization
iteration=0
MAX_ITER=40

exec > >(tee -a ../output/log.jnl) 2>&1

cd ../dir
echo "iteration,C" > ../output/complexity_log.dat
echo "iteration,NbrVer" > ../output/NbrVer_log.dat
echo "iteration,PSNR,delta" > ../output/PSNR_log.dat

while true; do
    echo ">>> Iteration $iteration"

    # run main_adaptation, output: metric.sol
    echo "Running main_adaptation..."
    if [ $iteration -eq 0 ]; then
        ./../build/adaptation maillage.mesh maillage.niveaugris.sol
    else
        ./../build/adaptation maillage.adapte.mesh maillage.niveaugris.itp.sol
    fi

    # Extract complexity value C from log file
    C=$(grep "C =" ../output/log.jnl | tail -n 1 | awk '{print $3}')
    echo "$iteration,$C" >> ../output/complexity_log.dat
    NbrVer=$(grep "NbrVer =" ../output/log.jnl | tail -n 1 | awk '{print $4}')
    echo "$iteration,$NbrVer" >> ../output/NbrVer_log.dat

    # run feflo
    echo "Running feflo..."
    if [ $iteration -eq 0 ]; then
        ./../bin/linux/feflo.a_2d -in maillage.mesh -itp maillage.niveaugris.sol -met ../output/maillage.met.sol \
            -hgrad 1.5 -out maillage.adapte.mesh -noref
    else
        ./../bin/linux/feflo.a_2d -in maillage.adapte.mesh -itp maillage.niveaugris.itp.sol -met ../output/maillage.met.sol \
            -hgrad 1.5 -out maillage.adapte.adapte.mesh -noref
    fi

    # run main_error，output : error
    echo "Running main_error..."
    if [ $iteration -eq 0 ]; then
        ./../build/error maillage.mesh maillage.niveaugris.sol maillage.adapte.mesh maillage.niveaugris.itp.sol
        PSNR=$(grep "PSNR =" ../output/log.jnl | tail -n 1 | awk '{print $3}')
        PSNR_prev=$PSNR
        echo "$iteration,$PSNR,0.0" >> ../output/PSNR_log.dat
    else 
        ./../build/error maillage.adapte.mesh maillage.niveaugris.itp.sol maillage.adapte.adapte.mesh maillage.niveaugris.itp.itp.sol
        PSNR=$(grep "PSNR =" ../output/log.jnl | tail -n 1 | awk '{print $3}')
        delta=$(echo "$PSNR - $PSNR_prev" | bc -l)
        echo "$iteration,$PSNR,$delta" >> ../output/PSNR_log.dat
        PSNR_prev=$PSNR
    fi 

    mv maillage.niveaugris.itp.itp.sol maillage.niveaugris.itp.sol
    mv maillage.adapte.adapte.mesh maillage.adapte.mesh
    
    ((iteration++))
    if [ "$iteration" -ge "$MAX_ITER" ]; then
        echo "Reached max iterations ($MAX_ITER). Stopping."
        break
    fi
done

echo "Adaptation finished."

# Visualization
./../vizir4/vizir4.exe -in ../dir/maillage.adapte.mesh -sol ../dir/maillage.niveaugris.itp.sol

echo "Plotting complexity curve with Gnuplot..."
gnuplot -persist << EOF
set datafile separator ","
set terminal pngcairo size 800,600 enhanced font 'Arial,10'
set output "../output/complexity_vs_iteration.png"
set xlabel "Iteration"
set ylabel "Complexity C"
set grid
plot "../output/complexity_log.dat" using 1:2 with linespoints title "C vs iteration"
EOF

echo "Plotting NbrVer curve with Gnuplot..."
gnuplot -persist << EOF
set datafile separator ","
set terminal pngcairo size 800,600 enhanced font 'Arial,10'
set output "../output/NbrVer_vs_iteration.png"
set xlabel "Iteration"
set ylabel "NbrVer"
set grid
plot "../output/NbrVer_log.dat" using 1:2 with linespoints title "NbrVer vs iteration"
EOF

echo "Plotting PSNR curve with Gnuplot..."
gnuplot -persist << EOF
set datafile separator ","
set terminal pngcairo size 800,600 enhanced font 'Arial,10'
set output "../output/PSNR_vs_iteration.png"
set xlabel "Iteration"
set ylabel "PSNR"
set grid
plot "../output/PSNR_log.dat" using 1:2 with linespoints title "PSNR vs iteration"
EOF

echo "Plotting ΔPSNR curve with Gnuplot..."
gnuplot -persist << EOF
set datafile separator ","
set terminal pngcairo size 800,600 enhanced font 'Arial,10'
set output "../output/delta_PSNR_vs_iteration.png"
set xlabel "Iteration"
set ylabel "Δ PSNR"
set grid
plot "../output/PSNR_log.dat" using 1:3 with linespoints title "ΔPSNR vs iteration"
EOF