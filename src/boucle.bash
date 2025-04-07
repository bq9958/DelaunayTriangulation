#!/bin/bash

THRESHOLD=0.01

# compilation
make -f makeadaptation.make
make -f makeerror.make

# Initialization
error=1.0
iteration=0
MAX_ITER=10


while (( $(echo "$error > $THRESHOLD" | bc -l) )); do
    echo ">>> Iteration $iteration"

    # run main_adaptation, output: metric.sol
    echo "Running main_adaptation..."
    ./main_adaptation

    # run feflo
    echo "Running feflo..."
    feflo.a_2d -in maillage.mesh -itp maillage.niveaugris.sol -met maillage.met.sol \
        -hgrad 1.5 -out maillage.adapte.mesh -noref

    # run main_error，output : error
    echo "Running main_error..."
    error=$(./main_error)
    echo "Current error: $error"
    
    ((iteration++))
    if [ "$iteration" -ge "$MAX_ITER" ]; then
        echo "Reached max iterations ($MAX_ITER). Stopping."
        break
    fi
done

echo "Adaptation finished. Final error: $error"