#!/bin/bash

# Valores posibles para cada opción
n_values=("8")
mode_values=("YKR" "YSR" "YOR")
tolerance_values=("f" "ftolerance")

# Ejecutar todas las combinaciones
for n in "${n_values[@]}"; do
  for mode in "${mode_values[@]}"; do
    for tol in "${tolerance_values[@]}"; do
      # Verificar que no se use 'r' con 'YOR'
      if [ "$mode" == "YOR" ] && [ "$tol" == "r" ]; then
        continue
      fi
      mkdir -p "graphs/$mode/$n/$tol"
      # Ejecutar el comando
      echo "Ejecutando: python3 generate_graphics.py --n $n --mode $mode --tol $tol"
      python3 generate_graphics.py --n $n --mode $mode --tol $tol
    done
  done
done
