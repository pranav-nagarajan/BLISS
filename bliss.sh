#!/bin/bash
# Complete BLISS analysis pipeline

python bliss.py $1 --background $2 --background $4 --background $6 --cutoff 2  --multi --compress --output bliss_output_one.txt
echo "Processing: ON File #1 Complete."
python bliss.py $3 --background $2 --background $4 --background $6 --cutoff 2  --multi --compress --output bliss_output_two.txt
echo "Processing: ON File #2 Complete."
python bliss.py $5 --background $2 --background $4 --background $6 --cutoff 2  --multi --compress --output bliss_output_three.txt
echo "Processing: ON File #3 Complete."

python plotter.py --inputs bliss_output_one.txt --inputs bliss_output_two.txt --inputs bliss_output_three.txt --signal $1 --signal $3 --signal $5 --background $2 --background $4 --background $6
