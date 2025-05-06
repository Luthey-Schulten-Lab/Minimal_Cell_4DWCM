source ~/anaconda3/etc/profile.d/conda.sh
conda activate LM_Cell

date="April2025"

# start_time and end_time are in minutes
python output_concs.py -date $date -start_time 1 -end_time 94