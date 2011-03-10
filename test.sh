#$ -S /bin/bash
#$ -t 1-2290
#$ -o /home/zzj/temp
#$ -j y
#$ -q all.q

cd ~/Research/gamli/src/gex
php  test.php config.json $SGE_TASK_ID
