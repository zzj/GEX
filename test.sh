#$ -S /bin/bash
#$ -t 1-18998
#$ -o /home/zzj/temp
#$ -j y
#$ -q all.q

cd ~/Research/gex
php  test.php gse13870.json $SGE_TASK_ID
