#$ -S /bin/bash
#$ -t 1-2856
#$ -o /home/zzj/temp
#$ -j y
#$ -q all.q
#$ -cwd
php  test.php config.json $SGE_TASK_ID
