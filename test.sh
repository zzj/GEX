#$ -S /bin/bash
#$ -t 1-952
#$ -o /home/zzj/temp
#$ -j y
#$ -q all.q
#$ -cwd
php  test.php config.json $SGE_TASK_ID
