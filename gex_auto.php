<?php
function wait_grid(){
     while( ($ret=system('qstat')) != ""){
          sleep(60);
     }
}
$config=$argv[1];
system('php gex.php '.$config." kinship_analysis");
system('php gex.php '.$config." variance_global_opt_analysis");
system('qsub test.sh');
wait_grid();
system('php gex.php '.$config." emma_variance_analysis");
system('qsub test.sh');
wait_grid();
?>