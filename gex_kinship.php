<?php

include_once('gex_info.php');

function gex_kinship_analysis($config){
	$command=("R CMD BATCH --no-save --no-restore '--args a=1 b=19 genotypefolder=\"".gex_get_genotype_folder($config)."\" datafolder=\"".gex_get_kinship_folder($config).'"' ."' kinship.analysis.R kinship.analysis.Rout ");
    printf($command."\n");
    system($command);
}


?>

	 