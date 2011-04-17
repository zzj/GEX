<?php

include_once('gex_info.php');
include_once('gex_init.php');


function gex_summary($config, $type){
	 $summaryfolder=gex_get_summary_folder($config);
	 $geneexpfolder=gex_get_gene_expression_folder($config);
     $string=gex_r_command($config);
     $string=gex_r_append_parameter($string,'resultfolder',$summaryfolder);
     $string=gex_r_append_parameter($string,'geneexpfolder',$geneexpfolder);
     $string=gex_r_append_parameter($string,'fun',$type);
     $string=gex_r_append_parameter($string,'step', $config->kinship_region_size);
     $string=gex_r_append_parameter($string,'rough.scan', $config->rough_scan);
     $string=gex_r_append_parameter($string,'root',gex_get_project_result_folder($config));
     $string=gex_r_command_end($string,'summary.R','/dev/tty');
     printf($string."\n");
     system($string);
}

?>