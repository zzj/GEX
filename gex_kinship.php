<?php

include_once('gex_info.php');

function gex_kinship_analysis($config){
     $string=gex_r_command($config);
     $string=gex_r_append_parameter($string,'datafolder',gex_get_kinship_folder($config));
     $string=gex_r_append_parameter($string,'root',gex_get_project_result_folder($config));
     $string=gex_r_append_parameter($string,'step',$config->kinship_region_size);
     $string=gex_r_append_parameter($string,'overlap',$config->overlap);
     $string=gex_r_command_end($string,'kinship.analysis.R');
     printf($string."\n");
     system($string);
}


?>

	 