<?php

include_once('gex_info.php');
include_once('gex_init.php');


function gex_std_analysis($config){
	 gex_r_analysis($config, gex_get_std_folder($config), "std.analysis.R", "analysis");
}

function gex_std_plot($config){
	 gex_r_analysis($config, gex_get_std_folder($config), "std.analysis.R", "plot");
}


?>