<?php

include_once('gex_info.php');
include_once('gex_init.php');


function gex_lasso_analysis($config){
	 gex_r_analysis($config, gex_get_lasso_folder($config), "lasso.analysis.R", "analysis");
}

function gex_lasso_plot($config){
	 gex_r_analysis($config, gex_get_lasso_folder($config), "lasso.analysis.R", "plot");
}


?>