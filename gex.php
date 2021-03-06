<?php

include_once("gex_init.php");
include_once("gex_kinship.php");	 
include_once("gex_variance.php");	 
include_once("gex_emma.php");	 
include_once("gex_std.php");	 
include_once("gex_lasso.php");	 
include_once("gex_summary.php");	 
if (!isset($argv[1])){
	 $config_file="config.json";
}
else {
	 $config_file=$argv[1];
}


$config=json_decode(file_get_contents($config_file));

if ($config==NULL) die("Can not parse ". $config_file);
echo " Start initializing Project $config->project_name ...\n" ;
$config->emma_result_folder="emma";
$config->variance_result_folder="variance";
$config->std_result_folder="std";
$config->lasso_result_folder="lasso";
$config->summary_result_folder="summary";
$config->kinship_result_folder="kinship";


$actions=split(',', $config->actions);
if (isset($argv[2])){ $actions=split(',',$argv[2]);}
system("rm ".$config->command_list);
foreach($actions as $action){
	 if ($action=='init')
		  gex_init($config);
	 else if ($action=='kinship_analysis')
		  gex_kinship_analysis($config);
	 else if ($action=='emma_analysis')
		  gex_emma_analysis($config);
	 else if ($action=='variance_analysis')
		  gex_variance_analysis($config);
	 else if ($action=='variance_local_opt_analysis')
		  gex_variance_local_opt_analysis($config);
	 else if ($action=='variance_global_opt_analysis')
		  gex_variance_global_opt_analysis($config);
	 else if ($action=='std_analysis')
		  gex_std_analysis($config);
	 else if ($action=='std_plot')
		  gex_std_plot($config);
	 else if ($action=='lasso_analysis')
		  gex_lasso_analysis($config);
	 else if ($action=='lasso_plot')
		  gex_lasso_plot($config);
	 else if (substr($action,-7)=='summary')
		  gex_summary($config, substr($action, 0, -8));
	 else die('unknown action '.$action."\n");
}

?>