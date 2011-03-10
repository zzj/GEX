<?php

include_once('gex_info.php');
include_once('gex_init.php');


function gex_variance_analysis($config){
	 $variancefolder=gex_get_variance_folder($config);
	 $geneexpfolder=gex_get_gene_expression_folder($config);
	 $meta_idx=gex_load_gene_meta($config);
	 $gene_idx=gex_load_gene_expression_file_info($config);
	 $fd=fopen("command_list","w+");
	 foreach ($gene_idx as $key => $valuep){
		  if (!isset($meta_idx[$key]))
			   continue;
		  $value=$meta_idx[$key];
		  $string=gex_r_command($config);
		  $string=gex_r_append_parameter($string,'chrid',(int)$value[0]);
		  $string=gex_r_append_parameter($string,'genestart',(int)$value[1]);
		  $string=gex_r_append_parameter($string,'geneend',(int)$value[2]);
		  $string=gex_r_append_parameter($string,'rangestart',$value[1]-5000000);
		  $string=gex_r_append_parameter($string,'rangeend',$value[2]+5000000);
		  $string=gex_r_append_parameter($string,'phenotypename',$key);
		  $string=gex_r_append_parameter($string,'phenotypefile', $geneexpfolder.$gene_idx[$key].'/'.$key);
		  $resultfolder=$variancefolder.$gene_idx[$key].'/';
		  gex_check_folder($resultfolder);
		  $string=gex_r_append_parameter($string,'datafolder',$resultfolder);
		  $string=gex_r_command_end($string,'variance.1M.analysis.R');
		  fprintf($fd,$string."\n");
	 }
	 fclose($fd);
}

function gex_variance_local_opt_analysis($config){
	 $variancefolder=gex_get_variance_folder($config);
	 $geneexpfolder=gex_get_gene_expression_folder($config);
	 $meta_idx=gex_load_gene_meta($config);
	 $gene_idx=gex_load_gene_expression_file_info($config);
	 $fd=fopen("command_list","w+");
	 foreach ($gene_idx as $key => $valuep){
		  if (!isset($meta_idx[$key]))
			   continue;
		  $value=$meta_idx[$key];
		  $string=gex_r_command($config);
		  $string=gex_r_append_parameter($string,'chrid',(int)$value[0]);
		  $string=gex_r_append_parameter($string,'genestart',(int)$value[1]);
		  $string=gex_r_append_parameter($string,'geneend',(int)$value[2]);
		  $string=gex_r_append_parameter($string,'rangestart',$value[1]-5000000);
		  $string=gex_r_append_parameter($string,'rangeend',$value[2]+5000000);
		  $string=gex_r_append_parameter($string,'phenotypename',$key);
		  $string=gex_r_append_parameter($string,'phenotypefile', $geneexpfolder.$gene_idx[$key].'/'.$key);
		  $resultfolder=$variancefolder.$gene_idx[$key].'/';
		  gex_check_folder($resultfolder);
		  $string=gex_r_append_parameter($string,'datafolder',$resultfolder);
		  $string=gex_r_command_end($string,'variance.local.opt.analysis.R');
		  fprintf($fd,$string."\n");
	 }
	 fclose($fd);
}


?>

