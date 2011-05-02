<?php

include_once('gex_info.php');
include_once('gex_init.php');


function gex_emma_analysis($config,$fun){
	 $emmafolder=gex_get_emma_folder($config);
	 $geneexpfolder=gex_get_gene_expression_folder($config);
	 $variancefolder=gex_get_variance_folder($config);
	 $meta_idx=gex_load_gene_meta($config);
	 $gene_idx=gex_load_gene_expression_file_info($config);
	 $fd=fopen("command_list","w+");
     $skip=$config->rough_scan;
	 foreach ($gene_idx as $key => $valuep){
		  if (!isset($meta_idx[$key]))
			   continue;
          if ($skip==0) {
               $skip=$config->rough_scan;
          }
          else {
               $skip--;
               continue;
          }
		  $value=$meta_idx[$key];
		  $string=gex_r_command($config);
		  $string=gex_r_append_parameter($string,'chrid',(int)$value[0]);
		  $string=gex_r_append_parameter($string,'genestart',(int)$value[1]);
		  $string=gex_r_append_parameter($string,'geneend',(int)$value[2]);
		  $string=gex_r_append_parameter($string,'rangestart',$value[1]-5000000);
		  $string=gex_r_append_parameter($string,'rangeend',$value[2]+5000000);
		  $string=gex_r_append_parameter($string,'phenotypename',$key);
		  $string=gex_r_append_parameter($string,'fun',$fun);
		  $string=gex_r_append_parameter($string,'phenotypefile', $geneexpfolder.$gene_idx[$key].'/'.$key);
		  $string=gex_r_append_parameter($string,'variancefile', $variancefolder.$gene_idx[$key].'/'.$key."_global_nlminb.Rdata");
		  $resultfolder=$emmafolder.$gene_idx[$key].'/';
		  gex_check_folder($resultfolder);
		  $string=gex_r_append_parameter($string,'datafolder',$resultfolder);
		  $string=gex_r_command_end($string,'emma.analysis.R');
		  fprintf($fd,$string."\n");
	 }
	 fclose($fd);
}


?>

