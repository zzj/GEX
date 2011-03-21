<?php


function gex_check_folder($folder){
	 if (!file_exists($folder)){
		  mkdir($folder);
	 }
	 return;
  }


function gex_get_emma_folder($config){

	 $folder=gex_get_project_result_folder($config).$config->emma_result_folder."/";
	 if (!file_exists($folder)){
		  mkdir($folder);
	 }
	 return $folder;
}

function gex_get_summary_folder($config){
	 $folder=gex_get_project_result_folder($config).$config->summary_result_folder."/";
	 if (!file_exists($folder)){
		  mkdir($folder);
	 }
	 return $folder;
}

function gex_get_variance_folder($config){

	 $folder=gex_get_project_result_folder($config).$config->variance_result_folder."/";
	 if (!file_exists($folder)){
		  mkdir($folder);
	 }
	 return $folder;
}


function gex_get_std_folder($config){

	 $folder=gex_get_project_result_folder($config).$config->std_result_folder."/";
	 if (!file_exists($folder)){
		  mkdir($folder);
	 }
	 return $folder;
}


function gex_get_lasso_folder($config){

	 $folder=gex_get_project_result_folder($config).$config->lasso_result_folder."/";
	 if (!file_exists($folder)){
		  mkdir($folder);
	 }
	 return $folder;
}



function gex_get_kinship_folder($config){

	 $folder=gex_get_project_result_folder($config).$config->kinship_result_folder."/";
	 if (!file_exists($folder)){
		  mkdir($folder);
	 }
	 return $folder;
  }

function gex_get_project_folder($config){
	 if (!file_exists($config->project_folder)){
          mkdir($config->project_folder);
     }
	 return $config->project_folder;
}



function gex_get_project_result_folder($config){
     
     $resultfolder=gex_get_project_folder($config).$config->project_alias;
     if (!file_exists($resultfolder)){
          mkdir($resultfolder);
     }

	 return $resultfolder.'/';
}


function gex_get_genotype_folder($config){
	 
     $resultfolder=gex_get_project_result_folder($config)."genotype";
     if (!file_exists($resultfolder)){
          mkdir($resultfolder);
     }

	 return $resultfolder.'/';

}

function gex_get_gene_expression_folder($config){
	 
     $resultfolder=gex_get_project_result_folder($config)."gene_expression";
     if (!file_exists($resultfolder)){
          mkdir($resultfolder);
     }
	 return $resultfolder.'/';

}


function gex_r_append_parameter($string, $name, $value){
	 if (is_string($value)){
		  $quate='"';
	 }
	 else $quate="";

	 $string=$string.' '.$name."=".$quate.$value.$quate.' ';
	 return $string;
}

function gex_r_command($config){
	 $command= "R CMD BATCH --no-save --no-restore '--args a=1 b=100 genotypefolder=\"".gex_get_genotype_folder($config)."\" kinshipfolder=\"".gex_get_kinship_folder($config).'"';
	 return $command;
}

function gex_r_command_end($string,$file,$output=NULL){
	 if ($output==NULL)
          $string=$string."'  $file "." $file"."out";
     else $string=$string."' $file $output";
 	 return $string;
}

function gex_load_gene_expression_file_info($config){
	 if ($config->verbose)
		  printf("loading gene expression information file\n");
	 $info=file(gex_get_gene_expression_folder($config)."gene_list", $flag=FILE_IGNORE_NEW_LINES);
	 $index=array();
	 foreach($info as $value){
		  $content=split("\t",$value);
//		  print_r($content);
		  $index[$content[0]]=$content[1];
	 }
	 return $index;
}


function gex_r_analysis($config, $resultfolder, $r_file, $r_fun=NULL){
	 $geneexpfolder=gex_get_gene_expression_folder($config);
	 $meta_idx=gex_load_gene_meta($config);
	 $gene_idx=gex_load_gene_expression_file_info($config);
	 $fd=fopen($config->command_list,"a");
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
          if ((int)$value[0]==0) continue;
		  $string=gex_r_command($config);
		  $string=gex_r_append_parameter($string,'chrid',(int)$value[0]);
		  $string=gex_r_append_parameter($string,'genestart',(int)$value[1]);
		  $string=gex_r_append_parameter($string,'geneend',(int)$value[2]);
		  $string=gex_r_append_parameter($string,'rangestart',$value[1]-5000000);
		  $string=gex_r_append_parameter($string,'rangeend',$value[2]+5000000);
		  $string=gex_r_append_parameter($string,'phenotypename',$key);
		  $string=gex_r_append_parameter($string,'phenotypefile', $geneexpfolder.$gene_idx[$key].'/'.$key);
		  if ($r_fun!=NULL){
			   $string=gex_r_append_parameter($string, 'fun', $r_fun);
		  }
		  $resultfolder2=$resultfolder.$gene_idx[$key].'/';
		  gex_check_folder($resultfolder2);
		  $string=gex_r_append_parameter($string,'datafolder',$resultfolder2);
          $string=gex_r_append_parameter($string,'kinshipfolder',gex_get_kinship_folder($config));
          $string=gex_r_append_parameter($string,'step',$config->kinship_region_size);
		  $string=gex_r_command_end($string,$r_file);
		  fprintf($fd,$string."\n");
	 }
	 fclose($fd);
}

?>