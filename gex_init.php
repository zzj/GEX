<?php
  /** this function create a folder named project_alias in project_folder

	  then create a folder named gene_expression_folder,
	  it includes
	  (1) gene_meta.txt and gene_expression.txt
	  (2) a list of folders, which contains a file named by gene id, and encoded as json file,
	  the json file contains two fields, one is the meta info, and the other is gene expression data.

	  if the genotype file is standard, then it has only original genotype data. 

	  if the genotype file is hmm,
	  in this folder, it has two different type of genotype data,
	  one is original genotype data, and the other is probability of founder pairs 
  */

include_once('lib/getCSVvalues.php');
/* Load gene meta info file, create a dictionary with gene's position information */

function gex_load_gene_meta($config){

	 if ($config->verbose){
		  printf("Loading gene meta info\n");
	 }
	 
	 $meta=file($config->gene_meta_file, $flag=FILE_IGNORE_NEW_LINES);
	 $meta_idx=array();
	 foreach ($meta as $key => $value){
		  if ($key==0) continue;
		  $content=getCSVvalues($value);
		  $meta_idx[$content[0]]=array($content[$config->gene_meta_chr],$content[$config->gene_meta_start_position],$content[$config->gene_meta_end_position],$value);
	 }
	 return $meta_idx;
}


function gex_init($config){
     
     if (!file_exists($config->project_folder)){
          mkdir($config->project_folder);
     }
     
     $resultfolder=$config->project_folder.$config->project_alias;
     if (!file_exists($resultfolder)){
          mkdir($resultfolder);
     }
     $infofile=$config->project_alias.".info";

	 $subject_list=gex_init_joint_strains($config);
	 $subject_list=array_values($subject_list);
	 printf("There are total ".count($subject_list)." strains founded in both datasets\n");
	 $subject_list=gex_init_gene_expression($config,$subject_list);
	 gex_dump_strain_list($config, $resultfolder, $subject_list);
	 gex_init_genotype($config,$subject_list);
}

function gex_dump_strain_list($config, $resultfolder, $subject_list){
     $fd=fopen($resultfolder."/strain_list","w+");
     foreach ($subject_list as $value){
		  fprintf($fd, "%s\n", $value);
	 }
	 fclose($fd);
}

function gex_init_joint_strains($config){
	 list($gene_list, $subject_list, $data_list)=gex_load_gene_expression($config,NULL, false);
	 $hmm_strain_name=file($config->genotype_folder."strain_list", $flag=FILE_IGNORE_NEW_LINES);
	 print_r($hmm_strain_name);
	 return array_intersect($hmm_strain_name, $subject_list);
}

function gex_init_genotype($config, $subject_list){
	 $chr_list=split(',',$config->genotype_chr_folder);
	 $folder=$config->genotype_folder;
	 $strain_name=$subject_list;
	 $target_folder=$config->project_folder.$config->project_alias.'/'."genotype/";
	 if (!file_exists($target_folder)){
		  mkdir($target_folder);
	 }
	 $finfo=fopen($target_folder."marker_list","w+");
	 $id=1;
	 if ($config->verbose)
		  printf("Start generating genotype folder\n");
	 
	 for($chr_id=0;$chr_id<count($chr_list);$chr_id++){
		  if ($config->verbose)
			   printf("Working on Chromosome ".$chr_list[$chr_id]."\n");

		  $current_folder=$folder.$chr_list[$chr_id]."/";
		  if (!file_exists($current_folder)){
			   mkdir($current_folder);
		  }
		  $files=array();
		  $genotypes=array();
		  $markers=file($current_folder."Extra.txt");
		  for ($i=0;$i < count($strain_name); $i++){
			   if ($config->have_probability) $files[$i]=file($current_folder.$strain_name[$i].".txt.prob",$flag=FILE_IGNORE_NEW_LINES);
			   $genotypes[$i]=file($current_folder.$strain_name[$i].".txt",$flag=FILE_IGNORE_NEW_LINES);
		  }
          /*
		  if ($config->verbose){
			   printf("dumping position information ...\n");
		  }
		  $positions=array();
		  $fp=fopen($target_folder.$chr_list[$chr_id].".position","w+");
		  foreach ($markers as $value){
			   $info=split(',',$value);
			   array_push($positions, $info[2]);
		  }
          */
		  if ($config->verbose)
			   printf("dumping genotype data .. \n");
		  $fg=fopen($target_folder.$chr_list[$chr_id].".genotype","w+");
		  for ($i=1;$i <count($genotypes[0]); $i++){
			   for ($j=0; $j<count($genotypes);$j++){
					fprintf($fg,"%s\t", $genotypes[$j][$i]);
			   }
			   fprintf($fg,"\n");
		  }
		  fclose($fg);
		  if ($config->have_probability) {
			   if ($config->verbose)
					printf("dumping probability table .. \n");
			   for ($i=1;$i < count($files[0]);$i++){
					$sd="";
					if (($id-1) % 1000==0){
						 if (!file_exists($target_folder.floor(($id-1)/1000)))
							  mkdir($target_folder.floor(($id-1)/1000));
					}
					for ($j=0; $j < count($files); $j++){
						 $temp=trim($files[$j][$i]);
						 $t=split(',',$temp);
						 #$sd=$sd.sprintf("%d\t", $j+1);
						 for($k=1;$k<count($t);$k++){
							  $sd=$sd.sprintf( "%s\t", $t[$k]);
						 }
						 $sd=$sd.sprintf("\n");
					}
					$fd=fopen($target_folder.floor(($id-1)/1000)."/".$id.".txt","w+");
					fprintf($fd, $sd);
					fclose($fd);
					fprintf($finfo,"%d\t%s\t%d\t%d\n", $id, $chr_list[$chr_id],$positions[$i-1], $i);
					$id++;
			   }
		  }
	 }
	 if ($config->verbose)
		  printf("Finish generating genotype folder\n");

}

function gex_init_gene_expression($config, $subject_list){
	 $gene_folder=$config->project_folder.$config->project_alias.'/'."gene_expression/";
	 if (!file_exists($gene_folder)){
		  mkdir($gene_folder);
	 }
	 
	 $gene_meta_idx=gex_load_gene_meta($config);

	 //copy gene file to gene folder
	 copy($config->gene_meta_file, $gene_folder."gene_info");
	 copy($config->gene_expression_file, $gene_folder."gene_expression");
	 list($gene_list, $subject_list, $data_list)=gex_load_gene_expression($config, $subject_list);
	 $i=0;
	 $finfo=fopen($gene_folder."gene_list","w+");
	 foreach ($gene_list as $key => $value){
		  $currfolder=$gene_folder.(int)($i/1000);
		  if ($i%1000==0) {
			   if (!file_exists($currfolder))
					mkdir($currfolder);
		  }
		  $fd=fopen($currfolder."/".$value, 'w+');
		  fprintf($finfo,"%s\t%s\n", $value, (int)($i/1000));
		  foreach ($data_list as $k => $value){
			   fprintf($fd, "%s\n",$value[$key]);
		  }
		  fclose($fd);
		  $i++;
	 }

	 fclose($finfo);
	 return $subject_list;
	 //return object lists
	 
}

/**
   If the function has arguments named subject_list, the result will be the same order in $subject_list,
   Otherwise, it will be loaded by the order in the data file. 
*/
   
function gex_load_gene_expression($config, $specified_subject_list=NULL, $require_data=true){
	 if ($config->verbose){
		  printf("Loading gene expression info\n");
	 }

	 $meta=file($config->gene_expression_file, $flag=FILE_IGNORE_NEW_LINES);
	 $subject_idx=array();
	 $subject_list=array();
	 $data_list=array();
	 $idx=1;
 	 foreach ($meta as $key => $value){
		  if ($key==0) {
			   $gene_name_list=split("\t", $value);
			   unset($gene_name_list[0]);
			   /* build meta index */
			   
			   $gender=file($config->gene_expression_gender, $flag=FILE_IGNORE_NEW_LINES);
			   $gene_gender_list=split("\t", $gender[0]);
		  }
		  else {
			   if ($require_data) $gene_data=split("\t", $value);
			   else $gene_data=split("\t",$value,2);
			   if ($specified_subject_list==NULL || in_array($gene_data[0], $specified_subject_list)){
					print($gene_gender_list[$idx]);
					print($config->gene_gender);
					if  ((!$config->gene_filter_gender || $config->gene_gender==$gene_gender_list[$idx])){
						 array_push($subject_list, $gene_data[0]);
						 if ($require_data){
							  printf("individual ".$gene_data[0]." \n");
							  for ($i=0;$i<count($gene_data);$i++) if ($gene_data[$i]=="") $gene_data[$i]="NA";
							  $subject_idx[$gene_data[0]]=count($subject_list)-1;
							  unset($gene_data[0]);
							  array_push($data_list, $gene_data);
						 }
					}
			   }
		  }
		  $idx++;
	 }
	 if ($config->verbose){
		  printf("There are %d genes loaded \n", count($gene_name_list));
	 }
	 if ($specified_subject_list==NULL)
		  return array(  $gene_name_list, $subject_list, $data_list);
	 else {
		  $new_data_list=array();
		  $i=0;
		  print_r($specified_subject_list);
		  print_r($subject_list);
		  
		  foreach ($specified_subject_list as $value){
			   if (isset($subject_idx[$value])){
					continue;
			   }
			   else {
					die("The function of gex_load_gene_expression detects that the specified_subject_list contains unknown subject name $value");
			   }
		  }
		  return array( $gene_name_list, $subject_list, $data_list);
	 }
			   
}

?>