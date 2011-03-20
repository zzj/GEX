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

function gene_position_cmp($a,$b){
     if (is_numeric($a[0]) && is_numeric($b[0])){
          if ($a[0]==$b[0]){
               return $a[1]-$b[1];
          }
          return $a[0]-$b[0];
     }
     else {
          if (is_numeric($a[0])) return -1;
          else if (is_numeric($b[0])) return 1;
          else  if  ($a[0]==$b[0]){
               return $a[1]-$b[1];
          }
          return $a[0]-$b[0];
     }
}

/* Load gene meta info file, create a dictionary with gene's position information */

/* meta_idx is an array, and key is the probe id, and the value is position information
   0 chr id
   1 start position
   2 end position
*/
function gex_load_gene_meta($config){

	 if ($config->verbose){
		  printf("Loading gene meta info\n");
	 }
	 
	 $meta_idx=array();

     $meta=fopen($config->gene_meta_file, "r");
     if ($meta){
          $key=0;
          while(($value=(fgets($meta,1000000))) !== false) {
               $value=trim($value);
               $key++;
               if ($key==0) continue;
               $content=getCSVvalues($value);
               
               $meta_idx[$content[0]]=array($content[$config->gene_meta_chr],$content[$config->gene_meta_start_position],$content[$config->gene_meta_end_position]);
          }
          if (!feof($meta)) {
               debug_print_backtrace();
               die("Unexpected fgets fail\n");
          }
     }
     else {
          debug_print_backtrace();
          die("Can not open file ".$config->gene_meta_file);
     }
     $result=uasort($meta_idx,'gene_position_cmp');
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
     /* find joint strain set from both genotype data and gene expression data*/

	 $subject_list=gex_init_joint_strains($config);
	 $subject_list=array_values($subject_list);
	 printf("There are total ".count($subject_list)." strains founded in both datasets\n");

     /* init gene expression data */
	 $subject_list=gex_init_gene_expression($config,$subject_list);

     /* init genotype data */
	 gex_init_genotype($config,$subject_list);

     /*The strain list must be dumped here, because the subject_list is right at this time*/
	 gex_dump_strain_list($config, $resultfolder, $subject_list);
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
	 $gid=1;
	 if ($config->verbose)
		  printf("Start generating genotype folder\n");
	 
	 for($chr_id=0;$chr_id<count($chr_list);$chr_id++){
		  if ($config->verbose)
			   printf("Working on Chromosome ".$chr_list[$chr_id]."\n");

		  //create a folder for probability data of chromosome genotype
          $current_folder=$folder.$chr_list[$chr_id]."/";
		  if (!file_exists($current_folder)){
			   mkdir($current_folder);
		  }

          $files=array();
		  $genotypes=array();
          for ($i=0;$i < count($strain_name); $i++){
			   if ($config->have_probability)
                    $files[$i]=file($current_folder.$strain_name[$i].".txt.prob",$flag=FILE_IGNORE_NEW_LINES);
			   $genotypes[$i]=fopen($current_folder.$strain_name[$i].".txt","r");
               if (!$genotypes[$i]){
                    debug_print_backtrace();
                    die( " can not open file ". $strain_name[$i].".txt");
               }
		  }
          if ($config->verbose){
               printf("Loading position information ...\n");
          }
		  $markers=file($current_folder."Extra.txt");
          $positions=array();
          foreach ($markers as $value){
               $info=split(',',$value);
               array_push($positions, $info[2]);
          }
		  if ($config->verbose)
			   printf("dumping genotype data .. \n");
		  $fg=fopen($target_folder.$chr_list[$chr_id].".genotype","w+");
          $id=1;
          $pos_idx=0;
          while(!feof($genotypes[0])){
               $out="";
               $ignore=false;
               $same=true;
               $last="";
               for ($j=0; $j<count($genotypes);$j++){
                    if (fscanf($genotypes[$j],"%s",$temp)==0) break;
					$out.=$temp."\t";
                    if ($last==""){
                         $last=$temp;
                    }
                    else {
                         if ($last!=$temp) $same=false;
                    }
                    # must use continue, because other genotype file handles does not reach the same position
                    if ($temp=="N") {$ignore=true;}
			   }
               if ($same) $ignore=true;
               
               $pos_idx++;
               if (feof($genotypes[0])) break;
               if ( $ignore) {continue;}
               fprintf($fg,$out."\n");
               if (!$config->have_probability){
                    //printf("%d\t%d\n", $positions[$id-1], $id);
                    fprintf($finfo,"%d\t%s\t%d\t%d\n", $gid, $chr_list[$chr_id],$positions[$pos_idx-1], $id);
                    $id++; $gid++;
               }
          }
		  fclose($fg);
          for ($i=0;$i<count($genotypes);$i++){
               fclose($genotypes[$i]);
          }
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

	 //copy gene file to gene folder
	 copy($config->gene_meta_file, $gene_folder."gene_info");
	 copy($config->gene_expression_file, $gene_folder."gene_expression");
	 list($gene_list, $subject_list, $data_list)=gex_load_gene_expression($config, $subject_list);
     $sorted_meta_idx=gex_load_gene_meta($config);
	 $i=0;
	 $finfo=fopen($gene_folder."gene_list","w+");
     $gene_list=array_flip($gene_list);
	 foreach ($sorted_meta_idx as $key => $value){
          if (!isset($gene_list[$key])) continue;
		  $currfolder=$gene_folder.(int)($i/1000);
		  if ($i%1000==0) {
			   if (!file_exists($currfolder))
					mkdir($currfolder);
		  }
          $fd=fopen($currfolder."/".$key, 'w+');
		  fprintf($finfo,"%s\t%s\t%s\t%s\t%s\n", $key, (int)($i/1000), $value[0], $value[1],$value[2]);
		  foreach ($data_list as $k => $value2){
               if (!isset($value2[$gene_list[$key]])){
                    print($k);
                    print($key);
                    print($gene_list[$key]);
               }
			   fprintf($fd, "%s\n",$value2[$gene_list[$key]]);
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

	 $subject_idx=array();
	 $subject_list=array();
	 $data_list=array();
	 $idx=1;
	 //$meta=file($config->gene_expression_file, $flag=FILE_IGNORE_NEW_LINES);
 	 //foreach ($meta as $key => $value){
     $gene_expression=fopen($config->gene_expression_file, "r");
     $key=0;
     if ($gene_expression){
          while(($value=(fgets($gene_expression,10000000))) !== false) {
               if ($key==0) {
                    $gene_name_list=split("\t", $value);
                    unset($gene_name_list[0]);
                    /* build gene_expression index */
                    if ($config->gene_filter_gender) {
                         $gender=file($config->gene_expression_gender, $flag=FILE_IGNORE_NEW_LINES);
                         $gene_gender_list=split("\t", $gender[0]);
                    }
               }
               else {
                    if ($require_data) $gene_data=split("\t", $value);
                    else $gene_data=split("\t",$value,2);
                    // Check whether it is specified in list, otherwise, should be ignored
                    if ($specified_subject_list==NULL || in_array($gene_data[0], $specified_subject_list)){
                         // check whether match the gender information
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
               $idx++;$key++;
          }
          if (!feof($gene_expression) ) {
               debug_print_backtrace();
               die("Unexpected fgets fail\n");
          }

	 }
     else {
          debug_print_backtrace();
          die("Can not open file ".$config->gene_expression_file);
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