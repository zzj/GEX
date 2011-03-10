<?php
if (!isset($argv[1])){
	 $config_file="config.json";
}
else {
	 $config_file=$argv[1];
}
$config=json_decode(file_get_contents($config_file));
$f=file($config->command_list);
if (isset($argv[2]))
	 system($f[(int)($argv[2])]);
?>