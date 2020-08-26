<?php
	ini_set('display_errors', '1');

	$win=200; $step=10;	
	if (isset($_POST['json'])) {
		$json=true;
		$data=json_decode($_POST['json']);
		$ref=$data->{'ref'};
		$region=$data->{'region'};
	} elseif (isset($_REQUEST['ref'])) {
		$json=false;
		$ref=$_REQUEST['ref'];
		$region=$_REQUEST['region'];
	} else {
		system("Hello, world!");
		exit(-1);
	}

	$method="gRNA";
	$samtools="bin/samtools";
	$se="bin/se";
	
	if ($method=='gRNA') {
		$script="$samtools view bam/$ref.bam $region";
		$cols=array('chr', 'start', 'end', 'seq', 'strand', 'efficiency');
		$header=array('chr'=>'Chr', 'start'=>'Start', 'end'=>'End', 'seq'=>'Seq', 'strand'=>'Strand', 'efficiency'=>'Efficiency');
		$dotdot=array('chr'=>'...', 'start'=>'', 'end'=>'', 'seq'=>'', 'strand'=>'', 'efficiency'=>'');
		$parser=function($line) {
			$data = preg_split("/[\s:]+/", $line); $strand='+'; if ($data[2]==16) $strand='-'; $ef=$data[14];
			return array('chr'=>$data[3], 'start'=>$data[4], 'end'=>$data[4]+20, 'seq'=>$data[1], 'strand'=>$strand, 'efficiency'=>$ef);
		};
		$log="buf/$method.$ref.$region.log";
		$output="buf/$method.$ref.$region.txt";
	}
	
	if (!file_exists($output)) {
		system("$script 2> $log > $output");
	}
	
	if ($json) {
		$lines = file($output);
		$result=array(); $result[]=$header; $n=0; $more='';
		foreach ($lines as $line) {
			if ($n>=100) {
				$result[]=$dotdot;
				$more='+';
				break;
			}
			$result[]=$parser($line);
			$n++;
		}
		
		$count=count($result);
		$msg="Found $n$more candidates at $region from $ref.";
		$ret=array('good'=>true, 'message'=>"$msg", 'script'=>$script, 'cols'=>$cols, 'table'=>$result);
		print json_encode($ret);
	} else {
		header("Content-Type: text/plain");
		include($output);
	}
?>
