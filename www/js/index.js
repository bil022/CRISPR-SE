(function () {
  'use strict';
  function revComp(seq) {
     var s=seq.toUpperCase();
     var o=s.split();for(var i=0;i<s.length;i++)o[i]=s.charAt(i).replace(/[ACGT]/g,function(c){return {'A':'T','C':'G','G':'C','T':'A'}[c];});
     return o.reverse().join('');
  }

  function countGuideRNA(seq) {
    seq=seq.replace(/\s+/gm,'');
    var match=null;
    if (!seq.match(/^[ACTGN]+$/i)) {
      return { msg: "Error: input does not match /^[ACTGN]+$/i" };
    }
    var fwd=0, rev=0;
    var nggRE=/[ACTG]GG+/ig;
    while ((match=nggRE.exec(seq))!=null) {
      var pam=match[0];
      for (var off=0; off<=pam.length-3; off++) {
        var idx=match.index+off-20;
        if (idx<0) continue;
        var gRNA=seq.substr(idx, 20);
        if (gRNA.length<20) continue;
        fwd++;
      }
    }

    var ccnRE=/CC+[ACTG]/ig;
    while ((match=ccnRE.exec(seq))!=null) {
      var pam=match[0];
      for (var off=0; off<=pam.length-3; off++) {
        var idx=match.index+off+3;
        var gRNA=seq.substr(idx, 20);
        if (gRNA.length<20) continue;
        rev++;
      }
    }

    return {fwd: fwd, rev: rev, total: fwd+rev, msg: ""};
  } 

  function scanInput(reg) {
    console.log('Input: '+reg+" "+reg.length);
    var regionRE=/^(\S+):([\d,]+)-([\d,]+)$/;
    var match=regionRE.exec(reg);
    if (match!=null) {
      var chr=match[1];
      var s=match[2].replace(/,/g, "");
      var e=match[3].replace(/,/g, "");
      s=parseInt(s); e=parseInt(e);
      if (s>e) {
        var t=s; s=e; e=t;
      }
      return {type: 'region', region: chr+":"+s+"-"+e, len: e-s };
    }

    var lst=[];
    var nggRE=/[ACTG]GG+/ig;
    while ((match=nggRE.exec(reg))!=null) {
      var pam=match[0];
      for (var off=0; off<=pam.length-3; off++) {
        var idx=match.index+off-20;
        if (idx<0) continue;
        var seq=reg.substr(idx, 20);
        if (seq.length<20) continue;
        lst.push({index: idx, gRNA: seq, strand: '+'});
      }
    }

    var ccnRE=/CC+[ACTG]/ig;
    var len=reg.length;
    while ((match=ccnRE.exec(reg))!=null) {
      var pam=match[0];
      for (var off=0; off<=pam.length-3; off++) {
        var idx=match.index+off+3;
        var seq=reg.substr(idx, 20);
        if (seq.length<20) continue;
        seq=revComp(seq);
        lst.push({index: idx, gRNA: seq, strand: '-'});
      }
    }

    if (!lst.length) {
      return "";
    } 
 
    lst=lst.sort(function(a,b){return a.index-b.index;});
    lst.forEach(function(a){ console.log(a); })
    return {type: 'gRNA', lst: lst};
  }

  angular
      .module('MyApp',['ngMaterial', 'ngMessages', 'material.svgAssetsCache'])
      .controller('AppCtrl', AppCtrl)
      .filter('regionInput', function() {
          return function(input, ref) {
 	    var parsed=scanInput(input);
	    console.log(parsed); 
            if (parsed.type=="region") {
               return parsed.region+" "+parsed.len+" bp";
            }
            if (parsed.type=="gRNA") {
               return "Found "+parsed.lst.length+" potential gRNA(s)";
            }
	    return "N/A";
          }
        })
      .filter('scanSequence', function() {
        return function(input) {
          var counted=countGuideRNA(input);
          if (counted.msg.length>0) {
            return counted.msg;
          }
          var msg = "No gRNA found";
          if (counted.total==0) {
            return msg;
          }
          if (counted.fwd>0 && counted.rev>0) {
            msg=counted.total + " gRNA(s)";
            msg+=", "+counted.fwd+ " on forward strand";
            msg+=", "+counted.rev+ " on reverse strand.";
          } else if (counted.fwd>0 ) {
            msg=counted.total + " gRNA(s) on forward strand.";
          } else {
            msg=counted.total + " gRNA(s) on reverse strand.";
          }
          return msg;
        }
      });

  function AppCtrl ($scope, $http, $log) {
    var tabs = [
          { title: 'HOME', html: 'include/home.html' },
          { title: 'View', html: 'include/view.html' },
          { title: 'BLAT', html: 'include/tools.html' },
          { title: 'UCSC', html: 'include/ucsc.html' },
          { title: 'Scan', html: 'include/scan.html' },
          { title: 'Downloads', html: 'include/downloads.html' },
          { title: 'Software', html: 'include/software.html' },
          { title: 'About', html: 'include/about.html' }
        ],
        selected = null,
        previous = null;
    $scope.data = "hello"; 
    $scope.tabs = tabs;
    $scope.selectedIndex = 0;
    $scope.$watch('selectedIndex', function(current, old){
      previous = selected;
      selected = tabs[current];
    });
    $scope.constant={
      refs: [ {id:'mm9'}, {id:'mm10'}, {id:'hg18'}, {id:'hg19'}, {id:'hg38'} ],
      hgChrs: [ { id: 'chr1'}, {id: 'chr2'}, {id: 'chr3'}, {id: 'chr4'}, {id: 'chr5'}, {id: 'chr6'}, {id: 'chr7'}, {id: 'chr8'}, {id: 'chr9'}, {id: 'chr10'}, {id: 'chr11'}, {id: 'chr12'}, {id: 'chr13'}, {id: 'chr14'}, {id: 'chr15'}, {id: 'chr16'}, {id: 'chr17'}, {id: 'chr18'}, {id: 'chr19'}, {id: 'chr20'}, {id: 'chr21'}, {id: 'chr22'}, {id: 'chrX'}, {id: 'chrY'} ],
      mmChrs: [ { id: 'chr1'}, {id: 'chr2'}, {id: 'chr3'}, {id: 'chr4'}, {id: 'chr5'}, {id: 'chr6'}, {id: 'chr7'}, {id: 'chr8'}, {id: 'chr9'}, {id: 'chr10'}, {id: 'chr11'}, {id: 'chr12'}, {id: 'chr13'}, {id: 'chr14'}, {id: 'chr15'}, {id: 'chr16'}, {id: 'chr17'}, {id: 'chr18'}, {id: 'chr19'}, {id: 'chrX'}, {id: 'chrY'} ],
      method: [ { id: 'sgRNA'}, { id: 'Tiling'}]
   };
   $scope.downloads=[ 
     {id:'mm9', info: 'Mus musculus (mouse) MGSCv37'},
     {id:'mm10', info: 'Mus musculus (mouse) GRCm38'},
     {id:'hg18', info: 'Homo sapiens (human) NCBI36'},
     {id:'hg19', info: 'Homo sapiens (human) GRCh37'},
     {id:'hg38', info: 'Homo sapiens (human) GRCh38'},
     {id:'aaegL5', info:'Aedes aegypti (yellow fever mosquito) aaegL5'},
     {id:'anoGam3', info:'Anopheles gambiae (mosquito) anoGam3'},
     {id:'danRer11', info:'Danio rerio (zebrafish) danRer11'},
     {id:'dm6', info:'Drosophila melanogaster (fruit fly) dm6'},
     {id:'dmel-r6.34', info:'Drosophila melanogaster (fly) dmel-r6.34'},
     {id:'ecoli', info:'E. coli K-12'},
     {id:'rn6', info:'Rattus norvegicus (Norway rat) rn6'},
     {id:'sacCer3', info:'Saccharomyces cerevisiae (yeast) S288C'},
     {id:'tair10', info:'Arabidopsis thaliana (Thale cress) tair10'},
     {id:'xenLae2', info:'Xenopus laevis (African clawed frog) xenLae2'},
    ];
  
  $scope.fastaUrl={
    'mm9': 'http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz',
    'mm10': 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz',
    'hg18': 'http://hgdownload.cse.ucsc.edu/goldenpath/hg18/chromosomes/',
    'hg19': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz',
    'hg38': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz',
    'aaegL5': 'ftp://ftp.ensemblgenomes.org/pub/metazoa/release-47/fasta/aedes_aegypti_lvpagwg/dna/Aedes_aegypti_lvpagwg.AaegL5.dna.toplevel.fa.gz',
    'anoGam3': 'https://hgdownload.soe.ucsc.edu/goldenPath/anoGam3/bigZips/anoGam3.fa.gz',
    'danRer11': 'https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz',
    'dm6': 'https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz',
    'dmel-r6.34': 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-chromosome-r6.34.fasta.gz',
    'ecoli': 'ftp://ftp.ensemblgenomes.org/pub/bacteria/release-48/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz',
    'rn6': 'https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz',
    'sacCer3': 'https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz',
    'tair10': 'https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas',
    'xenLae2': 'https://hgdownload.soe.ucsc.edu/goldenPath/xenLae2/bigZips/xenLae2.fa.gz' 
  }


   $scope.selectedRef='hg19';
    $scope.data={
      ref: 'hg19',
      region: 'chr6:31132114-31238451',
      scan: 'ACCGAACTTTAAAATCTGTGTGGC',
   //   chr: 'chr6',
   //   s: 31132114,
   //   e: 31238451,
   //   method: 'sgRNA',
   //   valid: 1
   };
   //$scope.help=0;
   //$scope.hgChr='chr1';
   //$scope.mmChr='chr1';
   $scope.php={};

   $scope.blat = function() {
     $scope.selectedIndex=1;  
   }
   $scope.scan = function() {
     $scope.selectedIndex=2;  
   }
   $scope.showIdx = function(idx) {
     $scope.selectedIndex=idx;  
   }
   $scope.region = function(seq) {
     $scope.data.region=seq;
   }
   $scope.scanDetail = function(seq) {
    seq=seq.replace(/\s+/gm,'');

     $scope.scanner={};
     $scope.scanner.input=seq;
     $scope.scanner.cols=['offset', 'gRNA', 'strand'];

    if (!seq.match(/^[ACTGN]+$/i)) {
      return;
    }

    $scope.scanner.revcmp=revComp(seq);

     var match=null;
     var lst=[];
     var nggRE=/[ACTG]GG+/ig;
     while ((match=nggRE.exec(seq))!=null) {
       var pam=match[0];
       for (var off=0; off<=pam.length-3; off++) {
         var idx=match.index+off-20;
         if (idx<0) continue;
         var s=seq.substr(idx, 20);
         if (s.length<20) continue;
         lst.push({offset: idx, gRNA: s, strand: '+'});
       }
     }
 
     var ccnRE=/CC+[ACTG]/ig;
     //var len=seq.length;
     while ((match=ccnRE.exec(seq))!=null) {
       var pam=match[0];
       for (var off=0; off<=pam.length-3; off++) {
         var idx=match.index+off+3;
         var s=seq.substr(idx, 20);
         if (s.length<20) continue;
         lst.push({offset: idx, gRNA: revComp(s), strand: '-'});
       }
     }
 
     lst=lst.sort(function(a,b){return a.offset-b.offset;});
     $scope.scanner.table=lst;
   }

   $scope.valid = function(seq) {
     var re=/^(\S+):([\d,]+)-([\d,]+)$/;
     if (re.exec(seq)) return true;
     var gRE=/([ACTG]{20}).GG/i;
     if (gRE.exec(seq)) return true;
     var rRE=/CC.([ACTG]{20})/i;
     if (rRE.exec(seq)) return true;
     return false;
   }
   $scope.go4table = function() {
      $scope.help++;
      var strJson=JSON.stringify($scope.data);
      var postData = 'json='+strJson;
      // $log.debug(postData);
      console.log(strJson);
      $scope.php.table="";
      $scope.php.message='Checking '+strJson;
      $scope.php.err="";
      $http({
        method : 'POST',
        url    : 'table.php',
        data: postData,
        headers : {'Content-Type': 'application/x-www-form-urlencoded'} 
      }).success(function(data) {
        $log.debug(data);
        if (data.good) {
          $scope.php.cols=data.cols;
          $scope.php.table=data.table;
          $scope.php.message=data.message;
          $scope.php.script=data.script;
        } else {
          $scope.php.message="";
          $scope.php.err=data;
        }
        //$scope.debug=data;
        //$scope.php.table=data;
        //$scope.php.message=data;
     }).error(function(error) {
        $log.debug(error);
        //$scope.php.message=error;
     });
   }
  }
})();

/**
Copyright 2016 Google Inc. All Rights Reserved. 
Use of this source code is governed by an MIT-style license that can be in foundin the LICENSE file at http://material.angularjs.org/license.
**/
