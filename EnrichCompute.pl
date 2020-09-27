#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use PerlIO::gzip;

my $usage=<<"USAGE";

name:    $0
	 
	 It is designed for functional genomic and large-scale genetic
	 studies from which large number of gene lists (e.g. differentially
	 expressed gene sets, co-expressed gene sets, or differential epigenomic
	 modification gene sets etc) are continuously generated. EAGE incorporates
	 information from different public resources and provides an easy way for
	 biologists to make sense out of gene lists.

author:  luhanlin.hi\@qq.com

date:	 2015-08-06

usage:   perl $0

	 -L <str>  gene list, texting all gene symbols/name by one per line;

	 -S <str>  specify the species names, default is bbe;
		   optional: 
                        bbe : Branchiostoma belcheri
			hsa : human 
  			mmu : mouse
                        ..

	 -R <str>  specify the path of Rscript, default is "Rscript";

	 -T <str>  specify the type of analysis, such as GO, pathway, default is GO;
		   optional: GO, KEGG;
 
	 -B <str>  specify the background gene sets using gene sets list text, default is whole genome;

	 -A <str>  specify Multiple Test Adjustment method, such as BH, BY, none, default is none;
		   optional: BH, fdr, BY, holm, none;

	 -O <str>  prefix of output file, default is "./EC"
	
example:  perl $0 -L deg.list -S hsa -T GO -A none 

USAGE

my ($List, $Species, $Type, $Background, $Adjustment, $help, $Out, $Rscript);
GetOptions ( 
	"L=s" => \$List,
	"S:s" => \$Species,
	"T:s" => \$Type,
	"B:s" => \$Background,
	"A:s" => \$Adjustment,
	"O:s" => \$Out,
        "R:s" => \$Rscript,
	"h|?|help" => \$help
);  

die $usage if $help;
die $usage unless $List;

$Species ||= 'bbe';
$Type ||= 'GO';
$Background ||= 'genome';
$Adjustment ||= 'none';
$Out ||= "./EC.$Type";

if(defined $Rscript){
    unless(-e $Rscript){
	print STDERR "[err] $Rscript does't exists!\n";
	print STDERR "please specify the path of Rscript\n";
	exit 0;
    }
}else{
	print STDERR "please specify the path of Rscript\n";
        exit 0;
}

my $gene2annt = "$Bin/../data/$Species.$Type.map.gz";
my $id_map = "$Bin/../data/$Species.$Type.dec.gz";
my %gene_list;
my %ref;
my %genome_list;
my %background_list;
my $background_total;
my $gene_total;
my %id2name;
my $bak_gene_sets = "$Bin/../data/$Species.annotation.gz";

### step 01
open LI, $List or die "[err] can not open the gene sets list : $List\n";
while(<LI>){
	chomp;
	next if /^\s*$/;
	my $gene = (split)[0];
	$gene_list{ $gene } = 1;
}
close LI;
print STDERR "[ok] read-in $List\n";
$gene_total = scalar keys %gene_list;

## step 02
my @term_list;
open IN, "zcat $gene2annt | " or die "[err] can not open the file : $gene2annt\n";
while(<IN>){
	chomp;
	next if /^#/;
	my ($gene, $term) = (split(/\t/, $_))[0,1];
	$genome_list{$gene} = 1;
	$ref{$term}{$gene} = 1;
}
close IN;
print STDERR "[ok] read-in $gene2annt\n";

## step 03
if($Background eq "genome"){
	open IN, "zcat $bak_gene_sets | " or die "[err] can not open the background gene list : $bak_gene_sets\n";
	while(<IN>){
		chomp;
		next if /^#/;
		my @arr = split(/\t/, $_);
		my $gene = $arr[0];
		$background_list{ $gene } = 1;
	}
	close IN;
	print STDERR "[ok] read-in Background gene list : $bak_gene_sets\n";
}
else{
	open IN, $Background or die "[err] can not open the background gene sets: $Background\n";
	while(<IN>){
		chomp;
		next if /^\s*$/;
		next if /^#/;
		my @arr = split;
		$background_list{ $arr[0] } = 1;
	}
	close IN;
	print STDERR "[ok] read-in Background gene list : $Background\n";
}
$background_total = scalar keys %background_list;
if($background_total>1){
	print STDERR "[ok] calculated background gene sets number is : $background_total\n";
}
else{
	print STDERR "[err]  calculated background gene sets number is < 1\n";
	exit 1;
}

### step 04
open OUT, ">$Out.tmp1.xls" or die "[err] can not open the output: $Out.tmp1.xls\n";
print OUT "Term\tCategory\tObserved\tExpected\tFoldChange\tGenelist\n";
foreach my $term (keys %ref){
	my $num_in_C = 0;
	my $num_in_O = 0;
	my @glist = ();
	foreach my $gene (keys %{$ref{$term}}){
		if( $ref{$term}{$gene} == 1 ){
			if( exists $background_list{$gene} ){
				$num_in_C += 1;
			}
			if( exists $gene_list{$gene} ){
				$num_in_O += 1;
				push @glist, $gene;
			}
		}
	}
	if($num_in_O > 0){
		my $expected = $background_total ? sprintf("%.6f", $num_in_C / $background_total * $gene_total) : 0.1;
		my $Foldchange = $expected != 0 ? sprintf("%.6f", $num_in_O / $expected) : 1;
		my $list = join(";", @glist);
		print OUT join("\t", $term, $num_in_C, $num_in_O, $expected, $Foldchange, $list), "\n";
	}
}
print STDERR "[ok] stat against background gene sets is over\n";

## step 05
print STDERR "[run] $Rscript $Bin/../src/Htest.R $Out.tmp1.xls $background_total $gene_total $Adjustment $Out.tmp2.xls\n";
system("$Rscript $Bin/../src/Htest.R $Out.tmp1.xls $background_total $gene_total $Adjustment $Out.tmp2.xls");

## step 06
open MAP, "zcat $id_map | " or die "[err] can not open $id_map\n";
while(<MAP>){
	chomp;
	next if /^#/;
	my @arr = split(/\t/, $_);
	my $id = $arr[0];
	my $name = $arr[1];
	$id2name{$id} = $name;
}
close MAP;

open XLS, "<$Out.tmp2.xls" or die "[err] can not open $Out.tmp2.xls\n";
open OUT, ">$Out.xls" or die "[err] can not open $Out.xls\n";
print OUT "#This table is generate by $0 written by H.L. Lu\n";  # an explaination for output as the final results. 
print OUT "#organism / species: $Species\n";
print OUT "#input gene number: $gene_total\n";
print OUT "#background gene number: $background_total\n";
print OUT "#Col.1 -> Term: GO or KEGG pathway ID\n";
print OUT "#Col.2 -> Category: gene number in this Term according to background\n";
print OUT "#Col.3 -> Observed: gene number in this Term according to input gene list\n";
print OUT "#Col.4 -> Expected: gene number in this Term according to expectation\n"; 
print OUT "#Col.5 -> FoldChange: fold change / eriched factor for this Term\n";
print OUT "#Col.6 -> rawP: the raw p value of Hypergeometric test\n";
print OUT "#Col.7 -> adjP: the adjust p value by $Adjustment\n";
print OUT "#Col.8 -> Name: the description for this Term\n";
print OUT "#Col.9 -> Genelist: gene list in this Term  according to input\n";
while(<XLS>){
	chomp;
	my @arr = split(/\t/,$_);
	if($.==1){
		print OUT join("\t", @arr[0..6], "Name", $arr[7]) . "\n";
	}
	else{
		if(exists $id2name{$arr[0]}){
		    print OUT join("\t", @arr[0..6], $id2name{$arr[0]}, $arr[7]), "\n";
                }
	}
}
close XLS;
close OUT;	
print STDERR "[ok] Term id convert to name finished..\n";
system("rm -vf $Out.tmp1.xls $Out.tmp2.xls");
print STDERR "[ok] $0 is finished\n";

exit 0;
#
