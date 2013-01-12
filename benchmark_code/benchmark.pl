use Benchmark;
use warnings;
use strict;
use POSIX;

my @files = (
"/norstore/user/mironov/workspace/omcl/omcl2.0.3/sample/goodProteins.blast" ,
 	  "/norstore/user/mironov/workspace/omcl/omcl2.0.3/test/goodProteins.blast"
	     ,
 	  "/norstore/user/mironov/workspace/omcl/omcl2.0.3/2-set/goodProteins.blast",
 	  "/norstore/user/mironov/workspace/omcl/omcl2.0.3/3-set/goodProteins.blast"
,
 	  "/norstore/user/mironov/workspace/omcl/omcl2.0.3/4-set/goodProteins.blast",
 	  "/norstore/user/mironov/workspace/omcl/omcl2.0.3/5-set/goodProteins.blast",
 	  "/norstore/user/mironov/workspace/omcl/omcl2.0.3/6-set/goodProteins.blast",
 	  "/norstore/user/mironov/workspace/omcl/omcl2.0.3/7-set/goodProteins.blast",
 	  "/norstore/user/mironov/workspace/omcl/omcl2.0.3/sprot/goodProteins.blast -mp 100",
 	  "/norstore/user/mironov/workspace/omcl/omcl2.0.3/gexo/goodProteins.blast"
	  );
for my $a (0 .. $#files)
{    
    open(my $result, ">", "result_benchmark.txt") or die $!;
    my $cmd_basis = "./orthaGogue -i " . $files[$a];
    printf("-\tRuns file %s\n", $cmd_basis);
    printf($result "-\tRuns file %s\n", $cmd_basis);
    for my $cpu (1 .. 7) {
	my $cmd = $cmd_basis . " -c " . $cpu;
#    printf("Executes: %s\n", $cmd);
	my $t0 = Benchmark->new;
	system($cmd);
	my $t1 = Benchmark->new;
	my $td = timediff($t1, $t0);
	my $wc=  `wc -l all.abc`;
	chomp($wc);
	print "\t", $cpu , "\t ", $wc, " lines: ", timestr($td), "\n";
	printf($result "\t%d %s lines: %s\n", $cpu , $wc, timestr($td));
    }
    close($result) or die $!;
}


    # ... your code here ...
