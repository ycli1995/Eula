
use strict;
use warnings;

use Data::Dumper;
use FindBin qw/$Bin/;

use lib "../lib";
use Eula::Utils::Options;

print "TEST real_path(\$opts, \$defaults)\n";
my $opts = { A => 'AA', B => 2, C => undef, D => 'dd' };
my $defaults = { A => 'a', B => 'b', C => 'c', D => 'd'};
get_defaults($opts, $defaults);
print Dumper $opts;

print "TEST norm_option_array(\$opt_arr)\n";
my $opt_arr = [ 
	'geneAA', 
	'geneBB,geneCC', 
	"$Bin/my_genes1.list", 
	undef,
	"$Bin/my_genes1.list,$Bin/my_genes2.list" 
	];
$opt_arr = norm_option_array($opt_arr);
print Dumper $opt_arr;

print "TEST get_options_from_args(\@args)\n";
my @args;
my %opts_out;
@args = qw/argA argB -a AA -b -c --num_arg 10 --bool_arg/;
%opts_out = get_options_from_args(@args);
print Dumper \%opts_out;
