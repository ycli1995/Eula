
use strict;
use warnings;

use lib "../lib";
#use Eula::Utils::DEBUG;
use Eula;

sub foo {
    my ($opts) = @_;

    my $input  = require_opt($opts, 'input');
    my $output = require_opt($opts, 'output');

    print "input=$input, output=$output\n";
}

print "TEST require_opt()\n";
foo({input => "a.txt", output => "b.txt", });

foo({output => "a.txt"});
print "\n";

