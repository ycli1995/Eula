
use strict;
use warnings;

use Data::Dumper;

use lib "../lib";
use Eula::Utils::DEBUG;

print "TEST deref()\n";
print deref("A word") . "\n";
print join(" ", deref(['An', 'array'])) . "\n";
print deref({A => 1, B => 2});
print "\n";

