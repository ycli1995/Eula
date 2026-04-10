
use strict;
use warnings;

use lib "../lib";
#use Eula::Utils::DEBUG;
use Eula;

print "TEST time_log()\n";
time_log();
print "\n";

print "TEST excu_shell(\"echo 'Hello world.'\")\n";
excu_shell("echo 'Hello world.'");
print "\n";

print "TEST excu_shell(\"echo 'Hello world. print = 0'\", 1, 0)\n";
excu_shell("echo 'Hello world. print = 0'", 1, 0);
print "\n";

print "TEST excu_shell(\"echo 'Hello world. run = 0'\", 0, 0)\n";
excu_shell("echo 'Hello world. run = 0'", 0, 0);
print "\n";

print "TEST excu_shell(\"mkdir xxx/yyy\")\n";
excu_shell("mkdir xxx/yyy");
print "\n";

print "This should not be printed.\n";

