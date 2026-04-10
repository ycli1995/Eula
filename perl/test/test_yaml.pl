
use strict;
use warnings;

use FindBin qw/$Bin/;
use YAML;

use lib "../lib";
use Eula::Utils::ReadWrite;
use Eula::Utils::DEBUG qw/excu_shell/;
use Eula::Utils::PATH qw/rm_paths/;

print "TEST read_yaml('my_yaml.yaml')\n";
my $conf = read_yaml("$Bin/my_yaml.yaml");
print YAML::Dump($conf) . "\n";
print "\n";

