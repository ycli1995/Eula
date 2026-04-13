
use strict;
use warnings;

use FindBin qw/$Bin/;

use lib "../lib";
use Eula::Utils::ReadWrite;
use Eula::Utils::DEBUG qw/excu_shell/;
use Eula::Utils::PATH qw/rm_paths/;

my @lines;

print "TEST read_lines('my_lines.list')\n";
@lines = read_lines("$Bin/my_lines.list");
print join("\n", @lines) . "\n";
print "\n";

print "TEST read_lines('my_lines.list', {skip_blank => 0})\n";
@lines = read_lines("$Bin/my_lines.list", {skip_blank => 0});
print join("\n", @lines) . "\n";
print "\n";

print "TEST read_lines('my_lines.list', {trim => 0})\n";
@lines = read_lines("$Bin/my_lines.list", {trim => 0});
print join("\n", @lines) . "\n";
print "\n";

print "TEST write_lines('my_write_lines.list')\n";
write_lines("$Bin/my_write_lines.list");
excu_shell("cat $Bin/my_write_lines.list");
print "\n";

print "TEST write_lines('my_write_lines.list', \\\@lines)\n";
write_lines("$Bin/my_write_lines.list", @lines);
excu_shell("cat $Bin/my_write_lines.list");
print "\n";

rm_paths("$Bin/my_write_lines.list");

print "TEST read_shell('my_shell.sh')\n";
@lines = read_shell("$Bin/my_shell.sh");
print join("\n", @lines) . "\n";
print "\n";

print "TEST write_shell('my_shell2.sh')\n";
write_shell("$Bin/my_shell2.sh", \@lines, {skip_patterns => ["Hello"]});
excu_shell("cat $Bin/my_shell2.sh");
print "\n";

rm_paths("$Bin/my_shell2.sh");

