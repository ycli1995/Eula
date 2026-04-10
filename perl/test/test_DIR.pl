
use strict;
use warnings;

use Data::Dumper;

use lib "../lib";
use Eula::DIR;
use Eula::Utils::PATH;

my $dir = Eula::DIR->new({name => 'test1'});

$dir->add("step1", {ordered => 1});
$dir->add("step2", {ordered => 1});
$dir->add("shell", {ordered => 0});
$dir->add("step3", {ordered => 1});

$dir->get("step1")->add("sub");
$dir->fd("sub")->add("subA");
$dir->fd("sub")->add("subB");

$dir->get("step2")->add("subA", {ordered => 1});
$dir->get("step2")->add("subB", {ordered => 1});

print "TEST DIR object.\n";
print Dumper $dir;
print "\n";

print "TEST DIR->fd()\n";
my $subdir;

$subdir = $dir->fd("subA")->path();
print "BFS: $subdir\n";

$subdir = $dir->fd("subA", {method => "DFS"})->path();
print "DFS: $subdir\n";
print "\n";

rm_paths($dir->path());

