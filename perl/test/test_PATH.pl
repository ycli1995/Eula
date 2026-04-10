
use strict;
use warnings;

use File::Path qw/make_path/;

use lib "../lib";
use Eula::Utils::PATH;

print "TEST real_path('../lib')\n";
print real_path('../lib') . "\n";
print "\n";

print "TEST real_path('xxx')\n";
print real_path('xxx/yyy') . "\n";
print "\n";

print "TEST make_path\n";
make_path("test_make_path/xxx/yyy");
print real_path('test_make_path/xxx/yyy', {must_work => 1}) . "\n";
print "\n";

print "TEST real_path('xxx/yyy', {must_work => 1})\n";
make_path("xxx");
rm_paths("test_make_path", "xxx");
print real_path('xxx/yyy', {must_work => 1}) . "\n";
print "\n";

