use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Eula::Catalog;

subtest 'dump_catalog with single level' => sub {
    my $root = Eula::Catalog->new();
    $root->subtree('file1.txt', 'File One');
    $root->subtree('file2.txt', 'File Two');
    
    my $output = $root->dump_catalog();

    like($output, qr/^├── file1\.txt/, 'output contains file1.txt');
    like($output, qr/\n└── file2\.txt/, 'output contains file2.txt');
    like($output, qr/File One/, 'output contains File One');
    like($output, qr/File Two/, 'output contains File Two');
};

subtest 'dump_catalog with nested levels' => sub {
    my $root = Eula::Catalog->new();
    
    my $child1 = $root->subtree('child1.txt', 'Child1');
    $child1->subtree('grandchild1.txt', 'Grandchild1');
    $child1->subtree('grandchild2.txt', 'Grandchild2');
    
    my $child2 = $root->subtree('child2.txt', 'Child2');
    
    my $output = $root->dump_catalog();
    
    like($output, qr/^├── child1\.txt/, 'output contains child1.txt');
    like($output, qr/\n└── child2\.txt/, 'output contains child2.txt');
    like($output, qr/\n│   ├── grandchild1\.txt/, 'output contains grandchild1.txt');
    like($output, qr/\n│   └── grandchild2\.txt/, 'output contains grandchild2.txt');
    like($output, qr/Grandchild1/, 'output contains Grandchild1');
    like($output, qr/Grandchild2/, 'output contains Grandchild2');
};

subtest 'dump_catalog with empty catalog' => sub {
    my $root = Eula::Catalog->new();
    
    my $output = $root->dump_catalog();
    
    is($output, '', 'empty catalog returns empty string');
};

subtest 'dump_catalog with single child' => sub {
    my $root = Eula::Catalog->new();
    $root->subtree('only_child.txt', 'Only Child');
    
    my $output = $root->dump_catalog();
    
    like($output, qr/└── only_child\.txt/, 'single child uses ending marker');
    unlike($output, qr/├── /, 'single child does not use branch marker');
};

subtest 'complex tree structure' => sub {
    my $root = Eula::Catalog->new();
    
    my $branch1 = $root->subtree('branch1.txt', 'Branch 1');
    $branch1->subtree('leaf1.txt', 'Leaf 1');
    $branch1->subtree('leaf2.txt', 'Leaf 2');
    
    my $branch2 = $root->subtree('branch2.txt', 'Branch 2');
    my $subbranch = $branch2->subtree('subbranch.txt', 'Subbranch');
    $subbranch->subtree('deepleaf.txt', 'Deep Leaf');
    
    my $branch3 = $root->subtree('branch3.txt', 'Branch 3');
    
    my $output = $root->dump_catalog();
    print $output . "\n";
    
    like($output, qr/^├── branch1\.txt/, 'contains branch1');
    like($output, qr/\n│   ├── leaf1\.txt/, 'contains leaf1');
    like($output, qr/\n│   └── leaf2\.txt/, 'contains leaf2');
    like($output, qr/\n├── branch2\.txt/, 'contains branch2');
    like($output, qr/\n└── branch3\.txt/, 'contains branch3');
    like($output, qr/\n│   └── subbranch\.txt/, 'contains subbranch');
    like($output, qr/\n│       └── deepleaf\.txt/, 'contains deepleaf');
};

done_testing();
