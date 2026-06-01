use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Eula::DIR;
use File::Temp qw(tempdir);
use Cwd qw(getcwd);

# Test new and init
subtest 'new and init' => sub {
    my $tmp_dir = tempdir(CLEANUP => 1);
    
    my $dir = Eula::DIR->new({wkdir => $tmp_dir, name => 'test'});
    ok(defined $dir, 'DIR object created');
    is($dir->{_name}, 'test', 'name is set correctly');
    like($dir->{_path}, qr/test$/, 'path ends with name');
    is($dir->path(), "$tmp_dir/test", 'path is set correctly');
    is($dir->parent(), undef, 'parent is undef for new DIR');
};

# Test add
subtest 'add' => sub {
    my $tmp_dir = tempdir(CLEANUP => 1);
    
    my $dir = Eula::DIR->new({wkdir => $tmp_dir, name => 'root', mkdir => 0});
    $dir->mkdir();

    my $children_names = ["child1", "child2"];

    $dir->add($children_names, {mkdir => 1});
    is(scalar @{$dir->{_child_names}}, @{$children_names}, 'two children added');
    foreach (@$children_names) {
        ok(exists $dir->{_child}{$_}, "$_ exists");
        ok(-d "$tmp_dir/root/$_", "$_ directory created");
    }

    push @$children_names, 'child3';
    $dir->add('child3', {mkdir => 1});
    is(scalar @{$dir->{_child_names}}, @{$children_names}, 'third child added');
    
    $dir->add($children_names->[0], {mkdir => 0});
    is(scalar @{$dir->{_child_names}}, @{$children_names}, 'duplicate child not added');
};  

# Test add with ordered
subtest 'add with ordered' => sub {
    my $tmp_dir = tempdir(CLEANUP => 1);
    
    my $dir = Eula::DIR->new({wkdir => $tmp_dir, name => 'root', mkdir => 0});
    $dir->mkdir();
    
    $dir->add(['a', 'b', 'c'], {mkdir => 1, ordered => 1});
    
    my $child_a = $dir->get('a');
    is($child_a->{_index}, '1.', 'first child has index 1.');
    like($child_a->{_path}, qr/1\.a$/, 'first child path has index prefix');
    
    my $child_b = $dir->get('b');
    is($child_b->{_index}, '2.', 'second child has index 2.');
    
    my $child_c = $dir->get('c');
    is($child_c->{_index}, '3.', 'third child has index 3.');
};

# Test get
subtest 'get' => sub {
    my $tmp_dir = tempdir(CLEANUP => 1);
    
    my $dir = Eula::DIR->new({wkdir => $tmp_dir, name => 'root', mkdir => 0});
    $dir->mkdir();
    $dir->add(['x', 'y'], {mkdir => 0});
    
    my $child = $dir->get('x');
    ok(defined $child, 'get returns child');
    is($child->{_name}, 'x', 'child has correct name');
    
    my $nonexistent = $dir->get('nonexistent');
    is($nonexistent, undef, 'get returns undef for nonexistent child');
};

# Test path
subtest 'path' => sub {
    my $tmp_dir = tempdir(CLEANUP => 1);
    
    my $dir = Eula::DIR->new({wkdir => $tmp_dir, name => 'root', mkdir => 0});
    $dir->mkdir();
    $dir->add('subdir', {mkdir => 0});
    
    my $root_path = $dir->path();
    like($root_path, qr/root$/, 'path returns root path');
    
    my $subdir_path = $dir->path('subdir');
    like($subdir_path, qr/subdir$/, 'path returns child path');
};

# Test parent
subtest 'parent' => sub {
    my $tmp_dir = tempdir(CLEANUP => 1);
    
    my $dir = Eula::DIR->new({wkdir => $tmp_dir, name => 'root', mkdir => 0});
    $dir->mkdir();
    $dir->add('child', {mkdir => 0});
    
    my $child = $dir->get('child');
    my $parent = $child->parent();
    ok(defined $parent, 'child has parent');
    is($parent->{_name}, 'root', 'parent has correct name');
    is($parent->{_path}, "$tmp_dir/root", 'parent has correct path');
};

# Test fd with BFS
subtest 'fd with BFS' => sub {
    local $Eula::DIR::TraversalMethod = 'BFS';
    
    my $tmp_dir = tempdir(CLEANUP => 1);
    
    my $dir = Eula::DIR->new({wkdir => $tmp_dir, name => 'root', mkdir => 0});
    $dir->mkdir();
    $dir->add('level1', {mkdir => 0});
    
    my $level1 = $dir->get('level1');
    $level1->add('level2', {mkdir => 0});
    
    my $self = $dir->fd('.');
    is($self->{_name}, 'root', 'fd . returns self');
    
    my $found_level1 = $dir->fd('level1');
    is($found_level1->{_name}, 'level1', 'fd finds level1');
    
    my $found_level2 = $dir->fd('level2');
    is($found_level2->{_name}, 'level2', 'fd finds level2 in BFS');
};

# Test fd with DFS
subtest 'fd with DFS' => sub {
    local $Eula::DIR::TraversalMethod = 'DFS';
    
    my $tmp_dir = tempdir(CLEANUP => 1);
    
    my $dir = Eula::DIR->new({wkdir => $tmp_dir, name => 'root', mkdir => 0});
    $dir->mkdir();
    $dir->add('level1', {mkdir => 0});
    
    my $level1 = $dir->get('level1');
    $level1->add('level2', {mkdir => 0});
    
    my $found_level2 = $dir->fd('level2');
    is($found_level2->{_name}, 'level2', 'fd finds level2 in DFS');
};

# Test mkdir recursive
subtest 'mkdir recursive' => sub {
    my $tmp_dir = tempdir(CLEANUP => 1);
    
    my $dir = Eula::DIR->new({wkdir => $tmp_dir, name => 'root', mkdir => 0});
    $dir->mkdir();
    $dir->add('child1', {mkdir => 0});
    
    my $child1 = $dir->get('child1');
    $child1->add('grandchild', {mkdir => 0});
    
    $dir->mkdir({recursive => 1});
    
    ok(-d "$tmp_dir/root", 'root directory created');
    ok(-d "$tmp_dir/root/child1", 'child1 directory created');
    ok(-d "$tmp_dir/root/child1/grandchild", 'grandchild directory created');
};

done_testing();
