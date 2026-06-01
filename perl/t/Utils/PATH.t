use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use lib "$Bin/../../lib";
use Eula::Utils::PATH;
use File::Temp qw(tempdir);
use Cwd qw(abs_path);

# Test real_path
subtest 'real_path' => sub {
    plan tests => 6;

    my $abs_path = abs_path($Bin);
    is(
        real_path($Bin),
        $abs_path,
        'real_path returns absolute path for existing directory'
    );

    my $nonexistent = '/nonexistent/path/that/does/not/exist';
    is(
        real_path($nonexistent),
        $nonexistent,
        'real_path returns original path for non-existent path'
    );

    my $path_with_trailing = $Bin . '/';
    is(
        real_path($path_with_trailing),
        $abs_path,
        'real_path removes trailing slashes'
    );

    eval { real_path($nonexistent, {must_work => 1}); };
    like(
        $@,
        qr/No such file or directory/,
        'real_path dies when must_work=1 and path does not exist'
    );

    is(
        real_path($Bin, {must_work => 1}),
        $abs_path,
        'real_path works with must_work=1 for existing path'
    );

    my $tmp_dir = tempdir(CLEANUP => 1);
    my $tmp_file = "$tmp_dir/test_file.txt";
    open my $fh, '>', $tmp_file or die;
    close $fh;
    is(
        real_path($tmp_file),
        abs_path($tmp_file),
        'real_path returns absolute path for existing file'
    );
};

# Test rm_paths
subtest 'rm_paths' => sub {
    my $tmp_dir = tempdir(CLEANUP => 1);

    my $tmp_file = "$tmp_dir/test_to_remove.txt";
    open my $fh, '>', $tmp_file or die;
    close $fh;
    ok(-e $tmp_file, 'test file exists before rm_paths');

    rm_paths($tmp_file);
    ok(!-e $tmp_file, 'rm_paths removes file');

    my $tmp_dir2 = tempdir(CLEANUP => 1);
    my @files;
    for my $i (1 .. 3) {
        my $file = "$tmp_dir2/file_$i.txt";
        open my $fh, '>', $file or die;
        close $fh;
        push @files, $file;
    }
    my $count_before = scalar(grep { -e $_ } @files);
    is($count_before, 3, '3 test files exist before rm_paths');

    rm_paths("$tmp_dir2/*");
    my $count_after = scalar(grep { -e $_ } @files);
    is($count_after, 0, 'rm_paths with glob removes all files in dir');

    my $tmp_dir3 = tempdir(CLEANUP => 1);
    ok(-e $tmp_dir3, 'temp directory exists before rm_paths');

    rm_paths($tmp_dir3);
    ok(!-e $tmp_dir3, 'rm_paths removes directory');
};

# Test cp_paths
subtest 'cp_paths' => sub {
    plan tests => 6;

    my $tmp_dir = tempdir(CLEANUP => 1);

    my $src_file = "$tmp_dir/source.txt";
    open my $fh, '>', $src_file or die;
    print $fh "test content\n";
    close $fh;

    my $dst_file = "$tmp_dir/destination.txt";
    cp_paths($src_file, $dst_file);
    ok(-e $dst_file, 'cp_paths copies file');

    open $fh, '<', $dst_file or die;
    my $content = <$fh>;
    close $fh;
    is($content, "test content\n", 'cp_paths preserves file content');

    my $src_dir = "$tmp_dir/source_dir";
    mkdir $src_dir;
    open $fh, '>', "$src_dir/file1.txt" or die;
    print $fh "file1\n";
    close $fh;
    open $fh, '>', "$src_dir/file2.txt" or die;
    print $fh "file2\n";
    close $fh;

    my $dst_dir = "$tmp_dir/destination_dir";
    cp_paths($src_dir, $dst_dir);
    ok(-d $dst_dir, 'cp_paths copies directory');
    ok(-e "$dst_dir/file1.txt", 'cp_paths copies files in directory');
    ok(-e "$dst_dir/file2.txt", 'cp_paths copies all files in directory');

    my $new_dir = "$tmp_dir/new_dir";
    cp_paths($src_file, "$new_dir/new_file.txt");
    ok(-e "$new_dir/new_file.txt", 'cp_paths creates parent directories');
};

# Test empty_dir
subtest 'empty_dir' => sub {
    plan tests => 4;

    my $tmp_dir = tempdir(CLEANUP => 1);

    my $empty_dir = "$tmp_dir/empty";
    mkdir $empty_dir;
    is(empty_dir($empty_dir), 1, 'empty_dir returns 1 for empty directory');

    my $non_empty_dir = "$tmp_dir/non_empty";
    mkdir $non_empty_dir;
    open my $fh, '>', "$non_empty_dir/file.txt" or die;
    close $fh;
    is(empty_dir($non_empty_dir), 0, 'empty_dir returns 0 for non-empty directory');

    my $has_subdir = "$tmp_dir/has_subdir";
    mkdir $has_subdir;
    mkdir "$has_subdir/sub";
    is(empty_dir($has_subdir), 0, 'empty_dir returns 0 for directory with subdir');

    my $nonexistent = "$tmp_dir/nonexistent";
    eval { empty_dir($nonexistent); };
    like($@, qr/Cannot open directory/, 'empty_dir dies for non-existent directory');
};

done_testing();
