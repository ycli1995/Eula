use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use lib "$Bin/../../lib";
use Eula::Utils::ReadWrite;
use File::Temp qw(tempdir);
use YAML;

# Test write_lines and read_lines
subtest 'write_lines and read_lines' => sub {
    my $tmp_dir = tempdir(CLEANUP => 1);
    my $tmp_file = "$tmp_dir/test_lines.txt";

    my $lines = ['line1', 'line2', 'line3'];
    write_lines($tmp_file, $lines);
    ok(-e $tmp_file, 'write_lines creates file');

    my $read_lines = read_lines($tmp_file);
    is_deeply($read_lines, $lines, 'read_lines reads lines correctly');

    my $lines_with_blanks = ['line1', '', 'line2', ' ', 'line3'];
    write_lines($tmp_file, $lines_with_blanks);
    $read_lines = read_lines($tmp_file);
    is_deeply($read_lines, $lines, 'read_lines skips blanks');

    $read_lines = read_lines($tmp_file, {skip_blank => 0, trim => 0});
    is_deeply(
        $read_lines,
        $lines_with_blanks,
        'read_lines with skip_blank=0 and trim=0 includes blanks'
    );
    
    my $lines_with_spaces = ['  line1  ', ' line2 '];
    write_lines($tmp_file, $lines_with_spaces);
    $read_lines = read_lines($tmp_file);
    is_deeply($read_lines, ['line1', 'line2'], 'read_lines trims whitespace');

    $read_lines = read_lines($tmp_file, {trim => 0});
    is_deeply(
        $read_lines,
        $lines_with_spaces,
        'read_lines with trim=0 preserves whitespace'
    );
};

# Test write_yaml and read_yaml
subtest 'write_yaml and read_yaml' => sub {
    my $tmp_dir = tempdir(CLEANUP => 1);
    my $tmp_file = "$tmp_dir/test.yaml";

    my $config = {
        key1 => 'value1',
        key2 => ['a', 'b', 'c'],
        key3 => {a => 1, b => 2, c => 3},
    };
    write_yaml($tmp_file, $config);
    ok(-e $tmp_file, 'write_yaml creates file');

    my $read_config = read_yaml($tmp_file);
    is($read_config->{key1}, $config->{key1}, 'read_yaml reads key1 correctly');
    is_deeply($read_config->{key2}, $config->{key2}, 'read_yaml reads array correctly');
    foreach (keys %{$config->{key3}}) {
        is_deeply($read_config->{key3}{$_}, $config->{key3}{$_}, 'read_yaml reads hash correctly');
    }

    my $simple_hash = {name => 'test', value => 42};
    write_yaml($tmp_file, $simple_hash);
    $read_config = read_yaml($tmp_file);
    is($read_config->{name}, $simple_hash->{name}, 'read_yaml reads name correctly');
    is($read_config->{value}, $simple_hash->{value}, 'read_yaml reads value correctly');

    my $array_config = [1, 2, 3, 4];
    write_yaml($tmp_file, $array_config);
    $read_config = read_yaml($tmp_file);
    is_deeply($read_config, $array_config, 'read_yaml reads array yaml');
};

# Test read_table
subtest 'read_table' => sub {
    my $tmp_dir = tempdir(CLEANUP => 1);
    my $tmp_file = "$tmp_dir/test_table.txt";

    my $arr = [['a', 'b', 'c'], ['1', '2', '3'], ['4', '5', '6']];

    my $table_content = join("\n", map { join("\t", @$_) } @{$arr});
    open my $fh, '>', $tmp_file or die;
    print $fh "$table_content\n";
    close $fh;

    my $table = read_table($tmp_file);
    is_deeply($table, $arr, 'read_table reads tab-separated table');

    my $csv_content = join("\n", map { join(',', @$_) } @{$arr});
    open $fh, '>', $tmp_file or die;
    print $fh "$csv_content\n";
    close $fh;

    $table = read_table($tmp_file, {sep => ','});
    is_deeply($table, $arr, 'read_table reads comma-separated table');

    my $skip_arr = [['skip1'], ['skip2'], ['keep1'], ['keep2'], ['keep3']];
    
    my $skip_content = join("\n", map { join('', @$_) } @{$skip_arr});
    open $fh, '>', $tmp_file or die;
    print $fh "$skip_content\n";
    close $fh;

    $table = read_table($tmp_file, {skip => '0-2', sep => "\n"});
    is_deeply(
        $table,
        [['keep1'], ['keep2'], ['keep3']],
        'read_table skips lines with range'
    );

    $table = read_table($tmp_file, {skip => '1,3', sep => "\n"});
    is_deeply(
        $table,
        [['skip1'], ['keep1'], ['keep3']],
        'read_table skips lines with comma-separated indices'
    );
};

done_testing();
