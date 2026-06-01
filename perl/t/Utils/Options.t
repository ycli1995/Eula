use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use lib "$Bin/../../lib";
use Eula::Utils::Options;
use File::Temp qw(tempdir);

# Test get_defaults
subtest 'get_defaults' => sub {
    plan tests => 6;

    my $opts = {};
    my $defaults = {key1 => 'value1', key2 => 'value2'};
    get_defaults($opts, $defaults);
    foreach (keys %$defaults) {
        is($opts->{$_}, $defaults->{$_}, "get_defaults sets missing key $_");
    }

    $opts = {key1 => 'existing'};
    get_defaults($opts, $defaults);
    is($opts->{key1}, 'existing', 'get_defaults preserves existing key');
    is($opts->{key2}, 'value2', 'get_defaults sets missing key2');

    $opts = {key1 => '', key2 => 'existing'};
    get_defaults($opts, $defaults);
    is($opts->{key1}, '', 'get_defaults preserves empty string');
    is($opts->{key2}, 'existing', 'get_defaults preserves non-empty value');
};

# Test get_defaults with undef defaults
subtest 'get_defaults with undef defaults' => sub {
    plan tests => 1;

    my $opts = {key1 => 'value1'};
    get_defaults($opts, undef);
    is($opts->{key1}, 'value1', 'get_defaults handles undef defaults');
};

# Test norm_option_array
subtest 'norm_option_array' => sub {
    # Test undef input
    my $result = norm_option_array(undef);
    is_deeply($result, [], 'norm_option_array with undef returns empty array');

    # Test empty string
    $result = norm_option_array('');
    is_deeply($result, [], 'norm_option_array with empty string');

    # Test simple scalar
    $result = norm_option_array('single_value');
    is_deeply(
        $result,
        ['single_value'],
        'norm_option_array with scalar returns array'
    );

    # Test comma-separated string
    $result = norm_option_array('a,b,c');
    is_deeply(
        $result,
        ['a', 'b', 'c'],
        'norm_option_array splits comma-separated string'
    );

    # Test array reference
    $result = norm_option_array(['x', 'y', 'z']);
    is_deeply(
        $result,
        ['x', 'y', 'z'],
        'norm_option_array with array ref returns same array'
    );

    # Test nested array reference
    $result = norm_option_array(['a,b', ['c', 'd,e']]);
    is_deeply(
        $result,
        ['a', 'b', 'c', 'd', 'e'],
        'norm_option_array handles nested arrays recursively'
    );

    # Test file input
    my $opts = ['line1', 'line2', 'line3'];
    my $tmp_dir = tempdir(CLEANUP => 1);
    my $tmp_file = "$tmp_dir/options.txt";
    open my $fh, '>', $tmp_file or die;
    print $fh join("\n", @{$opts}) . "\n";
    close $fh;

    $result = norm_option_array($tmp_file);
    is_deeply($result, $opts, 'norm_option_array reads options from file');

    # Test non-existent file (should be treated as scalar)
    $result = norm_option_array('/non/existent/path.txt');
    is_deeply(
        $result,
        ['/non/existent/path.txt'],
        'norm_option_array treats non-existent file as scalar'
    );
};

done_testing();
