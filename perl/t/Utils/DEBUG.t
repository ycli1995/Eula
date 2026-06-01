use strict;
use warnings;
use Test::More tests => 12;
use FindBin qw($Bin);
use lib "$Bin/../../lib";
use Eula::Utils::DEBUG;

# Test ftime
subtest 'ftime' => sub {
    plan tests => 3;

    my $result = ftime();
    like(
        $result,
        qr/^\d{4}\-\d{2}\-\d{2}\-\d{2}\-\d{2}\-\d{2}$/,
        'ftime() returns correct format'
    );

    my $result_no_space = ftime(1);
    like(
        $result_no_space,
        qr/^\d{4}\-\d{2}\-\d{2} \d{2}:\d{2}:\d{2}$/,
        'ftime(1) returns correct format with space'
    );

    my $result_default = ftime(0);
    like(
        $result_default,
        qr/^\d{4}\-\d{2}\-\d{2}\-\d{2}\-\d{2}\-\d{2}$/,
        'ftime(0) returns correct format'
    );
};

# Test check_ref
subtest 'check_ref' => sub {
    plan tests => 6;

    ok(check_ref({}, 'HASH'), 'check_ref with HASH ref returns true');
    ok(check_ref([], 'ARRAY'), 'check_ref with ARRAY ref returns true');
    ok(
        !check_ref({}, 'ARRAY', {stop => 0}),
        'check_ref with wrong ref type returns false'
    );
    ok(
        !check_ref('scalar', 'HASH', {stop => 0}),
        'check_ref with non-ref returns false'
    );

    ok(
        check_ref({}, 'HASH', {stop => 0}),
        'check_ref with stop => 0 and correct type'
    );
    ok(
        !check_ref([], 'HASH', {stop => 0}),
        'check_ref with stop => 0 and wrong type'
    );
};

# Test deref
subtest 'deref' => sub {
    plan tests => 8;

    my $scalar = 'test';
    is(deref($scalar), $scalar, 'deref scalar returns value');

    my $scalar_ref = \$scalar;
    is(deref($scalar_ref), $scalar, 'deref SCALAR ref returns value');

    my $array = [1, 2, 3];
    my @array_result = deref($array);
    is_deeply(\@array_result, $array, 'deref ARRAY ref returns array');

    my $hash = {a => 1, b => 2};
    my %hash_result = deref($hash);
    is_deeply(\%hash_result, {a => 1, b => 2}, 'deref HASH ref returns hash');

    my $ref_ref = \$scalar_ref;
    is(deref($ref_ref), 'test', 'deref REF ref returns value');

    my $code = sub { return 'code'; };
    is(deref($code), 'code', 'deref CODE ref calls and returns result');

    my $glob = \*STDOUT;
    ok(deref($glob), 'deref GLOB ref returns glob');

    my $undef;
    is(deref($undef), undef, 'deref undef returns undef');
};

# Test get_omits
subtest 'get_omits' => sub {
    plan tests => 5;

    my $result = get_omits('1,3-5,7');
    is_deeply(
        $result,
        [[1,2], [3,5], [7,8]],
        'get_omits parses comma-separated values'
    );

    $result = get_omits('');
    is_deeply($result, [], 'get_omits with empty string returns empty array');

    $result = get_omits(undef);
    is_deeply($result, [], 'get_omits with undef returns empty array');

    $result = get_omits('10-15');
    is_deeply($result, [[10,15]], 'get_omits parses single range');

    $result = get_omits('2,5:7,9-12');
    is_deeply(
        $result,
        [[2,3], [5,7], [9,12]],
        'get_omits parses complex string'
    );
};

# Test omit_array
subtest 'omit_array' => sub {
    plan tests => 3;

    my $array = [1, 2, 3, 4, 5];
    omit_array($array, [1, 3]);
    is_deeply(
        $array,
        [1, 3, 5],
        'omit_array removes elements at specified indices'
    );

    $array = [1, 2, 3];
    omit_array($array, []);
    is_deeply($array, [1, 2, 3], 'omit_array with empty omits does nothing');

    $array = [1, 2, 3];
    omit_array($array, undef);
    is_deeply($array, [1, 2, 3], 'omit_array with undef omits does nothing');
};

# Test omit_array_by_ranges
subtest 'omit_array_by_ranges' => sub {
    plan tests => 3;

    my $array = [1, 2, 3, 4, 5, 6, 7];
    omit_array_by_ranges($array, [[1, 3], [5, 7]]);
    is_deeply(
        $array,
        [1, 4, 5],
        'omit_array_by_ranges removes elements in ranges'
    );

    $array = [1, 2, 3];
    omit_array_by_ranges($array, []);
    is_deeply(
        $array,
        [1, 2, 3],
        'omit_array_by_ranges with empty omits does nothing'
    );

    $array = [1, 2, 3];
    omit_array_by_ranges($array, undef);
    is_deeply(
        $array,
        [1, 2, 3],
        'omit_array_by_ranges with undef omits does nothing'
    );
};

# Test require_opt
subtest 'require_opt' => sub {
    plan tests => 3;

    my $opts = {key1 => 'value1', key2 => 'value2'};
    is(
        require_opt($opts, 'key1'),
        'value1',
        'require_opt returns value for existing key'
    );

    eval { require_opt($opts, 'nonexistent'); };
    like($@, qr/is required in function/, 'require_opt dies for missing key');

    eval { require_opt('not_hash', 'key'); };
    like(
        $@,
        qr/Expected a reference of type 'HASH'/,
        'require_opt dies for non-hashref'
    );
};

# Test stop
subtest 'stop' => sub {
    plan tests => 1;

    eval { stop('test error'); };
    like($@, qr/ERROR: test error/, 'stop throws exception with message');
};

# Test warning
subtest 'warning' => sub {
    plan tests => 1;

    my $warning;
    local $SIG{__WARN__} = sub { $warning = shift; };
    warning('test warning');
    like($warning, qr/WARN: test warning/, 'warning outputs warning message');
};

# Test excu_shell
subtest 'excu_shell' => sub {
    plan tests => 5;

    my $output = '';
    open my $fh, '>', \$output or die;
    my $old_fh = select($fh);

    excu_shell('echo "hello"', {run => 0});
    like($output, qr/echo "hello"/, 'excu_shell prints command when run=0');

    $output = '';
    excu_shell('echo "hello"', {run => 1, print => 0});
    is($output, '', 'excu_shell does not print when print=0');

    $output = '';
    excu_shell('echo "hello"', {run => 1, print => 1});
    like($output, qr/echo "hello"/, 'excu_shell prints command when print=1');

    select($old_fh);
    close $fh;

    my $result = excu_shell('echo "test"');
    is($result, undef, 'excu_shell returns undef on success');

    eval { excu_shell('false'); };
    like($@, qr/Error when runing shell command/, 'excu_shell throws on failure');
};

# Test time_log
subtest 'time_log' => sub {
    plan tests => 1;

    my $output = '';
    open my $fh, '>', \$output or die;
    my $old_fh = select($fh);

    time_log('test message');

    select($old_fh);
    close $fh;

    like(
        $output,
        qr/\[\d{4}\-\d{2}\-\d{2}\-\d{2}\-\d{2}\-\d{2}\] test message/,
        'time_log outputs formatted message'
    );
};

# Test message
subtest 'message' => sub {
    plan tests => 3;

    my $output = '';
    open my $fh, '>', \$output or die;
    my $old_fh = select($fh);

    message('Info', ' test info');
    like(
        $output,
        qr/\[.*\]\[.*Info.*\] test info/,
        'message outputs info message'
    );

    $output = '';
    message('Error', ' test error');
    like(
        $output,
        qr/\[.*\]\[.*Error.*\] test error/,
        'message outputs error message'
    );

    $output = '';
    message('Warning', ' test warning');
    like(
        $output,
        qr/\[.*\]\[.*Warning.*\] test warning/,
        'message outputs warning message'
    );

    select($old_fh);
    close $fh;
};

done_testing();
