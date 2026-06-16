package Eula::Utils::DEBUG;

use strict;
use warnings;

use Scalar::Util qw(blessed);
use Carp qw(carp confess);
use Term::ANSIColor qw(:constants);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
    check_keys_required
    check_ref
    deref
    excu_shell
    ftime
    get_omits
    is_instance_of
    message
    omit_array
    omit_array_by_ranges
    require_class
    require_opt
    stop
    time_log
    warning
);

sub time_log {
    my $time = "[" . ftime() . "]";
    print "$time @_\n";
}

sub warning {
    carp "WARN: @_\n";
}

sub stop {
    confess "ERROR: @_\n";
}

sub stop2 {
    my $ename = shift;
    confess "[$ename] ERROR:\n @_\n";
}

sub ftime {
    my ($no_space) = @_;
    $no_space //= 0;

    my (
        $sec, $min, $hour, $day, $mon,
        $year, $wday, $yday, $isdst
    ) = localtime();
    $year += 1900;
    $mon += 1;

    my $pattern = $no_space
        ? "%d-%02d-%02d %02d:%02d:%02d"
        : "%d-%02d-%02d-%02d-%02d-%02d";

    my $ftime = sprintf($pattern, $year, $mon, $day, $hour, $min, $sec);
    return $ftime;
}

sub excu_shell {
    my ($cmd, $opts) = @_;

    $opts //= {};
    my $print = $opts->{print} // 1;
    my $run = $opts->{run} // 1;
    $print = 1 unless ($run);

    print $cmd . "\n" if ($print);
    return unless $run;

    my $out = system($cmd);
    return unless $out;

    stop("Error when runing shell command: $cmd");
}

sub check_keys_required {
    my $opts = shift;
    foreach (@_) {
        next unless (!defined $opts->{$_} || $opts->{$_} eq '');
        stop("Required key '$_' not found in the input hash.\n");
    }
}

sub require_opt {
    my ($opts, $key) = @_;

    check_ref($opts, "HASH");

    stop("<\$opts> must be a hashref") unless (ref($opts) eq 'HASH');

    return $opts->{$key} if (exists $opts->{$key});

    my $func = (split /::/, (caller(1))[3])[-1];

    stop("<\$opts->{$key}> is required in function <$func>");
}

sub is_instance_of {
    my ($obj, $class) = @_;
    return blessed($obj) && $obj->isa($class);
}

sub require_class {
    my ($obj, $class) = @_;
    return if (is_instance_of($obj, $class));

    my $func = (split /::/, (caller(1))[3])[-1];
    stop("Expected a reference of class '$class' in function <$func>");
}

sub message {
    my $word = shift;
    my $color = GREEN;
    if ($word eq 'Error') {
        $color = RED;
    } elsif ($word eq 'Warning') {
        $color = YELLOW;
    }
    print(
        "[", ftime(), "]",
        "[", $color, $word, RESET, "]",
        @_
    );
}

sub deref {
    my $things = shift;
    return $things unless (ref $things);
    return deref($$things) if (ref $things eq 'REF');
    return $$things if (ref $things eq 'SCALAR');
    return @{$things} if (ref $things eq 'ARRAY');
    return %{$things} if (ref $things eq 'HASH');
    return &$things if (ref $things eq 'CODE');
    return *$things if (ref $things eq 'GLOB');
    return $things;
}

sub get_omits {
    my $omits = shift;

    $omits //= "";

    return [] unless ($omits);

    my @all_omits;
    my @tmp = split /,/, $omits;
    foreach (@tmp) {
        if (/^\d+$/) {
            push @all_omits, [$_, $_ + 1];
            next;
        }
        if (/(\d+)\s*[:-]\s*(\d+)/) {
            my ($start, $end) = ($1, $2);
            if ($end < $start) {
                stop("Invalid omits values: [$start > $end]");
            }
            push @all_omits, [$start, $end];
        }
    }
    return \@all_omits;
}

sub omit_array {
    my ($arr, $omits) = @_;
    check_ref($arr, "ARRAY");

    $omits //= [];
    check_ref($omits, "ARRAY");

    foreach (sort { $b <=> $a } @{$omits}) {
        splice(@{$arr}, $_, 1);
    }
}

sub omit_array_by_ranges {
    my ($arr, $omits) = @_;
    check_ref($arr, "ARRAY");

    $omits //= [];
    check_ref($omits, "ARRAY");
    foreach (sort { $b->[0] <=> $a->[0] } @{$omits}) {
        my ($start, $end) = @$_;
        splice(@{$arr}, $start, $end - $start);
    }
}

sub check_ref {
    my ($var, $expected_ref, $opts) = @_;

    $opts //= {};
    my $stop = $opts->{stop} // 1;

    my $actual_ref = ref($var);
    if (!$actual_ref) {
        return 0 unless ($stop);
        stop(
            "Expected a reference of type '$expected_ref', "
            . "got a non-reference value"
        );
    }
    if ($actual_ref ne $expected_ref) {
        return 0 unless ($stop);
        stop(
            "Expected a reference of type '$expected_ref', "
            . "got '$actual_ref'"
        );
    }
    return 1;
}

1;
