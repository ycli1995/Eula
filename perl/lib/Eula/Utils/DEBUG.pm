package Eula::Utils::DEBUG;

use strict;
use warnings;

use Carp qw(carp confess);
use Term::ANSIColor qw(:constants);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw/
	check_keys_required 
	deref 
	excu_shell 
	ftime 
	message 
	require_opt 
	stop 
	time_log 
	warning
/;

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

	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime();
	$year += 1900;
	$mon += 1;

	my $pattern = $no_space ? "%d\-%02d\-%02d %02d:%02d:%02d" : "%d\-%02d\-%02d-%02d-%02d-%02d";

	my $ftime = sprintf($pattern, $year, $mon, $day, $hour, $min, $sec);
	return $ftime;
}

sub excu_shell {
	my ($cmd, $run, $print) = @_;

	$run //= 1;
	$print //= 1;
	$print = 1 unless $run;

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

    stop("<\$opts> must be a hashref") unless (ref($opts) eq 'HASH');

    return $opts->{$key} if (exists $opts->{$key});

    my $func = (split /::/, (caller(1))[3])[-1];

    stop("<\$opts->{$key}> is required in function <$func>");
}

sub message {
	my $word = shift;
	my $color = GREEN;
	if ($word eq 'Error') {
		$color = RED;
	} elsif ($word eq 'Warning') {
		$color = YELLOW;
	}
	print ("[", ftime(), "]", "[", $color, $word, RESET, "]", @_);
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

1;
