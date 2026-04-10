package Eula::Utils::Options;

use strict;
use warnings;

use Eula::Utils::ReadWrite qw/read_lines/;
use Eula::Utils::DEBUG;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_defaults get_options_from_args norm_option_array);

sub get_defaults {
	my ($opts, $defaults) = @_;
	$defaults //= {};
	foreach (keys %{$defaults}) {
		$opts->{$_} //= $defaults->{$_};
	}
}

sub norm_option_array {
	my ($opt) = @_;

	# An empty array	
	return [] unless defined $opt;

	# Fetch an array from a list file
	if (-e $opt && -s $opt) {
		my @arr = read_lines($opt);
		return \@arr;
	}
	# Fetch an array from a string, separated by ',' 
	if ($opt =~ /\S,\S/) {
		my @all_opts = split /,/, $opt;
		return norm_option_array(\@all_opts);
	}
	# Fetch a scalar
	return [ $opt ] unless (ref($opt));
	# Fetch from an array recursively
	if (ref($opt) eq 'ARRAY') {
		my @new_opts = map { @{norm_option_array($_)} } @{$opt};
		return \@new_opts;
	}
	return $opt;
}

sub get_options_from_args {
	my %opts;
	my @args;
	my $args_end = 0;
	my $key = undef;
	foreach (@_) {
		$args_end = 1 if (/^-/);
		# A value
		if (/^[^-]/) {
			if (defined $key) {
				$opts{$key} = $_;
				$key = undef;
				next
			}
			next if ($args_end);
			push @args, $_;
			next;
		}
		if (defined $key) {
			$opts{$key} = 'TRUE';
			$key = undef;
		}
		# A key
		if (/^-[^-]/) {
			$key = $_;
			$key =~ s/^-//;
			print $key . "\n";
			next;
		}
		print STDERR "Ignore invalied option '$_'.\n";
		$key = undef;
	}
	$opts{$key} = 'TRUE' if (defined $key);
	$opts{ARGV} = \@args;
	return %opts;
}

sub make_option_help {
	my ($opt) = @_;

	my @keys = qw/name type default help/;
	foreach (@keys) {
		next if exists $opt->{$_};
		stop("'$_' is not defined for an option.");
	}
	my $default = defined $opt->{$_} ? $opt->{$_} : "NULL";

	my $help = <<HELP;
\t-$opt->{name}\t$opt->{help} default: $default
HELP

	return $help;
}

1;
