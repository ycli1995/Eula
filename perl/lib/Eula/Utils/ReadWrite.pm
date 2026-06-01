package Eula::Utils::ReadWrite;

use strict;
use warnings;

use YAML;

use Eula::Utils::DEBUG qw/check_ref get_omits stop/;
use Eula::Utils::PATH qw/real_path/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw/
    read_lines
    read_shell
    read_table
    read_yaml
    write_lines
    write_shell
    write_yaml
/;

$YAML::Preserve = 1;

sub read_lines {
	my ($file, $opts) = @_;
	$file = real_path($file, {must_work => 1});

	$opts //= {};
	$opts->{skip_blank} //= 1;
	$opts->{trim} //= 1;

	open my $fh, "<", $file or stop($!);
	my @lines;
	while (<$fh>) {
		chomp;
		next if ($opts->{skip_blank} && !/\S/);
		if ($opts->{trim}) {
			s/^\s+|\s+$//g;
		}
		push @lines, $_;
	}
	return \@lines;
}

sub write_lines {
	my ($file, $arr) = @_;

	open my $fh, ">", $file or stop($!);
	print $fh "$_\n" foreach (@{$arr});
	close $fh;
}

sub read_yaml {
	my ($file) = @_;
	$file = real_path($file, {must_work => 1});

	my $config = YAML::LoadFile($file);
	return $config;
}

sub write_yaml {
	my ($file, $config) = @_;

	open my $fh, ">", $file or stop($!);
	print $fh YAML::Dump $config;
	close $fh;
}

sub read_table {
	my ($file, $opts) = @_;

	$file = real_path($file, {must_work => 1});

	$opts //= {};
	my $sep = $opts->{sep} // "\t";
	my $skip = $opts->{skip} // 0;
	if ($skip) {
		$skip = get_omits($skip);
	}

	my @arr;
	my $ncol;
	my $i = 0;
	open my $fh, "<", $file or die $!;
	while (<$fh>) {
		chomp;
		if (_is_to_skip($skip, $i)) {
			$i++;
			next;
		}
		my @tmp = split $sep;
		push @arr, \@tmp;
		$i++;
	}
	close $fh;
	return \@arr;
}

sub _is_to_skip {
	my ($omits, $i) = @_;
	return 0 unless (ref($omits) eq 'ARRAY');

	foreach (@$omits) {
		return 1 if ($i >= $_->[0] && $i < $_->[1]);
	}
	return 0;
}

sub read_shell {
	my ($file, $opts) = @_;

	$file = real_path($file, {must_work => 1});

	my @cmds;
	open my $fh, "<", $file or die $!;
	while (<$fh>) {
		chomp;
		my $add = filter_cmd($_);
		next unless (defined $add);
		push @cmds, $add;
	}
	close $fh;
	return \@cmds;
}

sub write_shell {
	my ($file, $cmds, $opts) = @_;

	check_ref($cmds, "ARRAY");

	$opts //= {};

	open my $fh, ">", $file or stop($!);
	print $fh "#!/bin/bash\n";
	print $fh "set -e\n";
	foreach (@{$cmds}) {
		my $line = filter_cmd($_, $opts);
		print $fh "$line\n" if (defined $line);
	}
	close $fh;
}

sub filter_cmd {
	my ($cmd, $opts) = @_;

	$opts //= {};
	$opts->{skip_patterns} //= [];
	check_ref($opts->{skip_patterns}, "ARRAY");

	return undef unless $cmd;
	return undef if ($cmd =~ /^set -e$/);
	return undef if ($cmd =~ /^#/);
	foreach (@{$opts->{skip_patterns}}) {
		return undef if ($cmd =~ /\Q$_\E/);
	}
	return $cmd;
}

1;
