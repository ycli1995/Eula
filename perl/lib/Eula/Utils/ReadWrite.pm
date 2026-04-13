package Eula::Utils::ReadWrite;

use strict;
use warnings;

use YAML;

use Eula::Utils::DEBUG qw/stop/;
use Eula::Utils::PATH qw/real_path/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(read_lines read_shell read_yaml write_lines write_shell write_yaml);

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
	return @lines;
}

sub write_lines {
	my $file = shift;

	open my $fh, ">", $file or stop($!);
	print $fh "$_\n" foreach (@_);
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

sub read_shell {
	my ($file, $opts) = @_;

	my @cmds;
	open my $fh, "<", $file or die $!;
	while (<$fh>) {
		chomp;
		my $add = filter_cmd($_);
		next unless (defined $add);
		push @cmds, $add;
	}
	close $fh;
	return @cmds;
}

sub write_shell {
	my ($file, $cmds, $opts) = @_;

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
	return undef unless $cmd;
	return undef if ($cmd =~ /^set -e$/);
	return undef if ($cmd =~ /^#/);
	foreach (@{$opts->{skip_patterns}}) {
		return undef if ($cmd =~ /\Q$_\E/);
	}
	return $cmd;
}

1;
