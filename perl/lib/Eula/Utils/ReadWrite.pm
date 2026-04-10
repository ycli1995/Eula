package Eula::Utils::ReadWrite;

use strict;
use warnings;

use YAML;

use Eula::Utils::DEBUG qw/stop/;
use Eula::Utils::PATH qw/real_path/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(read_lines read_yaml write_lines write_yaml);

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

1;
