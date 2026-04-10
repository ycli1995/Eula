package Eula::Utils::PATH;

use strict;
use warnings;

use Cwd qw/abs_path/;
use File::Remove qw/rm/;

use Eula::Utils::DEBUG qw/stop/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(real_path rm_paths);

sub real_path {
	my ($path, $opts) = @_;

	$opts //= {};
	$opts->{must_work} //= 0;

	my $norm_path = abs_path($path);
	if (defined $norm_path) {
		$norm_path = undef unless (-e $norm_path);
	}
	return $norm_path if defined $norm_path;

	$path =~ s/\/+$//g;
	return $path unless $opts->{must_work};

	stop('No such file or directory: ', $path);
}

sub rm_paths {
	rm(\1, $_) foreach (@_);
}

1;
