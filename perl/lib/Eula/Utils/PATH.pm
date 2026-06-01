package Eula::Utils::PATH;

use strict;
use warnings;

use Cwd qw/abs_path/;
use File::Basename qw/dirname/;

use Eula::Utils::DEBUG qw/stop/;

use lib dirname(__FILE__)."/../../vendor/lib/perl5";
use File::Copy::Recursive qw/rcopy/;
use File::Remove qw/rm/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(cp_paths empty_dir real_path rm_paths);

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

sub cp_paths {
	my ($src, $dst) = @_;
	rcopy($src, $dst);
}

sub empty_dir {
	my ($dir) = @_;
    opendir(my $dh, $dir) or die "Cannot open directory $dir: $!";

    while (my $entry = readdir($dh)) {
        next if ($entry eq '.' || $entry eq '..');
        closedir($dh);
        return 0;
    }
    closedir($dh);
    return 1;
}

1;
