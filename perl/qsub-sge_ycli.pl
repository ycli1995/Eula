#!/usr/bin/env perl

use strict;
use warnings;

my $USAGE =<<USAGE;
Usage:

  perl $0 [options] <shell_file>

  options:
	--qsub          <str> Path of `qsub` command. Default: the system `which qsub`
	--queue	        <str> Specify the queue to use, Default: all.q
	--resource      <str> The required resource used in `qsub -l` option. Default: mf=10G
	--maxjob        <int> The maximum number of jobs to throw out. Default: 4
	--ppn           <int> The slot number for a job. Default: 1
	--interval      <int> Interval time (in seconds) of checking by `qstat`. Default: 3 seconds
	--lines         <int> Number of lines to form a job. Default: 999999
	--prefix        <str> The prefix tag for qsubed jobs. Default: 'work'
	--nsub          <int> Max qsub times for a job. Default: 1
	--verbose       Show the progress.
	--help

Developed by: Yuchen Li <ycli1995\@outlook.com>

USAGE

use Carp qw(confess);
use Data::Dumper;
use File::Basename;
use File::Path qw/make_path mkpath/;
use Cwd qw/abs_path/;

use Getopt::Long qw/GetOptionsFromArray/;
Getopt::Long::Configure("no_ignore_case", "no_auto_abbrev", "bundling_override", "prefix_pattern=(--)");

my $COMPLETE_STR = "This-Work-is-Completed!";

my $RUNNING_JOBS = {};  # job_id => shell_name
my $QSUB_INFO = {};

END { 
	&kill_all_running_jobs(); 
}

for my $sig (qw/TERM INT KILL/) {
	$SIG{$sig} = sub {
		print "\nCaught SIG$sig.\nWait for killing all running jobs.\n";
		&kill_all_running_jobs();
		exit 1;
	};
}

&qsub_main(@ARGV);

# Sub-routines ############################
sub default_qsub_opts {
	return {
		qsub => 'qsub',
		queue => "all.q",
		resource => "mf=10G",
		ppn => 1,
		maxjob => 4,
		interval => 3,
		lines => 999999,
		prefix => "work",
		nsub => 1,
	};
}

sub qsub_main {
	die $USAGE unless @_;

	my $opts = {};

	GetOptionsFromArray(
		\@_,
		$opts,
		'qsub=s',
		'queue=s',
		'resource=s',
		'ppn=s',
		'lines=i',
		'maxjob=i',
		'prefix=s',
		'nsub=i',
		'verbose',
		'help',
	) or die $USAGE;

	my @invalid_opts = grep { /^-/ } @_;
	die "Invalid options: " . join(", ", @invalid_opts) . "\n\n" . $USAGE if (@invalid_opts);

	die $USAGE if ($opts->{help});

	my $shell_file = shift;
	die $USAGE unless $shell_file;

	## Get default options
	my $default_opts = default_qsub_opts();
	get_default_options($opts, $default_opts);

	## Init qsub
	time_log("Initiate qsub") if $opts->{verbose};
	$QSUB_INFO->{shell_file} = real_path($shell_file, {must_work => 1});
	$QSUB_INFO->{outdir} = "$QSUB_INFO->{shell_file}.$$.qsub";
	$QSUB_INFO->{shells} = split_shell($QSUB_INFO->{shell_file}, $QSUB_INFO->{outdir}, $opts);
	$QSUB_INFO->{ALL_JOBS} = { map { $_ => [] } @{$QSUB_INFO->{shells}} };  # shell_name => [ job_id1, job_id2, ... ]
	$QSUB_INFO->{COUNT_JOBS} = { map { $_ => 0 } @{$QSUB_INFO->{shells}} };  # shell_name => nsub
	$QSUB_INFO->{OK_JOBS} = {};

	normal_qsub_opts($QSUB_INFO, $opts);

	## Run qsub
	time_log("Run qsub jobs parallely") if $opts->{verbose};
	run_qsub_parallel($QSUB_INFO, $RUNNING_JOBS, $opts);

	## Summary
	time_log("Summary qsub") if $opts->{verbose};
	summary_qsub($QSUB_INFO);
}

sub normal_qsub_opts {
	my ($info, $opts) = @_;

	$opts->{nsub} = 1 if ($opts->{nsub} < 1);
	$opts->{maxjob} = 1 if ($opts->{maxjob} < 1);

	my $n_shells = scalar @{$info->{shells}};
	$opts->{maxjob} = $n_shells if ($opts->{maxjob} > $n_shells);
}

sub kill_all_running_jobs {
	my @jobs = keys %{$RUNNING_JOBS};
	return unless @jobs;
	
	my $cmd = "qdel " . join(" ", @jobs);
	print $cmd . "\n";
	system($cmd);

	delete @{$RUNNING_JOBS}{@jobs};
}

sub summary_qsub {
	my $QSUB_INFO = shift;

	my @fail_cmds;
	foreach (@{$QSUB_INFO->{shells}}) {
		next if (exists $QSUB_INFO->{OK_JOBS}->{$_});
		my $cmds = read_shell("$QSUB_INFO->{outdir}/$_");
		push @fail_cmds, @{$cmds};
	}

	return unless (@fail_cmds);

	my $fail_shell = "$QSUB_INFO->{shell_file}.$$.fail.sh";
	write_shell(\@fail_cmds, $fail_shell, {for_qsub => 0});

	print "Some jobs failed.\n";
	print "qsub: " . $QSUB_INFO->{outdir} . "/\n";
	print "failed shell: " . $fail_shell . "\n";
	exit 1;
}

sub run_qsub_parallel {
	my ($QSUB_INFO, $RUNNING_JOBS, $opts) = @_;

	while (1) {
		sleep $opts->{interval};

		update_running_qsub($QSUB_INFO, $RUNNING_JOBS);
		my $n_remain = scalar(keys %{$QSUB_INFO->{COUNT_JOBS}});
		last unless ($n_remain);

		$opts->{maxjob} = $n_remain if ($opts->{maxjob} > $n_remain);
		my $n_running = scalar keys %{$RUNNING_JOBS};
		next unless ($n_running < $opts->{maxjob});

		fill_running_qsub($QSUB_INFO, $RUNNING_JOBS, $opts);
	}	
}

sub update_running_qsub {
	my ($QSUB_INFO, $RUNNING_JOBS) = @_;  # job_id => shell_name

	my @running_jobs = keys %{$RUNNING_JOBS};
	foreach (@running_jobs) {
		my $info = { job_id => $_ };
		update_qstat_info($info);

		# Ignore running job.
		next if ($info->{exists});

		# Delete job that doesn't exist.
		my $shell_name = $RUNNING_JOBS->{$_};
		delete $RUNNING_JOBS->{$_};

		# Only update shell state for this job id (OK or fail?)
		updata_shell_result($QSUB_INFO, $shell_name);
	}
}

sub fill_running_qsub {
	my ($QSUB_INFO, $RUNNING_JOBS, $opts) = @_;  # %RUNNING_JOBS: job_id => shell_name

	my @shells = sort keys %{$QSUB_INFO->{COUNT_JOBS}};
	foreach (@shells) {
		last unless (keys %{$RUNNING_JOBS} < $opts->{maxjob});

		# Check if this shell is completed or failed.
		updata_shell_result($QSUB_INFO, $_);
		next if (exists $QSUB_INFO->{OK_JOBS}->{$_});  # Already ok
		next unless (exists $QSUB_INFO->{COUNT_JOBS}->{$_});  # Cannot run anymore

		if ($QSUB_INFO->{COUNT_JOBS}->{$_} < $opts->{nsub}) {
			next if (is_running_shell($QSUB_INFO, $RUNNING_JOBS, $_));  # Already running

			time_log("Submit $_") if ($opts->{verbose});
			qsub_a_shell($QSUB_INFO, $RUNNING_JOBS, $_, $opts);
			next;
		}
		delete $QSUB_INFO->{COUNT_JOBS}->{$_};
	}
}

sub is_running_shell {
	my ($QSUB_INFO, $RUNNING_JOBS, $shell) = @_;

	update_running_qsub($QSUB_INFO, $RUNNING_JOBS);

	foreach (keys %{$RUNNING_JOBS}) {
		return 1 if ($shell eq $RUNNING_JOBS->{$_});
	}
	return 0;
}

sub qsub_a_shell {
	my ($QSUB_INFO, $RUNNING_JOBS, $shell, $opts) = @_;

	my $shfile = $QSUB_INFO->{outdir} . "/" . basename($shell);
	my $job_id = qsub_job($shfile, $opts);

	$RUNNING_JOBS->{$job_id} = $shell;
	push @{$QSUB_INFO->{ALL_JOBS}->{$shell}}, $job_id;
	$QSUB_INFO->{COUNT_JOBS}->{$shell} += 1;
}

sub updata_shell_result {
	my ($QSUB_INFO, $shell) = @_;

	my $ok = check_shell_complete($QSUB_INFO->{outdir}, $shell);
	return unless $ok;

	# $shell is already ok.
	$QSUB_INFO->{OK_JOBS}->{$shell} = $QSUB_INFO->{ALL_JOBS}->{$shell};
	delete $QSUB_INFO->{COUNT_JOBS}->{$shell};
}

# shell output
sub check_shell_complete {
	my ($indir, $shell) = @_;

	my $oe = find_shell_output($indir, $shell);
	foreach (@{$oe->{o}}) {
		my $ofile = $indir . "/" . basename($_);
		return 1 if check_o_complete($ofile);
	}
	return 0;
}

sub check_o_complete {
	my $ofile = shift;
	$ofile = real_path($ofile, {must_work => 1});

	my $out = 0;
	open my $fh, "<", $ofile or die $!;
	while (<$fh>) {
		chomp;
		next unless ($_ eq $COMPLETE_STR);
		$out = 1;
		last;
	}
	close $fh;
	return $out;
}

sub find_shell_output {
	my ($indir, $shell_name) = @_;

	my $target = basename($shell_name);

	my (@e, @o);

	opendir my $dh, $indir or die $!;
	my @all_files = grep { /\.sh\.o[0-9]+$/ || /\.sh\.e[0-9]+$/ } readdir($dh);

	foreach (@all_files) {
		my $file = $_;
		$file =~ s/[0-9]+$//g;
		push(@o, $indir . "/" . $_) if ($file eq "$target.o");
		push(@e, $indir . "/" . $_) if ($file eq "$target.e");
	}
	return { o => \@o, e => \@e };
}

# qsub
sub update_qstat_info {
	my $info = shift;
	check_ref($info, "HASH");

	my $job_id = $info->{job_id};
	die "No 'job_id'." unless defined $job_id;

	my $stat = `qstat -j $job_id 2>&1`;
	if ($stat =~ /do not exist/) {
		$info->{exists} = 0;
		return;
	}
	$info->{exists} = 1;
}

sub qsub_job {
	my ($shfile, $opts) = @_;

	$opts //= {};

	$shfile = real_path($shfile, {must_work => 1});

	my $shdir = dirname($shfile);
	my $qcmd = qsub_cmd($opts);

	my $job_out = `$qcmd -e $shdir -o $shdir $shfile`;
	my $job_id = $1 if ($job_out =~ /Your job (\d+)/);
	stop("Fetching job id failed: \n$job_out") if ($job_id eq '0' || !defined $job_id);

	return $job_id;
}

sub qsub_cmd {
	my $opts = shift;

	$opts //= {};

	my $default_opts = default_qsub_opts();
	get_default_options($opts, $default_opts);

	&check_keys_required($opts, qw/qsub queue resource ppn/);

	my $cmd = "$opts->{qsub} -q $opts->{queue} -l $opts->{resource} -cwd -S /bin/bash -V -sync no -pe smp $opts->{ppn} ";
	return $cmd;
}

# Shell
sub read_shell {
	my ($file, $opts) = @_;

	$opts //= {};
	$opts->{lines} //= 999999;

	$file = real_path($file, {must_work => 1});

	my $cmds = [];
	my $cmd = [];
	open IN, "<", $file or die $!;
	while (<IN>) {
		chomp;
		my $add = filter_cmd($_);
		next unless (defined $add);

		push @{$cmd}, $add;
		next unless (@{$cmd} == $opts->{lines});

		push @{$cmds}, join("\n", @{$cmd});
		$cmd = [];
	}
	close IN;
	push @{$cmds}, join("\n", @{$cmd}) unless (@{$cmd} == 0);
	return $cmds;
}

sub split_shell {
	my ($file, $outdir, $opts) = @_;

	$opts //= {};
	$opts->{lines} //= 999999;
	$opts->{pfx} //= "work";

	my $cmds = read_shell($file, $opts);

	my $shells = [];
	my $job_id = "0001";
	my $cmd;
	foreach (@{$cmds}) {
		my $outfile = "$outdir/$opts->{pfx}_${job_id}.sh";
		make_path($outdir);
		write_shell([$_], $outfile);
		push @{$shells}, basename($outfile);
		$job_id++;
	}
	return $shells;
}

sub write_shell {
	my ($cmds, $outfile, $opts) = @_;

	$opts //= {};
	$opts->{for_qsub} //= 1;

	my $script_pfx = "#!/bin/bash\nset -e";
	my $script_sfx = "";

	if ($opts->{for_qsub}) {
		$script_pfx .= "\necho \$HOSTNAME";
		$script_sfx = "echo  $COMPLETE_STR";
	}

	open OUT, '>', $outfile or stop($!);
	print OUT "$script_pfx\n";
	foreach (@{$cmds}) {
		print OUT "$_\n";
	}
	print OUT "$script_sfx\n" unless ($script_sfx eq "");
	close OUT;

	return basename($outfile);
}

sub filter_cmd {
	my $cmd = shift;

	return undef unless $cmd;
	return undef if ($cmd =~ /^set -e$/);
	return undef if ($cmd =~ /^#/);
	return undef if ($cmd =~ /\Q$COMPLETE_STR\E/);
	return undef if ($cmd =~ /^echo \$HOSTNAME$/);

	return $cmd;
}

# PATH
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

# Options
sub check_keys_required {
	my $opts = shift;
	foreach (@_) {
		next unless (!defined $opts->{$_} || $opts->{$_} eq '');
		stop("Required key '$_' not found in the input hash.");
	}
}

sub get_default_options {
	my ($opts, $default_opts) = @_;
	foreach (keys %{$default_opts}) {
		$opts->{$_} //= $default_opts->{$_};
	}
}

# DEBUG
sub check_ref {
	my ($x, $ref) = @_;
	return if (ref $x eq $ref);
	stop("Require a ref of '$ref'.");
}

sub stop {
	confess "ERROR:\n @_\n";
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

sub time_log {
	my $time = "[" . ftime(1) . "]";
	print "$time @_\n";
}

