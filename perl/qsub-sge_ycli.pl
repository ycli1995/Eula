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
	--nsub          <int>
	
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

my %RUNNING_JOBS;  # job_id => shell_name
my %QSUB_INFO;

END { 
	&kill_all_running_jobs(); 
}

$SIG{INT} = sub {
	print "\nInterupted. \nWait for killing all running jobs.\n";
	&kill_all_running_jobs();
	exit 1;
};

$SIG{TERM} = sub {
	print "\nTerminated. \nWait for killing all running jobs.\n";
	&kill_all_running_jobs();
    exit 1;
};

&qsub_main(@ARGV);

# Sub-routines ############################
sub default_options {
	my %opts = (
		qsub => 'qsub',
		queue => "all.q",
		resource => "mf=10G",
		ppn => 1,
		maxjob => 4,
		interval => 3,
		lines => 999999,
		prefix => "work",
		nsub => 1,
	);
	return %opts;
}

sub qsub_main {
	die $USAGE unless @_;

	my %opts;
	GetOptionsFromArray(
		\@_,
		\%opts,
		'qsub=s',
		'queue=s',
		'resource=s',
		'ppn=s',
		'lines=i',
		'maxjob=i',
		'prefix=s',
		'nsub=i',
		'help',
	) or die $USAGE;

	my @invalid_opts = grep { /^-/ } @_;
	die "Invalid options: " . join(", ", @invalid_opts) . "\n\n" . $USAGE if (@invalid_opts);
	die $USAGE if (@_ == 0 || $opts{help});

	my $shell_file = shift;
	die $USAGE unless $shell_file;

	## Get default options
	my %default_opts = default_options();
	&get_default_options(\%opts, \%default_opts);

	$opts{nsub} = 1 if ($opts{nsub} < 1);
	$opts{maxjob} = 1 if ($opts{maxjob} < 1);

	## Get shell commands
	$QSUB_INFO{shell_file} = real_path($shell_file, 1);
	$QSUB_INFO{outdir} = $QSUB_INFO{shell_file} . ".$$.qsub";
	my @shells = split_shell($QSUB_INFO{shell_file}, $QSUB_INFO{outdir}, $opts{lines}, $opts{prefix});

	my $n_shells = scalar @shells;
	$opts{maxjob} = $n_shells if ($opts{maxjob} > $n_shells);

	## Init qsub
	$QSUB_INFO{shells} = \@shells;
	$QSUB_INFO{opts} = \%opts;
	$QSUB_INFO{ALL_JOBS} = { map { $_ => [] } @{$QSUB_INFO{shells}} };  # shell_name => [ job_id ]
	$QSUB_INFO{COUNT_JOBS} = { map { $_ => 0 } @{$QSUB_INFO{shells}} };  # shell_name => nsub
	$QSUB_INFO{OK_JOBS} = {};

	## Run qsub
	run_qsub_parallel(\%QSUB_INFO, \%RUNNING_JOBS);

	## Summary
	summary_qsub(\%QSUB_INFO);
}

sub kill_all_running_jobs {
	my @jobs = keys %RUNNING_JOBS;
	return unless @jobs;
	
	my $cmd = "qdel " . join(" ", @jobs);
	print $cmd . "\n";
	system($cmd);

	delete @RUNNING_JOBS{@jobs};
}

sub summary_qsub {
	my $QSUB_INFO = shift;

	my @fail_cmds;
	foreach (@{$QSUB_INFO->{shells}}) {
		next if (exists $QSUB_INFO->{OK_JOBS}->{$_});

		my $shfile = $QSUB_INFO->{outdir} . "/" . $_;
		my @cmds = read_shell($shfile);
		@fail_cmds = (@fail_cmds, @cmds);
	}

	unless (@fail_cmds) {
		print "All jobs finished!\n";
		return;
	}

	my $fail_shell = $QSUB_INFO->{shell_file} . ".$$.fail.sh";

	write_shell($fail_shell, join("\n", @fail_cmds));

	print "Some jobs failed.\n";
	print "qsub: " . $QSUB_INFO->{outdir} . "/\n";
	print "failed shell: " . $fail_shell . "\n";
	exit 1;
}

sub run_qsub_parallel {
	my $QSUB_INFO = shift;
	my $RUNNING_JOBS_HASH = shift;

	my $interval = $QSUB_INFO->{opts}->{interval};
	my $maxjob = $QSUB_INFO->{opts}->{maxjob};

	while (1) {
		sleep $interval;

		update_running_qsub($QSUB_INFO, $RUNNING_JOBS_HASH);
		my $n_remain = scalar(keys %{$QSUB_INFO->{COUNT_JOBS}});

		last unless ($n_remain);

		$maxjob = $n_remain if ($maxjob > $n_remain);

		my $n_running = scalar keys %{$RUNNING_JOBS_HASH};
		next unless ($n_running < $maxjob);

		fill_running_qsub($QSUB_INFO, $RUNNING_JOBS_HASH, $maxjob);
	}	
}

sub update_running_qsub {
	my $QSUB_INFO = shift;
	my $RUNNING_JOBS_HASH = shift;  # job_id => shell_name

	my @running_jobs = keys %{$RUNNING_JOBS_HASH};
	foreach (@running_jobs) {
		my %info = (job_id => $_);
		update_qstat_info(\%info);

		# Ignore running job.
		next if ($info{exists});

		# Delete job that doesn't exist.
		my $shell_name = $RUNNING_JOBS_HASH->{$_};
		delete $RUNNING_JOBS_HASH->{$_};

		# Only update shell state for this job id (OK or fail?)
		updata_shell_result($QSUB_INFO, $shell_name);
	}
}

sub fill_running_qsub {
	my $QSUB_INFO = shift;
	my $RUNNING_JOBS_HASH = shift;  # job_id => shell_name
	my $maxjob = shift;

	my $nsub = $QSUB_INFO->{opts}->{nsub};

	my @shells = sort keys %{$QSUB_INFO->{COUNT_JOBS}};
	foreach (@shells) {
		last unless (keys %{$RUNNING_JOBS_HASH} < $maxjob);

		# Check if this shell is completed or failed.
		updata_shell_result($QSUB_INFO, $_);
		next if (exists $QSUB_INFO->{OK_JOBS}->{$_});  # Already ok
		next unless (exists $QSUB_INFO->{COUNT_JOBS}->{$_});  # Cannot run anymore

		if ($QSUB_INFO->{COUNT_JOBS}->{$_} < $nsub) {
			next if (is_running_shell($QSUB_INFO, $RUNNING_JOBS_HASH, $_));  # Already running

			qsub_a_shell($QSUB_INFO, $RUNNING_JOBS_HASH, $_);
			next;
		}	
		delete $QSUB_INFO->{COUNT_JOBS}->{$_};
	}
}

sub is_running_shell {
	my $QSUB_INFO = shift;
	my $RUNNING_JOBS_HASH = shift;
	my $shell = shift;

	update_running_qsub($QSUB_INFO, $RUNNING_JOBS_HASH);

	my @running_jobs = keys %{$RUNNING_JOBS_HASH};
	foreach (@running_jobs) {
		return 1 if ($shell eq $RUNNING_JOBS_HASH->{$_});
	}
	return 0;
}

sub qsub_a_shell {
	my $QSUB_INFO = shift;
	my $RUNNING_JOBS_HASH = shift;
	my $shell = shift;

	my $shfile = $QSUB_INFO->{outdir} . "/" . basename($shell);
	my %info = qsub_job($shfile, %{$QSUB_INFO->{opts}});
	my $job_id = $info{job_id};

	$RUNNING_JOBS_HASH->{$job_id} = $shell;

	push @{$QSUB_INFO->{ALL_JOBS}->{$shell}}, $job_id;
	$QSUB_INFO->{COUNT_JOBS}->{$shell} += 1;
}

sub updata_shell_result {
	my $QSUB_INFO = shift;
	my $shell_name = shift;

	my $ok = &check_shell_complete($QSUB_INFO->{outdir}, $shell_name);
	return unless $ok;

	#print $shell_name . " ok.\n";
	$QSUB_INFO->{OK_JOBS}->{$shell_name} = $QSUB_INFO->{ALL_JOBS}->{$shell_name};
	delete $QSUB_INFO->{COUNT_JOBS}->{$shell_name};
}

# shell output
sub check_shell_complete {
	my $indir = shift;
	my $shell_name = shift;

	my $oe = &find_shell_output($indir, $shell_name);
	my @o = @{$oe->{o}};
	foreach (@o) {
		my $ofile = $indir . "/" . basename($_);
		return 1 if &check_o_complete($ofile);
	}
	return 0;
}

sub check_o_complete {
	my $ofile = shift;

	$ofile = real_path($ofile, 1);

	my $out = 0;
	open my $fh, "<", $ofile or die $!;
	while (<$fh>) {
		chomp;
		#print $_ . "\n";
		if ($_ eq $COMPLETE_STR) {
			$out = 1;
			last;
		}
	}
	return $out;
}

sub find_shell_output {
	my $indir = shift;
	my $shell_name = shift;

	my $target = basename($shell_name);

	my (@e, @o);

	opendir my $dh, $indir or die $!;
	my @all_files = grep { /\.sh\.o[0-9]+$/ || /\.sh\.e[0-9]+$/ } readdir($dh);

	foreach (@all_files) {
		my $file = $_;
		$file =~ s/[0-9]+$//g;
		push(@o, $indir . "/" . $_) if ($file eq $target . ".o");
		push(@e, $indir . "/" . $_) if ($file eq $target . ".e");
	}
	return { o => \@o, e => \@e };
}

# qsub
sub update_qstat_info {
	my $info_hash = shift;

	my $job_id = $info_hash->{job_id};
	die "No 'job_id'." unless defined $job_id;

	my $stat = `qstat -j $job_id 2>&1`;
	if ($stat =~ /do not exist/) {
		$info_hash->{exists} = 0;
		return;
	}
	$info_hash->{exists} = 1;
	my @lines = split /\n/, $stat;
	my $id;
	foreach my $line (@lines) {
		$id = $1 if ($line =~ /^([^:]+):(.*)/);
		next unless defined $id;

		$id =~ s/\s+/_/g;
		$line =~ s/^[^:]+://;
		$line =~ s/^\s*//;
		$info_hash->{$id} .= $line;
	}
}

sub qsub_job {
	my $shfile = shift;
	my %opts = @_;

	$shfile = real_path($shfile, 1);
	my $shdir = dirname($shfile);

	my $qcmd = qsub_cmd(%opts);

	my $job_out = `$qcmd -e $shdir -o $shdir $shfile`;
	my $job_id = $1 if ($job_out =~ /Your job (\d+)/);
	die "Fetching job id failed: \n$job_out" if ($job_id eq '0' || !defined $job_id);

	my %info = (job_id => $job_id);
	return %info;
}

sub qsub_cmd {
	my %opts = @_;

	my %default_opts = &default_options();
	&get_default_options(\%opts, \%default_opts);

	my @keys_required = qw/qsub queue resource ppn/;
	&check_keys_required(\%opts, @keys_required);

	my $cmd = "$opts{qsub} -q $opts{queue} -l $opts{resource} -cwd -S /bin/bash -V -sync no -pe smp $opts{ppn} ";
	return $cmd;
}

# Shell
sub read_shell {
	my $file = shift;
	my $lines = shift;

	$lines //= 999999;

	$file = real_path($file, 1);

	my @cmds;
	my $n_line = 0;
	my $cmd = "";

	open IN, "<", $file or die $!;
	while (<IN>) {
		chomp;
		my $add = &filter_cmd($_);
		next unless (defined $add);

		$cmd .= $add . "\n";
		$n_line++;
		next unless ($n_line == $lines);

		push @cmds, $cmd;
		$n_line = 0;
		$cmd = "";
	}
	close IN;
	push @cmds, $cmd unless $cmd eq '';
	return @cmds;
}

sub split_shell {
	my $file = shift;
	my $outdir = shift;
	my $lines = shift;
	my $pfx = shift;

	$lines //= 999999;
	$pfx //= "work";

	$file = real_path($file, 1);
	
	my $script_pfx = join("\n", "#!/bin/bash", "set -e", "echo \$HOSTNAME");
	my $script_sfx = "echo " . $COMPLETE_STR;

	my @shells;

	my $n_line = 0;
	my $job_id = "0001";
	my $cmd = "";
	my $pattern = "%s/%s_%s.sh";

	open IN, "<", $file or die $!;
	while (<IN>) {
		chomp;
		my $add = &filter_cmd($_);
		next unless (defined $add);
	
		$cmd .= $add . "\n";
		$n_line++;
		next unless ($n_line == $lines);

		# Write $cmd to target
		my $outfile = sprintf($pattern, $outdir, $pfx, $job_id);
		make_path($outdir);
		write_shell($cmd, $outfile);
		
		# Record target shell
		push @shells, basename($outfile);

		# Refresh
		$job_id++;
		$n_line = 0;
		$cmd = "";
	}
	close IN;
	return @shells if ($cmd eq "");

	my $outfile = sprintf($pattern, $outdir, $pfx, $job_id);
	make_path($outdir);
	write_shell($cmd, $outfile);
	push @shells, basename($outfile);

	return @shells;
}

sub write_shell {
	my $cmd = shift;
	my $outfile = shift;
	my $for_qsub = shift;

	$for_qsub //= 1;

	my $script_pfx = join("\n", "#!/bin/bash", "set -e", "echo \$HOSTNAME");
	my $script_sfx = "";

	$script_pfx .= "\necho \$HOSTNAME" if $for_qsub;
	$script_sfx = "echo " . $COMPLETE_STR if $for_qsub;
	
	open OUT, '>', $outfile or stop($!);
	print OUT join("\n", $script_pfx, $cmd, $script_sfx) . "\n";
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
	my $path = shift;
	my $must_work = shift;

	$must_work //= 0;

	my $norm_path = abs_path($path);
	
	return $norm_path if defined $norm_path;

	$path =~ s/\/+$//g;
	return $path unless $must_work;

	stop('No such file or directory: ', $path);
}

# Options
sub check_keys_required {
	my $opts = shift;
	my @keys_arr = @_;
	foreach (@keys_arr) {
		next unless (!defined $opts->{$_} || $opts->{$_} eq '');
		stop("Required key '$_' not found in the input hash.");
	}
}

sub get_default_options {
	my $opts = shift;
	my $default_opts = shift;
	foreach (keys %{$default_opts}) {
		$opts->{$_} //= $default_opts->{$_};
	}
}

# DEBUG
sub stop {
	confess "ERROR:\n @_\n";
}

