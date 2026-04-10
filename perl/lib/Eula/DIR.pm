package Eula::DIR;

use strict;
use warnings;

use Cwd qw/getcwd/;
use File::Basename qw/basename/;
use File::Path qw/make_path/;

use Eula::Utils::DEBUG qw/stop/;

our $Ordered = 0;
our $Verbose  = 1;
our $Mkdir = 1;
our $TraversalMethod = 'BFS';

my $StartFrom = 1;

sub new {
	my ($class, $opts) = @_;
	my $self = {};
	bless $self, $class;
	$self->init($opts);
	return $self;
}

sub init {
	my ($self, $opts) = @_;

	$opts->{wkdir} //= getcwd();
	$opts->{name} //= '';

	$opts->{wkdir} = '/' if ($opts->{name} =~ /^\//);

	my $path = "$opts->{wkdir}/$opts->{name}";

	$self->{_name} = basename($path);
	$self->{_path} = $path;
	$self->{_parent} = undef;
	$self->{_child_names} = [];
	$self->{_child} = {};
	$self->{_index} = '';
	$self->{_next_index} = $StartFrom;
}

sub parent {
	my $self = shift;
	return $self->{_parent};
}

sub mkdir {
	my ($self, $opts) = @_;

	$opts->{recursive} //= 0;

	make_path($self->{_path});
	return unless ($opts->{recursive});
	
	foreach ($self->{_child_names}) {
		$self->{_child}->{$_}->mkdir({recursive => 1});
	}
}

sub add {
	my ($self, $names, $opts) = @_;

	$opts->{ordered} //= $Ordered;
	$opts->{mkdir} //= $Mkdir;

	$names = [ $names ] unless (ref $names);
	stop('$names must be a ref for ARRAY.') unless (ref $names eq 'ARRAY');

	foreach (@{$names}) {
		next if (exists $self->{_child}->{$_});

		my $next_idx = '';
		if ($opts->{ordered}) {
			$next_idx = $self->_next_index() . ".";
			$self->{_next_index} += 1;
		}

		my $child = {
			_name => $_,
			_parent => $self,
			_child_names => [],
			_child => {},
			_index => $next_idx,
			_next_index => $StartFrom,
		};
		$child->{_path} = "$self->{_path}/$child->{_index}$child->{_name}";
		$child = bless $child, ref $self;

		$child->mkdir() if ($opts->{mkdir});

		$self->{_child}->{$_} = $child;
		push @{$self->{_child_names}}, $_;
	}
}

sub _next_index {
	my $self = shift;
	return $self->{_next_index};
}

sub get {
	my ($self, $name) = @_;
	return $self->{_child}->{$name};
}

sub path {
	my ($self, $name, $opts) = @_;
	return $self->{_path} unless (defined $name);

	my $child = $self->fd($name, $opts);
	return $child->{_path};
}

sub fd {
	my ($self, $name, $opts) = @_;

	return $self if ($name eq '.');

	if ($name eq '..') {
		my $parent = $self->{_parent};
		return $parent if (defined $parent);
		stop("No parent for searching '..'. This is the root DIR.");
	}

	$opts->{method} //= $TraversalMethod;

	my $child = undef;
	if ($opts->{method} eq 'BFS') {
		$child = $self->_fd_BFS($name);
	} else {
		$child = $self->_fd_DFS($name);
	}
	stop("'$name' is not found in DIR: $self->{_path}") unless (defined $child);
	return $child;
}

sub _fd_BFS {
	my ($self, $name, $queue) = @_;

	$queue //= [];
	my $out = undef;
	foreach (@{$self->{_child_names}}) {
		my $child = $self->{_child}->{$_};
		if ($_ eq $name) {
			$out = $child;
			last;
		}
		next unless (keys %{$child->{_child}});
		push @{$queue}, $child;
	}
	return $out if (defined $out);

	while (@{$queue}) {
		my $child = shift @{$queue};
		my $new_child = $child->_fd_BFS($name, $queue);
		return $new_child if (defined $new_child);
	}
	return $out
}

sub _fd_DFS {
	my ($self, $name) = @_;
	
	my $out = undef;
	foreach (@{$self->{_child_names}}) {
		$out = $self->{_child}->{$_};
		return $out if ($_ eq $name);
		
		$out = $out->_fd_DFS($name);
		return $out if (defined $out);
	}
	return $out;
}

1;
