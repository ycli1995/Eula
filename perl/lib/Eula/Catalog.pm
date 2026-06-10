package Eula::Catalog;

use strict;
use warnings;

use Eula::Utils::DEBUG qw/require_opt/;

require Exporter;
our @ISA = qw(Exporter);

my %scaffold = (
    empty => '    ',
    long => '│   ',
    branch => '├── ',
    ending => '└── ',
);

sub new {
    my ($class, $parent, $opts) = @_;

    my $self = {};
    $self->{parent} = $parent || undef;

    $self->{file} = "";
    $self->{name} = "";
    if (defined $self->{parent}) {
        require_opt($opts, 'file');
        require_opt($opts, 'name');
        $self->{file} = $opts->{file};
        $self->{name} = $opts->{name};
    }
    $self->{subtrees} = $opts->{subtrees} || [];
    $self->{branch} = "";  # add to self
    $self->{pref} = "";  # add to subtree
    $self->{is_last} = 1;

    bless $self, $class;
    return $self;
}

sub subtree {
    my ($self, $file, $name) = @_;
    my $tree = __PACKAGE__->new($self, { file => $file, name => $name });
    
    if (@{$self->{subtrees}}) {
        my $last_subtree = $#{$self->{subtrees}};
        my $last_tree = $self->{subtrees}->[$last_subtree];
        $last_tree->{is_last} = 0;
    }
    push @{$self->{subtrees}}, $tree;
    
    return $tree;
}

sub _extend_branch {
    my ($self) = @_;
    
    if (defined $self->{parent}) {
        $self->{branch} = $self->{is_last} ? $scaffold{ending} : $scaffold{branch};
    }
    for my $tree (@{$self->{subtrees}}) {
        $tree->_extend_branch();
        $tree->{pref} = $self->{pref};
        if (@{$tree->{subtrees}}) {
            $tree->{pref} .= $tree->{is_last} ? $scaffold{empty} : $scaffold{long};
        }
    }
}

sub dump_catalog {
    my ($self) = @_;
    $self->_extend_branch();
    my $out = "";
    for my $i (0 .. $#{$self->{subtrees}}) {
        my $tree = $self->{subtrees}->[$i];
        $out .= $tree->{parent}->{pref} . $tree->{branch} . $tree->{file} . "    " . $tree->{name} . "\n";
        $out .= $tree->dump_catalog();
    }
    return $out;
}

1;
