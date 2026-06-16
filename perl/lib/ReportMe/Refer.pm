package ReportMe::Refer;

use strict;
use warnings;

use Eula::Utils qw/ 
    require_opt 
    stop 
/;

require Exporter;
our @ISA = qw(Exporter);

sub new {
    my ($class, %opts) = @_;
    my $self = {};
    bless $self, $class;

    $self->{names} = [];
    $self->{refers} = {};

    return $self;
}

sub add_refer {
    my ($self, $id, $name, $refer) = @_;
    return if (exists($self->{refers}->{$id}));

    push @{$self->{names}}, $id;
    $self->{refers}->{$id}->{name} = $name;
    $self->{refers}->{$id}->{refer} = $refer;
}

sub get_refer {
    my ($self, $id) = @_;
    $self->_check_ref($id);
    return $self->{refers}->{$id}->{refer};
}

sub get_refer_name {
    my ($self, $id) = @_;
    $self->_check_ref($id);
    return $self->{refers}->{$id}->{name};
}

sub _check_ref {
    my ($self, $id) = @_;
    stop("id is not found") unless (exists($self->{refers}->{$id}));
}

sub refer_order {
    my ($self, $id) = @_;
    $self->_check_ref($id);

    my $i = 0;
    for my $n (@{$self->{names}}) {
        $i += 1;
        return $i if ($n eq $id);
    }
}

sub get_refer_html {
    my ($self, $id) = @_;
    $self->_check_ref($id);
    my $order = $self->refer_order($id);
    return qq(<sup>[<a href="#cite:$id" id="recite:$id">$order</a>]</sup>);
}

sub dump {
    my ($self) = @_;
    
    my $html = "<p><ul>";
    for my $id (@{$self->{names}}) {
        my $refer = $self->get_refer($id);
        my $order = $self->refer_order($id);
        $html .= qq([<a href="#recite:$id" id="recite:$id">$order</a>] $refer\n);
    }
    $html .= "</ul></p>";
    $html .= '<script type="text/javascript">
			function close_help(){
				var f = document.getElementById("show_help");
				var b = document.getElementById("bgbox");
				f.style.display = "none";
				document.body.removeChild(b);
			}</script>';
	return $html;
}

1;
