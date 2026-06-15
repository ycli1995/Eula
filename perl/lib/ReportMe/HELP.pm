package ReportMe::HELP;

use strict;
use warnings;

use File::Basename qw/basename/;

use Eula::Utils::DEBUG qw/check_ref/;

use ReportMe::HTML;

require Exporter;
our @ISA = qw(Exporter ReportMe::HTML);

sub new {
    my ($class, $parent, %opts) = @_;

    stop("<parent> is required") unless ($parent);

    my $self = ReportMe::HTML->new($parent, %opts);
    bless $self, $class;

    $self->{path}->{outdir} = $self->{parent}->{path}->{outdir} . "/src/doc/";
    $self->{path}->{name} = "help";
    $self->{nonlazy} = 1;
    return $self;
}

sub _add_help_html {
    my ($self, $id, $title, $html) = @_;
    my $add_html = <<HTML;
    <div name="$id">
        <h2>$title</h2>
        $html
    </div>
HTML
    $self->add_html($add_html);
}

sub add_help_desc {
    my ($self, $id, $title, $desc) = @_;
    $self->_add_help_html($id, $title, qq(<p>$desc</p>));
}

sub add_help_tab {
    my ($self, $id, $title, $tab_arr) = @_;
    check_ref($tab_arr, "ARRAY");

    my $html = qq(<table style="width:900px" border="1" border-color="#b6ff00">\n);
    for my $row (@$tab_arr) {
        my @tmp_row;
        @tmp_row = ref($row) eq "ARRAY" ? @$row : split(/\t+/, $row);
        $row = join("</td><td>", @tmp_row);
        $html .= qq(<tr><td>$row</td></tr>\n);
    }
    $html .= "</table>\n";
    $self->_add_help_html($id, $title, $html);
}

sub main_html_head {
    my $self = shift;

    my $parent = $self->{parent};
    foreach (qw/report.css help_page.css help_document.css/) {
        $parent->require_src("css", $_);
    }
    my $title = $parent->require_content_config("HELP_TITLE");

    my $head = <<HTML;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<!-- Basic Information -->
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<title>$title</title>        
<!-- CSS Document -->
<link rel="stylesheet" type="text/css" href="../css/report.css" />
<link rel="stylesheet" type="text/css" href="../css/help_page.css" />
<link rel="stylesheet" type="text/css" href="../css/help_document.css" />
</head>
<body>
<a name="help"></a>
HTML
    return $head;
}

sub main_html_tail {
    my ($self, %opts) = @_;

    my $tail = <<HTML;
<br/><br/><br/><br/><br/><br/>
</div>
</body>
</html>
HTML

    return $tail;
}

1;
