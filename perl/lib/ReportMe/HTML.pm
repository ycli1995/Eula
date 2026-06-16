package ReportMe::HTML;

use strict;
use warnings;

use File::Basename qw/basename dirname/;

use lib dirname(__FILE__)."/../../vendor/lib/perl5";
use JSON;

use Eula::Utils qw/
    check_ref 
    cp_paths 
    real_path 
    require_opt 
    rm_paths 
    stop 
    time_log
    write_lines
/;

use ReportMe::HELP;

require Exporter;
our @ISA = qw(Exporter);

sub new {
    my ($class, $parent, %opts) = @_;

    $opts{'-class'} ||= "normal_content";
    $opts{'-tag'} ||= "div";

    my $self = {};
    bless $self, $class;

    $parent ? $self->{parent} = $parent : $self->_new_root_section(\%opts);

    my $attrs = opts2attrs(%opts);
    $self->{head} = "<$opts{'-tag'} $attrs>\n";
    $self->{tail} = "</$opts{'-tag'}>\n";
    $self->{main} = "";

    return $self;
}

sub _new_root_section {
    my ($self, $opts) = @_;

    time_log("Create new root section ... ");

    # path info
    my $outdir = $opts->{'-outdir'} || ".";
    delete $opts->{'-outdir'};
    $self->{path}->{outdir} = $outdir;

    my $name = $opts->{'-name'} || "index";
    delete $opts->{'-name'};
    $self->{path}->{name} = $name;

    # content config
    my $content_config = $opts->{'-content_config'} || {};
    delete $opts->{'-content_config'};
    
    $self->set_content_config($content_config);

    # the public vars
    $self->{resp_htabs_cnt} = 1;  # for parent horizontal tabs
    $self->{resp_vtabs_cnt} = 1;  # for parent vertical tabs
    $self->{child_vtabs_cnt} = 1; # for child tabs in parent horizontal tab
    
    $self->{menu_cnt} = 0;
    $self->{submenu_cnt} = 0;
    $self->{ssubmenu_cnt} = 0;

    $self->{section} = [];
    $self->{split_line} = "";

    $self->{nonlazy} = $opts->{'-nonlazy'} || 0;
    $self->{getPic} = "";
    $self->{json} = "";

    $self->_set_default_src();
}

sub show_info {
    my ($self) = @_;

    my $me = $self->{parent} ? $self->{parent} : $self;

    # path
    print "path:\n";
    for my $p (keys %{$me->{path}}) {
        print " $p: $me->{path}->{$p}\n";
    }

    # src
    print "src:\n";
    for my $src (keys %{$me->{src}}) {
        for my $file (keys %{$me->{src}->{$src}}) {
            my $infile = $me->require_src($src, $file);
            print " $src/$file: $infile\n";
        }
    }

    # content config
    print "content config:\n";
    for my $key (keys %{$me->get_content_config()}) {
        my $val = $me->require_content_config($key);
        print " $key: $val\n";
    }

    # counts
    for my $key (qw/menu_cnt submenu_cnt ssubmenu_cnt/) {
        print "$key: $me->{$key}\n";
    }
}

sub get_content_config {
    my ($self) = @_;
    my $content_config = $self->{content_config};
    return $content_config if ($content_config);
    return {} unless ($self->{parent});
    return $self->{parent}->get_content_config();
}

sub set_content_config {
    my ($self, $content_config) = @_;
    $content_config ||= {};
    
    $content_config->{TITLE} ||= "Analysis";
    $content_config->{CONTENT_TITLE} ||= "Contents";
    $content_config->{CLOSE} ||= "Close";
    $content_config->{OPEN} ||= "Open";
    $content_config->{GO_TO_TOP} ||= "Back to top";
    $content_config->{HELP_TITLE} ||= "Help";
    $content_config->{URL} ||= "";
    $content_config->{SEARCH_PLACEHOLDER} ||= "Search";

    $self->{content_config} = $content_config;
}

sub require_content_config {
    my ($self, $key) = @_;
    my $content_config = $self->get_content_config();
    if (!exists $content_config->{$key}) {
        stop("content config $key is missing.");
    }
    return $content_config->{$key};
}

sub get_src {
    my ($self, $type, $key) = @_;
    return $self->{src}->{$type}->{$key} || "";
}

sub set_src {
    my ($self, $type, $key, $file) = @_;
    $file = real_path($file, {must_work => 1});
    $self->{src}->{$type}->{$key} = $file;
}

sub require_src {
    my ($self, $type, $key) = @_;
    if (!exists $self->{src}->{$type}->{$key}) {
        stop("src $type/$key is missing.");
    }
    my $file = $self->{src}->{$type}->{$key};
    $file = real_path($file, {must_work => 1});
    return $file;
}

my %JS = (
    HEAD => [ qw/
        jquery-1.9.1-min.js
        modernizr-min.js
        jquery.jumpto.js
        toggle.js
        jquery.nicescroll-min.js
        easyResponsiveTabs-min.js
        show_help-min.js
        jquery.dataTables.min.js
    / ],
    TAIL => [ qw/report_init.js/ ],
    LAZY_LOAD => [ qw/lazy_load_pre.js lazy_load_main.js/ ],
);

my %CSS = (
    HEAD => [ qw/
        index.css
        jquery.dataTables.min.css
        report.css
        jumpto.css
        easy-responsive-tabs.css
        toggle.css
    / ],
    HELP => [ qw/
        help_page.css
        help_document.css
    / ],
);

my %IMAGE = (
    TAIL => [ qw/goTop.jpg/ ],
    CHECK => [ qw/no.png yes.png warn.png/ ],
    HELP => [ qw/help.png/ ],
);

sub get_src_keys {
    my ($self, $type) = @_;
    return %IMAGE if ($type eq "image");
    return %CSS if ($type eq "css");
    return %JS if ($type eq "js");
    stop("src type $type is invalid. Must be image, css, or js.");
}

sub _set_default_src {
    my ($self) = @_;
    my $root = _get_root_dir();
    for my $key (keys %IMAGE) {
        for my $file (@{$IMAGE{$key}}) {
            $self->set_src("image", $file, "$root/image/$file");
        }
    }
    for my $key (keys %CSS) {
        for my $file (@{$CSS{$key}}) {
            $self->set_src("css", $file, "$root/css/$file");
        }
    }
    for my $key (keys %JS) {
        for my $file (@{$JS{$key}}) {
            $self->set_src("js", $file, "$root/js/$file");
        }
    }
}

sub _get_root_dir {
    my $root = dirname(__FILE__) . "/src";
    $root = real_path($root, {must_work => 1});
    return $root;
}

sub _cp_src {
    my ($self) = @_;
    my $outdir = $self->{path}->{outdir};

    time_log("Fetch the 'src' start ... ");

    mkdir $outdir unless (-d $outdir);
    $outdir = real_path($outdir, {must_work => 1});
    
    mkdir "$outdir/src" unless (-d "$outdir/src");
    rm_paths("$outdir/src/*");

    for my $src (keys %{$self->{src}}) {
        for my $file (keys %{$self->{src}->{$src}}) {
            my $infile = $self->require_src($src, $file);
            print " $file: $infile\n";
            cp_paths($infile, "$outdir/src/$src/$file");
        }
    }
    time_log("Fetch the 'src' done ... ");
}

sub get_split_line {
    my ($self) = @_;
    my $split_line = $self->{split_line};
    return $split_line if ($split_line);
    return "" unless ($self->{parent});
    return $self->{parent}->get_split_line();
}

sub set_split_line {
    my ($self, $split_line) = @_;
    $self->{split_line} = $split_line;
}

# return the object of HTML with tag 'section'
sub _section {
    my ($class, $parent, %opts) = @_;

    $opts{'class'} ||= "normal_content";
    $opts{'-tag'} = "section";

    my $self = __PACKAGE__->new($parent, %opts);

    $self->{nonlazy} = $parent->{nonlazy};
    $self->{src} = $parent->{src};

    my $split_line = $self->get_split_line();
    if ($opts{'-break'} || $opts{'-split'}) {
        $self->{head} = qq(<div style="page-break-after:always;"></div>\n$split_line\n$self->{head}\n);
    } elsif ($opts{'-page_head'}) {
        $self->{head} = qq($split_line\n$self->{head}\n);
    }
    return $self;
}

sub section {
    my ($self, %opts) = @_;
    my $section = __PACKAGE__->_section($self, %opts);
    return $section;
}

sub add_section {
    my ($self, $section) = @_;
    push @{$self->{section}}, $section;
    return $section;
}

sub help {
    my ($self, %opts) = @_;
    my $help = ReportMe::HELP->new($self);
    $self->{help} = $help;
    return $help;
}

# return the inner code of HTML object
sub innerHTML {
    my $self = shift;
    return $self->{head} . $self->{main} . $self->{tail};
}

# add the html code to HTML object
sub add_html {
    my ($self, $str) = @_;
    $self->{main} .= $str;
}

sub main_html_head {
    my ($self) = @_;

    my $title = $self->require_content_config("TITLE");
    my $close = $self->require_content_config("CLOSE");
    my $open = $self->require_content_config("OPEN");

    # fetch used files
    my @css = map {
        $self->require_src("css", $_);
        qq(<link rel="stylesheet" type="text/css" href="src/css/$_" />);
    } @{$CSS{HEAD}};
    my $css_html = join "\n", @css;

    my @js = map {
        $self->require_src("js", $_);
        qq(<script src="src/js/$_"></script>);
    } @{$JS{HEAD}};
    my $js_html = join "\n", @js;

    my $head = <<HTML;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<!-- Basic Information -->
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">        
<!-- CSS Document -->
$css_html
<!-- JS Script -->
$js_html
</head>
<body>
    <section>
        <div id="header_banner">
        <div id="banner_logo"></div>
        <div id="banner_title">$title</div>
        <div id="banner_bg_image"></div>
        </div>
    </section>
    <div class="toggleNav">
        <span class="fold1">$close</span>
        <span class="fold2 close">$open</span>
    </div>
    <div id="report_body">
HTML
    return $head;
}

sub main_html_tail {
    my ($self) = @_;

    my $content_config = $self->get_content_config();

    my $url = $self->require_content_config("URL");
    my $go_to_top = $self->require_content_config("GO_TO_TOP");
    my $help_page_title = $self->require_content_config("HELP_TITLE");
    my $content_title = $self->require_content_config("CONTENT_TITLE");

    $self->require_src("image", "goTop.jpg");
    my @js = map {
        $self->require_src("js", $_);
        qq(<script src="src/js/$_"></script>);
    } @{$JS{TAIL}};
    my $js_html = join "\n", @js;

    my $lazy_load = $self->lazy_load();
    my $config = encode_json({
        content_title => $content_title,
        resp_vtabs_count => $self->{resp_vtabs_cnt},
        resp_htabs_count => $self->{resp_htabs_cnt},
        child_vtabs_count => $self->{child_vtabs_cnt},
    });
    
    my $tail = <<HTML;
</div>    
<div id="goTop" style="display:none;">
    <a title="${go_to_top}" class="backtotop">
        <img class="back-top" src="src/image/goTop.jpg">
    </a>
</div>    
<!-- Help Window -->
<div id="show_help">
    <h3>${help_page_title}</h3>
    <iframe id="help_page" name="help_page" src="$url"></iframe>
</div>
<!-- JS Plugin Initialization -->
<script>window.REPORT_CONFIG = $config;</script>
$js_html
$lazy_load
</body>
</html>
HTML

    return $tail;
}

sub _new_menu {
    my ($self, $str, $order, $type, %opts) = @_;

    if (defined $opts{'-check'}) {
        my $check_icon = $IMAGE{CHECK}[$opts{'-check'}] || "";
        if (!$check_icon eq "") {
            $self->require_src("image", $check_icon);
            $str .= qq(<img class="check_icon" src="src/image/$check_icon"> );
        }
    }

    my $attrs = opts2attrs(%opts);
    my $help = $self->opts2help(%opts);
    return qq(<$type $attrs>$order $str$help</$type>\n);
}

# add the <h3> to HTML object
sub menu {
    my ($self, $str, %opts) = @_;
    
    my $parent = $self->{parent};
    $parent->{menu_cnt} ++;
    $parent->{submenu_cnt} = 1;
    $parent->{img_cnt} = 0;
    $parent->{tab_cnt} = 0;

    my $order = $parent->{menu_cnt};
    my $html = $self->_new_menu($str, $order, "h3", %opts);
    $opts{'-return'} ? return $html : $self->add_html($html);
}

# add the <h5> to HTML object
sub submenu {
    my ($self, $str, %opts) = @_;

    my $parent = $self->{parent};
    $parent->{img_cnt} = 0;
    $parent->{tab_cnt} = 0;

    my $order = join ".", @{$parent}{qw/menu_cnt submenu_cnt/};
    $parent->{submenu_cnt} ++;
    $parent->{ssubmenu_cnt} = 1;
    
    my $html = $self->_new_menu($str, $order, "h5", %opts);
    $opts{'-return'} ? return $html : $self->add_html($html);
}

# add the <h4> to HTML object
sub ssubmenu {
    my ($self, $str, %opts) = @_;

    my $parent = $self->{parent};
    my $submenu_cnt = $parent->{submenu_cnt} - 1;
    my $menu_cnt = $parent->{menu_cnt};
    my $ssubmenu_cnt = $parent->{ssubmenu_cnt};
    my $order = join ".", ($menu_cnt, $submenu_cnt, $ssubmenu_cnt);

    $parent->{ssubmenu_cnt}++;
    
    my $html = $self->_new_menu($str, $order, "h4", %opts);
    $opts{'-return'} ? return $html : $self->add_html($html);
}

sub which_menu {
    my $self = shift;
    my $parent = $self->{parent};
    return 0 if ($parent->{menu_cnt} == 0 || $self->{main} eq ''); ## pre menu
    return 1 if (exists $parent->{submenu_cnt} && $parent->{submenu_cnt} == 1); ## menu
    return 2 if (exists $parent->{ssubmenu_cnt} && $parent->{ssubmenu_cnt} == 1); ## submenu
    return 3; ## ssubmenu
}

## add menu to next level 
sub next_level_menu {
    my $self = shift;
    my $level = $self->which_menu();
    stop("<ssubmenu> is the last level menu already!") if ($level == 3);
    $self->menu(@_) if ($level == 0);
    $self->submenu(@_) if ($level == 1);
    $self->ssubmenu(@_) if ($level == 2);
}

sub write {
    my ($self, %opts) = @_;
    
    my $outdir = $self->{path}->{outdir};
    my $name = $self->{path}->{name} || "index";

    time_log("The html '$name.html' will be create :)");
    
    my $head = $self->main_html_head();
    my $tail = $self->main_html_tail();
    my $html = $self->{main};
    foreach my $section (@{$self->{section}}) {
        $html .= $section->innerHTML;
    }

    delete $self->{parent};
    $self->_cp_src();
    write_lines("$outdir/$name.html", [$head, $html, $tail]);

    $self->{help}->write() if ($self->{help});
    unless ($self->{nonlazy}) {
        open JSON , ">" , "$outdir/src/js/pic.js" or die $!;
        print JSON "getPic ({$self->{json}})";
        close JSON;
    }

    time_log("The html report '$name.html' was create done :)");
}

# pack the report folder, this function will use the pack cmd of system, so just for Linux now.
# support pack format: gz, bz2, zip, rar
sub pack {
    my ($self, %opts) = @_;
    
    my $outdir = $self->{path}->{outdir};
    $outdir =~ s/\/$//;

    my $out_dir = dirname($outdir);
    my $out_name = basename($outdir);

    my $format = $opts{'-format'} || "bz2";
    
    my $cmd = "";
    if ($format eq "bz2") {
        $cmd = "cd $out_dir; tar -h -cjf $out_name.tar.bz2 $out_name";
    } elsif ($format eq "gz") {
        $cmd = "cd $out_dir; tar -h -czf $out_name.tar.gz $out_name";
    } elsif ($format eq "zip") {
        $cmd = "cd $out_dir; zip -r $out_name.zip $out_name";
    } elsif ($format eq "rar") {
        $cmd = "cd $out_dir; rar a -r $out_name.rar $out_name";
    } else {
        stop("pack format stop, now just support format: gz, bz2, zip, rar.",$format);
    }
    system("$cmd 1>/dev/null");
    time_log("The html report has been packed :)");
}

# add a page break
sub break {
    my $self = shift;
    my $split_line = $self->get_split_line();
    my $html = qq(<div style="page-break-after:always;"></div>\n$split_line\n);
    $self->add_html($html);
}

# add <p> to HTML object
sub desc {
    my ($self, $str, %opts) = @_;

    my $attrs = opts2attrs(%opts);
    my $help = $self->opts2help(%opts);
    
    $str = "<pre>$str</pre>" if ($opts{'-pre'});

    my $html = qq(<p $attrs>$str$help</p>\n);
    $opts{'-return'} ? return $html : $self->add_html($html);
}

sub img_order {
    my $self = shift;
    $self->{parent}->{img_cnt}++;

    my $submenu_cnt = $self->{parent}->{submenu_cnt} - 1;
    my $menu_cnt = $self->{parent}->{menu_cnt};
    my $img_cnt = $self->{parent}->{img_cnt};
    my $order = join ".", ($menu_cnt, $submenu_cnt, $img_cnt);
    return qq(Fig <span style="color:red">$order</span>);
}

sub tab_order {
    my $self = shift;
    $self->{parent}->{tab_cnt}++;
    my $submenu_cnt = $self->{parent}->{submenu_cnt} - 1;
    my $menu_cnt = $self->{parent}->{menu_cnt};
    my $tab_cnt = $self->{parent}->{tab_cnt};
    my $order = join ".", ($menu_cnt, $submenu_cnt, $tab_cnt);
    return qq(Tab <span style="color:red">$order</span>);
}

# generate the title of figure or table
sub title {
    my ($self, %opts) = @_;
    
    # defined the img/table title and its description
    my $name = $opts{'-name'} || "";
    return $name unless ($opts{'-type'});

    my $note = $opts{'-note'} || "";
    my $help = $self->opts2help(%opts);
    my $order = $opts{'-type'} eq "image" ? $self->img_order() : $self->tab_order();

    my $note_html = "";
    if ($note) {
        $note_html = qq(<p class="img_note">$note</p><br /><br />);
    }

    $name = <<TEMP;
<div>
    <p class="img_title">$order $name$help</p>
    $note_html
</div>
TEMP

    return $name;
}

# turn files list to <ul> and <li>, with file description (must be defined)
sub files2list {
    my ($self, %opts) = @_;
    
    my $files = require_opt(\%opts, '-files');
    my $desc = require_opt(\%opts, '-desc');
    
    my $short_dir = $opts{'-short_dir'} || 0 ;

    @$files = map { s/^\.\./\./r  } @$files;
    my @fnames = map { $short_dir ? basename($_) : s/^\.\.//r; } @$files;

    my @lies = map { 
        qq(<li>$desc->[$_]: <a href="$files->[$_]" target="_blank">$fnames[$_]</a></li>) 
    } 0 .. $#$files;
    my $li = join "\n" , @lies;
    
    my $attrs = opts2attrs(%opts);
    my $html = "<ul $attrs>\n$li\n</ul>\n";

    $opts{'-return'} ? return $html : $self->add_html($html);
}

sub matrix2html {
    my ($self, %opts) = @_;

    my $matrix = $opts{'-matrix'} or stop("-matrix must be defined in function <matrix2html>");

    my $top = $opts{top} ? $opts{top} : $opts{'-top'} ? $opts{'-top'} : 10;
    my $header = defined $opts{'-header'} ? $opts{'-header'} : 1;
    my $class_name = $opts{'-class'} || "func_table";
    my $width = $opts{'-width'} || "100%";
    my $max_chars = $opts{'-max_chars'} || 0;
    my $skip = $opts{'-skip'} || 0;
    
    my $caption = $self->_tsv_caption(%opts);

    my $tab = "";
    map { shift @$matrix } 1 .. $skip if ($skip);
    if ($header) {
        my $header_row = shift @$matrix;
        $tab .= _format_data_row(1, $max_chars, $opts{'-omits'}, @$header_row);
    }
    
    $tab .= "<tbody>\n";
    my $i = 0;
    foreach (@$matrix) {
        $i ++;
        last if $i > $top;
        $tab .= _format_data_row(0, $max_chars, $opts{'-omits'}, @$_);
    }
    $tab .= "</tbody>\n";

    my $html = _format_table_html($class_name, $width, $caption, $tab);
    
    $opts{'-return'} ? return $html : $self->add_html($html);
}

# fetch the top lines of table, and turn it to html table
sub tsv2html {
    my ($self, %opts) = @_;

    my $file = require_opt(\%opts, '-file');
    
    my $top = $opts{top} ? $opts{top} : $opts{'-top'} ? $opts{'-top'} : 10;
    my $header = defined $opts{'-header'} ? $opts{'-header'} : 1;
    my $class_name = $opts{'-class'} || "func_table";
    my $width = $opts{'-width'} || "100%";
    my $max_chars = $opts{'-max_chars'} || 0;
    my $skip = $opts{'-skip'} || 0;
    
    my $caption = $self->_tsv_caption(%opts);

    my $tab = "";
    my $ncol = 0;
    open my $fh, '<', $file or die "$file $!";
    if ($skip) {
        map { my $tmp = <$fh>; }  1 .. $skip;
    }
    if ($header) {
        my $header_line = <$fh>;
        stop("$file is empty, please check it !") unless ($header_line);
        my @values = split /\t/, $header_line;
        chomp $values[-1];
        $ncol = scalar @values;
        $tab .= _format_data_row(1, $max_chars, $opts{'-omits'}, @values);
    }
    
    $tab .= "<tbody>\n";
    my $i = 0;
    while (<$fh>) {
        next if (/^#/);

        $i ++;
        last if $i > $top;

        my @values = $opts{'-force'} ? split /\t/, $_, $ncol : split /\t/;
        chomp $values[$#values];
        
        if ($ncol && $ncol != $#values + 1) {
            stop("the number of columns is not equals to the header in LINE $. ", $file);
        } elsif (!$ncol) {
            $ncol = scalar @values;
        }
        $tab .= _format_data_row(0, $max_chars, $opts{'-omits'}, @values);
    }
    close $fh;
    $tab .= "</tbody>\n";

    my $html = _format_table_html($class_name, $width, $caption, $tab);
    
    $opts{'-return'} ? return $html : $self->add_html($html);
}

sub _tsv_caption {
    my ($self, %opts) = @_;
    return "" if ($opts{'-no_order'});

    my $name = require_opt(\%opts, '-name');
    my $help = $self->opts2help(%opts);
    my $order = $self->tab_order();
    return "<caption>$order $name$help</caption>";
}

sub _format_data_row {
    my ($is_header, $max_chars, $omits, @values) = @_;
    
    my $tag1 = $is_header ? "th" : "td";
    my $tag2 = $is_header ? "<thead><tr>" : "<tr>";
    my $tag2_end = $is_header ? "</tr></thead>" : "</tr>";

    @values = &omits_some_columns($omits, @values);
    @values = map {
        my $len = length $_;
        if ($max_chars && $len - 4 > $max_chars) {
            qq(<$tag1 class="abbrTab" data=$max_chars>$_</$tag1>);
        } else {
            "<$tag1>$_</$tag1>";
        }
    } @values;
    
    return $tag2 . join("", @values) . $tag2_end . "\n";
}

sub _format_table_html {
    my ($class_name, $width, $caption, $tab) = @_;
    
    return <<HTML;
<table class="$class_name nowrap" width="$width">
    $caption
    $tab
</table>
<br />
HTML
}

sub img2html {
    my ($self, %opts) = @_;

    my $file = require_opt(\%opts, '-file');
    my $name = require_opt(\%opts, '-name');

    $opts{'-type'} = "image";
    my $title = $self->title(%opts);

    my $dir = _clean_image_path($file);
    my $width = $opts{'-width'} || "auto";
    my $src = $self->_get_src_attribute();
    my $img_link = _generate_img_link($dir, $src, $width eq "auto" ? undef : $width);

    my $img_html = qq(<td>$img_link</td>);
    if ($opts{'-desc'}) {
        $img_html = <<IMG;
<td style="width: $width">$img_link</td>
<td class="pic_table_desc" style="width: 50%"><p>$opts{'-desc'}</p></td>
IMG
    }

    my $html = <<HTML;
<table class="pic_table">
    <tr>
        $img_html
    </tr>
    <tr>
        <td class="img_title">$title</td>
        <td></td>
    </tr>
</table>
HTML
    
    $opts{'-return'} ? return $html : $self->add_html($html);
}

# add two images to the html
sub img2html2 {
    my ($self, %opts) = @_;

    my $file1 = require_opt(\%opts, '-file1');
    my $name1 = require_opt(\%opts, '-name1');
    my $file2 = require_opt(\%opts, '-file2');
    my $name2 = require_opt(\%opts, '-name2');

    $file1 = _clean_image_path($file1);
    $file2 = _clean_image_path($file2);
    
    my $desc1 = $opts{'-desc1'} || "";
    my $desc2 = $opts{'-desc2'} || "";
    
    my $help1 = $opts{'-help1'} ? $self->opts2help('-help'=>$opts{'-help1'}) : "";
    my $help2 = $opts{'-help2'} ? $self->opts2help('-help'=>$opts{'-help2'}) : "";
    
    my $space = $opts{'-space'} || 10;
    my $width = (100 - $space) / 2;
    
    my $src = $self->_get_src_attribute();
    my $title1 = $self->_generate_image_title($name1, $help1);
    my $title2 = $self->_generate_image_title($name2, $help2);
    my $img_link1 = _generate_img_link($file1, $src);
    my $img_link2 = _generate_img_link($file2, $src);

    my $html = <<HTML;
<table class="pic_table">
    <tr>
        <td style="width: $width%">$img_link1</td>
        <td style="width: $space%"></td>
        <td style="width: $width%">$img_link2</td>
    </tr>
    <tr>
        <td class="img_title">$title1</td>
        <td style="width: $space%"></td>
        <td class="img_title">$title2</td>
    </tr>
    <tr>
        <td align="left">$desc1</td>
        <td style="width: $space%"></td>
        <td align="left">$desc2</td>
    </tr>
</table>
<br />
HTML

    $opts{'-return'} ? return $html : $self->add_html($html);
}

sub imgs2html {
    my ($self, %opts) = @_;

    my $placeholder = $self->require_content_config("SEARCH_PLACEHOLDER");

    my $images = require_opt(\%opts, '-files');
    my $names = require_opt(\%opts, '-names');
    
    $opts{'-type'} = "image";
    my $title_html = $self->title(%opts);

    my $resp_tabs_cnt = $self->{parent}->{resp_vtabs_cnt};
    $self->{parent}->{resp_vtabs_cnt} ++;

    stop("the number of names is not equal to the number of images") unless ($#$names == $#$images);

    my $names_li = join("\n", map { qq(        <li>$_</li>) } @$names);

    my $images_div = $self->_get_img_div_html($images, %opts);

    my $html = _generate_vtab_structure(
        $resp_tabs_cnt, 
        $names_li, 
        $images_div, 
        $title_html, 
        $placeholder
    );

    $self->imgs_into_json(
        -files => $images,
        -order => $resp_tabs_cnt,
        -type => "parentVerticalTab"
    ) unless $self->{'nonlazy'};
    $opts{'-return'} ? return $html : $self->add_html($html);
}

sub plot_grid {
    my ($self, %opts) = @_;

    my $files = require_opt(\%opts, '-files');
    my $nfiles = @$files;

    $opts{'-byrow'} //= 1;

    use POSIX qw/ceil/;
    if (exists $opts{'-ncol'} && exists $opts{'-nrow'}) {
        unless ($opts{'-ncol'} * $opts{'-nrow'} >= $nfiles) {
            stop("function <plot_grid> : -ncol x -nrow must be larger or equal than the number of -files.");
        }
    } elsif (exists $opts{'-ncol'}) {
        stop("function <plot_grid> : -ncol must be larger than 0.") unless ($opts{'-ncol'} > 0);
        $opts{'-nrow'} = POSIX::ceil($nfiles / $opts{'-ncol'});
    } elsif (exists $opts{'-nrow'}) {
        stop("function <plot_grid> : -nrow must be larger than 0.") unless ($opts{'-nrow'} > 0);
        $opts{'-ncol'} = POSIX::ceil($nfiles / $opts{'-nrow'});
    } else {
        if ( $opts{'-byrow'} ) {
            $opts{'-ncol'} = POSIX::ceil(sqrt($nfiles));
            $opts{'-nrow'} = POSIX::ceil($nfiles / $opts{'-ncol'});
        } else {
            $opts{'-nrow'} = POSIX::ceil(sqrt($nfiles));
            $opts{'-ncol'} = POSIX::ceil($nfiles / $opts{'-nrow'});
        }
    }
        
    my $colspan = 1;
    my $out = "<div><table class='pic_table'><tbody>";
    if ($opts{'-byrow'}) {
        my $max_col = $nfiles < $opts{'-ncol'} ? $nfiles : $opts{'-ncol'};
        while (@$files && $opts{'-nrow'} --) {
            my @files_in_a_row = map { $self->_img2html(-file => $_) } splice @$files, 0, $max_col;
            ## in case files number is less than max_col
            push @files_in_a_row, ('') x ($max_col - @files_in_a_row); 
            $out .= join '', "<tr>", (map { "<td>$_</td>" } @files_in_a_row), "</tr>";
        }
        $colspan = $max_col;
    } else {
        my $max_row = $nfiles < $opts{'-nrow'} ? $nfiles : $opts{'-nrow'};
        my @t; ## t for transpose
        while (@$files && $opts{'-ncol'} --) {
            my @files_in_a_col = map { $self->_img2html(-file => $_) } splice @$files, 0, $max_row;
            push @files_in_a_col, ( '' ) x ( $max_row - @files_in_a_col );
            push @t, \@files_in_a_col;
        }
        for my $i ( 0 .. $max_row - 1 ) {
            $out .= join '', "<tr>", (map { "<td>$_</td>" } map { $_->[$i] } @t), "</tr>";
        }
        $colspan = @t;
    }
    if (exists $opts{'-name'}) {
        my $order = $self->img_order();
        $out .= "<tr><td colspan='$colspan' class='img_title'>$order $opts{'-name'}</td></tr>";
    }
    $out .= "</tbody></table><p></p></div>";

    return $out;
}

# add one image to the html 
sub _clean_image_path {
    my ($path) = @_;
    $path =~ s/^\.\.\///;
    $path =~ s/^image/src\/image/;
    return $path;
}

sub _get_src_attribute {
    my ($self) = @_;
    return $self->{'nonlazy'} ? "src" : "data-src";
}

sub _generate_img_link {
    my ($path, $src_attr, $width) = @_;
    $width = defined $width ? "width=\"$width\"" : "";
    return qq(<a href="$path" target="_blank"><img $src_attr="$path" $width/></a>);
}

sub _generate_image_title {
    my ($self, $name, $help) = @_;
    my $order = $self->img_order();
    return "$order $name$help";
}

sub _generate_img_with_desc {
    my ($path, $src_attr, $desc) = @_;
    return qq(<a href="$path" target="_blank"><img $src_attr="$path" /></a>) unless ($desc);
    return <<TEMP;
    <div height="80%"><a href="$path" target="_blank"><img $src_attr="$path" /></a></div>
    <div height="20%"><p>$desc</p></div>
TEMP
}

sub _generate_vtab_structure {
    my ($resp_tabs_cnt, $names_li, $images_div, $title_html, $placeholder) = @_;
    $placeholder ||= "Search";
    $placeholder .= "...";
    
    return <<HTML;
<div id="parentVerticalTab$resp_tabs_cnt" class="VerticalTab">
    <div style="text-align:left;margin-bottom:10px;">
        <input class="iput" placeholder="${placeholder}" oninput="search(this.value,'#resp-vtabs-list$resp_tabs_cnt')">
    </div>
    <ul id="resp-vtabs-list$resp_tabs_cnt" class="resp-tabs-list hor_$resp_tabs_cnt">
$names_li
    </ul>
    <div id="resp-vtabs-container$resp_tabs_cnt" class="resp-tabs-container hor_$resp_tabs_cnt">
$images_div
    </div>
$title_html
    </div>
<br />
HTML
}

# add multi images to the html 
sub _get_img_div_html {
    my ($self, $images, %opts) = @_;

    my $images_div = "";
    return $images_div unless ($self->{'nonlazy'});

    my $i = 0;
    my $src = $self->_get_src_attribute();
    foreach (@$images) {
        my $img_path = _clean_image_path($_);
        my $img_content = $opts{'-desc'} 
            ? _generate_img_with_desc($img_path, $src, $opts{'-desc'}->[$i])
            : _generate_img_link($img_path, $src);
        $images_div .= qq(<div>\n$img_content</div>\n);
        $i ++;
    }
    return $images_div;
}

# add multi images to the html with two columns
sub imgs2html2 {
    my ($self, %opts) = @_;

    my $images1 = require_opt(\%opts, '-files1');
    my $images2 = require_opt(\%opts, '-files2');
    my $names = require_opt(\%opts, '-names');
    my $placeholder = $self->require_content_config("SEARCH_PLACEHOLDER");

    my $desc1 = $opts{'-desc1'};
    my $desc2 = $opts{'-desc2'};
    
    # defined the img title and its description
    $opts{'-type'} = "image";
    my $name = $self->title(%opts);

    my $resp_tabs_cnt = $self->{parent}->{resp_vtabs_cnt};
    $self->{parent}->{resp_vtabs_cnt} ++;

    stop("the number of names is not equal to the number of images1") unless ($#$names == $#$images1);
    stop("the number of names is not equal to the number of images2") unless ($#$names == $#$images2);

    my $src = $self->{'nonlazy'} ? "src" : "data-src";
    my $names_li = join("\n", map { qq(        <li>$_</li>) } @$names);
    
    my $images_div = "";
    my @tabs;
    foreach (0 .. $#$images1) {
        $images1->[$_] =~ s/^\.\.\///;
        $images2->[$_] =~ s/^\.\.\///;

        my @tds;
        if ($opts{'-desc1'} && $opts{'-desc2'}) {
            $images_div .= <<TEMP;
<div>
<table class="pic_table">
    <tr>
        <td>
            <a href="$images1->[$_]" target="_blank"><img $src="$images1->[$_]" /></a>
            <p>$desc1->[$_]</p>
        </td>
        <td>
            <a href="$images2->[$_]" target="_blank"><img $src="$images2->[$_]" /></a>
            <p>$desc2->[$_]</p>
        </td>
    </tr>
</table>
</div>
TEMP
            @tds = (
                qq("$images1->[$_]"), 
                qq("$desc1->[$_]"),
                qq("$images2->[$_]"),
                qq("$desc2->[$_]")
            );
        } else {
            $images_div .= <<TEMP;
<div>
<table class="pic_table">
    <tr>
        <td>
            <a href="$images1->[$_]" target="_blank"><img $src="$images1->[$_]" /></a>
        </td>
        <td>
            <a href="$images2->[$_]" target="_blank"><img $src="$images2->[$_]" /></a>
        </td>
    </tr>
    <tr>
        <td></td>
        <td></td>
    </tr>
</table>
</div>
TEMP
            @tds = (
                qq("$images1->[$_]"),
                qq(""),
                qq("$images2->[$_]"),
                qq("")
            );
        }
        my $tds = join ",\n" , @tds;
        push @tabs , qq("tb$_":[$tds]);
    }
    
    $images_div = "" unless $self->{'nonlazy'};

    my $html = <<HTML;
<div id="parentVerticalTab$resp_tabs_cnt" class="VerticalTab">
    <div style="text-align:left;margin-bottom:10px;">
        <input class="iput" placeholder="${placeholder}..." oninput="search(this.value,'#resp-vtabs-list$resp_tabs_cnt')">
    </div>
    <ul id="resp-vtabs-list$resp_tabs_cnt" class="resp-tabs-list hor_$resp_tabs_cnt">
$names_li
    </ul>
    <div id="resp-vtabs-container$resp_tabs_cnt" class="resp-tabs-container hor_$resp_tabs_cnt">
$images_div
    </div>
    $name
</div>
<br />
HTML
    
    $self->{parent}->{getPic} .= "\t\t\t\taddTabPic(res.container$resp_tabs_cnt, $resp_tabs_cnt)\n";
    my $tabs = join ",\n" , @tabs;
    $self->{parent}->{json} .= <<JSON;
"container$resp_tabs_cnt":{
    $tabs
},
JSON
    
    $opts{'-return'} ? return $html : $self->add_html($html);
}

# turn img to html
sub _img2html {
    my ($self, %opts) = @_;
    
    my $file = _clean_image_path(require_opt(\%opts, "-file"));
    my $src = _get_src_attribute($self);
    
    return _generate_img_with_desc($file, $src, $opts{'-desc'});
}

sub div_pack {
    my $self = shift;
    my @objs = @_;
    my @divs = map { "<div>$_</div>" } @objs;
    return ${[ join "" , @divs ]}[0];
}

#  create a vertical slides which can contains anything
sub any2vertab {
    my ($self, $main_div, %opts) = @_;
    
    my $names = require_opt(\%opts, "-names");
    my $title = $self->title(%opts);

    my $resp_tabs_cnt = $self->{parent}->{resp_vtabs_cnt};
    $self->{parent}->{resp_vtabs_cnt} ++;

    my $html = _vh_tab_html($main_div, $names, $resp_tabs_cnt, $title, "vtab");

    $opts{'-return'} ? return $html : $self->add_html($html);
}

# create a horizontal sildes which can contains anything
sub any2hontab {
    my ($self, $main_div, %opts) = @_;

    my $names = require_opt(\%opts, "-names");
    my $title = $self->title(%opts);

    my $resp_tabs_cnt = $self->{parent}->{resp_htabs_cnt};
    $self->{parent}->{resp_htabs_cnt}++;
    
    my $html = _vh_tab_html($main_div, $names, $resp_tabs_cnt, $title, "htab");

    $opts{'-return'} ? return $html : $self->add_html($html);
}

sub any2child_vtab {
    my ($self, $main_div, %opts) = @_;

    my $names = require_opt(\%opts, "-names");
    my $title  = $self->title(%opts);

    my $placeholder = $self->require_content_config("SEARCH_PLACEHOLDER");
    
    my @names_li = map { "<li>$_</li>" } @$names;
    my $names_li = join "\n" , @names_li;
    
    my $child_tabs_cnt = $self->{parent}->{child_vtabs_cnt};
    $self->{parent}->{child_vtabs_cnt} ++;

    my $html = <<HTML;
<div id="ChildVerticalTab$child_tabs_cnt" class="VerticalTab">
    <div style="text-align:left;margin-bottom:10px;">
        <input class="iput" placeholder="${placeholder}..."
            oninput="search(this.value,'#child-vtabs-list$child_tabs_cnt')">
    </div>
    <ul id="child-vtabs-list$child_tabs_cnt" class="resp-tabs-list ver_$child_tabs_cnt">
$names_li
    </ul>
    <div id="child-vtabs-container$child_tabs_cnt" class="resp-tabs-container ver_$child_tabs_cnt">
$main_div
    </div>
    $title
</div>
<br />
HTML
    
    return $html;
}

sub _vh_tab_html {
    my ($main_div, $names, $resp_tabs_cnt, $title, $type) = @_;

    my $class;
    my $short_class;
    if ($type eq "vtab") {
        $class = "VerticalTab";
        $short_class = "vtab";
    } elsif ($type eq "htab") {
        $class = "HorizontalTab";
        $short_class = "htab";
    } else {
        stop("type must be vtab or htab");
    }

    check_ref($names, "ARRAY");
    my $names_li = join("\n", map { "<li>$_</li>" } @$names);

    my $html = <<HTML;
<div id="parent${class}${resp_tabs_cnt}" class="$class">
    <ul id="resp-${short_class}-list${resp_tabs_cnt}" class="resp-tabs-list hor_${resp_tabs_cnt}">
$names_li
    </ul>
    <div id="resp-${short_class}-container${resp_tabs_cnt}" class="resp-tabs-container hor_${resp_tabs_cnt}">
$main_div
    </div>
    $title
</div>
<br />
HTML

    return $html;
}

sub tsvs2html {
    my ($self, %opts) = @_;
    my $files = require_opt(\%opts, '-files');
    
    my @tabs_div = map { $self->_tsv2html(-file => $_ ,  %opts) } @$files;

    if ($opts{'-orient'} && $opts{'-orient'} eq "vertical") {
        $self->any2vertab($self->div_pack(@tabs_div), -type => "table", %opts);
    } else {
        $self->any2hontab($self->div_pack(@tabs_div), -type => "table", %opts);
    }
}

# fetch the top lines of table, and turn it to html table
sub _tsv2html {
    my ($self, %opts) = @_;

    $opts{'-force'} = 0;
    $opts{'-no_order'} = 1;
    $opts{'-return'} = 1;

    return $self->tsv2html(%opts);
}

sub imgs2mtxhtml {
    my ($self, %opts) = @_;

    my $rownames = require_opt(\%opts, '-rownames');
    my $colnames = require_opt(\%opts, '-colnames');
    my $files = require_opt(\%opts, '-files');

    my $nrow = $#$rownames + 1;
    my $ncol = $#$colnames + 1;

	if (ref $files->[0] eq 'ARRAY') {
		stop("Length of '-files' isnot equal to '-colnames'") unless ($#$files == $ncol - 1);
	} else {
		$files = [map {[splice @$files, 0, $nrow]} 0 .. $ncol-1];
	}
	if ( ref $rownames->[0] eq 'ARRAY' ) {
		stop("Number of '-rownames' isnot equal to '-colnames'.") unless ($nrow == $ncol);
		for my $i ( 0 .. $ncol-1 ) {
			stop("Number of '-rownames'->[$i] isnot equal to '-files'->[$i].") unless ($#{$rownames->[$i]} == $#{$files->[$i]});
		}
	} else {
		stop("the number of image files is not enough!") unless ($nrow * $ncol == (map { @$_ } @$files));
		$rownames = [ map{ $rownames } 0 .. $ncol - 1 ];
	}

    my @child_imgs = map { 
		my $subfile = $files->[$_];
		my $sub_rownames = $rownames->[$_];

        # save images to json file for lazy load
        $self->imgs_into_json(
            -files => $subfile,
            -order => $self->{parent}->{child_vtabs_cnt},
            -type  => "ChildVerticalTab"
        ) unless ($self->{nonlazy});
        
        # create child vertical tabs 
        my @sub_img_div = map { $self->_img2html(-file=>$_) } @$subfile;
        my $child_div = $self->{nonlazy} ? $self->div_pack(@sub_img_div) : "";
        my $child_vtabs = $self->any2child_vtab($child_div, -names => $sub_rownames);

        $child_vtabs;
    } 0 .. $ncol-1;

    my $parent_img_div = $self->div_pack(@child_imgs);
    $self->any2hontab(
        $parent_img_div,
        -names => $colnames,
        -type => "image",
        -name => $opts{'-name'}, 
        -note => $opts{'-note'}, 
        -help => $opts{'-help'}
    );
}

# lazy load plots
sub lazy_load {
    my $self = shift;
    my $getPic = $self->{getPic};

    $self->require_src("js", "lazy_load_pre.js");
    $self->require_src("js", "lazy_load_main.js");

    my $js = <<JS;
        <script src="src/js/lazy_load_pre.js"></script>
        <script type="text/javascript">
            // fetch the data from json file
            function getPic(res){
                $getPic
            }
        </script>
        <script src="src/js/lazy_load_main.js"></script>
        <script src="src/js/pic.js?callback=getPic"></script>
JS

    return $js;
}

# push images to json
sub imgs_into_json {
    my ($self, %opts) = @_;

    my $images = require_opt(\%opts, '-files');
    my $resp_tabs_cnt = require_opt(\%opts, '-order');
    
    my $type = $opts{'-type'} || "parentHorizontalTab";
    
    my %funcs = (
        "parentHorizontalTab" => "addParentHorPic",
        "parentVerticalTab" => "addParentVerPic",
        "ChildHorticalTab" => "addChildHorPic" ,
        "ChildVerticalTab" => "addChildVerPic" 
    );
    my %containers = (
        "parentHorizontalTab" => "phcontainer",
        "parentVerticalTab" => "pvcontainer",
        "ChildHorticalTab" => "chcontainer",
        "ChildVerticalTab" => "cvcontainer"
    );
    my $container = $containers{$type} . $resp_tabs_cnt;
    
    $self->{parent}->{getPic} .= "$funcs{$type}(res.$container, $resp_tabs_cnt)\n";
    my $imgsrc = join ",\n" , (map { qq("$_") } @$images);

    $self->{parent}->{json} .= <<JSON;
"$container":[
    $imgsrc
],
JSON

}

# return the help code
sub opts2help {
    my ($self, %opts) = @_;

    if ($opts{'-help_id'} && !$opts{'-help'}) {
        $opts{'-help'} = "help.html#$opts{'-help_id'}";
    }
    return "" unless ($opts{'-help'});

    $self->require_src("image", "help.png");
    my $html = <<HTML;
<a href="src/doc/$opts{'-help'}" target="help_page" onclick="show_help();">
    <img src="src/image/help.png" class="help_logo">
</a>
HTML

    return $html;
}

# turn options to attributes of html tag
sub opts2attrs {
    my %opts = @_;

    my @besides_keys = qw/parent help files desc pre page_head check/;
    my %besides = map { $_ => 1 } @besides_keys;

    my @attrs = map {
        my $name = $_;
        my $val = $opts{$name};
        if (substr($name, 0, 1) eq "-") {
            $name = substr($name, 1);
        }
        if ($besides{$name}) {
            "";
        } elsif ($val =~ /^\d+$/) {
            "$name=$val";
        } else {
            qq($name="$val")
        }
    } keys %opts;
    return join " " , @attrs;
}

sub omits_some_columns {
    my ($omits, @values) = @_;
    
    return @values unless ($omits);

    my @tmp = split /,/, $omits;
    my @omits;
    foreach (@tmp) {
        if (/^\d+$/) {
            push @omits , [$_, $_];
        } elsif (/(\d+)\s*[:-]\s*(\d+)/) {
            my ($start, $end) = ($1, $2);
            if ($end < $start) {
                stop("the omits values defined is stop, [$start > $end] !");
            }
            push @omits , [$start, $end];
        }
    }
    foreach (sort {$b->[0] <=> $a->[0]} @omits) {
        my ($start, $end) = @$_;
        my $len = $end - $start + 1;
        splice(@values, $start - 1, $len, ("..."));
    }
    return @values;
}

1;
