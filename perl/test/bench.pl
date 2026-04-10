
use strict;
use warnings;
use Benchmark qw(cmpthese);

sub method_copy {
    my ($opts, @keys) = @_;
    foreach my $k (@keys) {
        next if defined $opts->{$k} && $opts->{$k} ne '';
    }
}

sub method_shift {
    my $opts = shift;
    foreach my $k (@_) {
        next if defined $opts->{$k} && $opts->{$k} ne '';
    }
}

# 构造测试数据
my %opts = map { ("k$_" => $_) } (1..100);

# 构造参数（模拟很多 key）
my @keys = map { "k$_" } (1..100);

cmpthese(-3, {
    copy  => sub { method_copy(\%opts, @keys) },
    shift => sub { method_shift(\%opts, @keys) },
});
