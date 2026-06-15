package ReportMe;

use strict;
use warnings;

use ReportMe::HTML;
use ReportMe::HELP;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = (
	@ReportMe::HTML::EXPORT,
	@ReportMe::HELP::EXPORT,
);

1;
