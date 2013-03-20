#!/usr/bin/perl
use lib qw (MediaWiki/);
use lib qw (util/MediaWiki/);
use lib qw (util/);
use MediaWiki::Bot;
use strict;

my $bot = MediaWiki::Bot->new();
$bot->set_wiki('www.snpedia.com','/');

$bot->{api}->{use_http_get} = 1;

print STDERR "Getting all genotypes...\n";
my @genotypes = $bot->get_pages_in_category('Category:Is a genotype', {max => 0}) ;
print STDERR "GTs: @genotypes\n";
printf STDERR "\nQuerying genotypes... #= " . scalar(@genotypes);
#@genotypes = grep {/^rs/i} @genotypes;
printf STDERR "\nFiltered genotypes... #= " . scalar(@genotypes);

foreach my $rs (@genotypes) {
    print STDERR "GT: $rs\n";
    my $text = $bot->get_text($rs);
    open(F, ">snpedia/$rs.wiki") || die($rs);
    print F $text;
    close(F);
    sleep(1);
}
exit 0;
