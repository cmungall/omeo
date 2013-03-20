#!/usr/bin/perl
use lib qw (MediaWiki/);
use lib qw (util/);
use MediaWiki::Bot;
my $bot = MediaWiki::Bot->new();
$bot->set_wiki('www.snpedia.com','/');

$bot->{api}->{use_http_get} = 1;
my $rs = shift @ARGV;
my $text = $bot->get_text($rs);
print $text;
