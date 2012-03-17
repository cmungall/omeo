#!/usr/bin/perl
while (<>) {
    my ($d,$syms,$gene) = split(/\|/,$_);
    $d =~ s/\s+$//;
    $d =~ s/\s*\(\d+\)//;
    if ($d =~ /\s(\d+)$/) {
	print "MIM:$1\tMIM:$gene\n" unless $1 < 100; # eg Alopecia areata 1 (2)  |AA1|104000|18p11.3-p11.2
    }
}
