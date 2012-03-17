#!/usr/bin/perl
my $SP = shift;
while(<>) {
    s/(\S+)/$1\t$SP chromosome $1/;
    print;
}
