#!/usr/bin/perl
my $t = shift;
while (<>) {
    chomp;
    split(/\t/);
    if ($_[15] eq $t && $_[16] eq $t) {
        printf "NCBIGene:$_[1]\tNCBIGene:$_[2]\n";
    }
}
