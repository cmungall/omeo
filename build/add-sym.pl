#!/usr/bin/perl
while(<>) {
    split;
    print "$_[0]\t$_[1]\n";
    print "$_[1]\t$_[0]\n";
}
