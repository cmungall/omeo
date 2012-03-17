#!/usr/bin/perl
while (<>) {
    next if /^\!/;
    @vs = split(/\t/,$_);
    my $id = "$vs[0]:$vs[1]";
    next if $done{$id};
    @syns = split(/\|/,$vs[10]);
    print "$id\t$_\n" foreach @syns;
    $done{$id}=1;
}
