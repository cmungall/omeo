#!/usr/bin/perl
print "ontology: _\n\n";
while (<>) {
    #next if /^\!/;
    @vs = split(/\t/,$_);
    my $id = "$vs[0]:$vs[1]";
    next if $done{$id};
    @syns = split(/\|/,$vs[10]);
    if (@syns) {
        print "[Term]\n";
        print "id: $id\n";
        foreach (@syns) {
            s/\"/\\\"/g;
            print "synonym: \"$_\" EXACT []\n";
        }
        print "\n";
    }
    $done{$id}=1;
}
