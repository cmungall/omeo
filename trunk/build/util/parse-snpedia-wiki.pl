#!/usr/bin/perl
use strict;
my @files = glob(shift @ARGV);

foreach my $f (@files) {
    if ($f =~ /(Rs\d+)\((\S*);(\S*)\)/) {
        my ($snp,$a1,$a2) = ($1,$2,$3);
        $snp = lc($snp);
        open(F,$f);
        my $content = join('',<F>);
        close(F);
        $content =~ s/\n/ /g;
        my @blocks = split(/}}/, $content);
        #$/ = '}}';
        foreach (@blocks) {
            if (m@(.*){{\s*(\S+)\s*(.*)@m) {
                my ($text,$type,$infobox) = ($1,$2,$3);
                #print STDERR $infobox;
                w($snp,$a1,$a2,'text',$text);
                my @infos = split(/\|/,$infobox);
                #printf STDERR " N=%d\n", scalar(@infos);
                foreach my $info (@infos) {
                    #print STDERR "  INFO:".$info."\n";
                    my ($t,$v) = split(/\s*=\s*/,$info);
                    $t =~ s/ //g;
                    my $p = $type . '_' . $t;
                    w($snp,$a1,$a2,$p,$v);
                }
            }
            else {
                w($snp,$a1,$a2,'text',$_);
            }
        }
        #close(F);
    }
}
exit 0;

sub w {
    my ($snp,$a1,$a2,$p,$v) = @_;
    return unless defined $v;
    $v =~ s/^\s+//;
    $v =~ s/\s+$//;
    return if $v eq '';
    print join("\t",@_)."\n";
    
}
