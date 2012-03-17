#!/usr/bin/perl


print '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',"\n";
print '<catalog prefer="public" xmlns="urn:oasis:names:tc:entity:xmlns:xml:catalog">', "\n";

while (<>) {
    chomp;
    s@^\./@@;
    printf '    <uri id="User Entered Import Resolution" name="http://purl.obolibrary.org/obo/omeo/%s" uri="%s"/>%s', $_,$_,"\n";
}

print "</catalog>\n";
