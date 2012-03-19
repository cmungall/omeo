#!/usr/bin/perl


print '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',"\n";
print '<catalog prefer="public" xmlns="urn:oasis:names:tc:entity:xmlns:xml:catalog">', "\n";

my %ext = 
    (
     "http://purl.obolibrary.org/obo/so.owl" => "so.owl",
     "http://purl.obolibrary.org/obo/ro.owl" => "ro.owl",
     "http://purl.obolibrary.org/obo/iao/ontology-metadata.owl" => "ontology-metadata.owl",
    );

while (<>) {
    chomp;
    s@^\./@@;
    printf '    <uri id="User Entered Import Resolution" name="http://purl.obolibrary.org/obo/omeo/%s" uri="%s"/>%s', $_,$_,"\n";
}

foreach my $uri (keys %ext) {
    printf '    <uri id="User Entered Import Resolution" name="%s" uri="../external/%s"/>%s', $uri,$ext{$uri},"\n";

}

print "</catalog>\n";
