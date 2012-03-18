

% ----------------------------------------
% MART
% ----------------------------------------

tax_db_species(7955,'drerio_gene_ensembl','Danio_rerio').
tax_db_species(9606,'hsapiens_gene_ensembl','Homo_sapiens').
tax_db_species(10090,'mmusculus_gene_ensembl','Mus_musculus').
tax_db_species(10116,'rnorvegicus_gene_ensembl','Rattus_norvegicus').
tax_db_species(6239,'celegans_gene_ensembl','Caenorhabditis_elegans').
tax_db_species(7227,'dmelanogaster_gene_ensembl','Drosophila_melanogaster').
tax_db_species(4932,'scerevisiae_gene_ensembl','Saccharomyces_cerevisiae').
tax_db_species(9031,'ggallus_gene_ensembl','Gallus_gallus').
tax_db_species(7719,'cintestinalis_gene_ensembl','Ciona_intestinalis').
tax_db_species(4896,'spombe_eg_gene','Schizosaccharomyces_pombe').
%tax_db_species(7757,'pmarinus_gene_ensembl','Petromyzon_marinus').

non_refg(7719).
non_refg(7757).


anc_combo('Opisthokont',['Danio_rerio',
                         'Homo_sapiens',
                         'Mus_musculus',
                         'Rattus_norvegicus',
                         'Caenorhabditis_elegans',
                         'Drosophila_melanogaster',
                         'Saccharomyces_cerevisiae']).

anc_combo('Fungi',['Schizosaccharomyces_pombe',
                   'Saccharomyces_cerevisiae']).

anc_combo('Mammal',['Homo_sapiens',
                    'Mus_musculus',
                    'Rattus_norvegicus']).

anc_combo('Mammal-hm',['Homo_sapiens',
                       'Mus_musculus']).

anc_combo('Vertebrate',['Homo_sapiens',
                        'Mus_musculus',
                        'Rattus_norvegicus',
                        'Gallus_gallus',
                        'Danio_rerio']).

anc_combo('Vertebrate-hmz',['Homo_sapiens',
                            'Mus_musculus',
                            'Danio_rerio']).


/*
tax_idprop_mod(9606,mim_gene_accession,'MIM').
%tax_idprop_mod(9606,mim_gene_accession,'MIM').
%tax_idprop_mod(9606,mim_morbid_accession,'MIM').
tax_idprop_mod(7955,'zfin_id','ZFIN').
tax_idprop_mod(10090,'mgi_id','MGI').
tax_idprop_mod(10116,'rgd','RGD').
tax_idprop_mod(6239,'wormbase_gene','WB').
tax_idprop_mod(7227,'flybase_gene_id','FB').
tax_idprop_mod(4932,'sgd','SGD').
*/

tax_idprop(9606,mim_gene_accession).
tax_idprop(9606,mim_morbid_accession).
tax_idprop(7955,'zfin_id').
tax_idprop(10090,'mgi_id').
tax_idprop(10116,'rgd').
tax_idprop(6239,'wormbase_gene').
tax_idprop(7227,'flybase_gene_id').
tax_idprop(4932,'sgd').
tax_idprop(9031,'entrezgene').
tax_idprop(7719,'entrezgene').
tax_idprop(7757,'hgnc_id').
tax_idprop(4896,pombase_gene).

tax_mod(9606,'NCBIGene').
tax_mod(7955,'ZFIN').
tax_mod(10090,'MGI').
tax_mod(10116,'RGD').
tax_mod(6239,'WB').
tax_mod(7227,'FB').
tax_mod(4932,'SGD').
tax_idprop(4896,'SPombe').
tax_mod(9031,'NCBIGene').
tax_mod(7719,'NCBIGene').
tax_mod(7757,'ENSEMBL').
tax_mod(4896,'PomBase').

sp_prefix(SP,Pre) :- tax_db_species(Tax,_,SP),tax_mod(Tax,Pre).
sp(SP) :- tax_db_species(_,_,SP).

attr_pair(ensembl_gene_id,uniprot_swissprot_accession).
attr_pair(ensembl_gene_id,external_gene_id).
% attr_pair(ensembl_gene_id,hgnc_symbol). human only
attr_pair(ensembl_gene_id,description).
attr_pair(ensembl_gene_id,gene_biotype).
attr_pair(ensembl_gene_id,chromosome_name).
attr_pair(ensembl_gene_id,entrezgene).
attr_pair(ensembl_gene_id,interpro).
attr_pair(ensembl_gene_id,external_id). % snp rs#, etc. Only for some species; not currently used in any owl
attr_pair(ensembl_gene_id,interpro).
attr_pair(ensembl_gene_id,family).
attr_pair(family, family_description).

% this determines what goes in owl:imports
subont_tax('label').
subont_tax('comment').
subont_tax('taxon').
subont_tax('equiv-mod').
subont_tax('gene-chrom').
subont_tax('chrom').
subont_tax('chrom2label').
subont_tax('panther').
subont_tax('entrez').
%subont_tax('panther-ids').
subont_tax('type').
subont_tax('uniprot').
%subont_tax('entrez',T) :- T\=7757. % no entrez for lamprey
%subont_tax('uniprot',10090).
%subont_tax('uniprot',9606).
subont_tax('uniprot-pro',10090).
subont_tax('uniprot-pro',9606).
subont_tax('ensfm',T) :- T\=4896.
subont_tax('family-description',T) :- T\=4896.
subont_tax(S,_) :- subont_tax(S).
subont_tax(S,_) :- subont_tax(S).
subont_tax(S,T,main) :- subont_tax(S,T).
subont_tax('goannotations',T,ann) :- \+ non_refg(T).  % don't include in imports
subont_tax('expr',10090,ann).
subont_tax('expr',9606,ann).
subont_tax('expr',7955,ann).

sp_attr_pair(SP,A1,A2) :-
        tax_db_species(_,_,SP),
        attr_pair(A1,A2).
sp_attr_pair(SP,mod_id,ensembl_gene_id) :-
        tax_db_species(Tax,_,SP).


all_mart <-- Deps,
       {findall(t([SP,-,A1,-,A2,'.mart']),
                sp_attr_pair(SP,A1,A2),
                Deps)}.

all_owl_subont <-- Deps,
       {findall(t([SP,-,Type,'.owl']),
                (   tax_db_species(Tax,_,SP),subont_tax(Type,Tax,_)),
                Deps)}.

all_owl <-- Deps,
       {findall(t([SP,'.owl']),
                tax_db_species(_,_,SP),
                Deps)}.

all_obo <-- Deps,
       {findall(t([SP,'.obo']),
                tax_db_species(_,_,SP),
                Deps)}.

all_ensembl_obo <-- Deps,
       {findall(t([SP,'-ensembl.obo']),
                tax_db_species(_,_,SP),
                Deps)}.

all_mod_obo <-- Deps,
       {findall(t([SP,'-MOD.obo']),
                tax_db_species(_,_,SP),
                Deps)}.

all <-- [all_owl,all_obo].

'owl-$SP' <-- Deps,
       {findall(t([SP,-,Type,'.owl']),
                (   tax_db_species(Tax,_,SP),subont_tax(Type,Tax)),
                Deps)},
       'touch $@'.

species_dependency(SP, File, t([SP,-,Type,'.owl'])) :-
        tax_db_species(Tax,_,SP),
        subont_tax(Type,Tax),
        atomic_list_concat([SP,-,Type,'.owl'],File).
species_dependency(SP, File, t([SP,'-gaf-syns.owl'])) :-
        tax_db_species(Tax,_,SP),
        \+ non_refg(Tax),
        atom_concat(SP,'-gaf-syns.owl',File).

species_dependencies(SP,Deps,Files) :-
        sp(SP),
        findall(Dep,
                species_dependency(SP,_,Dep),
                Deps),
        findall(File,
                species_dependency(SP,File,_),
                Files).

        

%'$SP.owl' <-- ['owl-$SP'],
'$SP.owl' <-- Deps,
  {sp(SP),
   species_dependencies(SP,Deps,Files),
   maplist(atom_concat('http://purl.obolibrary.org/obo/omeo/'),Files,URLs),
   atomic_list_concat(URLs,' ',FA)},
  'owltools --create-ontology omeo/$SP.owl --add-imports-declarations $FA http://purl.obolibrary.org/obo/so.owl http://purl.obolibrary.org/obo/ro.owl -o file://`pwd`/$@'.
%  'owltools --create-ontology omeo/$SP.owl $FA so.owl --merge-support-ontologies -o file://`pwd`/$@'.

% makes merged ontologies [DEPRECATED]
'combo/$Anc.owl' <-- [],
  {anc_combo(Anc,SpList),
   findall(F,(member(Sp,SpList),atom_concat(Sp,'.owl',F)),Fs),
   atomic_list_concat(Fs,' ',FA)
   },
  'owltools  --create-ontology omeo/$Anc.owl $FA --merge-support-ontologies -o file://`pwd`/$@'.


'$SP.obo' <-- Deps,
  {species_dependencies(SP,Deps,Files),atomic_list_concat(Files,' ',FilesAtom)},
  'owltools --create-ontology omeo/$SP.owl $FilesAtom --merge-support-ontologies -o -f obo $@'.

% map to MOD IDs
'$SP-MOD.obo' <-- '$SP.obo',
  {sp_prefix(SP,Prefix)},
  'obo-rewrite-to-equivalent-class.pl -t $Prefix $< > $@'.

%'$SP-MOD.owl' <-- '$SP-MOD.obo',
%   'obolib-owl2obo --allow-dangling -o $@ $<'.

'%-label.obo' <-- '%-label.owl',
  'owltools $< --translate-undeclared-to-classes -o -f obo $@'.

'%.obo' <-- '%.owl',
  {sp(SP)},
  'obolib-owl2obo -o $@ $<'.

'$Binomial-$F1-$F2.mart' <-- [],
  {F1\=mod_id,tax_db_species(_,DB,Binomial)},
  'biomart-fetch.pl -d $DB -a $F1 -a $F2 > $@.tmp && mv $@.tmp $@'.

'$Binomial-mod_id-$F1.mart' <-- [],
  {tax_idprop(Tax,F2),
   tax_db_species(Tax,DB,Binomial)},
    'biomart-fetch.pl -d $DB -a $F2 -a $F1 > $@.tmp && mv $@.tmp $@'.

'$SP-gene2tax.tbl' <-- ['$SP-ensembl_gene_id-external_gene_id.mart'],
  {tax_db_species(Tax,_,SP)},
  'perl -npe "s/\t.*/\tNCBITaxon:$Tax/" $< > $@'.

'$SP-chrom.tbl' <-- ['$SP-ensembl_gene_id-chromosome_name.mart'],
  'cut -f2 $< | sort -u > $@'.

'$SP-biotype.tbl' <-- ['$SP-ensembl_gene_id-gene_biotype.mart'],
  'blip-findall -consult biotype_cls.pro -i $< -f "tbl(g2t)" "g2t(G,TN),biotype_cls(TN,T)" -select G-T -no_pred > $@'.

'$SP-equiv-mod.owl' <-- ['$SP-mod_id-ensembl_gene_id.mart'],
  {sp_prefix(SP,MOD)},
  'owltools --create-ontology omeo/$@ --parse-tsv -a EquivalentClasses --iri-prefix 1 $MOD --iri-prefix 2 ENSEMBL $< -o file://`pwd`/$@'.

'$SP-entrez.owl' <-- ['$SP-ensembl_gene_id-entrezgene.mart'],
  'owltools --create-ontology omeo/$@ --parse-tsv -a EquivalentClasses --iri-prefix 1 ENSEMBL --iri-prefix 2 NCBIGene $< -o file://`pwd`/$@'.

'$SP-label.owl' <-- ['$SP-ensembl_gene_id-external_gene_id.mart'],
  'owltools --create-ontology omeo/$@ --parse-tsv -l --iri-prefix 1 ENSEMBL $< -o file://`pwd`/$@'.

'$SP-comment.owl' <-- ['$SP-ensembl_gene_id-description.mart'],
  'owltools --create-ontology omeo/$@ --parse-tsv --comment --iri-prefix 1 ENSEMBL $< -o file://`pwd`/$@'.

'$SP-taxon.owl' <-- ['$SP-gene2tax.tbl', 'ro-subset.owl'],
  'owltools --create-ontology omeo/$@ ro-subset.owl --merge-support-ontologies --parse-tsv -p RO:0002162 -a SubClassOf --iri-prefix 1 ENSEMBL $< -o file://`pwd`/$@'.

% todo - infer this?
'$SP-panther.owl' <-- ['$SP-pthr.tbl'],
  'owltools --create-ontology omeo/$@ --parse-tsv -a SubClassOf --iri-prefix 2 PANTHER $< -o file://`pwd`/$@'.

'$SP-ensfm.owl' <-- ['$SP-ensembl_gene_id-family.mart'],
  'owltools --create-ontology omeo/$@ --parse-tsv -a SubClassOf --iri-prefix 1 ENSEMBL --iri-prefix 2 ENSEMBL $< -o file://`pwd`/$@'.

'$SP-family-description.owl' <-- ['$SP-family-family_description.mart'],
  'owltools --create-ontology omeo/$@ --parse-tsv -l --iri-prefix 1 ENSEMBL $< -o file://`pwd`/$@'.


'$SP-uniprot.owl' <-- ['$SP-ensembl_gene_id-uniprot_swissprot_accession.mart', 'ro-subset.owl'],
  'owltools --create-ontology omeo/$@ ro-subset.owl --merge-support-ontologies --parse-tsv -p RO:0002205 -a SubClassOf --iri-prefix 1 ENSEMBL --iri-prefix 2 UniProtKB $< -o file://`pwd`/$@'.

%'$SP-uniprot-pro.owl' <-- ['$SP-uniprot.owl', 'uniprot2pro.owl'],
%  'owltools $< uniprot2pro.owl --mcat -o file://`pwd`/$@'.

'$SP-uniprot-pro.owl' <-- ['pro_uniprotmapping-$SP.owl', '$SP-uniprot.owl', 'pro-$SP.obo'],
  'owltools $< $SP-uniprot.owl pro-$SP.obo --merge-support-ontologies -o file://`pwd`/$@'.


'$SP-gene-chrom.owl' <-- ['$SP-ensembl_gene_id-chromosome_name.mart'],
  {tax_db_species(Tax,_,SP)},
  'owltools --create-ontology omeo/$@ --parse-tsv -p BFO:0000050 -a SubClassOf --iri-prefix 1 ENSEMBL --iri-prefix 2 CHROMOSOME-$Tax  $< -o file://`pwd`/$@'.

'$SP-chrom.owl' <-- ['$SP-chrom.tbl'],
  {tax_db_species(Tax,_,SP)},
  'owltools --create-ontology omeo/$@ --parse-tsv -a SubClassOf --iri-prefix 1 CHROMOSOME-$Tax --default2 SO:0000340  $< -o file://`pwd`/$@'.

'$SP-chrom2label.txt' <-- ['$SP-chrom.tbl'],
  {tax_db_species(Tax,_,SP)},
  './util/add-chrom-label.pl $SP $< > $@'.

% TODO
'$SP-chrom2label.owl' <-- ['$SP-chrom2label.txt'],
  {tax_db_species(Tax,_,SP)},
  'owltools --create-ontology omeo/$@ --parse-tsv -l --iri-prefix 1 CHROMOSOME-$Tax  $< -o file://`pwd`/$@'.

'$SP-type.owl' <-- ['$SP-biotype.tbl'],
  'owltools --create-ontology omeo/$@ --parse-tsv -a SubClassOf --iri-prefix 1 ENSEMBL $< -o file://`pwd`/$@'.

'$SP-ensembl.obo' <-- ['$SP.obo'],
  'obo-rewrite-to-equivalent-class.pl -t ENSEMBL $< | obo-order-tags.pl  - > $@'.

'$SP-ensembl.owl' <-- ['$SP-ensembl.obo'],
  'obolib-obo2owl -o $@ $<'.


/*
% TOO SLOW!!
  '$SP-ensembl.owl' <-- ['$SP.owl'],
  'owltools $< --merge-equivalent-classes -t http://purl.obolibrary.org/obo/ensembl_ -o file://`pwd`/$@'.

'$SP-ensembl.obo' <-- ['$SP-ensembl.owl'],
  'obolib-owl2obo -o $@ $<'.
*/

% ----------------------------------------
% PANTHER
% ----------------------------------------

'RefGenomeOrthologs.tar' <-- 'RefGenomeOrthologs.tar.gz',
  'gzip -d $<'.

'RefGenomeOrthologs.tar.gz' <-- [],
  'wget ftp://ftp.pantherdb.org/ortholog/current/$@'.

'refg_ortho_family.txt' <-- 'RefGenomeOrthologs.txt',
      'process-RefGenomeOrthologs.pl $< > $@'.

'hom.tbl' <-- 'refg_ortho_family.txt',
     'perl -ne "print unless /\\tP\\t/" $< | cut -f1,3 > $@.tmp && sort -u $@.tmp > $@'.

% symmetric
'shom.tbl' <-- 'hom.tbl',
     './add-sym.pl $< | sort -u > $@'.

% orthologs from panther
%% TODO - do per species instead? Instance-level?
%'ortho.owl' <-- 'shom.tbl',
%     'owltools --create-ontology panther ro-subset.owl --merge-support-ontologies --parse-tsv -a SubClassOf -p RO:0002158  $< -o file://`pwd`/$@'.

'g2family.txt' <-- 'refg_ortho_family.txt',
     'cut -f1,2,7 $< > $@.tmp && sort -u $@.tmp > $@'.

% rename? gene2family
'ens2family.txt' <-- 'g2family.txt',
   'cut -f1,3 $< > $@'.

'pthr_ids.txt' <-- 'g2family.txt',
   'cut -f3 $< | sort -u > $@'.

'Danio_rerio-pthr.tbl' <-- 'ens2family.txt',   'egrep "^(ZFIN|ENSEMBL:ENSDARG)" $< > $@'.
'Homo_sapiens-pthr.tbl' <-- 'ens2family.txt',   'grep ^ENSEMBL:ENSG0 $< > $@'.
'Mus_musculus-pthr.tbl' <-- 'ens2family.txt',   'grep ^MGI $< > $@'.
'Drosophila_melanogaster-pthr.tbl' <-- 'ens2family.txt',   'grep ^FB $< > $@'.
'Caenorhabditis_elegans-pthr.tbl' <-- 'ens2family.txt',   'grep ^WB $< > $@'.
'Saccharomyces_cerevisiae-pthr.tbl' <-- 'ens2family.txt',   'grep ^SGD $< > $@'.
'Schizosaccharomyces_pombe-pthr.tbl' <-- 'ens2family.txt',   'grep ^GeneDB_Spombe $< > $@'.
'Rattus_norvegicus-pthr.tbl' <-- 'ens2family.txt',   'grep ^RGD $< > $@'.
'Gallus_gallus-pthr.tbl' <-- 'ens2family.txt',   'grep ^ENTREZ $< | perl -npe "s/ENTREZ/NCBIGene/" > $@'.  % TODO - check this
'$SP-pthr.tbl' <-- 'ens2family.txt',   'head -1 $< > $@'.  % TODO - more elegant solution required!

'$SP-pthr-ids.txt' <-- '$SP-pthr.tbl',
   'cut -f2 $< | sort -u > $@'.

'$SP-pthr-ids.owl' <-- '$SP-pthr-ids.txt',
  'owltools --create-ontology panther pthr_root.obo --merge-support-ontologies --parse-tsv -a SubClassOf --iri-prefix 2 PANTHER --default2 PTHR:00000  $< -o file://`pwd`/$@'.

'pro_uniprotmapping.txt' <-- [],
  'wget ftp://ftp.pir.georgetown.edu/databases/ontology/pro_obo/PRO_mappings/uniprotmapping.txt -O $@'.

'uniprot2pro.owl' <-- ['pro_uniprotmapping.txt'],
  'owltools --create-ontology omeo/$@ --parse-tsv -s -a SubClassOf  $< -o file://`pwd`/$@'.

'$SP-g-p.pro' <-- ['$SP-ensembl_gene_id-uniprot_swissprot_accession.mart'],
  'tbl2p -p x -prefix 2=UniProtKB: $< > $@'.

'pro_uniprotmapping-$SP.txt' <-- '$SP-g-p.pro',
  'blip-findall -index "user:x(-,1)" -i $< -i pro_uniprotmapping.txt "pro_uniprotmapping(P,U),x(G,U)" -select "pro_uniprotmapping(P,U)" -no_pred > $@.tmp && sort -u $@.tmp > $@'.

'pro_uniprotmapping-$SP.owl' <-- 'pro_uniprotmapping-$SP.txt',
  'owltools --create-ontology omeo/$@ --parse-tsv -s -a SubClassOf $< -o file://`pwd`/$@'.

'%-ids.pro' <-- '%.txt',
  'cut -f1 $< | tbl2p -p id > $@'.

'pro-$SP.obo' <-- 'pro_uniprotmapping-$SP-ids.pro',
  'blip ontol-query -i only_in_taxon.obo -r protein -query "blipkit:id(ID)" -consult $< -to obo > $@'.


% GFF
% ftp://ftp.ensembl.org/pub/release-64/gtf/

% ----------------------------------------
% GO
% ----------------------------------------

sp_gaf('Homo_sapiens',goa_human).
sp_gaf('Gallus_gallus',goa_chicken).
sp_gaf('Schizosaccharomyces_pombe','GeneDB_Spombe').
sp_gaf(SP,DB) :-
        tax_db_species(Tax,_,SP),tax_mod(Tax,Mod),downcase_atom(Mod,DB).

% todo - uniprot
'$SP.gaf' <-- [],
  {sp_gaf(SP,DB)},
  'wget http://geneontology.org/gene-associations/gene_association.$DB.gz -O $@.gz && gzip -d $@.gz'.

'$SP-gaf-syns.obo' <-- ['$SP.gaf'],
  'util/gaf2obosyns.pl $< > $@'.
'$SP-gaf-syns.owl' <-- ['$SP-gaf-syns.obo'],
  'obolib-obo2owl $< -o  $@'.

'$SP-goannotations.owl' <-- ['$SP.gaf'],




'owltools http://purl.obolibrary.org/obo/go.owl --gaf $< --gaf2owl -n go/annotations/$SP-goannotations -o file://`pwd`/$@'.


%'$SP-gaf-syns.tbl' <-- ['$SP.gaf'],
%  'util/gaf2syns.pl $< > $@'.

% ----------------------------------------
% EXPRESSION
% ----------------------------------------

'zfin-wildtype-expression.txt' <-- [],
  'wget http://zfin.org/data_transfer/Downloads/wildtype-expression.txt -O $@'.

'geisha.txt' <-- [],
  'wget http://geisha.arizona.edu/geisha/expression.txt -O $@'.

% expression.tsv is from bgee; extract annotations for which there is at least one high quality observation. Simplify to gene-AOterm pairs
'expression_hq.tsv' <-- 'expression.tsv',
  'grep high $< | cut -f2,3 | sort -u > $@'.

'Homo_sapiens-expr.txt' <-- 'expression_hq.tsv',
   'grep ^ENSG $< | perl -npe "s/^ENS/ENSEMBL:ENS/"  > $@'.

'Danio_rerio-expr.txt' <-- 'expression_hq.tsv',
   'grep ^ENSDAR $< | perl -npe "s/^ENS/ENSEMBL:ENS/"  > $@'.

'Mus_musculus-expr.txt' <-- 'expression_hq.tsv',
   'grep ^ENSMUSG $< | perl -npe "s/^ENS/ENSEMBL:ENS/"  > $@'.

'Xenopus_tropicalis-expr.txt' <-- 'expression_hq.tsv',
   'grep ^ENSXETG $< | perl -npe "s/^ENS/ENSEMBL:ENS/"  > $@'.

'Drosophila_melanogaster-expr.txt' <-- 'expression_hq.tsv',
   'grep ^FB $< | perl -npe "s/^FB/FB:FB/"  > $@'.

'$SP-expr-ubr.txt' <-- '$SP-expr.txt',
   'obo-map-ids.pl --use-xref-inverse -k 2 --regex-filter UBERON composite-metazoan.obo  $< > $@'.

'$SP-expr.owl' <-- '$SP-expr-ubr.txt',
   'owltools --create-ontology omeo/$@ --parse-tsv -p RO:00022060 -a SubClassOf  $< -o file://`pwd`/$@'.


% ----------------------------------------
% VARIATION
% ----------------------------------------
'Homo_sapiens-varpheno.txt' <-- [],
   'biomart-fetch.pl -c -d hsapiens_snp  -p varpheno  > $@.tmp && mv $@.tmp $@'.

'Homo_sapiens-gene-snp.txt' <-- 'Homo_sapiens-varpheno.txt',
   'cut -f1,2 $< | grep ^ENS | perl -npe "s/ENSG/ENSEMBL:ENSG/;s/\\trs/\\tdbSNP:rs/" > $@'.

'Homo_sapiens-gene-snp.owl' <-- 'Homo_sapiens-gene-snp.txt',
  'owltools --create-ontology omeo/$@ --parse-tsv -p BFO:0000050 -a ClassAssertion -t http://purl.obolibrary.org/obo/SO_0001059  $< -o file://`pwd`/$@'.

% ----------------------------------------
% BIOGRID
% ----------------------------------------
'$SP-biogrid.tbl' <-- 'biogrid.txt',
   {tax_db_species(Tax,_,SP)},
  'util/filter_biogrid_by_taxon.pl $Tax $< > $@.tmp && sort -u $@.tmp > $@'.

% ----------------------------------------
% OMIM
% ----------------------------------------
'omim.txt' <-- [],
        'wget ftp://ftp.ncbi.nih.gov/repository/OMIM/$@.Z && gzip -d $@.Z'.

% MOVED!
'morbidmap.txt' <-- [],
        'wget ftp://ftp.ncbi.nih.gov/repository/OMIM/morbidmap -O $@'.

'disorder_mimgene.txt' <-- 'morbidmap.txt',
    './mm2tbl.pl $< > $@'.

'disorder_ensgene.txt' <-- ['disorder_mimgene.txt', 'Homo_sapiens-equiv-mod.obo'],
   'blip-findall -i disorder_mimgene.txt -i Homo_sapiens-equiv-mod.obo -index "ontol_db:equivalent_class(1,1)" "disorder_mimgene(D,MG),equivalent_class(EG,MG),equivalent_class(EG,NG)" -select D-EG -no_pred > $@'.

download('http://www.genome.gov/admin/gwascatalog.txt').
download('ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/organism_data/OmimVarLocusIdSNP.bcp.gz').  % appears more up to date than the ensembl table

'ext-$F' <-- [],
  {download(URL),atomic_list_concat(Toks,/,URL),reverse(Toks,[F|_])},
  'wget $URL -O ext-$F'.

'Saccharomyces_cerevisiae-pheno.txt' <-- [],
  'wget http://downloads.yeastgenome.org/curation/literature/phenotype_data.tab -O $@'.

vp(hsapiens_snp_som). % somatic
vp(hsapiens_snp). % 132998 / 41427246 
%%%   'biomart-fetch.pl -c -d $*  -p varpheno -i with_variation_annotation > $@'.
%%%%   'biomart-fetch.pl -d $*  -a refsnp_id -a associated_gene -a phenotype_name -a phenotype_description > $@'. UNIQUE

% SEQXML
% e.g. http://seqxml.org/download/ReferenceProteomes/current/10090_mus_musculus.xml.gz

%tax_db_species(

% ----------------------------------------
% PHENOTYPE
% ----------------------------------------

'Mus_musculus-gene-phenotype.txt' <-- [],
  'wget "http://obo.svn.sourceforge.net/viewvc/obo/phenotype-commons/annotations/MGI/gene_phenotype.txt" -O $@'.

'Mus_musculus-genotype-phenotype.txt' <-- [],
  'wget "http://obo.svn.sourceforge.net/viewvc/obo/phenotype-commons/annotations/MGI/genotype_phenotype.txt" -O $@'.

'Homo_sapiens-gene-phenotype.txt' <-- [],
  'wget "http://obo.svn.sourceforge.net/viewvc/obo/phenotype-commons/annotations/Human/gene_phenotype.txt" -O $@'.

% ----------------------------------------
% CATALOGS
% ----------------------------------------

'catalog-v001.xml' <-- [],
  'find . -name "*.owl" | ./make-catalog-xml.pl  > $@'.

% ----------------------------------------
% SUBSETS
% ----------------------------------------
'../subsets/panther-$Fam-$Tax.owl' <-- '$Tax.owl',
  'owltools --use-catalog $Tax.owl --merge-import-closure --reasoner-query -r elk -d -c http://purl.obolibrary.org/obo/omeo/subsets/panther-$Fam-$Tax.owl -l PANTHER:$Fam -o file://`pwd`/$@'.

% ----------------------------------------
% SYNC
% ----------------------------------------

deploy <-- [],
  './deploy.sh'.
%% 'rsync -avz *.{obo,owl,txt,tbl,mart} .lbl.gov:/var/www/ontologies/omeo/'.


