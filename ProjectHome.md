OMEO integrates various Omics resources into a collection of interoperable OWL2 ontologies.

Entities and sources:

  * Genes (ENSEMBL, NCBI Gene, various MODs)
  * Transcripts (ENSEMBL)
  * Protein (ENSEMBL, UniProtKB, PRO)
  * Families and gene/protein homology assertions (ENSFAM, PANTHER)
  * Chromosomes (ENSEMBL)

The ontology has the base http://purl.obolibrary.org/obo/omeo/

Omeo is actually a collection of modularized ontologies that are combined via importer ontologies. There is a core importer for each major species. Examples include:

  * http://purl.obolibrary.org/obo/omeo/Mus_musculus.owl
  * http://purl.obolibrary.org/obo/omeo/Homo_sapiens.owl
  * http://purl.obolibrary.org/obo/omeo/Danio_rerio.owl

Note that this google code repository does not store the ontologies themselves - it contains the code for building the ontologies, which are then deployed to the above URL.
To build the ontologies yourself, see the README in the build/ directory.

## OWL Modeling ##

We currently utilize class axioms only (Note: this has changed), but this is subject to change. One reason for the existing class representation is to facilitate reasoning with Elk. Once Elk supports instances, we may employ ABox axioms for some entities.

  * genes are subtypes of SO:gene
  * PRO is used for proteins (each protein is a class)
  * product\_of relations ( http://purl.obolibrary.org/obo/RO_0002205 ) connect genes to proteins (no RNA products as yet)
  * chromosomes are subtypes of SO:chromosome
  * equivalent classes axioms are used to connect gene classes between MOD, ENSEMBL and NCBIGene URI spaces