# http://biostar.stackexchange.com/questions/1289/disease-associated-snps
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A  -D hg18 -e '
select
   concat(left(title1,30),"..."),
   omimId,
   S.name,
   S.func,
   G.chrom,
   S.chromStart,
   S.chromEnd
from
   omimGene as G,
   omimGeneMap as M,
   snp130 as S
where
  G.name=M.omimId and
  G.chrom=S.chrom and
  S.chromStart>=G.chromStart and
  S.chromEnd <= G.chromEnd
limit 10;'
