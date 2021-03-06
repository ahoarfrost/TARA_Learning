#download all prokaryote-enriched metagenome assemblies from TARA
##note cewd01_env.dat.gz does not seem to exist, but it has the identical sample
#ID as ceom01. See AssembliesToDownload.csv for sample IDs if you want to compare

urls="ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceng01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cent01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cenn01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cess01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceok01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cenp01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cewd01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceom01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cenf01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cenj01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceop01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceof01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cenx01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cenu01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cenz01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cenq01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cenw01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceny01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceno01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceoj01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceoq01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceoo01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceoc01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceqe01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cewa01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cepw01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceps01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cept01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cepv01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cepu01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevz01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceqd01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cepk01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cepl01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cepj01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cepa01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cetw01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceoy01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceos01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cepb01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceow01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceov01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceoz01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cepd01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cepc01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevy01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cerq01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceox01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cerw01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cerk01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cerm01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cerl01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cewb01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cerg01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceqw01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cere01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceqt01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cerb01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceqs01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceqj01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cequ01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceqh01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesa01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceqo01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesg01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesu01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesi01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesj01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cest01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cese01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cewr01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesf01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceso01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesh01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesw01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesn01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesd01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesy01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesp01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesv01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesq01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesr01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cewg01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesl01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cesm01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cewh01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cete01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceta01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cetq01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cetd01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cett01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cetk01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceti01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cewi01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cetl01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cetb01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cetr01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceuw01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cets01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cetm01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceuv01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevc01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceug01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cewj01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cewk01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceua01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cetz01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceue01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cety01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cetv01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceuh01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceuc01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceud01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceuf01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_ceub01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cetx01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevg01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevh01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevn01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cetu01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevi01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cewo01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevm01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevl01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevj01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevr01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevq01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevw01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevk01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevp01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevt01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevs01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevu01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevx01_env.dat.gz
ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/wgs/ce/wgs_cevo01_env.dat.gz"

keys="ceng01
cent01
cenn01
cess01
ceok01
cenp01
cewd01
ceom01
cenf01
cenj01
ceop01
ceof01
cenx01
cenu01
cenz01
cenq01
cenw01
ceny01
ceno01
ceoj01
ceoq01
ceoo01
ceoc01
ceqe01
cewa01
cepw01
ceps01
cept01
cepv01
cepu01
cevz01
ceqd01
cepk01
cepl01
cepj01
cepa01
cetw01
ceoy01
ceos01
cepb01
ceow01
ceov01
ceoz01
cepd01
cepc01
cevy01
cerq01
ceox01
cerw01
cerk01
cerm01
cerl01
cewb01
cerg01
ceqw01
cere01
ceqt01
cerb01
ceqs01
ceqj01
cequ01
ceqh01
cesa01
ceqo01
cesg01
cesu01
cesi01
cesj01
cest01
cese01
cewr01
cesf01
ceso01
cesh01
cesw01
cesn01
cesd01
cesy01
cesp01
cesv01
cesq01
cesr01
cewg01
cesl01
cesm01
cewh01
cete01
ceta01
cetq01
cetd01
cett01
cetk01
ceti01
cewi01
cetl01
cetb01
cetr01
ceuw01
cets01
cetm01
ceuv01
cevc01
ceug01
cewj01
cewk01
ceua01
cetz01
ceue01
cety01
cetv01
ceuh01
ceuc01
ceud01
ceuf01
ceub01
cetx01
cevg01
cevh01
cevn01
cetu01
cevi01
cewo01
cevm01
cevl01
cevj01
cevr01
cevq01
cevw01
cevk01
cevp01
cevt01
cevs01
cevu01
cevx01
cevo01"

for u in $urls
do
    curl -O "$u"
done