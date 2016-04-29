#download OM-RGC and KO functional annotations from TARA website, and kegg modules

#Kegg KO annotations
curl -o TARA_KO_annotations.gz http://ocean-microbiome.embl.de/data/TARA243.KO.profile.release.gz

#Kegg Module annotations
curl -o TARA_KeggModule_annotations.gz http://ocean-microbiome.embl.de/data/TARA243.KO-module.profile.release.gz

#OM-RGC counts; link from S. Sunagawa
curl -o OMRGC_gene_profile.gz http://ocean-microbiome.embl.de/data/TARA243.gene.profile.release.gz