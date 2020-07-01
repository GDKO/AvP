## Create environment
- Conda
 - ```conda create --name avp python=3```
 - ```conda activate avp```
- Programs
  - ``` conda install -c bioconda mafft blast trimal fasttree iqtree```
- Python libraries
  - ```pip install numpy networkx pyyaml ete3 six biopython docopt pybedtools```


### Databases

#### Using NR
- ~~Alieness server route~~
- Alternate route
  - Run BLAST
    - ```blastp -query [proteins.fasta] -db [nr] -seg no -evalue 1e-3 -outfmt '6 std staxids' -out [blast.out]```

#### Make a diamond database from swissprot
- Install diamond from (https://github.com/bbuchfink/diamond/releases)
- Download swissport.gz database (ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/)
- Download taxdump.tar.gz from NCBI (ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)
  - ```mkdir taxdump```
  - ```tar xvf taxdump.tar.gz -C taxdump```
- Create taxid file with swissprot_to_acc2taxid.pl
  - ```swissprot_to_acc2taxid.py -i uniprot_sprot.fasta.gz > sp.taxids```
- Create Diamond db
  - ```diamond makedb --in uniprot_sprot.fasta.gz --taxonmap sp.taxids --db uniprot_sprot.fasta.dmnd --taxonnodes taxdump/nodes.dmp --taxonnames taxdump/names.dmp```


#### Using user-created databases
- swissprot
  - Run diamond
    - ```diamond blastp -d [database.dmnd] --max-target-seqs 500 -q [proteins.fasta] --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids > [diamond.out]```

#### Create AI features file with calculate_ai.py (for groups.yaml see next section)
- ```calculate_ai.py -i [blast.out|diamond.out] -g groups.yaml > [ai.out]```

### Run AvP

#### Create yaml files
- Create groups.yaml
  - ```TOI``` is for which is the target of the HGT
  - ```EGP``` is for which taxonomic groups to exclude from calculations
- In the following example we have proteins from the nematode _Meloidogyne incognita_ and we want to find HGTs from Non Metazoa species to our species. For that we set ```TOI``` to Metazoa and ```EGP``` to the suborder Tylenchida which our species belongs to, to allow for HGTs that may be present also in other Tylenchida species

```
---
TOI:
 33208: Metazoa
EGP:
 6300: Tylenchida
```
- The config yaml file should be filled with the database paths. The file also contains parameters that can be changed.

```
---
max_threads: 2

# DB path
nr_db_path: [nr]
sp_fasta_path: [swissprot.fasta]

## Algorithm options
# prepare
ai_cutoff: 0
percent_identity: 100
cutoffextend: 20    # when toi hit is found, we take this hit + n hits
threads: 2
trimal: false
min_num_hits: 4   # select queries with at least that many blast hits
percentage_similar_hits: 0.7  # group queries based on this
mode: nr    # use nr for nr database, use sp for swissprot database
# evaluate, clasify
fastml: true  # Use fasttree instead of IQTree
node_support: 0  # nodes below that number will collapse
complex_per_toi: 20   # if H/(H+T) smaller than this then node is considered T
complex_per_hgt: 80   # if H/(H+T) greater than this then node is considered H
complex_per_node: 90  # if node contains percent number of this category, it is assigned

# Program specific options
mafft_options: '--anysymbol --auto'
trimal_options: '-automated1'

#IQ-Tree
iqmodel: '-mset WAG,LG,JTT -AICc -mrate E,I,G,R'
ufbootstrap: 1000
iq_threads: 4
```

#### Prepare files for downstream analyses
- ```prepare_files.py -a [ai.out] -o [output_dir] -f [protein.fasta] -b [blast.out|diamond.out] -g groups.yaml -c config.yaml```
- The following folders will be created inside output dir
  - ```tmp``` temporary files
  - ```fastagroups``` fasta sequences for each group
  - ```mafftgroups``` aligned fasta sequences for each group
- The following file will be created inside output dir
  - ```tmp/taxonomy_nexus.txt``` taxonomic lineage for each protein
  - ```groups.tsv``` file specifying the groups

_**if fastml was set to false then substitute fasttree for iqtree**_

#### Create and evaluate trees for presence/absence of HGTs
- ```phylogeny_evaluation.py -i [output_dir]/mafftgroups/ -o [output_dir] -g [output_dir]/groups.tsv -t [output_dir]/tmp/taxonomy_nexus.txt -c config.yaml```
- The following folders will be created inside output dir
  - ```fasttree``` phylogenetic results in newick format for each group
  - ```fasttree_nexus``` phylogenetic results in nexus format for each group. These are best for visualisation with the a tree viewer program like FigTree
- The following file will be created inside output dir
  - ```fasttree_general_results.txt``` shows general results for the analyses
  - ```fasttree_tree_results.txt``` shows the result for each protein

### *(OPTIONAL)* Further analyses

#### Classify HGTs in taxonomic ranks
- Create classification file. The file ```sample.classificiation_toi_Metazoa.txt``` is provided within depot directory which will be useful for most analyses of Non Metazoa to Metazoa HGTs

```
#rank members
NM_Eukaryota	Eukaryota;Viridiplantae;Fungi
Prokaryota	Bacteria;Archaea
Viriods	Viroids
Viruses	Viruses
```
- ```classify_trees.py -i [output_dir]/fasttree_nexus/ -t [output_dir]/fasttree_tree_results.txt -f [classification.file] -c config.yaml -o [output_dir]```
- The following folders will be created inside output dir
  - ```classification``` Inside the folder are folders based on the classification file. Each folder contains the nexus for each gene classfied as such.
- The following file will be created inside output dir
  - ```classification_results.txt``` shows general results for the analyses
  - ```classification_tree_results.txt``` shows the result for each protein

#### Test alternative topology
- ```test_alt_topology.py -i [output_dir]/mafftgroups/ -t [output_dir]/fasttree_tree_results.txt -m fasttree -o [output_dir] -x 2```
- The following folders will be created inside output dir
  - ```alt_topology``` Contains all intermediate files for the phylogenetic analyses
- The following file will be created inside output dir
  - ```alt_topology_results.tsv``` results for the alternative topology analyses

#### Provide a score based on neighbouring genes in the genome
- TODO
