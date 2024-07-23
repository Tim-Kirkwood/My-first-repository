# SyntenyQC
## Motivation: 
Synteny plots are widely used for the comparison of genomic neighbourhoods, particularly in the context of genome mining.  Whilst synteny plots are typically included as part of larger software suites (e.g. antiSMASH ClusterBlast module), various low-code stand-alone tools are now available that allow users to source candidate neighbourhoods and build their own synteny plots.  However, a gap remains between: 

**(i) tools that source these candidate neighbourhoods** (e.g. `cblaster`, which can find hundreds of candidates), 

**(ii) tools that build the synteny plots** (e.g. 'clinker', which struggles as the number of neighbourhoods exceeds 30-50) and 

**(iii) the synteny plots themselves**, which become much harder to analyse/present as the number of neighbourhoods they include increases.

## Description: 
`SyntenyQC` is a python app for the curation of neighbourhoods immediately prior to synteny plot creation. `SyntenyQC` supports the systematic definition and annotation of candidate neighbourhoods based on a direct integration to `cblaster`, a popular tool for finding clustered homologs to a user-supplied query.  `SytenyQC` also offers a flexible method for objectively removing redundant neighbourhoods (sourced using `cblaster` or any other tool) prior to synteny plot creation.  This is in some cases an absolute requirement (e.g. `cblaster` called via the `CAGECAT` webserver places a limit of 50 neighbourhoods).  

## References
cblaster 
clinker
Cagecat
antiSMASH ClusterBlast

## Installation 
```
pip install SyntenyQC
```

## Usage
### General help:
```
>SyntenyQC -h
usage: SyntenyQC [-h] {collect,sieve} ...

options:
  -h, --help       show this help message and exit

subcommands:
  {collect,sieve}  Synteny quality control options
```
### Collect subcommand:
```
SyntenyQC collect -h
usage: SyntenyQC collect [-h] -bp -ns -em [-fn] [-sp] [-wg]

Write genbank files corresponding to cblaster neighbourhoods from a specified CSV-format binary file loacted at
BINARY_PATH.  For each cblaster hit accession in the binary file:

1) A record is downloaded from NCBI using the accession.  NCBI requires a user EMAIL to search for this record
   programatically.  If WRITE_GENOMES is specified, this record is written to a local file according to FILENAMES
   (see final bulletpoint).
2) A neighbourhood of size NEIGHBOURHOOD_SIZE bp is defined, centered on the cblaster hits defined in the binary
   file for the target accession.
3) (If STRICT_SPAN is specified:) If the accession's record is too small to contain a neighbourhood of the
   desired size, it is discarded.  For example, if an accession record is a 25kb contig and NEIGHBOURHOOD_SIZE
   is 50000, the record is discarded.
4) If FILENAMES is "organism", the nighbourhood is written to file called *organism*.gbk. If FILENAMES is
   "accession", the neighbourhood is written to *accession*.gbk. Synteny softwares such as clinker can use these
   filesnames to label synetny plot samples.

Once COLLECT has been run, a new folder with the same name as the binary file should be created in the directory
that holds the binary file (i.e. the file "path/to/binary/file.txt" will generate the folder "path/to/binary/file").
This folder will have a subdirectory called "neighbourhood", containing all of the neighbourhood genbank files
(i.e. "path/to/binary/file/neighbourhood"). If WRITE_GENOMES is specified, a second direcory ("genome") will also
be present, containing the entire record associated with each cblaster accession (i.e. "path/to/binary/file/genome").
Finally, a log file will be present in the folder "path/to/binary/file", containing a summary of accessions whose 
eighbourhoods were discarded.

options:
  -h, --help            show this help message and exit
  -bp, --binary_path
                        Full filepath to the CSV-format cblaster binary file containing neighbourhoods that should
                        be extracted
  -n, --neighbourhood_size
                        Size (basepairs) of neighbourhood to be extracted (centered on middle of CBLASTER-defined
                        neighbourhood)
  -em, --email          Email - required for NCBI entrez querying
  -fn, --filenames
                        If "organism", all collected files will be named according to organism. If "accession", all
                        files will be named by NCBI accession. (default:
                        organism)
  -sp, --strict_span
                        If set, will discard all neighbourhoods that are smaller than neighbourhood_size bp. For
                        example, if you set a neighbourhood_size of 50000, a 50kb neighbourhood will be extracted
                        from the NCBI record associateed with each cblaster hit. If the record is too small for this
                        to be done (i.e. the record is smaller then 50kb) it is discarded.
  -wg, --write_genomes
                        If set, will write entire NCBI record containing a cblaster hit to file (as well as just the
                        neighbourhood)
```
### Sieve subcommand:
```
>SyntenyQC sieve -h
usage: SyntenyQC sieve [-h] -gf [-ev] [-mi] [-mts] [-mev] -sf

Filter redundant genomic neighbourhoods based on neighbourhood similarity:
- First, an all-vs-all BLASTP is performed with user-specified BLASTP settings and the neighbourhoods in GENBANK_FOLDER.
- Secondly, these are parsed to define reciprocal best hits between every pair of neighbourhoods.
- Thirdly, these reciprocal best hits are used to derive a neighbourhood similarity network.  Nodes are neighbourhood
  filenames and edges indicate two neighbourhood nodes that have a similarity > SIMILARITY_FILTER.
  Similarity = Number of RBHs / Number of proteins encoded in smallest neighbourhood in pair.
- Finally, this network is pruned to remove neighbourhoods that exceed the user's SIMILARITY_FILTER threshold.
  Nodes that remain are copied to the newly created folder 'genbank_folder/ClusterSieve'.

options:
  -h, --help            show this help message and exit
  -g, --genbank_folder
                        Full path to folder containing neighbourhood genbank files requiring de-duplication
  -ev, --e_value    BLASTP evalue threshold. (default: 1e-05)
  -mi, --min_percent_identity
                        BLASTP percent identity threshold. (default: 50)
  -mts, --max_target_seqs
                        BLASTP -max_target_seqs. Maximum number of aligned sequences to keep. (default: 200)
  -mev, --min_edge_view
                        Minimum similarity between two neighbourhoods for an edge to be drawn betweeen them in the RBH
                        graph. Purely for visualisation of the graph HTML file - has no impact on the graph pruning
                        results. (default: None)
  -sf, --similarity_filter
                        Similarity threshold above which two neighbourhoods are considered redundant
```
#### Sieve pruning algorithm:
```
Data: RBH graph
Result: Nodes from pruned RBH graph
Procedure:
    while max(node degrees in RBH graph) > 0:
        delete nodes = []
        for node in RBH graph:
            if node degree == max(node degrees in RBH graph):
                delete nodes + node
        delete node = random node from delete nodes
        RBH graph = RBH graph - delete node 
    return nodes in RBH graph 
```

# Example use:
## (1) Preliminaries
Take BGC from MIBIG, and and analyse with CBLASTER via the command line or CAGECAT web server. ⚠️**WARNING: cblaster must be run with these settings - **⚠️  
## (2) Collect neighbourhoods  
### Starting directory structure:
```
folder/with/binary.csv
```
        
### Command (neighbourhood size kb = 2 x BGC length):
     
```
SyntenyQC collect -bp path/to/BGC0000194_binary.txt -ns 42566 -em my_email@domain.com -fn organism -sp -wg
```
   
### Finishing directory structure: 
           
```
folder/with/binary/neighbourhood/organism1.gbk, organism2.gbk...organism157.gbk
                  /genomes      /organism1.gbk, organism2.gbk...organism157.gbk         ###only if -wg!
                  /log.txt
```
#### 🔴 157 neighbourhoods is a lot for a synteny plot 🔴


## (3) Sieve neighbourhoods
### Starting directory structure 
From `SyntenyQC collect` in this example, but any folder with at least one genbank file can be used:
```
folder/with/binary/neighbourhood/organism1.gbk, ...
```
### Command:
```
SyntenyQC sieve -gf folder/with/binary/neighbourhood -sf 0.7
```
### Finishing directory structure: 
```
folder/with/binary/neighbourhood/organism1.gbk, ...
                                /ClusterSieve/organism1.gbk, ...organism38.gbk
                                             /log.txt
```
#### :green_heart: 38 neighbourhoods is OK for a synteny plot :green_heart:
(note, filenames are to show number of files - `neighbourhood/ClusterSieve/organism1.gbk` is in the `neighbourhood` folder, but `neighbourhood/organism1.gbk` may be different to `neighbourhood/ClusterSieve/organism1.gbk`)
