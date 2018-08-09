Tuturial on the use of VIRULIGN
================


This pages decribes the instructions for download and use of VIRULIGN for viral sequence data, with an explanation of all optional parameters. This repository also contains a PDF version of this tutorial.  In this document, three example datasets (Dengue virus, Zika virus and HIV-1) are described where VIRULIGN was used to generate research-relevant output formats of the constructed codon-correct multiple sequence alignments.  


Rationale
---------------------------------------


Virus sequence data are an essential resource for reconstructing
spatiotemporal dynamics of viral spread as well as treatment and
prevention strategies. However, the potential benefit of using sequence
data for these applications critically depends on the accuracy and the
correct annotation of these alignments of genetically diverse data. In
particular, coding sequences of viral pathogens should be analyzed in
their corresponding open reading frame (ORF) to fully utilize their
biological information. Therefore, while the construction of multiple
sequence alignments (MSAs) can be done with a range of sequence
alignment software, MSAs of virus coding sequences in the correct
reading frame and annotated with respect to the proteins encoded in the
genome are more difficult to achieve.\
VIRULIGN is an easy-to-use command line application to construct
codon-correct alignments of large virus sequence datasets. Additionally,
VIRULIGN has support for standardized genome annotation and implements
various alignment export formats that are useful for various research
applications. VIRULIGN is an open-source project written in the C++
programming language and available under the GPLv2 license.\
VIRULIGN operates by aligning each target sequence (i.e., $t \in T$) of
the input file codon-correctly against the reference sequence ($r$).
Subsequently a multiple sequence alignment $MSA(r,T)$ is constructed
based on all codon-correct (cc) pairwise aligned target sequences
$A_{cc}(r,t)$ (Figure below).


![alt text](https://github.com/rega-cev/virulign-tutorial/blob/master/figures/Figure1.png)

 

The most recent version and executable of VIRULIGN can be downloaded
from the [GitHub project web
page](https://github.com/rega-cev/virulign):

    https://github.com/rega-cev/virulign

Instructions are provided to build and install the software. VIRULIGN is
a cross-platform application and has been tested on GNU/Linux, MacOS and
Windows.

Use and features
================

### Basic command 

VIRULIGN minimally requires a FASTA file with target sequences and a
reference sequence in order to generate a codon-correct alignment in a
predefined output format (see below). The reference sequence can be
either provided in FASTA format or embedded in an XML file (see below).
The default command for VIRULIGN is as follows:

    $ virulign [reference.fasta orf-description.xml] sequences.fasta

### Optional arguments 

Additional parameters can be specified to configure the alignment
construction and to export the alignment to a variety of output formats.
In case that a parameter has not been explicitly specified, the first
value of this optional parameter is used as the default value. The
following parameters can be used to configure the alignment and its
representation [^1]:

    $ virulign [reference.fasta orf-description.xml] 
               sequences.fasta
        --exportKind [Mutations PairwiseAlignments 
                      GlobalAlignment  PositionTable 
                      MutationTable] 
        --exportAlphabet [AminoAcids Nucleotides]
        --exportWithInsertions [yes no]
        --exportReferenceSequence [no yes]
        --gapExtensionPenalty doubleValue=>3.3
        --gapOpenPenalty doubleValue=>10.0
        --maxFrameShifts intValue=>3
        --progress [no yes]
        --nt-debug [dir]

        Output: The alignment will be printed to standard out and any progress 
        or error messages will be printed to the standard error. 
        This output can be redirected to files,    e.g.: virulign ref.xml 
        sequence.fasta  > alignment.mutations 2> alignment.err

To print these options, invoke the virulign command without any arguments.

### Redirecting output to files 

The alignment will be printed to standard out and any progress or error
messages will be printed to the standard error. The output of VIRULIGN
can be redirected to files. For example, the output is redirected to a
file with the following command.

    $ virulign [reference.fasta orf-description.xml] sequences.fasta 
    > output.file 

In case any progress or error messages should be suppressed, make the
following extension to your command.

    $ virulign [reference.fasta orf-description.xml] sequences.fasta 
    > output.file 2> error.file

Description of the optional parameters 
---------------------------------------

### exportKind 

The parameter `--exportKind` defines the output type of the alignment,
either a FASTA sequence file or a CSV mutation file. To display the
different options, consider this example FASTA input file including
sequences of different length and a full-length reference sequence for
comparison.

    >Ref
    CCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAA
    >Seq1
    ATTGACACTGTACCAGTAACATTAAAGCCAGGAATGGATGGACCAAAG
    >Seq2
    CCTATGGAAACTGTGCCAGTAAAATTAAAGCCAGGAATGGAT
    >Seq3
    CTCATTAGTCCTATTAGTGTAAAATTAAAACCAGGAATGGATGGCCCAAGG
    >Seq4
    AGTCCCATTGAAACTGTACCAGTAAAAGGAGATGGCCCAAAG

The option `GlobalAlignment` will generate a FASTA file of the target
sequences aligned against a single reference sequence and formatted as a
MSA.

    >Ref
    CCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAA
    >Seq1
    ------------ATTGACACTGTACCAGTAACATTAAAGCCAGGAATGGATGGACCAAAG
    >Seq2
    CCTATG---------GAAACTGTGCCAGTAAAATTAAAGCCAGGAATGGAT---------
    >Seq3
    CTCATTAGTCCTATTAGT---------GTAAAATTAAAACCAGGAATGGATGGCCCAAGG
    >Seq4
    ------AGTCCCATTGAAACTGTACCAGTAAAA---------GGA---GATGGCCCAAAG

The option `PairwiseAlignments` will generate a FASTA file of the target
sequences, with each sequence aligned separately against the reference
sequence.

    >Ref
    CCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAA
    >Seq1
    ------------ATTGACACTGTACCAGTAACATTAAAGCCAGGAATGGATGGACCAAAG
    >Ref
    CCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAA
    >Seq2
    CCTATG---------GAAACTGTGCCAGTAAAATTAAAGCCAGGAATGGAT---------

The option `PositionTable` will create a comma-separated value (CSV)
file where each position of the alignment is given as a separate column.
The CSV file is annotated according to the numerical positions in the
protein (e.g., Table below).

|ID | Pos1 | Pos2 | Pos3 | …
| ------ |-------| -----| -----| -----|
|Ref |P | I | S | …
|Seq1 | -|- | - | - | …
|Seq2 | P | M | - | …
|Seq3 | L | I | S | …

The option `MutationTable` will create a CSV file where each mutation
present at a specific position is given as a separate column in Boolean
representation. The CSV file is annotated according to the numerical
position in the protein (e.g., Table below).

 
|ID | Mut1P | Mut1L | Mut2I | Mut2M | …
| ------ |-------| -----| -----| -----|-----|
|Ref | y | n | y | n | …
|Seq1 | n | n | n | n | …
|Seq2 | y | y | n | y | …
|Seq3 | n | y | y | n | …

-   the option `Mutations` will output, for each sequence, a list of
    amino acids changes compared to the reference sequence (e.g., Table below).

 
|ID | Mutations |
|Seq1 | -
|Seq2 | 2M |
|Seq3 | 1L |

### exportAlphabet

The parameter `--exportAlphabet` defines the alphabet in which the
alignment is generated.

-   the option `AminoAcids` will export an alignment with the
    translation of the nucleotide codons.

-   the option `Nucleotides` will export an alignment with nucleotides

### exportWithInsertions

The parameter `--exportWithInsertions` determines whether insertions can
be added to the reference sequence.

-   the option `yes` will insert gaps into the reference sequence to
    accommodate the identification of codon insertions in the target
    sequences.

-   the option `No` removes codon insertions in the reference sequence
    that were generated during the alignment procedure.

### exportReferenceSequence

The parameter `--exportReferenceSequence` controls whether the reference
sequence is to be added to the alignment (yes/no).

### gapExtensionPenalty

The parameter `--gapExtensionPenalty` defines the value of the penalty
to extend an existing gap.

### gapOpenPenalty

The parameter `--gapOpenPenalty` defines the value of the penalty to
start a new gap.

### maxFrameShifts

The parameter `--maxFrameShifts` defines the maximum number of
frame-shifts allowed.

### progress

The parameter `--progress` allows to monitor the estimated time until
completion of the alignment.

When the option `yes` is used, a progress message stating the number and
percentage of sequences aligned as well as the estimated time left to
finalize the pending alignment is shown.

### nt-debug

The parameter `--nt-debug` allows to visualise sequences that could not
be aligned by VIRULIGN. When the option `dir` is used, pairwise sequence
alignments of failed target sequences and the reference sequence are
stored in the directory with the name ’dir’. This directory needs to be
created before the execution of the VIRULIGN command. This feature
allows to inspect each failed target sequence individually to understand
why the target sequence did not pass the quality control of the
alignment. Subsequently, errors in the target sequence can be corrected.

A reference sequence can be either provided to VIRULIGN in FASTA format
or embedded in an XML file. In this XML file, also an annotation of the
different proteins, regions or other structures can be given by the
positions relative to the reference genome.

For illustrational purposes, we present a toy XML file to align virus
sequence data against a reference sequence and the annotation of the
proteins A, B and C in the genome of VIRUS. Later in this document, we
will consider and use some more realistic annotations (i.e., ZIKV,
HIV-1).

    <?xml version="1.0" encoding="UTF-8"?>
        <orf name="VIRUS" 
             referenceSequence="atgaaaaacccaaaaaagaaatccgga" >
           <protein abbreviation="A" 
                    startPosition="1" stopPosition="7" />
           <protein abbreviation="B" 
                    startPosition="7" stopPosition="13" />
           <protein abbreviation="C" 
                    startPosition="13" stopPosition="17" />
        </orf>

Table below shows an example of an alignment that is
constructed with this XML file as reference, and is exported using the
tabular format.


|ID | A_1 | A_2 | … | B_1 | …
| ------ |-------| -----| -----| -----|-----|
|Seq1 | X | Y | … | C | …
|Seq2 | X | T | … | F | …
|Seq3 | X | Y | … | D | …
|Seq4 | R | M | … | D | …

### Converting Genbank file to XML file 

Currently, the direct use of a Genbank XML file is not supported, as we
aim for well-curated annotations of reference genomes for all virus
pathogens. However, the creation of a XML file can be based on
information provided within a Genbank description of the respective
(reference) genome. To facilitate this transfer for the user community,
a python script has been made to extract relevant information from a
Genbank description file into the format of a VIRULIGN annotation file,
which can be found at the [GitHub repository of VIRULIGN
tools](https://github.com/rega-cev/virulign-tools).

    https://github.com/rega-cev/virulign-tools

The command to run the script file for the conversion of the Genbank XML
file is:

    $ python genbank_to_virulign.py genbank_insdseq.xml 
    orf-name seq-start seq-end
     

The parameter `orf-name` will set the name for the defined ORF, while
parameters `seq-start` and `seq-end` define the start and stop
nucleotide position of the respective ORF in the reference genome.

To demonstrate the working of this script with a relevant example, we
downloaded the Genbank INSDSEQ XML file for the reference genome
NC\_001477 ([link](www.ncbi.nlm.nih.gov/nuccore/NC\_001477)) of Dengue Serotype 1. The Genbank XML file is available
at the [tutorial web page](https://github.com/rega-cev/virulign-tutorial). The command needed
for conversion is:

    $ python genbank_to_virulign.py NC_001477.gbc.xml  DENV1 95 10273
     

The output of this command may need to be improved into a correct
annotation file for VIRULIGN. For example, the presence of annotations
for the precursor or polyprotein together with the separate proteins can
be conflicting as they overlap in genome positioning.

    < ... "membrane glycoprotein precursor M" startPosition="343" stopPosition="841" />
    < ... "protein pr" startPosition="343" stopPosition="616" />
    < ... "membrane glycoprotein M" startPosition="616" stopPosition="841" />
    < ... "envelope protein E" startPosition="841" stopPosition="2326" />

Examples
========

We demonstrate the use of VIRULIGN by constructing sequence alignments
for three viral pathogens that are the causative agents for major
epidemics: HIV-1, Dengue virus serotype 1 (DENV-1) and Zika virus
(ZIKV). For each pathogen, a different feature of VIRULIGN is
demonstrated.\
Virus sequence datasets were collected from public databases (i.e.,
Genbank and the Stanford HIV Drug Resistance Database) and passed to
different alignment software applications for evaluation. We used a
minimum number of processing steps as possible between dataset retrieval
and alignment input to clearly illustrate the strength of VIRULIGN. All
data files of these examples can be found on the [tutorial web
page](https://github.com/rega-cev/virulign-tutorial):

    https://github.com/rega-cev/virulign-tutorial

DENV-1 alignment
----------------

We compare the output of VIRULIGN against three popular alignment tools
(MAFFT, MUSCLE and Clustal Omega) in their ability to construct an
accurate codon-correct alignment from a genome sequence dataset.

### Sequence dataset

Genome sequence data of DENV-1 (i.e., Dengue Serotype 1) were collected
from the Dengue Virus Variation Database ([link](www.ncbi.nlm.nih.gov/genomes/VirusVariation)) [@hatcher2017]. Only
full-length nucleotide sequences originating from a human host were
retained and identical sequences were collapsed. From a total of 3539
genome sequences, the corresponding serotype information was used to
select a subset of 1432 DENV-1 isolates.\
The input FASTA file ’denv-1.fasta’ can be found in the [tutorial Dengue
folder](https://github.com/rega-cev/virulign-tutorial/examples-alignments/DENV):

    https://github.com/rega-cev/virulign-tutorial/examples-alignments/DENV/

### Reference sequence

The NCBI Reference Sequence for DENV-1 is NC\_001477 ([link](www.ncbi.nlm.nih.gov/nuccore/NC\_001477)). Find a FASTA
file ’NC\_001477.fasta’ that contains this reference sequence in
[tutorial Dengue
folder](https://github.com/rega-cev/virulign-tutorial/examples-alignments/DENV):

    https://github.com/rega-cev/virulign-tutorial/examples-alignments/DENV/

### Alignments by different tools

Alignments were constructed using the default or recommended parameters
for each tool. The following versions were downloaded; VIRULIGN (v1.0),
MAFFT (v7.313) [@katoh2014], MUSCLE (v3.8.31) [@edgar2004] and Clustal
Omega (v1.2.3) [@sievers2011]. MUSCLE was used with the additional
option `-diags`, which is intended for alignments of highly similar
sequences. No additional parameters were used for the other programs,
although for individual cases, the use of specific parameters could
affect the speed or accuracy of the alignment construction process.\
The genome sequence of NC\_001477 was added to the target dataset to
facilitate the trimming of constructed alignments to the boundaries of
the reference CDS, in order to remove the alignment of the 5’/3’
untranslated regions. The codon-correctness of the alignment was then
evaluated by visually inspecting the amino acid translation of the
respective alignments.\
The following commands were used to obtain alignments:

    $ mafft --auto denv-1.fasta > denv-1-mafft.fasta
     
    $ muscle -maxiters 1 -diags -in denv-1.fasta -out denv-1-muscle.fasta

    $ clustalo --auto -i denv-1.fasta -o denv-1-clustalo.fasta

    $ virulign NC_001477.fasta denv-1.fasta  
        --exportKind GlobalAlignment 
        --exportAlphabet Nucleotides > denv-1-virulign.fasta

Each alignment was then trimmed to the start and stop position of the
coding sequence of the reference sequence, the size of this coding
region is 10188 nucleotides. All trimmed alignments can be found in
[tutorial Dengue
folder](https://github.com/rega-cev/virulign-tutorial/examples-alignments/DENV):

    https://github.com/rega-cev/virulign-tutorial/examples-alignments/DENV

From this evaluation, it can be observed that VIRULIGN is able to handle
insertions and deletions without disrupting the reading frame and
resulting in the absence of stop codons within the alignment, while
maintaining quality of the alignment. Figure below visualises
a selected window from constructed alignments to illustrate the
codon-correctness of VIRULIGN. We recorded the time needed for each
alignment construction (Table below). (Performed on a 3.6 GHz Intel Core i7 CPU with 12 GB of RAM, where each application had access to 1 CPU core.) This evaluation
shows that VIRULIGN is able to obtain these results while still being
computationally competitive with MAFFT.


|Command | Alignment | Run time 
| ---- |---- |---- |
|mafft | denv-1-mafft.fasta | 2m33s|
|muscle | denv-1-muscle.fasta | 58m|
|clustalo | denv-1-clustalo.fasta | 760m|
|virulign | denv-1-virulign.fasta | 19m46s|

  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- --
   ![ MAFFT (top window) and VIRULIGN (bottom window) constructed alignments visualized with SeaView. Stop codons are denoted by the presence of an asterisk.[]{data-label="fig:compare"}](Figure2.png "fig:")  
  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- --

### Annotation file

While the reference sequence for DENV-1 was provided by means of a
simple FASTA file, VIRULIGN can also be used with an XML file containing
both the reference sequence and the protein annotation of the reference
genome. To construct the same alignment as above, with an annotated
reference sequence, use this command:

    $ virulign DENV1-NC001477.xml denv-1.fasta  
        --exportKind GlobalAlignment 
        --exportAlphabet Nucleotides  > denv-1-virulign.fasta

XML annotation files for each of the four DENV serotypes are available
in the [VIRULIGN references
folder](https://github.com/rega-cev/virulign/references):

    https://github.com/rega-cev/virulign/references

More information on this XML feature is available in the section
\[features\] Use and Features. In the next section, we demonstrate how
this XML annotation file can be used to obtain alignments directed
towards specific research applications.\

In 2015, ZIKV caused a worldwide public health emergency, resulting in
an intensive community effort to identify genomic correlates of disease
manifestations of microcephaly and other neurological complications. As
we have shown recently [@theys2017], the rapid advance in ZIKV genomics
resulted in inconsistencies that complicate the interpretation,
reproducibility and comparison of findings from and across studies,
particularly due to the lack of a consensus on the standardized and
representative reference annotation. ZIKV reference genomes did not
match virus strains sampled from the global epidemic or showed high
level of heterogeneity in reported peptide lengths across their genome
annotations.\
To mitigate these concerns, we provided a correction with respect to the
NCBI reference sequence NC\_012532 for ZIKV (Figure below).
More information on the corrected reference sequence can be found [at
the Rega ZIKV reference sequence
website](https://rega.kuleuven.be/cev/reference-sequences/rega-zikv) ([Link](rega.kuleuven.be/cev/reference-sequences/rega-zikv)).

  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- --
   ![Accurate ZIKV genome positions of the different proteins based on corrections made in the NCBI reference sequence NC\_012532.[]{data-label="fig:annotation"}](Figure3.png "fig:") 
  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- --

This example shows how the functionality of the XML configuration file,
describing the genome annotation for all proteins, can greatly simplify
the analysis to find associations of genomic features with clinically,
epidemiologically or evolutionary parameters. We show that VIRULIGN
makes this possible, while keeping manual processing limited to a
minimum.

In particular, we replicate the evidence that indicated the necessity to
correct the reference genome that was proposed by the NCBI. Extensive
variation of a N-glycosylation motif in the Envelope (E) protein can be
observed between different hosts, different viral lineages and even
within a set of virus genomes derived from the historical MR766 strain.

### Sequence dataset

A specific subset of 19 (near-)complete ZIKV genomes was collected from
Genbank, the resulting FASTA file ’zikv.fasta’ can be found in the
[tutorial Zika
folder](https://github.com/rega-cev/virulign-tutorial/examples-alignments/ZIKV):

    https://github.com/rega-cev/virulign-tutorial/examples-alignments/ZIKV/

### Reference sequence

The reference genome sequence and the corresponding protein annotation
for ZIKV can be found in the XML configuration file. This XML file
’ZIKV-rega.xml’ can be found in the [VIRULIGN references
folder](https://github.com/rega-cev/virulign/references):

    https://github.com/rega-cev/virulign/references

### Motif investigation

Previous literature has shown variability of the glycosylation motif
around positions 150 - 165 in the E protein, and this has been suggested
to result from excessive *in vitro* passaging [@theys2017]. We used
VIRULIGN to create a position table of the amino acids in the alignment,
annotated according to the respective protein. The following command was
used:

    virulign ZIKV-rega.xml zikv.fasta 
        --exportKind PositionTable 
        --exportReferenceSequence yes  > zikv-aligned.csv

The resulting CSV file contains a column for each position in the
genome, and can be found in the [tutorial Zika
folder](https://github.com/rega-cev/virulign-tutorial/examples-alignments/ZIKV):

    https://github.com/rega-cev/virulign-tutorial/examples-alignments/ZIKV/

A simple R script shows the variability of the glycosylation motif
across the different virus variants.

    # import the alignment file
    data<-read.csv('zika-aligned.csv',header=TRUE)
    # determine the positions of the motif
    reg<-match('E_151',names(data)):match('E_163',names(data))
    # show the motif sequences
    tidyr::unite(data[,c(1,reg)],Motif,-1,sep='') 

The relevant region in the CSV file looks like:

    id,E_151,E_152,E_153,E_154,E_155,E_156,E_157,E_158,E_159,E_160,E_161,E_162,E_163
    Ref,M,I,V,N,D,T,G,H,E,T,D,E,N
    KF268948,M,I,V,N,D,I,G,H,E,T,D,E,N
    KF268949,M,I,V,N,-,-,-,-,-,-,D,E,N
    KF268950,M,I,V,N,D,I,G,H,E,T,D,E,N
    KU955595,M,I,V,N,D,T,G,H,E,T,D,E,N
    KY014323,M,I,V,N,D,T,G,H,E,T,D,E,N
    KU963574,M,I,V,N,-,-,-,-,-,-,D,E,N
    HQ234500,M,I,V,N,-,-,-,-,-,-,D,E,N
    KX369547,M,I,V,N,D,T,G,H,E,T,D,E,N
    ....

When additional meta-data is given as well, this analysis clearly
illustrates the presence of a VNDT motif in viruses sampled from the
recent epidemic, and an independence of the deletion regarding the host,
year of collection and viral lineage (Figure below).

  ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- --
   ![Variation in the glycosylation motif (red) of the ZIKV E protein across different hosts, virus lineage as well as date and country of collection. []{data-label="lbl:motif"}](Figure4.png "fig:")  
  ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- --

### Scripts to select individual proteins

An alignment CSV table, that includes protein and position data in the
header, can be easily processed by external tools. One could easily
develop their own scripts to operate on this format, but many
interesting manipulations can be done using default command line tools
as well.

As an example, based on the whole-genome ZIKV alignment above, we can
easily select a particular protein (e.g., the NS3 protein), using
csvkit ([link](https://github.com/wireservice/csvkit)).

    # define index of the NS3 headers
    ns3_headers=`csvcut -n ZIKA-pos.csv  | grep NS3 | cut -d ":" -f 2`
    # use comma as separator between header columns
    ns3_headers=`echo $ns3_headers | tr ' ' ','`
    # extract the NS3 region from the whole-genome alignment
    csvcut -c "seqid,${ns3_headers}" ZIKA-pos.csv  > ZIKA-NS3-pos.csv

HIV-1 alignment: gag
--------------------

The genome structure of HIV-1 is characterized by three ORFs, where each
frame determines different genes encoding the viral proteins. The
gag/pol/env gene organization, as for other retroviruses, encodes for
important structural proteins and enzymes, which are first translated as
large poly-proteins. Gag and pol have overlapping ORFs, requiring a
ribosomal frame shift to reveal the pol ORF. HIV-1 gag encodes for
several structural proteins and is considered as a potential target for
antiretroviral treatment [@tedbury2015hiv].

### Sequence dataset

HIV sequences were obtained from a large-scale analysis of HIV-1
diversity [@li2015], where 2966 sequences from different HIV-1 subtypes
have been analyzed. The FASTA input file ’hiv.fasta’ can be found in the
[tutorial HIV-1
folder](https://github.com/rega-cev/virulign-tutorial/examples-alignments/HIV-1):

    https://github.com/rega-cev/virulign-tutorial/examples-alignments/HIV-1/

### Reference sequence

The HXB2 sequence
[NC\_001802](https://www.ncbi.nlm.nih.gov/nuccore/9629357)[link](www.ncbi.nlm.nih.gov/nuccore/NC\_001802) was used
as the reference genome. An XML file with the corresponding coding
sequence and annotation is available for the different polyproteins of
relevance. The respective files ‘HIV-HXB2-env.xml‘, ‘HIV-HXB2-gag.xml‘
and ‘HIV-HXB2-pol.xml‘ are available in the [VIRULIGN references
folder](https://github.com/rega-cev/virulign/references):

    https://github.com/rega-cev/virulign/references

For this example, we used the file ‘HIV-HXB2-gag.xml‘.

### Alignment

To align the gag sequences, we used the following command

    virulign HIV-HXB2-gag.xml HIV.fasta 
        --exportKind GlobalAlignment 
        --exportAlphabet Nucleotides 
        --exportReferenceSequence yes  > HIVgag.fasta

An option of VIRULIGN can be used to avoid insertions towards the
reference sequence from being exported. This feature can be used to
inspect the quality of the sequence dataset when insertions are sparse
throughout the dataset or when insertions are not expected.

    virulign HIV-HXB2-gag.xml HIV.fasta 
      --exportKind  GlobalAlignment  
      --exportAlphabet Nucleotides  
      --exportWithInsertions no  
      --exportReferenceSequence yes  > HIVgag-NoInsertions.fasta

HIV-1 alignment: pol
--------------------

A second example of HIV-1 alignment is directed towards drug resistance
detection, which is still a major need for successful treatment, in
particular in developing countries as a result of the up-scale of
antiretroviral treatment. Therefore, the identification, understanding
and interpretation of resistance mutations remains an important research
topic. The pol polyprotein is cleaved into three viral enzymatic
proteins (protease, reverse transcriptase, and integrase), each of which
is an important drug target.

### Sequence dataset

We downloaded a large set of reverse transcriptase sequences (N=111223)
from the Stanford University HIV Drug Resistance Database [@Rhee2003],
heterogeneous in length and mapping of the complete reverse
transcriptase region. The resulting FASTA file ’HIVdb.fasta’ can be
found in the [tutorial HIV-1
folder](https://github.com/rega-cev/virulign-tutorial/examples-alignments/HIV-1):

    https://github.com/rega-cev/virulign-tutorial/examples-alignments/HIV-1/

### Reference sequence

The HXB2 sequence
[NC\_001802](https://www.ncbi.nlm.nih.gov/nuccore/9629357) was used
a reference genome. An XML file with the corresponding coding sequence
and annotation is available for each specific ORF. The respective files
‘HIV-HXB2-env.xml‘, ‘HIV-HXB2-gag.xml‘ and ‘HIV-HXB2-pol.xml‘ are
available in the [VIRULIGN references
folder](https://github.com/rega-cev/virulign-tutorial/references):

    https://github.com/rega-cev/virulign/references

For this example, we use the file ‘HIV-HXB2-pol.xml‘.

### Alignment

    virulign HIV-HXB2-pol.xml HIVdb.fasta 
      --exportKind  GlobalAlignment  
      --exportAlphabet Nucleotides  
      --exportWithInsertions no  
      --exportReferenceSequence yes  > HIVrt.fasta

111189 sequences could be aligned by VIRULIGN, and we visually inspected
the quality of the alignment. The constructed alignment can be used as
input for different applications to investigate drug resistance
mutations identification and interpretation. Thanks to VIRULIGN’s
computational complexity, our new method is able to deal well with large
dataset what is reflected in favorable run-times for this particular
analysis. VIRULIGN performed this alignment in 49 minutes, while it took
MAFFT 10 hours and 50 minutes. (Performed on the same hardware configuration as declared earlier.
)

### Detection of frameshift errors

34 sequences could not be aligned by VIRULIGN, which identifiers and
sequences are respectively stored in the files
‘HIVdb-errorsequences.txt‘ and ‘HIVdb-errorsequences.fasta‘ in the
[tutorial HIV-1
folder](https://github.com/rega-cev/virulign-tutorial/examples-alignments/HIV-1).
The identifiers of these sequences can be obtained by using the unix
command ’diff’ on the headers of the input and output files, and
subsequently the sequences can be extracted from the input file using
these identifiers.

However, the VIRULIGN parameter `--nt-debug` can be used to
automatically redirect failed sequences to a folder, which should be
created prior to command execution. Adding this parameter to the command
used above gives the following:

    virulign HIV-HXB2-pol.xml HIVdb.fasta 
      --exportKind  GlobalAlignment  
      --exportAlphabet Nucleotides  
      --exportWithInsertions no  
      --exportReferenceSequence yes  
      --nt-debug Failed > HIVrt.fasta

The directory ‘Failed‘ in the [tutorial HIV-1
folder](https://github.com/rega-cev/virulign-tutorial/examples-alignments/HIV-1)
contains pairwise sequence alignments of each failed target sequence
with the reference sequence. To demonstrate here the debug feature of
VIRULIGN without having to re-align the entire dataset again, we aligned
the subset of 34 sequences with VIRULIGN using the following command

    virulign HIV-HXB2-pol.xml HIVdb-errorsequences.fasta 
      --exportKind  GlobalAlignment  
      --exportAlphabet Nucleotides  
      --nt-debug Failed

Each alignment in the folder ‘Failed‘ can be subjected to closer
inspection in order to investigate the reason for the exclusion of the
sequence from the alignment. Figure below shows an example of
sequence 64344 which failed to be included in the final MSA.

  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ --
   ![Example of a failed sequence due to the presence of a stop codon and additional single nucleotide deletions requiring frameshift corrections.[]{data-label="fig:failed"}](Figure5.png "fig:") 
  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ --

### Detection of drug resistance mutations

To further demonstrate that VIRULIGN provides accurate and fast
codon-aware sequence alignments, we downloaded the set of sequences from
Genbank, that were described in a recent publication by Chaplin et al.
(2018) [@chaplin2018]. This study analyzed the distinct patterns of
thymidine analogue mutations with K65R in HIV-1 patients failing
tenofovir-based antiretroviral therapy. We have added the sequence file
‘chaplin2018-sequences.fasta‘ to the [tutorial HIV-1
folder](https://github.com/rega-cev/virulign-tutorial/examples-alignments/HIV-1).\
We aligned the set of sequences with VIRULIGN using the following
command

    virulign  HIV-HXB2-pol.xml chaplin2018-sequences.fasta
        > chaplin2018-mutations.csv 

We counted the number of occurrences of NRTI mutations K65R and K103N in
this dataset: respectively 92 times the K65R mutation was detected and
145 times the K103N mutation. This number matches exactly the mutation
frequencies obtained using the Stanford HIVdb pipeline that was
used in the study of Chaplin et al., supporting the confidence of
VIRULIGN to retrieve accurate alignments.

For comparison purposes, we also constructed a MSA of this set of
sequences with VIRULIGN and with MAFFT, and subsequently trimmed to the
first position of RT. Figure below provides a visual
inspection of the two constructed alignments, and an illustration of the
codon-awareness of VIRULIGN. It can be seen that insertions cause
frameshifts and the inclusion of stop codons in the MAFFT alignment,
while VIRULIGN accommodates these insertions without disturbing the
reading frames of the alignment.

With MAFFT:

    mafft  --auto chaplin2018-sequences.fasta 
         > chaplin2018-sequences-mafft.fasta

With VIRULIGN:

    virulign  HIV-HXB2-pol.xml chaplin2018-sequences.fasta 
        --exportKind GlobalAlignment 
        --exportAlphabet Nucleotides 
        --exportWithInsertions no 
        --exportReferenceSequence yes  > chaplin2018-sequences-virulign.fasta

  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ --
   ![Comparison of MSA’s constructed with MAFFT (upper panel) and VIRULIGN (lower panel), translated into amino acids and visualised with SeaView.[]{data-label="fig:chaplin"}](Figure6.png "fig:")
  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ --

List of case study examples
---------------------------

VIRULIGN has been used for a large number of analyses with respect to
virus genomics. We provide a non-exhaustive list of examples:

-   Identification of HIV-1 drug resistance mutations and the pathways
    emerging under drug selective pressure, as well as modeling the
    different factors leading to treatment failure and trends over time
    [@ngcapu2017; @theys2013; @vercauteren2013; @vercauteren2008; @deforche2008a; @deforche2008b].

-   Large-scale analysis of HIV-1 and HCV sequence datasets to explore
    genetic diversity at population level and to map structural and
    functional factors that shape viral evolution
    [@abecasis2013; @li2015; @cuypers2015; @cuypers2016].

-   Evaluation of the annotation and representativeness of current
    reference genomes, and support of correcting the NCBI reference
    sequence for the Zika virus [@theys2017].

-   Web-application to support surveillance and tracing of viral
    outbreaks (e.g., HIV-1 and DENV), necessitating efficient analysis
    of large sequence databases and phylogenetic trees [@libin2017].

-   VIRULIGN was recently integrated into the RegaDB data management and
    analysis platform for the clinical follow-up of HIV-1 patients
    [@vercauteren2013; @libin2013].

-   Evaluation of an automated framework for the virus typing (HIV-1,
    Dengue, Zika and other viruses) and resistance interpretation
    algorithms [@pineda2013; @snoeck2006; @theys2015]


References
---------------------------

[1] A. B. Abecasis, A. M. Wensing, D. Paraskevis, J. Vercauteren, K. Theys, D. A. Van de Vijver, J. Albert, B. Asjo, C. Balotta, D. Beshkov, R. J. Camacho, B. Clotet, C. De Gascun, A. Griskevicius, Z. Gross- man, O. Hamouda, A. Horban, T. Kolupajeva, K. Korn, L. G. Kostrikis, C. Kucherer, K. Liitsola, M. Linka, C. Nielsen, D. Otelea, R. Pare- des, M. Poljak, E. Puchhammer-Stockl, J. C. Schmit, A. Sonnerborg, D. Stanekova, M. Stanojevic, D. Struck, C. A. Boucher, and A. M. Van- damme. HIV-1 subtype distribution and its demographic determinants in newly diagnosed patients in Europe suggest highly compartmentalized epi- demics. Retrovirology, 10:7, Jan 2013.

[2] B. Chaplin, G. Imade, C. Onwuamah, G. Odaibo, R. Audu, J. Okpokwu, D. Olaleye, S. Meloni, H. Rawizza, M. Muazu, A. Z. Musa, J. Samuel, O. Agbaji, O. Ezechi, E. Idigbe, and P. J. Kanki. Distinct Pattern of Thymidine Analogue Mutations with K65R in Patients Failing Tenofovir- Based Antiretroviral Therapy. AIDS Res. Hum. Retroviruses, 34(2):228– 233, Feb 2018.

[3] L. Cuypers, G. Li, P. Libin, S. Piampongsant, A. M. Vandamme, and K. Theys. Genetic Diversity and Selective Pressure in Hepatitis C Virus Genotypes 1-6: Significance for Direct-Acting Antiviral Treatment and Drug Resistance. Viruses, 7(9):5018–5039, Sep 2015.

[4] L. Cuypers, G. Li, C. Neumann-Haefelin, S. Piampongsant, P. Libin, K. Van Laethem, A. M. Vandamme, and K. Theys. Mapping the genomic diversity of HCV subtypes 1a and 1b: Implications of structural and im- munological constraints for vaccine and drug development. Virus Evol, 2(2):vew024, Jul 2016.

[5] K. Deforche, R. J. Camacho, Z. Grossman, M. A. Soares, K. Van Laethem, D. A. Katzenstein, P. R. Harrigan, R. Kantor, R. Shafer, A. M. Van- damme, R. Kantor, D. A. Katzenstein, R. W. Shafer, R. J. Camacho, A. P. Carvalho, B. Wynhoven, P. R. Harrigan, P. Cane, J. Clarke, J. Weber, S. Sirivichayakul, P. Phanuphak, M. A. Soares, A. Tanuri, J. Snoeck, A. M. Vandamme, L. Morris, H. Rudich, Z. Grossman, J. M. Schapiro, R. Ro- drigues, L. F. Brigido, A. Holguin, V. Soriano, K. Ariyoshi, W. Sugiura, M. B. Bouzas, P. Cahn, D. Pillay, T. L. Katzenstein, and L. B. J?rgensen. Bayesian network analyses of resistance pathways against efavirenz and nevirapine. AIDS, 22(16):2107–2115, Oct 2008.

[6] K. Deforche, A. Cozzi-Lepri, K. Theys, B. Clotet, R. J. Camacho, J. Kjaer, K. Van Laethem, A. Phillips, Y. Moreau, J. D. Lundgren, and A. M. Van- damme. Modelled in vivo HIV fitness under drug selective pressure and estimated genetic barrier towards resistance are predictive for virological response. Antivir. Ther. (Lond.), 13(3):399–407, 2008.


[7] R. C. Edgar. MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res., 32(5):1792–1797, 2004.

[8] E. L. Hatcher, S. A. Zhdanov, Y. Bao, O. Blinkova, E. P. Nawrocki, Y. Ostapchuck, A. A. Schaffer, and J. R. Brister. Virus Variation Re- source - improved response to emergent viral outbreaks. Nucleic Acids Res., 45(D1):D482–D490, Jan 2017.

[9] K. Katoh and D. M. Standley. MAFFT: iterative refinement and additional methods. Methods Mol. Biol., 1079:131–146, 2014.

[10] G. Li, S. Piampongsant, N. R. Faria, A. Voet, A. C. Pineda-Pena, R. Khouri, P. Lemey, A. M. Vandamme, and K. Theys. An integrated map of HIV genome-wide variation from a population perspective. Retro- virology, 12:18, Feb 2015.

[11] P. Libin, G. Beheydt, K. Deforche, S. Imbrechts, F. Ferreira, K. Van Laethem, K. Theys, A. P. Carvalho, J. Cavaco-Silva, G. Lapadula, C. Torti, M. Assel, S. Wesner, J. Snoeck, J. Ruelle, A. De Bel, P. La- cor, P. De Munter, E. Van Wijngaerden, M. Zazzi, R. Kaiser, A. Ayouba, M. Peeters, T. de Oliveira, L. C. Alcantara, Z. Grossman, P. Sloot, D. Ote- lea, S. Paraschiv, C. Boucher, R. J. Camacho, and A. M. Vandamme. RegaDB: community-driven data management and analysis for infectious diseases. Bioinformatics, 29(11):1477–1480, Jun 2013.

[12] P. Libin, E. Vanden Eynden, F. Incardona, A. Nowe, A. Bezenchek, A. Son- nerborg, A. M. Vandamme, K. Theys, and G. Baele. PhyloGeoTool: inter- actively exploring large phylogenies in an epidemiological context. Bioin- formatics, 33(24):3993–3995, Dec 2017.

[13] S. Ngcapu, K. Theys, P. Libin, V. C. Marconi, H. Sunpath, T. Ndung’u, and M. L. Gordon. Characterization of Nucleoside Reverse Transcriptase Inhibitor-Associated Mutations in the RNase H Region of HIV-1 Subtype C Infected Individuals. Viruses, 9(11), Nov 2017.

[14] A. C. Pineda-Pena, N. R. Faria, S. Imbrechts, P. Libin, A. B. Abecasis, K. Deforche, A. Gomez-Lopez, R. J. Camacho, T. de Oliveira, and A. M. Vandamme. Automated subtyping of HIV-1 genetic sequences for clini- cal and surveillance purposes: performance evaluation of the new REGA version 3 and seven other tools. Infect. Genet. Evol., 19:337–348, Oct 2013.

[15] S. Y. Rhee, M. J. Gonzales, R. Kantor, B. J. Betts, J. Ravela, and R. W. Shafer. Human immunodeficiency virus reverse transcriptase and protease sequence database. Nucleic Acids Res., 31(1):298–303, 2003.

[16] F. Sievers, A. Wilm, D. Dineen, T. J. Gibson, K. Karplus, W. Li, R. Lopez, H. McWilliam, M. Remmert, J. Soding, J. D. Thompson, and D. G. Hig- gins. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Mol. Syst. Biol., 7:539, Oct 2011.


[17] J. Snoeck, R. Kantor, R. W. Shafer, K. Van Laethem, K. Deforche, A. P. Carvalho, B. Wynhoven, M. A. Soares, P. Cane, J. Clarke, C. Pillay, S. Sirivichayakul, K. Ariyoshi, A. Holguin, H. Rudich, R. Rodrigues, M. B. Bouzas, F. Brun-Vezinet, C. Reid, P. Cahn, L. F. Brigido, Z. Grossman, V. Soriano, W. Sugiura, P. Phanuphak, L. Morris, J. Weber, D. Pillay, A. Tanuri, R. P. Harrigan, R. Camacho, J. M. Schapiro, D. Katzenstein, and A. M. Vandamme. Discordances between interpretation algorithms for genotypic resistance to protease and reverse transcriptase inhibitors of hu- man immunodeficiency virus are subtype dependent. Antimicrob. Agents Chemother., 50(2):694–701, Feb 2006.

[18] Philip R Tedbury and Eric O Freed. Hiv-1 gag: an emerging target for antiretroviral therapy. In The Future of HIV-1 Therapeutics, pages 171– 201. Springer, 2015.

[19] K. Theys, A. Abecasis, P. Libin, P. T. Gomes, J. Cabanas, R. J. Camacho, and K. Van Laethem. Discordant predictions of residual activity could impact dolutegravir prescription upon raltegravir failure. J. Clin. Virol., 70:120–127, Sep 2015.

[20] K. Theys, P. Libin, K. Dallmeier, A. C. Pineda-Pena, A. M. Vandamme, L. Cuypers, and A. B. Abecasis. Zika genomics urgently need standardized and curated reference sequences. PLoS Pathog., 13(9):e1006528, 09 2017.

[21] K. Theys, J. Vercauteren, J. Snoeck, M. Zazzi, R. J. Camacho, C. Torti, E. Schulter, B. Clotet, A. Sonnerborg, A. De Luca, Z. Grossman, D. Struck, A. M. Vandamme, and A. B. Abecasis. HIV-1 subtype is an independent predictor of reverse transcriptase mutation K65R in HIV-1 patients treated with combination antiretroviral therapy including tenofovir. Antimicrob. Agents Chemother., 57(2):1053–1056, Feb 2013.

[22] J. Vercauteren, I. Derdelinckx, A. Sasse, M. Bogaert, H. Ceunen, A. De Roo, S. De Wit, K. Deforche, F. Echahidi, K. Fransen, J. C. Gof- fard, P. Goubau, E. Goudeseune, J. C. Yombi, P. Lacor, C. Liesnard, M. Moutschen, D. Pierard, R. Rens, Y. Schrooten, D. Vaira, A. van den Heuvel, B. van der Gucht, M. van Ranst, E. van Wijngaerden, B. Vander- cam, M. Vekemans, C. Verhofstede, N. Clumeck, A. M. Vandamme, and K. van Laethem. Prevalence and epidemiology of HIV type 1 drug resis- tance among newly diagnosed therapy-naive patients in Belgium from 2003 to 2006. AIDS Res. Hum. Retroviruses, 24(3):355–362, Mar 2008.

[23] J. Vercauteren, K. Theys, A. P. Carvalho, E. Valadas, L. M. Duque, E. Te- ofilo, T. Faria, D. Faria, J. Vera, M. J. Aguas, S. Peres, K. Mansinho, A. M. Vandamme, R. J. Camacho, K. Mansinho, A. Claudia Miranda, I. Aldir, F. Ventura, J. Nina, F. Borges, E. Valadas, M. Doroana, F. Antunes, M. Joao Aleixo, M. Joao Aguas, J. Botas, T. Branco, J. Vera, I. Vaz Pinto, J. Pocas, J. Sa, L. Duque, A. Diniz, A. Mineiro, F. Gomes, C. San- tos, D. Faria, P. Fonseca, P. Proenca, L. Tavares, C. Guerreiro, J. Nar- ciso, T. Faria, E. Teofilo, S. Pinheiro, I. Germano, U. Caixas, N. Faria, 
A. Paula Reis, M. Bentes Jesus, G. Amaro, F. Roxo, R. Abreu, and I. Neves. The demise of multidrug-resistant HIV-1: the national time trend in Por- tugal. J. Antimicrob. Chemother., 68(4):911–914, Apr 2013.
