# Snakemake workflow: `phage-assembly-and-annotation`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/ICTV-VBEG/phage-assembly-and-annotation/workflows/Tests/badge.svg?branch=main)](https://github.com/ICTV-VBEG/phage-assembly-and-annotation/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for `A snakemake workflow for phage assembly and annotation, following Shen and Millard et al. (2021): https://doi.org/10.1089/phage.2021.0015`


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=ICTV-VBEG%2Fphage-assembly-and-annotation).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) phage-assembly-and-annotationsitory and its DOI (see above).

# TODO

This is taken from the supplementary material of: 
Anastasiya Shen and Andrew Millard.Phage Genome Annotation: Where to Begin and End.PHAGE.Dec 2021.183-193.http://doi.org/10.1089/phage.2021.0015 

We will be adapting it to a TODO-list to use in creating a snakemake workflow following this walk-through.
Our general aims are to:

1. Follow the best practices regarding snakemake workflow folder structures: https://snakemake.readthedocs.io/en/latest/snakefiles/deployment.html#distribution-and-reproducibility
2. Use conda environments for automated software installation: https://snakemake.readthedocs.io/en/latest/snakefiles/deployment.html#integrated-package-management
3. Use snakemake wrappers for any steps where they are available: https://snakemake-wrappers.readthedocs.io/
4. Enable automatic inclusion of the workflow in the snakemake workflow catalog (this should work out-of-the-box, standardized usage can be enabled later): https://snakemake.github.io/snakemake-workflow-catalog/?rules=true
5. Use the example data provided for automated testing (continuous integration) of workflow functionality in a `.test/` folder and using GitHub Actions.


## Phage annotation walk through overview

The following provides a guided example of how a phage genome can be annotated and assembled. It is by no means the only way to assemble a phage genome. We utilise open source tools and the size of phage genomes, allowing for their assembly on an average laptop. There is an investment in time required to install some of the software and databases initially. However, the time taken to produce an assembled and annotated genome once installed can be completed in less than two hrs with this guide. The timings included at each step are approximate and based on the time to install the software and complete the computation task. 

For this example we use the reads from the previously published genome of phage Erla https://journals.asm.org/doi/full/10.1128/MRA.01354-20. All intermediary files that were produced as part of this process are available from https://leicester.figshare.com/articles/dataset/Bacteriophage_Genome_Annotation/16896277  


## Software installation

(Almost) all of the software installation will be handled using conda.
You can follow the snakemake installation instructions to get a working version of conda (or its faster drop-in replacement tool mamba): https://snakemake.readthedocs.io/en/latest/getting_started/installation.html#installation-via-conda-mamba

### snakemake wrappers

The easiest example is when a snakemake wrapper exists, for example the [snakemake wrapper for `fastqc`]().
Here, we can simply copy-paste the example code from the snakemake wrapper documentation into our workflow and adapt `input`, `output` and maybe the `params` specifications.
Software installation is automatically handled by the wrapper.

### (bio)conda installation

For a number of tools, no snakemake wrapper exists, but the software is available via the [`bioconda`](https://bioconda.github.io/) or [`conda-forge`](https://conda-forge.org/) conda channel.
Then, we can simply create a [`YAML`](https://koesterlab.github.io/data-science-for-bioinfo/data_formats/yaml.html) environment definition file in the folder `worfklow/envs/`.
For example, for `bandage` this would be a file called `workflow/envs/bandage.yaml` with the following contents:
```
channels:
  - conda-forge
  - bioconda
dependencies:
  - bandage =0.8
```

You can then reference the environment (and thus use bandage) in any snakemake rule with:
```
rule some_bandage_rule:
    input: "..."
    output: "..."
    conda: "../envs/bandage.yaml"
    shell: "bandage ..."
```

### manual installation

The only software that we will have to manually install to start with (but should try to [make available on bioconda](https://bioconda.github.io/contributor/index.html) in the mid-term), seems to `PhageTerm`: 

PhageTerm release can be downloaded here:https://gitlab.pasteur.fr/vlegrand/ptv/-/releases; or it can be installed with the following code: 

```
wget https://gitlab.pasteur.fr/vlegrand/ptv/-/archive/py3_release_1_light/ptv-py3_release_1_light.tar.gz  
tar -xvzf ptv-py3_release_1_light.tar.gz  
rm ptv-py3_release_1_light.tar.gz 
ptv-py3_release_1_light/PhageTerm.py -h 
```
 
For paired end data it might be necessary to provide interleaved fastq files.

### complete overview of software used

This list is just for a quick overview and to keep all the descriptions and links of the original supplement around:
* Python (a programming language; https://www.python.org/) 
* Biopython (a set of python tools; https://biopython.org/) 
* FastQC (a quality control tool; https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
* BBtools (bioinformatics tools that allow analysis of DNA and RNA sequence data; https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/) 
* Seqtk (a toolkit for reads subsampling; https://github.com/lh3/seqtk) 
* SPAdes (genome assembler; https://github.com/ablab/spades) 
* Bandage (allows assembly graphs to be visualised; https://github.com/rrwick/Bandage) 
* Samtools (​​programs for interacting with high-throughput sequencing data; http://www.htslib.org/) 
* Pilon (automatic assembly autocorrection; https://github.com/broadinstitute/pilon/wiki)  
* PhageTerm (a tool for phage termini prediction; https://sourceforge.net/projects/phageterm/) 
* Prokka (a tool for whole genome annotation; https://github.com/tseemann/prokka) 
 


## Data download

### reference data

#### prokka: external databases

This protocol we will use several prebuilt external databases for annotation improvement and comparison. Those are Caudovirales database, Virus orthologous groups database (pVOG; https://vogdb.org/) and Prokaryotic Virus Remote Homologous Groups database (PHROGs; https://phrogs.lmge.uca.fr/). 

* [ ] TODO: create a rule to download Caudovirales database (genus database of the Caudovirales group of viruses) with the following command:
```
wget http://s3.climb.ac.uk/ADM_share/crap/Caudovirales.tar.gz 
```

* [ ] TODO: create a rule to download PHROGs HMM database with the following commands: 
```
wget http://s3.climb.ac.uk/ADM_share/all_phrogs.hmm.gz -o all_phrogs.hmm.gz 
gunzip all_phrogs.hmm.gz 
```

*  [ ] TODO: generalize the two above rules for general downloads from the `https://s3.climb.ac.uk/` resource

### example data from SRA

* [ ] TODO: implement a rule for data download using an SRA accession, using the `sra-tools` snakemake wrapper for `fasterq-dump`: https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/sra-tools/fasterq-dump.html
* [ ] TODO: allow for the specification of SRA accessions in a `samples.tsv` file that gets parsed automatically in the Snakefile, using `pandas.read_tsv()` (see e.g.: https://github.com/snakemake-workflows/dna-seq-varlociraptor/blob/520e3611bb26f5f0a31b512470df1ba565494df3/workflow/rules/common.smk#L13 and https://github.com/snakemake-workflows/dna-seq-varlociraptor/blob/520e3611bb26f5f0a31b512470df1ba565494df3/.test/config-sra/units.tsv?plain=1#L2)
* [ ] TODO: create a `.test/` folder for continuous integration tests via [GitHub Actions](https://docs.github.com/en/actions), create a subfolder `config/` and put a testing `samples.tsv` with the below accession numbers specified, there

Information about the reads used as example data for this protocol can be found at European Nucleotide Archive: https://www.ebi.ac.uk/ena/browser/view/SRR13108336?show=reads

Download reads:
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/036/SRR13108336/SRR13108336.fastq.gz
 
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/001/SRR11901901/SRR11901901_1.fastq.gz


## workflow steps

### genome assembly using Illumina reads 

1. [ ] TODO: rule for quality control of reads with FastQC (Time < 2 min) 

Snakemake wrapper to use: https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastqc.html

Example command:
```
fastqc -o fastqc SRR13108336.fastq.gz
```

FastQC output is used for the adapters/overrepresented sequences removal and the poor quality bases trimming. 

2. [ ] TODO: rule for read quality trimming, phiX removal and adapters removal with BBduk (Time < 2 min)

Snakemake wrapper to use: https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bbtools/bbduk.html

All of the commands below, and especially their adapters and paramters need to be adjusted to fit the snakemake wrapper syntax.

* [ ] TODO: create rule to download `phix174_ill.ref.fa.gz` from bbmap repository into a `resources/` folder: https://github.com/BioInfoTools/BBMap/blob/v36.20/resources/phix174_ill.ref.fa.gz

Example command:
```
bbduk.sh in=SRR13108336.fastq.gz out=results/nophix/SRR13108336.fq.gz ref=resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=phiX_stats.txt
```

A reference file containing the phiX genome is used and all reads that have a 31 kmer match to phiX are removed, leaving a sequence file free of phiX reads. The phiX_stats.txt provides statistics on how many reads were removed.

The snakemake wrapper should be able to handle both single-end and paired-end data, so also commands like:
```
bbduk.sh in1=read1.gz in2=read2.gz out1=results/nophix/read1.gz out2=results/nophix/read2.gz
```

Post-phiX read removal, sequencing adapters and quality trimming of reads is also required. This can can again be done with bbduk.sh: 

```
bbduk.sh in=nophix_SRR13108336.fq.gz out=results/clean/SRR13108336.fq.gz ref=resources/adapters.fa ktrim=r hdist=1 tpe tbo minlen=100 qtrim=rl trimq=28
```

The above parameters require a minimum sequence length of 100 post trimming. Parameters may need to be adjusted depending on the quality of sequences. Generally with an excess of sequence data that occurs with sequencing phage genomes it is possible to be stringent on the parameters used. The file adapters.fa is part of the bbduk package and contains common sequencing adapters. Non-standard adapters can be added to this file or specified with the “--literal” flag (consult bbduk instructions). 

After read trimming and adapters removal the quality of reads can be reassessed with FastQC.

* [ ] TODO: implement another instance of the fastqc snakemake wrapper, or generalize the rule from above to also allow running on this output (asking for pre- and post-trimming output via the `rule all:`). 

```
fastqc -o fastqc/clean/ clean_SRR13108336.fq.gz 
```

3. [ ] TODO: create rule for calculating coverage for an average 100 kb genome based on the amount of data produced by sequencing

The information about number of reads and read length can be found in FastQC output. It should be read from it programmatically.

Estimated coverage per genome can be calculated using the formula below:
```
Coverage = (Number reads x Read Length) / Expected Genome size  
```

The `Number of reads = 643307` can be found in the FastQC with:

```
zcat clean_SRR13108336.fq.gz | grep "@SRR" | wc -l
```

Similarly for the `Average read length = 150 bp` (from the FastQC result most of the sequences are 150 bp long).

**===============================================================================================**
**FROM THIS POINT ONWARDS, THE INSTRUCTIONS ARE AS IN THE ORIGINAL SUPPLEMENT. PROCEED AS BEFORE.**

4| Reads subsampling with Seqtk (Time < 2 min) 

To calculate the number of reads for the subsampling to give a ~100x coverage of the genome can be used the formula below: 

Number of reads = (Expected coverage x Expected Genome size) / Read Length  

To get 100x coverage for the 100 kb phage genome we need to have ~66667 150 bp long reads. Subsample 66667 reads from FASTQ file: 

```
$ <path_to_installation_directory>/seqtk/seqtk sample -s 100 clean_SRR13108336.fq.gz 66667 > sub_SRR13108336.fq 
```

If paired end data was being used each set of reads would have to be sub-sampled separately. It is essential that the seed value (-s 100) is the same for each read set:  

seqtk sample -s 100 read1.fq 1000 > subR1.fq 
seqtk sample -s 100 read2.fq 1000 > subR2.fq 

The following step just compresses the reads file to make it smaller: 

$ gzip sub_SRR13108336.fq 

5| Genome assembly with SPAdes (Time < 15 min) 

$ spades.py --only-assembler -s sub_SRR13108336.fq.gz -o spades_result 

The output result is in a directory named spades_result and the FASTA file with the assembled contigs is named contigs.fasta. If using paired end data the -1 and -2 options would be used to specify read. The timing will depend on the number of threads available and the power of the machine. 

On rare occasions, it may be necessary to adjust the kmer values to improve the assembly, extending the largest kmer value above the default set by SPAdes.  

 

6| Preliminary checking of assembled contigs with Bandage (Time < 10 mins) 

The spades_results directory contains a file called assembly_graph.fastg which is an assembly graph. It can be visualized with Bandage: 

$ <path_to_installation_directory>/Bandage load spades_result/assembly_graph.fastg --draw  

7| Assembly validation by extracting the coverage per contig with BBMap (<10 min)  

$ gunzip sub_SRR13108336.fq.gz 
$ <path_to_installation_directory>/bbmap/bbmap.sh ref=spades_result/contigs.fasta in=sub_SRR13108336.fq covstats=contig_covstats.txt out=contig_mapped.sam 
$ gzip sub_SRR13108336.fq 

If using paired end data then both read files would be specified with the in and in2 flags e.g. in1=read1_file.fq and in2=read2_file.fq   

The BBMap command above creates two output files in the bbmap directory. The file called contig_constats.txt contains information about the coverage of the assembled contigs ("Avg_fold") and contig_mapped.sam contains sequencing reads aligned to the input genome. The SAM file can also be used for assembly correction larter.  

The resulting SPAdes assembly contains a single contig with high coverage. It is not usually the case and a short low-coverage contigs may be present, that most likely come from contaminating host DNA or spurious reads. It is recommended to keep only the long contig as this most likely represents an assembly of a complete phage genome. It is expected that DNA from the same source will have approximately similar coverage. Thus, if a phage genome is broken into two contigs these are likely to have similar coverage. Whereas contaminating bacterial DNA will be of lower abundance (or on some occasions higher!).  

Host contaminating DNA can be removed in python by using Biopython tools.  

For extracting the assembly of a phage genome will be used the header of the contig of interest.  

Display headers of a FASTA file: 

$ cat spades_result/contigs.fasta | grep '>' 

In this example the header NODE_1_length_41615_cov_83.163922 is of interest (note that ‘>’ is not included). 

Activate python: 

$ python  

To extract the contig of interest we will create a python script as detailed below. In this example our output only contains one contig, so the step is unnecessary. However, the code will demonstrate how to do this if more than one contig is present which often occurs. The name of your contig that you obtain from assembly will probably be different. You will need to adjust the line “if name == 'NODE_1_length_41615_cov_83.163922':” in the script below to match your header in the previous step. 

First type the command $ nano extract_contig.py, which will open a text editor. In the text editor paste in the python code: 

  

from Bio.SeqIO.FastaIO import SimpleFastaParser  
import gzip  
import pandas as pd  
file = 'spades_result/contigs.fasta'  
handle = open(file)  
textfile = open("phage_contig.fasta", "w")  
entries = list(SimpleFastaParser(handle))  
for name, seq in entries: 	 
	if 'NODE_1_length_41615_cov_83.163922' in name: 		 
		print(f">{name}\n{seq}\n", file=textfile) 

After press Cntrl+O to save the file (here you can also change a name of the file), then Enter and Cntrl+X to exit.  

Now we can execute the python script by typing: 

$ python extract_contig.py 

The contig representing an assembly of a phage genome is now saved to a file named phage_contig.fasta. 

8| Error detection and correction with Pilon (Time <10 min) 

Cleaned from contamination assembly can be automatically checked for errors using Pilon software. Pilon requires a BAM file containing sequencing reads aligned to the input genome, which has to be sorted and indexed. This can be done with Samtools (http://www.htslib.org/doc/samtools-sort.html). 

Sort and index: 

 

$ samtools view -bS -F4 contig_mapped.sam | samtools sort - -o contig_mapped_sorted.bam 
 
$ samtools index contig_mapped_sorted.bam  

 

The samtools view -bS -F4 bbmap/contig_mapped.bam converts the SAM file into a compressed BAM file, only keeping reads that are mapped to the reference (“-F4” flag determines this). The next part | samtools sort - -o contig_mapped_sorted.bam uses the “|” also known as a pipe, to take the output of the first command and pass into samtools sort command. With -o specifying the output will go to a file called contig_mapped_sorted.bam. 

 

Assembly autocorrection: 

 

$ pilon --genome phage_contig.fasta --frags contig_mapped_sorted.bam --output pilon_SRR13108336 --verbose --changes  

Pilon will detail any assembly errors detected and automatically attempt to fix them, details of changes will be recorded in output files (e.g. corrected_SRR13108336.changes). It may be necessary to complete multiple rounds of pilon polishing, until no more errors are detected. This requires mapping of reads to the corrected assembly, creation of a BAM file and running pilon again in an iterative process.  

In this case the important in the run output information is:  

Confirmed 41615 of 41615 bases (100.00%) 
Corrected 0 snps; 0 ambiguous bases; corrected 0 small insertions totaling 0 bases, 0 small deletions totaling 0 bases  

The *changes file is empty as there were no changes to make or fix.  

Corrected assembly can be re-assessed for changes with Bandage and BBMap. 

9| Preliminary identification of closest relatives with Blastn (<15 min) 

The closest relative of the assembled putative phage genome can be found using BLASTN web interface (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch), where the corrected_SRR13108336.fasta can be uploaded. In the field “Database” of “Choose Search Set” should be selected “Standard databases (nr etc.): Nucleotide collection (nr/nt)”. In the field “Organism” it is recommended to use “viruses (taxid:10239)”. In the “Program Selection” choose “Highly similar sequences (megablast)”. 

The genome of the closest relative can be directly downloaded from the web page of the BLASTN output. Select a genome with the highest percentage identity (Per. Ident) and click “GeneBank”. The link should redirect to the page https://www.ncbi.nlm.nih.gov/nuccore/MW291026.1. Here in “Send to” the format can be selected for the download. Further, we will use the genome in GenBank format for the annotation with the use of a closely related phage(Send to > Complete Record > File > Format:GenBank), a txt file with the coding DNA sequences (Send to > Coding Sequences > Format:FASTA Nucleotide) for the further comparison of the annotation results and the FASTA file with complete genome sequence (Send to > Complete Record > File > Format:FASTA) for manual genome reorder.  

Downloaded files can be moved to the working directory and renamed: 

$ mv sequence.gb erla_phage.gbk 
$ mv sequence.txt erla_phage.ffn 
$ mv sequence.fasta erla_phage.fasta 

10| Preliminary identification of closest relatives with PhageClouds 

The web interface PhageClouds (http://130.226.24.116/phagecompass/index.html#/Search%20phage%20 clouds) can be used to compare the assembled contig against a database of viruses with a visual representation of genome size that helps to identify the closest relatives and to judge the assembly completeness based on similarity and size of other phages. 

The file with the assembled contig (corrected_SRR13108336.fasta) can be directly uploaded to the PhageClouds and compared against all available databases. In the result you will get a list of phages with the top match. Our assembled contig forms a cluster with Microbacterium phages and it has the shortest distance with Microbacterium phage Erla - as expected given the reads where from this phage.  

Genome reordering  

11| Genome termini identification with PhageTerm (Time < 15 min) 

 

It is good practice to store reads in a compressed format (gz), unfortunately PhageTerm.py does not currently support compressed files, so the files needs to be expanded prior to use (gunzip).  

 

$ gunzip SRR13108336.fastq.gz 
 
$ ptv-py3_release_1_light/PhageTerm.py -f SRR13108336.fastq -r pilon_SRR13108336.fasta --report_title putative_phage  
 
$ gzip SRR13108336.fastq 

 

In the output there is a detailed PDF report (putative_phage_PhageTerm_report.pdf), a detailed statistics table (putative_phage_statistics.csv) and the phage genome sequence reorganized according to termini positions (putative_phage_sequence.fasta). 

 

In this instance PhageTerm did not identify termini and suggests the genome is circularly permuted (this does not mean the genome is circular).  

 

To confirm if it is circularly permuted the script apc.pl can be used(https://github.com/jfass/apc/blob/master/apc.pl) that looks for circular contigs. It identified the contig is circular with 77 bp that is identical on each end of the contig. This repeat can also be identified with Bandage using the BLASTn function. It is an artifact of assembly and is created by the assembler, it will not always be 77 bp. The 77 bp needs to be removed from one end, this can be done automatically using apc.pl or manually. In this instance we remove the 77 bp manually.  

 

The instructions provided above are for the PhageTerm version linked above. Previous versions have slightly different syntax and it should be noted that for paired end data, interleaved data might be required.  

 

Remove the assembly artifact with python: 

 

$ nano remove77.py 

 

from Bio import SeqIO 
record = SeqIO.read("pilon_SRR13108336.fasta", "fasta") 
with open("corrected_SRR13108336.fasta", "w") as out: 
    SeqIO.write(record[77:], out, "fasta") 

 

$ python remove77.py 

 

As the contig was previously found to be similar to other phages and does not have defined termini, it needs to be reordered in respect to a close relative.  

 

12| Manual genome reordering (<15 min) 

 

If the start and end of the genome are not identifiable by PhageTerm (which is the case for the assembled genome in this protocol), manual genome reordering can be performed by using information about the closest relative (representative genome). For the manual reordering the genome of the closest relative was downloaded from the web page of the BLASTN output (step 9).  

 

Reordering a genome is not always a trivial task and requires the manipulation of large sequence files. In the example below we will demonstrate one way this can be done with biopython to manipulate the sequence files. However, similar results can be achieved by the use of a text editor (never edit sequence files in word processors such as Word) and the copy and paste functions. 

 

As the arbitrary start of the genome will be taken the first CDS sequence of the closest relative, which we previously identified as phage Erla (https://www.ncbi.nlm.nih.gov/nuccore/1939416749). Looking at the Genbank file of Erla the first gene starts at base 1:  

 

 

 

 

 

Therefore we know the first 525 bases are from the first gene.  

 

Using python this can be extracted from the file erla_phage.fna that we have previously downloaded. Alternatively this can be gained through interaction with the NCBI website via the link https://www.ncbi.nlm.nih.gov/nuccore/1939416749.  

 

Create a python script by opening a text editor: 

$ nano extract_startseq.py 

 

from Bio import SeqIO 
record = SeqIO.read("erla_phage.fasta", "fasta") 
with open("startseq.fasta", "w") as out: 
    SeqIO.write(record[0:525], out, "fasta") 

Execute the python script by typing: 

$ python extract_startseq.py 

 

The sequence contained in startseq.fasta can then be used to find the start position in our newly assembled phage. An easy solution to this is to use the BLASTN web interface to align two sequences (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=MegaBlast&PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&BLAST_SPEC=blast2seq&DATABASE=n/a&QUERY=&SUBJECTS=). 

 

Here we can align the first gene from the reference genome to our newly assembled phage genome. In this example the “Query Sequence” was the newly assembled genome (corrected_SRR13108336.fasta) and the “Subject Sequence” is the first gene in the genome of our closest relative (in the file startseq.fasta). In the “Program Selection'' select “Highly similar sequences''.  

 

From the alignment section of the BLASTN search result, you can view the aligned sequence as well as coordinates of the query and subject sequences. 

 

 

In our case it can be seen the start point can be identified as 7930, but on the reverse strand. To make it easier to determine the position of the forward strand, the simplest action is to reverse complement the file corrected_SRR13108336.fasta. 

 

$ nano reverse.py 

 

from Bio import SeqIO 
for record in SeqIO.parse(open('corrected_SRR13108336.fasta'), 'fasta'): 
    reverse_seq = str(record.seq.reverse_complement()) 
with open("reversed_SRR13108336.fasta", "w") as out: 
    out.write('>' + record.id + '\n' + reverse_seq) 
 

 

$ python reverse.py 

 

Repeating the above blast analysis with the reverse complement sequence results in: 

 

 

 

 

Once we have the coordinates of the gene that we will use as the start of our genome the last step is to do the actual genome reordering. Be aware you might need to take into account reverse complementation with some genomes.  

 

$ nano reorder.py 

 

from Bio import SeqIO 
 
for contig_record in SeqIO.parse(open('corrected_SRR13108336.fasta'), 'fasta'): 
    contig = str(contig_record.seq) 
str1= contig[33608:] 
str2= contig[:33608] 
str_final='>'+contig_record.id + '\n' + str1 + str2 
with open("reordered_SRR13108336.fasta", "w") as out: 
    out.write(str_final) 

 

$ python reorder.py 

 

Here we need to repeat a pilon polishing step after reordering to ensure no errors have been introduced in the process. At this point we will also use all of the reads that were included, not just the subsampled reads.  

 

 
$<path_to_installation_directory>/bbmap/bbmap.sh ref=reordered_SRR13108336.fasta in=clean_SRR13108336.fq.gz covstats=final_covstats.txt out=all_reads_contig_mapped.sam  
 

View the file final_covstats.txt which will provide detail on the depth of sequencing, which is needed for submission of the genome. Record this information (2281 in this example).   

$ head final_covstats.txt  

 

We repeat the sorting and indexing of the bam file with the following two commands:  

 

$ samtools view -bS -F4 all_reads_contig_mapped.sam | samtools sort - -o all_reads_contig_sorted_mapped.bam 
$ samtools index all_reads_contig_sorted_mapped.bam 
 

 

At this point we can also collect all reads that did not map to the assembly:  

 

$ samtools view -bS -f4 all_reads_contig_mapped.sam | samtools sort - -o non_mapped_reads.bam 

 

Using the flag “-f4” rather than -F4 results in unmapped reads. These reads can then be converted back to fastq format and used as input for assembly:  

 

$ samtools fastq non_mapped_reads.bam > non_mapped_reads.fq  

  

In this example this does not result in the identification of any other phage contigs. However, it is a useful check to do when processing a new unknown phage genome sample. 

 

$ pilon --genome reordered_SRR13108336.fasta --frags all_reads_contig_sorted_mapped.bam --output final_SRR13108336 --verbose --changes  

 

In the run output information we see:  

Confirmed 41538 of 41538 bases (100.00%) 
Corrected 0 snps; 0 ambiguous bases; corrected 0 small insertions totaling 0 bases, 0 small deletions totaling 0 bases  

Note that the number of bases is 77bp shorter compared to what we had before genome reordering. The *changes file is empty as there were no changes to make or fix. 

Annotation - gene calling (Time ~1 hr for steps 13-16) 

 

Here we use and install a variety of databases to demonstrate the effect they have on annotation that might be obtained. The time taken for actual annotation is < 2 mins once databases are installed.   

 

Prior to any gene calling we will make the fasta header of the final phage contig compliant with requirements of Prokka.  

 

 

13| Initial annotation with Prokka 

 

Prokka requires contig id to be <= 37 chars long. Since our contig id in the final_SRR13108336.fasta might be longer, we can change it by: 

$ sed "1s/.*/>contig/" final_SRR13108336.fasta > SRR13108336.fasta 

 

$ prokka --outdir prokka --prefix mygenome --kingdom viruses SRR13108336.fasta  

 

If you run the command above without the “--kingdom” option, the program by default would use the kingdom Bacteria for annotation.   

 

The output files are in the directory specified by the option “--outdir”.  

View the list of output files: 

 

$ ls prokka 

 

View a specific file in the directory: 

 

$ less prokka/mygenome.faa 

 

The output directory “prokka” contains several files including FASTA protein file with predicted protein sequences (mygenome.faa), FASTA nucleotide file with contig sequences (mygenome.fna), GenBank file (mygenome.gbk) and Gene Feature Format file (mygenome.gff). 

 

As a result of gene calling 63 genes were predicted and no tRNAs, which is consistent with the submitted genbank file and original publication. Comparison of gene coordinates between the submitted Genbank file and gene calling here, identified six genes that differed in length (Table S1). Comparison of the annotation between the submitted Genbank file and annotation with “--kingdom Viruses” highlights the poor annotation with this database.  

 

For the annotation improvement and comparison we will use external databases for annotation in the next steps. To be able to use it in prokka, the database has to be added to the /path/to/prokka/db/directory. In this protocol we will use one genus database (has to be added to /path/to/prokka/db/genus) and HMM databases (has to be added to /path/to/prokka/db/hmm). The path to the Prokka db directory can be found by typing $ locate /pkgs/prokka-1.13-2/db/ which will output a list of prokka directories. If it was installed with conda the path it most likely /opt/anaconda3/pkgs/prokka-1.13-2/db. It might be necessary to index installed databases by typing: 

 

 $ prokka --setupdb  

 

14| Annotation with the use of Caudovirales database.  

 

$ tar -xvzf Caudovirales.tar.gz 
$ rm Caudovirales.tar.gz 

 

Move Caudovirales database to /prokka/db/genus: 

$ mv Caudovirales* /path/to/prokka/db/genus 

 

$ prokka --outdir prokka_caud --prefix mygenome SRR13108336.fasta --usegenus --genus Caudovirales --kingdom viruses  

 

Using the Caudovirales database improves the number of genes with any functional annotation with the built in Viruses database from 1 to 27 genes with the Caudovirales database. A number of these annotations are also uninformative e.g. “gp86” which provides no meaningful annotation (Table S2).  

 

15| Annotation with the use of PHROGs. 

 

$ wget http://s3.climb.ac.uk/ADM_share/all_phrogs.hmm.gz -o all_phrogs.hmm.gz  
$ gunzip all_PHROGs.tar.gz 

 

Move PHROGs database to /prokka/db/hmm: 

 

$ mv all_PHROGs.hmm /path/to/prokka/db/hmm 
$ hmmpress /path/to/prokka/db/hmm/all_PHROGs.hmm 

 

prokka --outdir prokka_phrog --prefix mygenome --hmms /path/to/prokka/db/hmm/all_PHROGs.hmm SRR13108336.fasta  

 

16| Annotation with the use of a well annotated closely related phage 

 

$ prokka --outdir prokka_closephage --prefix mygenome --proteins erla_phage.gbk SRR13108336.fasta  

 

In this last case the closest related phage is the submitted Genbank record that these reads are linked to. It allows comparison of the automated annotation here with the annotation of the original submitters. It is useful to view annotations and genes visually and this can be done with tools such as Artemis or Ugene.  

 

Initial rapid analysis of the number of genes predicted to have a function can be determined using some simple linux commands.  

  

For example the number of annotated protein coding sequences can be obtained with: 

 

$ grep -c '>' mygenome.faa  

 

Annotations can be displayed by:  

 

$ grep '>' mygenome.faa |less  

 

Or all proteins that are not hypothetical can be viewed by:  

 

$ grep '>' mygenome.faa | grep -v 'hypothetical protein'| less 

 

Each header in the FASTA file starts with ‘>’, a protein id followed by a number that corresponds to the sequence number, and the product name assigned to the protein sequence. 

 

The GenBank file contains a complete information about annotation of each coding region of a gene (CDS) that can be also examined and compared between results. Such information includes the amino acid sequence of the annotated protein, its coordinates and product description. The field “/inference” contains the information on the database used for annotation, which can be used to determine which external database the annotation was derived from.To view Genbank files in the command line, the following examples can be used: 

 

$ less mygenome.gbk 

or 

$ grep 'product' mygenome.gbk | less 

 

17| Comparison of annotation with different databases 

 

TableS2 contains the lists of protein products that were annotated by the use of different databases within this protocol.The use of different databases clearly highlights the quality of annotation that can be gained from automatic annotation. We recommend the use of PHROGs as a default database or using a closely related phage if there is one available.  

 

The output of automated annotation is merely a starting point and not the end point, prior to genome submission. The genomes should be checked prior to genome submission, the use of tools such as Ugene and Artemis allow visual interaction with genomes. However, as can be clearly seen the automated annotation here is largely congruent with manual annotation described in the original publication of phage Erla. The details included within this text allow getting to the point where annotations can be manually checked as quickly as possible.  

 

18| Genome submission 

 

For submitting the annotated genome to the European Nucleotide Archive (ENA), we will convert the *gff file to an EMBL file for submission. It may be necessary to manual edit the *gff file if further refinement of annotations is required. The conversion of the *gff to embl format will be done with a script available from https://github.com/sanger-pathogens/gff3toembl (See link for full instructions are provided). Briefly to install:  

 

$ git clone git@github.com:sanger-pathogens/gff3toembl.git  
$ python setup.py install 
 
$ gff3_to_embl -i "SubmitterName, A." --genome_type "linear" --classification "PHG" --locus_tag "EXAMPLE_LOCUS" --translation_table 11  "Microbacterium phage Erla" 2792036 'PRJ1234XX' 'Genus species subspecies strain of organism' mygenome.gff 

 

By default this produces a file called output.embl. 

 

The important options to note are:  

 

--genome_type "linear" all phages have linear genomes  

--classification "PHG" PHG is the descriptor for phage  

--locus_tag "EXAMPLE_LOCUS". We suggest the locus_tag is registered in advance  

 

"Microbacterium phage Erla" is the Organism  

 

2792036 is the Taxonomy ID. This has to be registered in advance of submission.  

 

Reorder.gff is the file we are converting  

 

PRJ1234XX is the project number that needs to be registered prior to submission. We recommend prior to sequencing.  

 

The resultant file still needs to be checked prior to submission. We strongly recommend you consult https://ena-docs.readthedocs.io/en/latest/submit/fileprep/flat-file-example.html 

 

Example of the lines from output.embl file:  

ID   XXX; XXX; linear; genomic DNA; STD; PHG; 41538 BP. 
XX 
AC   XXX; 
XX 
AC * _contig 
XX 
PR   Project:PRJ1234XXX; 
XX 
DE   XXX; 
XX 
RN   [1] 
RA   SubmitterName, A.; 
RT   "Draft assembly annotated with Prokka"; 
RL   Submitted (02-Nov-2021) to the INSDC. 
XX 
FH   Key             Location/Qualifiers 
FH 
FT   source          1..41538 
FT                   /organism="Microbacterium phage Erla" 
FT                   /mol_type="genomic DNA" 
FT                   /db_xref="taxon:2792036" 
FT                   /note="reorder" 
FT   CDS             7..525 
FT                   /product="Hypothetical protein" 
FT                   /inference="ab initio prediction:Prodigal:002006" 
FT                   /inference="protein motif:all_phrogs_ann:phrog_3678" 
FT                   /locus_tag="EXAMPLE_LOCUS_00001" 
FT                   /transl_table=11 
FT   CDS             522..1922 

 

 

 

We have highlighted some lines that can be edited. The text “all_phrogs_ann:” is not required and all occurrences can be removed prior to submission. Removal of the text can be achieved with simple linux or find/replace with a text editor.  

 

$ sed -i 's/all_phrogs_ann://gi' output.embl 

 

This uses sed to search for “all_phrogs_ann:” and replace with “” i. e. nothing. This processes, searches and replaces all instances of this string  in the output.embl file. 

 

When all changes have been made to the output.embl file it needs to be gzipped prior to submission:  

 

$ gzip output.embl  

 

This results in the file output.embl.gz. 

 

In addition to the embl file, a manifest file is also required see https://ena-docs.readthedocs.io/en/latest/submit/assembly/genome.html for full details. There are a number of required and optional fields. An example is listed below.  

STUDY           PRJ1234XXX 
SAMPLE          ERS1111111 
ASSEMBLYNAME     PhageAssembly1 
ASSEMBLY_TYPE   isolate 
COVERAGE        3196 
PROGRAM SPAdes v.3.13.0 
PLATFORM        Illumina MiSeq 
FLATFILE        output.embl.gz   

 

 

Again as study number is required, which needs to be registered in advance (https://ena-docs.readthedocs.io/en/latest/submit/assembly/genome.html#stage-1-pre-register-study-and-sample). In our example it is “PRJ1234XXX”. 

 

A sample number, which is also registered in advance. Our example is  ERS1111111. 

 

An assemblyname, which in our example is PhageAssembly1. 

 

Coverage was determined previously as 3196. 

 

Program is the assembly program, in our case SPades v.3.13.0 .  

 

Platform is the sequencing platform used, which in this case is Illumina MiSeq. This information should be supplied by the sequencing provider. 

 

This is the current minimum information that is required for submission in a manifest file, consult https://ena-docs.readthedocs.io/en/latest/submit/assembly/genome.html for full details. 

 

 

In addition to submission of an annotated genome it is good practice to submit raw reads. Again for this a manifest file is required see https://ena-docs.readthedocs.io/en/latest/submit/reads.html for full details. 

 

An example is shown below that would apply for sequencing of a phage with single reads. 

STUDY Requires registration 
SAMPLE Requires registration 
NAME User Defined 
INSTRUMENT Illumina MiSeq 
LIBRARY_SOURCE GENOMIC 
LIBRARY_SELECTION RANDOM 
LIBRARY_STRATEGY WGS 
FASTQ R1.fq.gz 

 

The STUDY requires registration and in our example would be PRJ1234XXX. Notice both reads and genome are part of the same study.  

 

Sample accession requires registration through the website.  

 

NAME - user defined value for experiment name.  

 

Library_Source, Library_Selection and Library_Strategy are controlled vocabulary values. The examples above are for genomic DNA that is randomly selected to produce whole genome sequencing.   

 

As in this example only single end sequencing was used, only one fastq is specified. If paired end sequencing was used then two fastq files would specified e.g.:  

 

FASTQ R1.fq.gz 
FASTQ R2.fq.gz 

 

EMBL flat files and reads can be submitted via the Webin Command Line interface  

 

See https://ena-docs.readthedocs.io/en/latest/submit/general-guide/webin-cli.html for full details  

 

19| Long read Assembly (< 15 min) 

 

Here we will provide a brief example of long read assembly with nanopore reads. A detailed description with all the possible variations is beyond the scope of this work.  

 

Install  FLYE (https://ena-docs.readthedocs.io/en/latest/submit/general-guide/webin-cli.html):  

 

$ git clone https://github.com/fenderglass/Flye  
$ cd Flye 
$ make 
 

 

We will use the minION read set from phage Esa https://journals.asm.org/doi/10.1128/MRA.00730-20. Details of reads are here  

https://www.ebi.ac.uk/ena/browser/view/SRR11901901?show=reads. 

 

Make directory to store the data:  

$mkdir nanopore  
$cd nanopore  

 

$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/001/SRR11901901/SRR11901901_1.fastq.gz 
 

Running flye without any parameters will bring up a list of options that can be used with flye:  

 

$ python <path_to_installation_directory>/bin/flye  

 

We will use:  

 

$ python <path_to_installation_directory>/bin/flye --nano-raw SRR11901901_1.fastq.gz --out-dir ./ --threads 10 --iterations 5 

 

We have selected to save the contents in the current directory (./), use 10 threads and use 5 iterations of polishing.  

 

A number of files and directories are produced on completion: 

 

00-assembly/  

10-consensus/ 

20-repeat/ 

30-contigger/ 

40-polishing/  

assembly.fasta  

assembly_graph.gfa  

assembly_graph.gv  

assembly_info.txt 

 

The file assembly.fasta contains the polished final assembly.  

 

The file assembly_info.txt ($ head assembly_info.txt) contains the following information on the assembly:  

 

#seq_name       length  cov.    circ.   repeat  mult.   alt_group       graph_path 
contig_1           142199  11      Y       N           1             *               1 

 

While a single contig of 143 kb has been assembled it only has 11x coverage.  

 

Using Prokka (Step 15) it is now possible to rapidly annotate this genome:  

 

$ prokka assembly.fasta --prefix nanopore 

 

Viewing the annotation file there are 339 predicted genes. This number is a very large number for a phage of 143 kb. Close inspection of the predicted functions and gene size, quickly identifies lots of short genes. This is a result of the low coverage nanopore sequencing. It is important to note that the published phage genome of phage Esa used a combination of nanopore and illumina reads for genome assembly. The illumina reads are also available from here:  


ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/000/SRR11901900/SRR11901900_1.fastq.gz 
 
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/000/SRR11901900/SRR11901900_2.fastq.gz 


These could be used to polish the nanopore genome pilon or used in assembly, adapting the steps described in this workflow.   

The workflow described in this document contains a number of steps that may seem daunting to a beginner. However, once the software is installed it is possible to go from raw sequence reads to an assembled and initial annotation of the subsequent genome in under a few hours.  