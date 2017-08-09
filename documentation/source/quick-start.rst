
.. _quickstart:

Canu Quick Start
================

Canu specializes in assembling PacBio or Oxford Nanopore sequences.  Canu operates in three phases:
correction, trimming and assembly.  The correction phase will improve the accuracy of bases in
reads.  The trimming phase will trim reads to the portion that appears to be high-quality sequence,
removing suspicious regions such as remaining SMRTbell adapter.  The assembly phase will order the
reads into contigs, generate consensus sequences and create graphs of alternate paths.

For eukaryotic genomes, coverage more than 20x is enough to outperform current hybrid methods,
however, between 30x and 60x coverage is the recommended minimum.  More coverage will let Canu use
longer reads for assembly, which will result in better assemblies.

Input sequences can be FASTA or FASTQ format, uncompressed or compressed with gzip (.gz), bzip2
(.bz2) or xz (.xz).  Note that zip files (.zip) are not supported.

Canu can resume incomplete assemblies, allowing for recovery from system outages or other abnormal
terminations.

Canu will auto-detect computational resources and scale itself to fit, using all of the resources
available and are reasonable for the size of your assembly.  Memory and processors can be explicitly
limited with with parameters :ref:`maxMemory` and :ref:`maxThreads`.  See section :ref:`execution`
for more details.

Canu will automaticall take full advantage of any LSF/PBS/PBSPro/Torque/Slrum/SGE grid available,
even submitting itself for execution.  Canu makes heavy use of array jobs and requires job
submission from compute nodes, which are sometimes not available or allowed.  Canu option
``useGrid=false`` will restrict Canu to using only the current machine, while option
``useGrid=remote`` will configure Canu for grid execution but not submit jobs to the grid.
See section :ref:`execution` for more details.

The :ref:`tutorial` has more background, and the :ref:`faq` has a wealth of practical advice.


Assembling PacBio or Nanopore data
----------------------

Pacific Biosciences released P6-C4 chemistry reads for Escherichia coli K12.  You can `download them
from their original release
<https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly>`_, but note that you
must have the `SMRTpipe software <http://www.pacb.com/support/software-downloads/>`_ installed to
extract the reads as FASTQ.  Instead, use a `FASTQ format 25X subset
<http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p6_25x.filtered.fastq>`_ (223MB).  Download from the command line
with::

 curl -L -o pacbio.fastq http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p6_25x.filtered.fastq

There doesn't appear to be any "official" Oxford Nanopore sample data, but the `Loman Lab
<http://lab.loman.net/>`_ released a `set of runs
<http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/>`_, also for Escherichia coli K12.
This is early data, from September 2015.  Any of the four runs will work; we picked `MAP-006-1
<http://nanopore.s3.climb.ac.uk/MAP006-PCR-1_2D_pass.fasta>`_ (243 MB).  Download from the command
line with::

 curl -L -o oxford.fasta http://nanopore.s3.climb.ac.uk/MAP006-PCR-1_2D_pass.fasta

By default, Canu will correct the reads, then trim the reads, then assemble the reads to unitigs.
Canu needs to know the approximate genome size (so it can determine coverage in the input reads)
and the technology used to generate the reads.

For PacBio::

 canu \
  -p ecoli -d ecoli-pacbio \
  genomeSize=4.8m \
  -pacbio-raw pacbio.fastq

For Nanopore::

 canu \
  -p ecoli -d ecoli-oxford \
  genomeSize=4.8m \
  -nanopore-raw oxford.fasta


Output and intermediate files will be in directories 'ecoli-pacbio' and 'ecoli-nanopore',
respectively.  Intermeditate files are written in directories 'correction', 'trimming' and
'unitigging' for the respective stages.  Output files are named using the '-p' prefix, such as
'ecoli.contigs.fasta', 'ecoli.contigs.gfa', etc.  See section :ref:`outputs` for more details on
outputs (intermediate files aren't documented).


Assembling With Multiple Technologies and Multiple Files
-------------------------------------------

Canu can use reads from any number of input files, which can be a mix of formats and technologies.
We'll assemble a mix of 10X PacBio reads in two FASTQ files and 10X of Nanopore reads in one FASTA
file::

 curl -L -o mix.tar.gz http://gembox.cbcb.umd.edu/mhap/raw/ecoliP6Oxford.tar.gz
 tar xvzf mix.tar.gz
 
 canu \
  -p ecoli -d ecoli-mix \
  genomeSize=4.8m \
  -pacbio-raw pacbio.part?.fastq.gz \
  -nanopore-raw oxford.fasta.gz


Correct, Trim and Assemble, Manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes, however, it makes sense to do the three top-level tasks by hand.  This would allow trying
multiple unitig construction parameters on the same set of corrected and trimmed reads, or skipping
trimming and assembly if you only want correced reads.

We'll use the PacBio reads from above.  First, correct the raw reads::

 canu -correct \
   -p ecoli -d ecoli \
   genomeSize=4.8m \
   -pacbio-raw  pacbio.fastq

Then, trim the output of the correction::

 canu -trim \
   -p ecoli -d ecoli \
   genomeSize=4.8m \
   -pacbio-corrected ecoli/ecoli.correctedReads.fasta.gz

And finally, assemble the output of trimming, twice, with different stringency on which overlaps to
use (see :ref:`correctedErrorRate <correctedErrorRate>`)::

 canu -assemble \
   -p ecoli -d ecoli-erate-0.039 \
   genomeSize=4.8m \
   correctedErrorRate=0.039 \
   -pacbio-corrected ecoli/ecoli.trimmedReads.fasta.gz

 canu -assemble \
   -p ecoli -d ecoli-erate-0.075 \
   genomeSize=4.8m \
   correctedErrorRate=0.075 \
   -pacbio-corrected ecoli/ecoli.trimmedReads.fasta.gz

Note that the assembly stages use different '-d' directories.  It is not possible to run multiple
copies of canu with the same work directory.


Assembling Low Coverage Datasets
----------------------------------

We claimed Canu works down to 20X coverage, and we will now assemble `a 20X subset of S. cerevisae
<http://gembox.cbcb.umd.edu/mhap/raw/yeast_filtered.20x.fastq.gz>`_ (215 MB).  When assembling, we
adjust :ref:`correctedErrorRate <correctedErrorRate>` to accomodate the slightly lower
quality corrected reads::

 curl -L -o yeast.20x.fastq.gz http://gembox.cbcb.umd.edu/mhap/raw/yeast_filtered.20x.fastq.gz

 canu \
  -p asm -d yeast \
  genomeSize=12.1m \
  correctedErrorRate=0.075 \
  -pacbio-raw yeast.20x.fastq.gz

Consensus Accuracy
-------------------

Canu consensus sequences are typically well above 99% identity.  Accuracy can be improved by
polishing the contigs with tools developed specifically for that task.  We recommend `Quiver
<http://github.com/PacificBiosciences/GenomicConsensus>`_ for PacBio and `Nanopolish
<http://github.com/jts/nanopolish>`_ for Oxford Nanpore data.
When Illumina reads are available, `Pilon <http://www.broadinstitute.org/software/pilon/>`_
can be used to polish either PacBio or Oxford Nanopore assemblies.
