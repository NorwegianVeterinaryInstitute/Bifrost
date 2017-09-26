# Third party software

The Bifrost pipeline depends on several third party packages. 
These have to be made available to the pipeline in some way.
The way that these are made available to the pipeline depends
on which system the pipeline is being run on.

Please note: not all of the software is used for all tracks. 
The track(s) that each software is used in is noted below.


## Currently used software

* FastQC - Track One and Three
* Ariba - Track Two
* Trimmomatic - Track Three
* SPAdes - Track Three
* QUAST - Track Three


## Tracks

There are currently three tracks:

* Track One: FastQC on input reads
* Track Two: Ariba MLST, virulence and AMR analysis
* Track Three: Trimming with trimmomatic followed by assembly
    with SPAdes. Trimming results are evaluated with MultiQC, and
    assemblies with QUAST
    
## Profiles

We currently have two profiles set up, standard and slurm.

### standard.config
This profile is used when running on a normal stand-alone
computer. This assumes that all software is available on
the command line, unless otherwise noted with a full path in
the standard.config file. 

### slurm.config
This profile is used when running on a system that uses the
slurm queue management system. At present, this also depends
heavily on the module system. Any software not in the module
system needs to either be available on the command line, or 
should be specified using the full path. 





