# dna_context
Starting with a database of genomes, pull out ORFs, use HMMER to identify homologs to a given protein, and identify the DNA context of that protein



#Pre steps
A) generate pHMM for sequence of interest. This can be done in multiple ways. One is to make a muscle alignment using a web interface (a .clw file). 
With a .clw file in an environment with hmm, run
>hmmbuild <output_file_name.hmm> <input_file_name.clw>

Then, put this .hmm file in dna_context/hmmer/proteins
and make a manifest of the proteins to be searched for called "GOI_manifest.txt" that has <output_file_name.hmm>

B) Acquire database of genomes of interest. Put the DNA fasta files in dna_context/database/genomes

##Step1##: Run orfipy

When running orfipy, you could change the desired translation table by adding "--table 1" flag to the orfipy command
from scripts/ copy the 'mass_orfipy.sh' script into dna_context/orfs/orfipy/ 

The script must be in the folder where you want the output to land. You should delete the script after the command has completed

In a conda environment with orfipy (https://github.com/urmi-21/orfipy), run the mass orfipy script
> bash mass_orfipy.sh

This generates folders for each genome include a .faa file with the genome name which is all its orfs

##Step2##: Collect Orfs

From scripts, if you open 'collect_orfipy_files.txt' you'll find the command that should be executed at dna_context/scripts/ that will pull
the orfs from their subfolders into the extracted_orfs/ folder

> find ../orfs/orfipy/ -name '*.faa' -exec cp -t ../orfs/extracted_orfs/ {} +

##Step3##: Make manifest

For the scripts, there must be a perfect match between file names. Navigate to /dna_context/orfs/extracted_orfs/ and execute
> ls >> bacteroides_manifest.txt
> mv bacteroides_manifest.txt ../../hmmer/

##Step 4##: Run kofamscan loop

In conda environment with kofamscan, in dna_context/hmmer/ execute the loop that iterates over kofamscan
https://github.com/takaram/kofam_scan

> bash hmmer_loop_gen4.sh

##Step5##: Capture e values, amino acid sequences, and DNA context

in dna_context/scripts/ with conda with python 3, execute the capture_relevant_results.py script
You could change the output file names at the top, or the length of the DNA that is captured
> python capture_relevant_results.py





