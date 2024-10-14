This VCF contains 15,225 genotypes for 7,155 wheat samples from the AGENT project and 8,070 IPK samples from the GB2.0 project. There are 151,488 markers included that have a minor allele frequency >= 0.0001. The file is intended for internal use in the AGENT consortium.

The sample designation is in the form of BioSample IDs and information can be retrieved from the BioSamples database (i.e. for sample SAMEA104431476 via https://wwwdev.ebi.ac.uk/biosamples/samples/SAMEA104431476).

To access passport information for batches of samples several options exist:

(1) The file IDs_231012.csv includes AGENT IDs and IPK ACCENUMBs alongside the BioSample IDs. Information thus can be retrieved from the AGENT Portal (https://agent.ipk-gatersleben.de/apex/agent/r/agent/home, username: agent_public, password: agent2020) and from EURISCO (https://eurisco.ipk-gatersleben.de/apex/eurisco_ws/r/eurisco/eurisco-download-by-species) and assigned to genotyped samples.

(2) The python script download_biosamples_agent.py (included in the repo \'93AGENT Barley VCF Assay\'94) can be used to download information from the BioSamples database for a list of BioSample IDs such as BIOSAMPLE_IDs_testlist.txt. You can execute the script in the following way:

python download_biosamples_agent.py --list BIOSAMPLE_IDs_testlist.txt --tsv

This will download one json file for each BioSample ID in --list into the current directory, and option --tsv will save the downloaded information in a table named accession_full_df.tsv.

(3) The IPK BioSamples Query Tool (https://divbrowse.ipk-gatersleben.de/biosamples-tool/#) can be used to download information from the BioSamples database for a list of BioSample IDs via a graphical user interface.

Please note that currently information for AGENT samples are not yet available for download from the BioSamples database, but option (1) - the AGENT Portal is available.