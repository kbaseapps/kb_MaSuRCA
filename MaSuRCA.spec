/*
   Name of module: kb_MaSuRCA

   This KBase module wraps the genome assembly software MaSuRCA(Maryland Super-Read Celera Assembler).
   MaSuRCA 3.2.3


   References:
   https://academic.oup.com/bioinformatics/article/29/21/2669/195975/The-MaSuRCA-genome-assembler
   https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btt476
   ftp://ftp.genome.umd.edu/pub/MaSuRCA/latest/

*/

module kb_MaSuRCA {
    typedef string assembly_ref;

    /* A boolean - 0 for false, 1 for true.
        @range (0, 1)
    */
                    
    typedef int bool;

    /* An X/Y/Z style KBase object reference
    */
    typedef string obj_ref;

   /* The workspace object name of a SingleEndLibrary or PairedEndLibrary file, whether of the
       KBaseAssembly or KBaseFile type.
    */
    typedef string read_lib;

    /* 
        Arguments for run_masurca

        *******for creating the sr_config.txt file*******
        1. DATA
        consisting of 5 fields: 1)two_letter_prefix 2)mean 3)stdev 4)fastq(.gz)_fwd_reads 5)fastq(.gz)_rev_reads.
        e.g.,
                PE= pe 180 20  /FULL_PATH/frag_1.fastq  /FULL_PATH/frag_2.fastq
                JUMP= sh 3600 200  /FULL_PATH/short_1.fastq  /FULL_PATH/short_2.fastq
                #pacbio reads must be in a single fasta file! make sure you provide absolute path
                PACBIO=/FULL_PATH/pacbio.fa
                OTHER=/FULL_PATH/file.frg

        2. PARAMETERS
        string GRAPH_KMER_SIZE - the k-mer size for deBruijn graph values between 25 and 127 are supported, 'auto' will compute the optimal size based on the read data and GC content
        bool USE_LINKING_MATES - set this to 1 for all Illumina-only assemblies; set this to 1 if you have less than 20x long reads (454, Sanger, Pacbio) and less than 50x CLONE coverage by Illumina, Sanger or 454 mate pairs; otherwise keep at 0
        int LIMIT_JUMP_COVERAGE - this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms
        CA_PARAMETERS: these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
        float cgwErrorRate=0.15 - set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
        int KMER_COUNT_THRESHOLD - minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if Illumina coverage >100
        bool CLOSE_GAPS - whether to attempt to close gaps in scaffolds with Illumina data (1) or not (0)
        int NUM_THREADS - auto-detected number of cpus to use
        int JF_SIZE  - this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*estimated_coverage (e.g., 200000000)
        bool SOAP_ASSEMBLY - set this to 1 to use SOAPdenovo contigging/scaffolding module.  Assembly will be worse but will run faster. Useful for very large (>5Gbp) genomes

        string workspace_name - the name of the workspace from which to take input and store output.
        int hash_length - an odd integer (if even, it will be decremented) <= 31
        string output_contigset_name - the name of the output contigset
        list<paired_end_lib> read_libraries - Illumina PairedEndLibrary files to assemble
        min_contig_length - integer to filter out contigs with length < min_contig_length
                     from the MaSuRCA output. Default value is 500 (where 0 implies no filter).

        @optional min_contig_length
        @optional cov_cutoff
        @optional ins_length
        @optional read_trkg
        @optional amos_file
        @optional exp_cov
        @optional long_cov_cutoff
     */

    typedef structure {
        string workspace_name;
        int hash_length;
        list<read_lib> read_libraries; 
        string output_contigset_name;
 
        int min_contig_length;
        float cov_cutoff;
        int ins_length;
        bool read_trkg;
        bool amos_file;
        float exp_cov;
        float long_cov_cutoff;
    } masurcaParams;
    
    /* Output parameter items for run_masurca

    report_name - the name of the KBaseReport.Report workspace object.
    report_ref - the workspace reference of the report.

    */
    typedef structure {
        string report_name;
        string report_ref;
    } masurcaResults;
    
    /* 
        Definition of run_masurca
     */
    funcdef run_masurca(masurcaParams params) returns (masurcaResults output) authentication required;
};
