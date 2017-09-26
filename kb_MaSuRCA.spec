/*
   Name of module: MaSuRCA

   This KBase module wraps the genome assembly software MaSuRCA(Maryland Super-Read Celera Assembler).
   MaSuRCA 3.2.3


   References:
   https://academic.oup.com/bioinformatics/article/29/21/2669/195975/The-MaSuRCA-genome-assembler
   https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btt476
   ftp://ftp.genome.umd.edu/pub/MaSuRCA/latest/

*/

module MaSuRCA {
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
        string graph_kmer_size - the k-mer size for deBruijn graph values between 25 and 127 are supported, 'auto' will compute the optimal size based on the read data and GC content
        bool use_linking_mates - set this to 1 for all Illumina-only assemblies; set this to 1 if you have less than 20x long reads (454, Sanger, Pacbio) and less than 50x CLONE coverage by Illumina, Sanger or 454 mate pairs; otherwise keep at 0
        int limit_jump_coverage - this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms
        CA_PARAMETERS: these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
        float cgwErrorRate=0.15 - set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
        int kmer_count_threshold - minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if Illumina coverage >100
        bool close_gaps - whether to attempt to close gaps in scaffolds with Illumina data (1) or not (0)
        int num_threads - auto-detected number of cpus to use, mandatory
        int jf_size  - this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*estimated_coverage (e.g., 2000000000)
        bool SOAP_ASSEMBLY - set this to 1 to use SOAPdenovo contigging/scaffolding module.  Assembly will be worse but will run faster. Useful for very large (>5Gbp) genomes
        bool do_homopolymer_trim - specifies if we do (1) or do not (0) want to trim long runs of homopolymers 

        string workspace_name - the name of the workspace from which to take input and store output.
        string output_contigset_name - the name of the output contigset
        list<paired_end_lib> read_libraries - Illumina PairedEndLibrary files to assemble

        @optional jump_libraries
        @optional jp_mean
        @optional jp_stdv
        @optional pacbio_reads
        @optional other_frg_file
        @optional graph_kmer_size
        @optional use_linking_mates
        @optional limit_jump_coverage
        @optional cgwErrorRate
        @optional kmer_count_threshold
        @optional close_gaps
        @optional soap_assembly
        @optional do_homopolymer_trim
     */

    typedef structure {
        string workspace_name;
        int num_threads;
        int jf_size;
        list<read_lib> read_libraries; 
        string pe_prefix;
        int pe_mean;
        int pe_stdv;
       
        list<read_lib> jump_libraries;
        string jp_prefix;
        int jp_mean;
        int jp_stdv;
        read_lib pacbio_reads;
        string other_frg_file;   
        string graph_kmer_size;
        bool use_linking_mates;
        int limit_jump_coverage;
        float cgwErrorRate;
        int kmer_count_threshold;
        bool close_gaps;
        bool soap_assembly;
        bool do_homopolymer_trim;

        string out_folder; 
        string output_contigset_name;
        bool create_report;
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
