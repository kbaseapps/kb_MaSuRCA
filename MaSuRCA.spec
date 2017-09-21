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
