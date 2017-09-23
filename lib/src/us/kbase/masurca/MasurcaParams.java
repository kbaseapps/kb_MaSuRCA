
package us.kbase.masurca;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: masurcaParams</p>
 * <pre>
 * Arguments for run_masurca
 * *******for creating the sr_config.txt file*******
 *  1. DATA
 * consisting of 5 fields: 1)two_letter_prefix 2)mean 3)stdev 4)fastq(.gz)_fwd_reads 5)fastq(.gz)_rev_reads.
 * e.g.,
 *         PE= pe 180 20  /FULL_PATH/frag_1.fastq  /FULL_PATH/frag_2.fastq
 *         JUMP= sh 3600 200  /FULL_PATH/short_1.fastq  /FULL_PATH/short_2.fastq
 *         #pacbio reads must be in a single fasta file! make sure you provide absolute path
 *         PACBIO=/FULL_PATH/pacbio.fa
 *         OTHER=/FULL_PATH/file.frg
 *  2. PARAMETERS
 * string graph_kmer_size - the k-mer size for deBruijn graph values between 25 and 127 are supported, 'auto' will compute the optimal size based on the read data and GC content
 * bool use_linking_mates - set this to 1 for all Illumina-only assemblies; set this to 1 if you have less than 20x long reads (454, Sanger, Pacbio) and less than 50x CLONE coverage by Illumina, Sanger or 454 mate pairs; otherwise keep at 0
 * int limit_jump_coverage - this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms
 * CA_PARAMETERS: these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
 * float cgwErrorRate=0.15 - set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
 * int kmer_count_threshold - minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if Illumina coverage >100
 * bool close_gaps - whether to attempt to close gaps in scaffolds with Illumina data (1) or not (0)
 * int num_threads - auto-detected number of cpus to use, mandatory
 * int jf_size  - this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*estimated_coverage (e.g., 2000000000)
 * bool SOAP_ASSEMBLY - set this to 1 to use SOAPdenovo contigging/scaffolding module.  Assembly will be worse but will run faster. Useful for very large (>5Gbp) genomes
 * bool do_homopolymer_trim - specifies if we do (1) or do not (0) want to trim long runs of homopolymers 
 * string workspace_name - the name of the workspace from which to take input and store output.
 * string output_contigset_name - the name of the output contigset
 * list<paired_end_lib> read_libraries - Illumina PairedEndLibrary files to assemble
 * @optional jump_libraries
 * @optional pacbio_reads
 * @optional other_frg_file
 * @optional graph_kmer_size
 * @optional use_linking_mates
 * @optional limit_jump_coverage
 * @optional cgwErrorRate
 * @optional kmer_count_threshold
 * @optional close_gaps
 * @optional soap_assembly
 * @optional do_homopolymer_trim
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "num_threads",
    "jf_size",
    "read_libraries",
    "jump_libraries",
    "pacbio_reads",
    "other_frg_file",
    "graph_kmer_size",
    "use_linking_mates",
    "limit_jump_coverage",
    "cgwErrorRate",
    "kmer_count_threshold",
    "close_gaps",
    "soap_assembly",
    "do_homopolymer_trim",
    "output_contigset_name",
    "create_report"
})
public class MasurcaParams {

    @JsonProperty("workspace_name")
    private java.lang.String workspaceName;
    @JsonProperty("num_threads")
    private Long numThreads;
    @JsonProperty("jf_size")
    private Long jfSize;
    @JsonProperty("read_libraries")
    private List<String> readLibraries;
    @JsonProperty("jump_libraries")
    private List<String> jumpLibraries;
    @JsonProperty("pacbio_reads")
    private java.lang.String pacbioReads;
    @JsonProperty("other_frg_file")
    private java.lang.String otherFrgFile;
    @JsonProperty("graph_kmer_size")
    private java.lang.String graphKmerSize;
    @JsonProperty("use_linking_mates")
    private Long useLinkingMates;
    @JsonProperty("limit_jump_coverage")
    private Long limitJumpCoverage;
    @JsonProperty("cgwErrorRate")
    private Double cgwErrorRate;
    @JsonProperty("kmer_count_threshold")
    private Long kmerCountThreshold;
    @JsonProperty("close_gaps")
    private Long closeGaps;
    @JsonProperty("soap_assembly")
    private Long soapAssembly;
    @JsonProperty("do_homopolymer_trim")
    private Long doHomopolymerTrim;
    @JsonProperty("output_contigset_name")
    private java.lang.String outputContigsetName;
    @JsonProperty("create_report")
    private Long createReport;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("workspace_name")
    public java.lang.String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(java.lang.String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public MasurcaParams withWorkspaceName(java.lang.String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("num_threads")
    public Long getNumThreads() {
        return numThreads;
    }

    @JsonProperty("num_threads")
    public void setNumThreads(Long numThreads) {
        this.numThreads = numThreads;
    }

    public MasurcaParams withNumThreads(Long numThreads) {
        this.numThreads = numThreads;
        return this;
    }

    @JsonProperty("jf_size")
    public Long getJfSize() {
        return jfSize;
    }

    @JsonProperty("jf_size")
    public void setJfSize(Long jfSize) {
        this.jfSize = jfSize;
    }

    public MasurcaParams withJfSize(Long jfSize) {
        this.jfSize = jfSize;
        return this;
    }

    @JsonProperty("read_libraries")
    public List<String> getReadLibraries() {
        return readLibraries;
    }

    @JsonProperty("read_libraries")
    public void setReadLibraries(List<String> readLibraries) {
        this.readLibraries = readLibraries;
    }

    public MasurcaParams withReadLibraries(List<String> readLibraries) {
        this.readLibraries = readLibraries;
        return this;
    }

    @JsonProperty("jump_libraries")
    public List<String> getJumpLibraries() {
        return jumpLibraries;
    }

    @JsonProperty("jump_libraries")
    public void setJumpLibraries(List<String> jumpLibraries) {
        this.jumpLibraries = jumpLibraries;
    }

    public MasurcaParams withJumpLibraries(List<String> jumpLibraries) {
        this.jumpLibraries = jumpLibraries;
        return this;
    }

    @JsonProperty("pacbio_reads")
    public java.lang.String getPacbioReads() {
        return pacbioReads;
    }

    @JsonProperty("pacbio_reads")
    public void setPacbioReads(java.lang.String pacbioReads) {
        this.pacbioReads = pacbioReads;
    }

    public MasurcaParams withPacbioReads(java.lang.String pacbioReads) {
        this.pacbioReads = pacbioReads;
        return this;
    }

    @JsonProperty("other_frg_file")
    public java.lang.String getOtherFrgFile() {
        return otherFrgFile;
    }

    @JsonProperty("other_frg_file")
    public void setOtherFrgFile(java.lang.String otherFrgFile) {
        this.otherFrgFile = otherFrgFile;
    }

    public MasurcaParams withOtherFrgFile(java.lang.String otherFrgFile) {
        this.otherFrgFile = otherFrgFile;
        return this;
    }

    @JsonProperty("graph_kmer_size")
    public java.lang.String getGraphKmerSize() {
        return graphKmerSize;
    }

    @JsonProperty("graph_kmer_size")
    public void setGraphKmerSize(java.lang.String graphKmerSize) {
        this.graphKmerSize = graphKmerSize;
    }

    public MasurcaParams withGraphKmerSize(java.lang.String graphKmerSize) {
        this.graphKmerSize = graphKmerSize;
        return this;
    }

    @JsonProperty("use_linking_mates")
    public Long getUseLinkingMates() {
        return useLinkingMates;
    }

    @JsonProperty("use_linking_mates")
    public void setUseLinkingMates(Long useLinkingMates) {
        this.useLinkingMates = useLinkingMates;
    }

    public MasurcaParams withUseLinkingMates(Long useLinkingMates) {
        this.useLinkingMates = useLinkingMates;
        return this;
    }

    @JsonProperty("limit_jump_coverage")
    public Long getLimitJumpCoverage() {
        return limitJumpCoverage;
    }

    @JsonProperty("limit_jump_coverage")
    public void setLimitJumpCoverage(Long limitJumpCoverage) {
        this.limitJumpCoverage = limitJumpCoverage;
    }

    public MasurcaParams withLimitJumpCoverage(Long limitJumpCoverage) {
        this.limitJumpCoverage = limitJumpCoverage;
        return this;
    }

    @JsonProperty("cgwErrorRate")
    public Double getCgwErrorRate() {
        return cgwErrorRate;
    }

    @JsonProperty("cgwErrorRate")
    public void setCgwErrorRate(Double cgwErrorRate) {
        this.cgwErrorRate = cgwErrorRate;
    }

    public MasurcaParams withCgwErrorRate(Double cgwErrorRate) {
        this.cgwErrorRate = cgwErrorRate;
        return this;
    }

    @JsonProperty("kmer_count_threshold")
    public Long getKmerCountThreshold() {
        return kmerCountThreshold;
    }

    @JsonProperty("kmer_count_threshold")
    public void setKmerCountThreshold(Long kmerCountThreshold) {
        this.kmerCountThreshold = kmerCountThreshold;
    }

    public MasurcaParams withKmerCountThreshold(Long kmerCountThreshold) {
        this.kmerCountThreshold = kmerCountThreshold;
        return this;
    }

    @JsonProperty("close_gaps")
    public Long getCloseGaps() {
        return closeGaps;
    }

    @JsonProperty("close_gaps")
    public void setCloseGaps(Long closeGaps) {
        this.closeGaps = closeGaps;
    }

    public MasurcaParams withCloseGaps(Long closeGaps) {
        this.closeGaps = closeGaps;
        return this;
    }

    @JsonProperty("soap_assembly")
    public Long getSoapAssembly() {
        return soapAssembly;
    }

    @JsonProperty("soap_assembly")
    public void setSoapAssembly(Long soapAssembly) {
        this.soapAssembly = soapAssembly;
    }

    public MasurcaParams withSoapAssembly(Long soapAssembly) {
        this.soapAssembly = soapAssembly;
        return this;
    }

    @JsonProperty("do_homopolymer_trim")
    public Long getDoHomopolymerTrim() {
        return doHomopolymerTrim;
    }

    @JsonProperty("do_homopolymer_trim")
    public void setDoHomopolymerTrim(Long doHomopolymerTrim) {
        this.doHomopolymerTrim = doHomopolymerTrim;
    }

    public MasurcaParams withDoHomopolymerTrim(Long doHomopolymerTrim) {
        this.doHomopolymerTrim = doHomopolymerTrim;
        return this;
    }

    @JsonProperty("output_contigset_name")
    public java.lang.String getOutputContigsetName() {
        return outputContigsetName;
    }

    @JsonProperty("output_contigset_name")
    public void setOutputContigsetName(java.lang.String outputContigsetName) {
        this.outputContigsetName = outputContigsetName;
    }

    public MasurcaParams withOutputContigsetName(java.lang.String outputContigsetName) {
        this.outputContigsetName = outputContigsetName;
        return this;
    }

    @JsonProperty("create_report")
    public Long getCreateReport() {
        return createReport;
    }

    @JsonProperty("create_report")
    public void setCreateReport(Long createReport) {
        this.createReport = createReport;
    }

    public MasurcaParams withCreateReport(Long createReport) {
        this.createReport = createReport;
        return this;
    }

    @JsonAnyGetter
    public Map<java.lang.String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(java.lang.String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public java.lang.String toString() {
        return ((((((((((((((((((((((((((((((((((((("MasurcaParams"+" [workspaceName=")+ workspaceName)+", numThreads=")+ numThreads)+", jfSize=")+ jfSize)+", readLibraries=")+ readLibraries)+", jumpLibraries=")+ jumpLibraries)+", pacbioReads=")+ pacbioReads)+", otherFrgFile=")+ otherFrgFile)+", graphKmerSize=")+ graphKmerSize)+", useLinkingMates=")+ useLinkingMates)+", limitJumpCoverage=")+ limitJumpCoverage)+", cgwErrorRate=")+ cgwErrorRate)+", kmerCountThreshold=")+ kmerCountThreshold)+", closeGaps=")+ closeGaps)+", soapAssembly=")+ soapAssembly)+", doHomopolymerTrim=")+ doHomopolymerTrim)+", outputContigsetName=")+ outputContigsetName)+", createReport=")+ createReport)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
