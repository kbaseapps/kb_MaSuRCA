
package us.kbase.kbmasurca;

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
 * string workspace_name - the name of the workspace from which to take input and store output.
 * int hash_length - an odd integer (if even, it will be decremented) <= 31
 * string output_contigset_name - the name of the output contigset
 * list<paired_end_lib> read_libraries - Illumina PairedEndLibrary files to assemble
 * min_contig_length - integer to filter out contigs with length < min_contig_length
 *              from the MaSuRCA output. Default value is 500 (where 0 implies no filter).
 * @optional min_contig_length
 * @optional cov_cutoff
 * @optional ins_length
 * @optional read_trkg
 * @optional amos_file
 * @optional exp_cov
 * @optional long_cov_cutoff
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "hash_length",
    "read_libraries",
    "output_contigset_name",
    "min_contig_length",
    "cov_cutoff",
    "ins_length",
    "read_trkg",
    "amos_file",
    "exp_cov",
    "long_cov_cutoff"
})
public class MasurcaParams {

    @JsonProperty("workspace_name")
    private java.lang.String workspaceName;
    @JsonProperty("hash_length")
    private Long hashLength;
    @JsonProperty("read_libraries")
    private List<String> readLibraries;
    @JsonProperty("output_contigset_name")
    private java.lang.String outputContigsetName;
    @JsonProperty("min_contig_length")
    private Long minContigLength;
    @JsonProperty("cov_cutoff")
    private Double covCutoff;
    @JsonProperty("ins_length")
    private Long insLength;
    @JsonProperty("read_trkg")
    private Long readTrkg;
    @JsonProperty("amos_file")
    private Long amosFile;
    @JsonProperty("exp_cov")
    private Double expCov;
    @JsonProperty("long_cov_cutoff")
    private Double longCovCutoff;
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

    @JsonProperty("hash_length")
    public Long getHashLength() {
        return hashLength;
    }

    @JsonProperty("hash_length")
    public void setHashLength(Long hashLength) {
        this.hashLength = hashLength;
    }

    public MasurcaParams withHashLength(Long hashLength) {
        this.hashLength = hashLength;
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

    @JsonProperty("min_contig_length")
    public Long getMinContigLength() {
        return minContigLength;
    }

    @JsonProperty("min_contig_length")
    public void setMinContigLength(Long minContigLength) {
        this.minContigLength = minContigLength;
    }

    public MasurcaParams withMinContigLength(Long minContigLength) {
        this.minContigLength = minContigLength;
        return this;
    }

    @JsonProperty("cov_cutoff")
    public Double getCovCutoff() {
        return covCutoff;
    }

    @JsonProperty("cov_cutoff")
    public void setCovCutoff(Double covCutoff) {
        this.covCutoff = covCutoff;
    }

    public MasurcaParams withCovCutoff(Double covCutoff) {
        this.covCutoff = covCutoff;
        return this;
    }

    @JsonProperty("ins_length")
    public Long getInsLength() {
        return insLength;
    }

    @JsonProperty("ins_length")
    public void setInsLength(Long insLength) {
        this.insLength = insLength;
    }

    public MasurcaParams withInsLength(Long insLength) {
        this.insLength = insLength;
        return this;
    }

    @JsonProperty("read_trkg")
    public Long getReadTrkg() {
        return readTrkg;
    }

    @JsonProperty("read_trkg")
    public void setReadTrkg(Long readTrkg) {
        this.readTrkg = readTrkg;
    }

    public MasurcaParams withReadTrkg(Long readTrkg) {
        this.readTrkg = readTrkg;
        return this;
    }

    @JsonProperty("amos_file")
    public Long getAmosFile() {
        return amosFile;
    }

    @JsonProperty("amos_file")
    public void setAmosFile(Long amosFile) {
        this.amosFile = amosFile;
    }

    public MasurcaParams withAmosFile(Long amosFile) {
        this.amosFile = amosFile;
        return this;
    }

    @JsonProperty("exp_cov")
    public Double getExpCov() {
        return expCov;
    }

    @JsonProperty("exp_cov")
    public void setExpCov(Double expCov) {
        this.expCov = expCov;
    }

    public MasurcaParams withExpCov(Double expCov) {
        this.expCov = expCov;
        return this;
    }

    @JsonProperty("long_cov_cutoff")
    public Double getLongCovCutoff() {
        return longCovCutoff;
    }

    @JsonProperty("long_cov_cutoff")
    public void setLongCovCutoff(Double longCovCutoff) {
        this.longCovCutoff = longCovCutoff;
    }

    public MasurcaParams withLongCovCutoff(Double longCovCutoff) {
        this.longCovCutoff = longCovCutoff;
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
        return ((((((((((((((((((((((((("MasurcaParams"+" [workspaceName=")+ workspaceName)+", hashLength=")+ hashLength)+", readLibraries=")+ readLibraries)+", outputContigsetName=")+ outputContigsetName)+", minContigLength=")+ minContigLength)+", covCutoff=")+ covCutoff)+", insLength=")+ insLength)+", readTrkg=")+ readTrkg)+", amosFile=")+ amosFile)+", expCov=")+ expCov)+", longCovCutoff=")+ longCovCutoff)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
