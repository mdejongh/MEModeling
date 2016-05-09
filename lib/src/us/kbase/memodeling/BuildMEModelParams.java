
package us.kbase.memodeling;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: BuildMEModelParams</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "model_ref",
    "workspace",
    "output_id"
})
public class BuildMEModelParams {

    @JsonProperty("model_ref")
    private String modelRef;
    @JsonProperty("workspace")
    private String workspace;
    @JsonProperty("output_id")
    private String outputId;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("model_ref")
    public String getModelRef() {
        return modelRef;
    }

    @JsonProperty("model_ref")
    public void setModelRef(String modelRef) {
        this.modelRef = modelRef;
    }

    public BuildMEModelParams withModelRef(String modelRef) {
        this.modelRef = modelRef;
        return this;
    }

    @JsonProperty("workspace")
    public String getWorkspace() {
        return workspace;
    }

    @JsonProperty("workspace")
    public void setWorkspace(String workspace) {
        this.workspace = workspace;
    }

    public BuildMEModelParams withWorkspace(String workspace) {
        this.workspace = workspace;
        return this;
    }

    @JsonProperty("output_id")
    public String getOutputId() {
        return outputId;
    }

    @JsonProperty("output_id")
    public void setOutputId(String outputId) {
        this.outputId = outputId;
    }

    public BuildMEModelParams withOutputId(String outputId) {
        this.outputId = outputId;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((("BuildMEModelParams"+" [modelRef=")+ modelRef)+", workspace=")+ workspace)+", outputId=")+ outputId)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
