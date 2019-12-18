package ca.corefacility.bioinformatics.irida.plugin.staramr;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import ca.corefacility.bioinformatics.irida.exceptions.IridaWorkflowNotFoundException;
import ca.corefacility.bioinformatics.irida.exceptions.PostProcessingException;
import ca.corefacility.bioinformatics.irida.model.sample.MetadataTemplateField;
import ca.corefacility.bioinformatics.irida.model.sample.Sample;
import ca.corefacility.bioinformatics.irida.model.sample.metadata.MetadataEntry;
import ca.corefacility.bioinformatics.irida.model.sample.metadata.PipelineProvidedMetadataEntry;
import ca.corefacility.bioinformatics.irida.model.workflow.IridaWorkflow;
import ca.corefacility.bioinformatics.irida.model.workflow.analysis.AnalysisOutputFile;
import ca.corefacility.bioinformatics.irida.model.workflow.analysis.type.AnalysisType;
import ca.corefacility.bioinformatics.irida.model.workflow.submission.AnalysisSubmission;
import ca.corefacility.bioinformatics.irida.pipeline.results.updater.AnalysisSampleUpdater;
import ca.corefacility.bioinformatics.irida.service.sample.MetadataTemplateService;
import ca.corefacility.bioinformatics.irida.service.sample.SampleService;
import ca.corefacility.bioinformatics.irida.service.workflow.IridaWorkflowsService;

/**
 * {@link AnalysisSampleUpdater} for AMR detection results to be written to
 * metadata of {@link Sample}s.
 */
public class StarAMRUpdater implements AnalysisSampleUpdater {
	private static final Logger logger = LoggerFactory.getLogger(StarAMRUpdater.class);
	private static final String STARAMR_SUMMARY = "staramr-summary.tsv";

	private static final Splitter SPLITTER = Splitter.on('\t');

	private static final String METADATA_STARAMR_GENE = "staramr/gene";
	private static final String METADATA_STARAMR_DRUG_CLASS = "staramr/drug-class";
	private static final String METADATA_STARAMR_QUALITY_MODULE = "staramr/quality-module";
	private static final String METADATA_STARAMR_PLASMID = "staramr/plasmid";
	private static final String METADATA_STARAMR_SCHEME = "staramr/scheme";
	private static final String METADATA_STARAMR_SEQUENCE_TYPE = "staramr/sequence-type";
	private static final String METADATA_STARAMR_GENOME_LENGTH = "staramr/genome-length";
	private static final String METADATA_STARAMR_N50 = "staramr/N50";
	private static final String METADATA_STARAMR_NUM_CONTIGS = "staramr/num-contigs";
	private static final String METADATA_STARAMR_QUALITY_MODULE_FEEDBACK = "staramr/quality-module-feedback";

	private MetadataTemplateService metadataTemplateService;
	private SampleService sampleService;
	private IridaWorkflowsService iridaWorkflowsService;

	/**
	 * Builds a new {@link StarAMRUpdater} with the following information.
	 * 
	 * @param metadataTemplateService The {@link MetadatTemplateService}.
	 * @param sampleService           The {@link SampleService}.
	 * @param iridaWorkflowsService   The {@link IridaWorkflowsService}.
	 */
	public StarAMRUpdater(MetadataTemplateService metadataTemplateService, SampleService sampleService,
			IridaWorkflowsService iridaWorkflowsService) {
		this.metadataTemplateService = metadataTemplateService;
		this.sampleService = sampleService;
		this.iridaWorkflowsService = iridaWorkflowsService;
	}

	/**
	 * Gets the staramr results from the given output file.
	 * 
	 * @param staramrFilePath The staramr output file containing the results.
	 * @return A {@link AMRResult} storing the results from staramr as
	 *         {@link String}s.
	 * @throws IOException             If there was an issue reading the file.
	 * @throws PostProcessingException If there was an issue parsing the file.
	 */
	private AMRResult getStarAMRResults(Path staramrFilePath) throws IOException, PostProcessingException {
		final int GENOTYPE_INDEX = 2;
		final int DRUG_INDEX = 3;
		final int QUALITY_MODULE_INDEX = 1;
		final int PLASMID_INDEX = 4;
		final int SCHEME_INDEX = 5;
		final int SEQUENCE_TYPE_INDEX = 6;
		final int GENOME_LENGTH_INDEX = 7;
		final int N50_INDEX = 8;
		final int NUM_CONTIGS_INDEX = 9;
		final int QUALITY_MODULE_FEEDBACK_INDEX = 10;

		final int MAX_TOKENS = 11;

		@SuppressWarnings("resource")
		BufferedReader reader = new BufferedReader(new FileReader(staramrFilePath.toFile()));
		String line = reader.readLine();
		List<String> tokens = SPLITTER.splitToList(line);
		if (tokens.size() != MAX_TOKENS) {
			throw new PostProcessingException("Invalid number of columns in staramr results file [" + staramrFilePath
					+ "], expected [" + MAX_TOKENS + "] got [" + tokens.size() + "]");
		}

		line = reader.readLine();
		tokens = SPLITTER.splitToList(line);
		String genotype = tokens.get(GENOTYPE_INDEX);
		String drug = tokens.get(DRUG_INDEX);
		String quality_module = tokens.get(QUALITY_MODULE_INDEX);
		String plasmid = tokens.get(PLASMID_INDEX);
		String scheme = tokens.get(SCHEME_INDEX);
		String sequence_type = tokens.get(SEQUENCE_TYPE_INDEX);
		String genome_length = tokens.get(GENOME_LENGTH_INDEX);
		String N50 = tokens.get(N50_INDEX);
		String num_contigs = tokens.get(NUM_CONTIGS_INDEX);
		String quality_module_feedback = tokens.get(QUALITY_MODULE_FEEDBACK_INDEX);

		line = reader.readLine();

		if (line == null) {
			return new AMRResult(genotype, drug, quality_module, plasmid, scheme, sequence_type, genome_length, N50, num_contigs, quality_module_feedback);
		} else {
			throw new PostProcessingException("Invalid number of results in staramr results file [" + staramrFilePath
					+ "], expected only one line of results but got multiple lines");
		}
	}


	public void update(Collection<Sample> samples, AnalysisSubmission analysis) throws PostProcessingException {
		if (samples.size() != 1) {
			throw new PostProcessingException(
					"Expected one sample; got '" + samples.size() + "' for analysis [id=" + analysis.getId() + "]");
		}

		final Sample sample = samples.iterator().next();

		AnalysisOutputFile staramrAof = analysis.getAnalysis().getAnalysisOutputFile(STARAMR_SUMMARY);

		Path staramrFilePath = staramrAof.getFile();


		Map<String, MetadataEntry> stringEntries = new HashMap<>();
		try {
			IridaWorkflow iridaWorkflow = iridaWorkflowsService.getIridaWorkflow(analysis.getWorkflowId());
			String workflowVersion = iridaWorkflow.getWorkflowDescription().getVersion();

			AMRResult staramrResult = getStarAMRResults(staramrFilePath);

			PipelineProvidedMetadataEntry staramrGenotypeEntry = new PipelineProvidedMetadataEntry(
					staramrResult.getGenotype(), "text", analysis);
			PipelineProvidedMetadataEntry staramrDrugEntry = new PipelineProvidedMetadataEntry(
					staramrResult.getDrugClass(), "text", analysis);
			PipelineProvidedMetadataEntry staramrQualityModuleEntry = new PipelineProvidedMetadataEntry(
					staramrResult.getQualityModule(), "text", analysis);
			PipelineProvidedMetadataEntry staramrPlasmidEntry = new PipelineProvidedMetadataEntry(
					staramrResult.getPlasmid(), "text", analysis);
			PipelineProvidedMetadataEntry staramrSchemeEntry = new PipelineProvidedMetadataEntry(
					staramrResult.getScheme(), "text", analysis);
			PipelineProvidedMetadataEntry staramrSequenceTypeEntry = new PipelineProvidedMetadataEntry(
					staramrResult.getSequenceType(), "text", analysis);
			PipelineProvidedMetadataEntry staramrGenomeLengthEntry = new PipelineProvidedMetadataEntry(
					staramrResult.getGenomeLength(), "text", analysis);
			PipelineProvidedMetadataEntry staramrN50Entry = new PipelineProvidedMetadataEntry(
					staramrResult.getN50(), "text", analysis);
			PipelineProvidedMetadataEntry staramrNumContigsEntry = new PipelineProvidedMetadataEntry(
					staramrResult.getNumContigs(), "text", analysis);
			PipelineProvidedMetadataEntry staramrQualityModuleFeedbackEntry = new PipelineProvidedMetadataEntry(
					staramrResult.getQualityModuleFeedback(), "text", analysis);

			stringEntries.put(appendVersion(METADATA_STARAMR_GENE, workflowVersion), staramrGenotypeEntry);
			stringEntries.put(appendVersion(METADATA_STARAMR_DRUG_CLASS, workflowVersion), staramrDrugEntry);
			stringEntries.put(appendVersion(METADATA_STARAMR_QUALITY_MODULE, workflowVersion), staramrQualityModuleEntry);
			stringEntries.put(appendVersion(METADATA_STARAMR_PLASMID, workflowVersion), staramrPlasmidEntry);
			stringEntries.put(appendVersion(METADATA_STARAMR_SCHEME, workflowVersion), staramrSchemeEntry);
			stringEntries.put(appendVersion(METADATA_STARAMR_SEQUENCE_TYPE, workflowVersion), staramrSequenceTypeEntry);
			stringEntries.put(appendVersion(METADATA_STARAMR_GENOME_LENGTH, workflowVersion), staramrGenomeLengthEntry);
			stringEntries.put(appendVersion(METADATA_STARAMR_N50, workflowVersion), staramrN50Entry);
			stringEntries.put(appendVersion(METADATA_STARAMR_NUM_CONTIGS, workflowVersion), staramrNumContigsEntry);
			stringEntries.put(appendVersion(METADATA_STARAMR_QUALITY_MODULE_FEEDBACK, workflowVersion), staramrQualityModuleFeedbackEntry);

			Map<MetadataTemplateField, MetadataEntry> metadataMap = metadataTemplateService
					.getMetadataMap(stringEntries);

			sample.mergeMetadata(metadataMap);
			sampleService.updateFields(sample.getId(), ImmutableMap.of("metadata", sample.getMetadata()));
		} catch (IOException e) {
			logger.error("Got IOException", e);
			throw new PostProcessingException("Error parsing staramr results", e);
		} catch (IridaWorkflowNotFoundException e) {
			logger.error("Got IridaWorkflowNotFoundException", e);
			throw new PostProcessingException("Workflow is not found", e);
		} catch (Exception e) {
			logger.error("Got Exception", e);
			throw e;
		}
	}

	/**
	 * Appends the name and version together for a metadata field name.
	 * 
	 * @param name    The name.
	 * @param version The version.
	 * @return The appended name and version.
	 */
	private String appendVersion(String name, String version) {
		return name + "/" + version;
	}

	/**
	 * Class used to store together data extracted from the staramr tables.
	 */
	private class AMRResult {
		private String genotype;
		private String drugClass;
		private String qualitymodule;
		private String plasmid;
		private String scheme;
		private String sequencetype;
		private String genomelength;
		private String N50;
		private String numcontigs;
		private String qualitymodulefeedback;

		public AMRResult(String genotype, String drugClass, String qualitymodule, String plasmid, 
        String scheme, String sequencetype, String genomelength, String N50, String numcontigs, String qualitymodulefeedback) {
			this.genotype = genotype;
			this.drugClass = drugClass;
			this.qualitymodule = qualitymodule;
			this.plasmid = plasmid;
			this.scheme = scheme;
			this.sequencetype = sequencetype;
			this.genomelength = genomelength;
			this.N50 = N50;
			this.numcontigs = numcontigs;
			this.qualitymodulefeedback = qualitymodulefeedback;
		}

		public String getGenotype() {
			return genotype;
		}

		public String getDrugClass() {
			return drugClass;
		}

		public String getQualityModule() {
			return qualitymodule;
		}
		public String getPlasmid() {
			return plasmid;
		}
		public String getScheme() {
			return scheme;
		}
		public String getSequenceType() {
			return sequencetype;
		}
		public String getGenomeLength() {
			return genomelength;
		}
		public String getN50() {
			return N50;
		}
		public String getNumContigs() {
			return numcontigs;
		}
		public String getQualityModuleFeedback() {
			return qualitymodulefeedback;
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public AnalysisType getAnalysisType() {
		return StarAMRPlugin.STAR_AMR;
	}
}
