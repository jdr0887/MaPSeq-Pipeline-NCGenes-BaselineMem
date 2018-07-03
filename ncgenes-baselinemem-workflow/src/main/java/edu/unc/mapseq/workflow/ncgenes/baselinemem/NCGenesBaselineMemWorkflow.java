package edu.unc.mapseq.workflow.ncgenes.baselinemem;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.lang.StringUtils;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.renci.jlrm.condor.CondorJob;
import org.renci.jlrm.condor.CondorJobBuilder;
import org.renci.jlrm.condor.CondorJobEdge;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.commons.ncgenes.baselinemem.RegisterToIRODSRunnable;
import edu.unc.mapseq.commons.ncgenes.baselinemem.SaveDepthOfCoverageAttributesRunnable;
import edu.unc.mapseq.commons.ncgenes.baselinemem.SaveFlagstatAttributesRunnable;
import edu.unc.mapseq.commons.ncgenes.baselinemem.SaveMarkDuplicatesAttributesRunnable;
import edu.unc.mapseq.dao.model.Attribute;
import edu.unc.mapseq.dao.model.Flowcell;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.dao.model.WorkflowRun;
import edu.unc.mapseq.dao.model.WorkflowRunAttempt;
import edu.unc.mapseq.module.sequencing.WriteVCFHeaderCLI;
import edu.unc.mapseq.module.sequencing.bwa.BWAMEMCLI;
import edu.unc.mapseq.module.sequencing.fastqc.FastQCCLI;
import edu.unc.mapseq.module.sequencing.fastqc.IgnoreLevelType;
import edu.unc.mapseq.module.sequencing.filter.FilterVariantCLI;
import edu.unc.mapseq.module.sequencing.gatk.GATKApplyRecalibrationCLI;
import edu.unc.mapseq.module.sequencing.gatk.GATKCountCovariatesCLI;
import edu.unc.mapseq.module.sequencing.gatk.GATKDepthOfCoverageCLI;
import edu.unc.mapseq.module.sequencing.gatk.GATKDownsamplingType;
import edu.unc.mapseq.module.sequencing.gatk.GATKFlagStatCLI;
import edu.unc.mapseq.module.sequencing.gatk.GATKIndelRealignerCLI;
import edu.unc.mapseq.module.sequencing.gatk.GATKPhoneHomeType;
import edu.unc.mapseq.module.sequencing.gatk.GATKRealignerTargetCreatorCLI;
import edu.unc.mapseq.module.sequencing.gatk.GATKTableRecalibrationCLI;
import edu.unc.mapseq.module.sequencing.gatk.GATKUnifiedGenotyperCLI;
import edu.unc.mapseq.module.sequencing.gatk.GATKVariantRecalibratorCLI;
import edu.unc.mapseq.module.sequencing.picard.PicardAddOrReplaceReadGroupsCLI;
import edu.unc.mapseq.module.sequencing.picard.PicardFixMateCLI;
import edu.unc.mapseq.module.sequencing.picard.PicardMarkDuplicatesCLI;
import edu.unc.mapseq.module.sequencing.picard.PicardSortOrderType;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsFlagstatCLI;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsIndexCLI;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsViewCLI;
import edu.unc.mapseq.workflow.WorkflowException;
import edu.unc.mapseq.workflow.sequencing.AbstractSequencingWorkflow;
import edu.unc.mapseq.workflow.sequencing.SequencingWorkflowJobFactory;
import edu.unc.mapseq.workflow.sequencing.SequencingWorkflowUtil;

public class NCGenesBaselineMemWorkflow extends AbstractSequencingWorkflow {

    private static final Logger logger = LoggerFactory.getLogger(NCGenesBaselineMemWorkflow.class);

    public NCGenesBaselineMemWorkflow() {
        super();
    }

    @Override
    public Graph<CondorJob, CondorJobEdge> createGraph() throws WorkflowException {
        logger.info("ENTERING createGraph()");

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(CondorJobEdge.class);

        int count = 0;
        String percentBadVariants = "0.05";

        Set<Sample> sampleSet = SequencingWorkflowUtil.getAggregatedSamples(getWorkflowBeanService().getMaPSeqDAOBeanService(),
                getWorkflowRunAttempt());
        logger.info("sampleSet.size(): {}", sampleSet.size());

        String sselProbe = getWorkflowBeanService().getAttributes().get("sselProbe");
        String siteName = getWorkflowBeanService().getAttributes().get("siteName");
        String knownVCF = getWorkflowBeanService().getAttributes().get("knownVCF");
        String referenceSequence = getWorkflowBeanService().getAttributes().get("referenceSequence");
        String icSNPIntervalList = getWorkflowBeanService().getAttributes().get("icSNPIntervalList");
        String flagstatIntervalList = getWorkflowBeanService().getAttributes().get("flagstatIntervalList");
        String depthOfCoverageIntervalList = getWorkflowBeanService().getAttributes().get("depthOfCoverageIntervalList");
        String unifiedGenotyperIntervalList = getWorkflowBeanService().getAttributes().get("unifiedGenotyperIntervalList");
        String readGroupPlatform = getWorkflowBeanService().getAttributes().get("readGroupPlatform");
        String readGroupPlatformUnit = getWorkflowBeanService().getAttributes().get("readGroupPlatformUnit");

        WorkflowRunAttempt attempt = getWorkflowRunAttempt();
        WorkflowRun workflowRun = attempt.getWorkflowRun();

        for (Sample sample : sampleSet) {

            if ("Undetermined".equals(sample.getBarcode())) {
                continue;
            }

            logger.debug(sample.toString());

            Flowcell flowcell = sample.getFlowcell();
            File outputDirectory = SequencingWorkflowUtil.createOutputDirectory(sample, workflowRun.getWorkflow());
            File tmpDirectory = new File(outputDirectory, "tmp");
            tmpDirectory.mkdirs();

            Set<Attribute> attributeSet = workflowRun.getAttributes();
            if (attributeSet != null && !attributeSet.isEmpty()) {
                Iterator<Attribute> attributeIter = attributeSet.iterator();
                while (attributeIter.hasNext()) {
                    Attribute attribute = attributeIter.next();
                    if ("GATKVariantRecalibrator.percentBadVariants".equals(attribute.getName())) {
                        percentBadVariants = attribute.getValue();
                    }
                    if ("sselProbe".equals(attribute.getName())) {
                        sselProbe = attribute.getValue();

                        // read probe version to determine interval list file(s)
                        switch (sselProbe) {
                            case "5":
                                flagstatIntervalList = "$NCGENES_RESOURCES_DIRECTORY/intervals/agilent_v5_capture_region_pm_100.shortid.interval_list";
                                depthOfCoverageIntervalList = "$NCGENES_RESOURCES_DIRECTORY/intervals/agilent_v5_capture_region_pm_100.shortid.interval_list";
                                unifiedGenotyperIntervalList = "$NCGENES_RESOURCES_DIRECTORY/intervals/agilent_v5_capture_region_pm_100.shortid.interval_list";
                                break;
                            case "5 + EGL":
                                flagstatIntervalList = "$NCGENES_RESOURCES_DIRECTORY/intervals/agilent_v5_egl_capture_region_pm_100.shortid.interval_list";
                                depthOfCoverageIntervalList = "$NCGENES_RESOURCES_DIRECTORY/intervals/agilent_v5_egl_capture_region_pm_100.shortid.interval_list";
                                unifiedGenotyperIntervalList = "$NCGENES_RESOURCES_DIRECTORY/intervals/agilent_v5_egl_capture_region_pm_100.shortid.interval_list";
                                break;
                            case "6":
                                flagstatIntervalList = "$NCGENES_RESOURCES_DIRECTORY/intervals/agilent_v6_capture_region_pm_100.shortid.interval_list";
                                depthOfCoverageIntervalList = "$NCGENES_RESOURCES_DIRECTORY/intervals/agilent_v6_capture_region_pm_100.shortid.interval_list";
                                unifiedGenotyperIntervalList = "$NCGENES_RESOURCES_DIRECTORY/intervals/agilent_v6_capture_region_pm_100.shortid.interval_list";
                                break;
                        }
                    }
                }
            }

            List<File> readPairList = SequencingWorkflowUtil.getReadPairList(sample);
            logger.info("fileList = {}", readPairList.size());

            // assumption: a dash is used as a delimiter between a participantId
            // and the external code
            // int idx = sample.getName().lastIndexOf("-");
            // String participantId = idx != -1 ? sample.getName().substring(0, idx) : sample.getName();
            // no longer expecting composite for participantId...can just use sample name
            String participantId = sample.getName();

            if (readPairList.size() == 2) {

                File r1FastqFile = readPairList.get(0);
                String r1FastqRootName = SequencingWorkflowUtil.getRootFastqName(r1FastqFile.getName());

                File r2FastqFile = readPairList.get(1);
                String r2FastqRootName = SequencingWorkflowUtil.getRootFastqName(r2FastqFile.getName());

                String rootFileName = String.format("%s_%s_L%03d", sample.getFlowcell().getName(), sample.getBarcode(),
                        sample.getLaneIndex());

                try {

                    // start site specific jobs

                    // new job
                    CondorJobBuilder builder = SequencingWorkflowJobFactory
                            .createJob(++count, WriteVCFHeaderCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
                    String flowcellProper = flowcell.getName().substring(flowcell.getName().length() - 9, flowcell.getName().length());
                    File writeVCFHeaderOut = new File(outputDirectory, String.format("%s.vcf.hdr", rootFileName));
                    builder.addArgument(WriteVCFHeaderCLI.BARCODE, sample.getBarcode())
                            .addArgument(WriteVCFHeaderCLI.RUN, flowcell.getName())
                            .addArgument(WriteVCFHeaderCLI.PARTICIPANTID, participantId)
                            .addArgument(WriteVCFHeaderCLI.STUDYNAME, sample.getStudy().getName())
                            .addArgument(WriteVCFHeaderCLI.LANE, sample.getLaneIndex().toString())
                            .addArgument(WriteVCFHeaderCLI.LABNAME, "Jonathan_Berg").addArgument(WriteVCFHeaderCLI.FLOWCELL, flowcellProper)
                            .addArgument(WriteVCFHeaderCLI.OUTPUT, writeVCFHeaderOut.getAbsolutePath());
                    CondorJob writeVCFHeaderJob = builder.build();
                    logger.info(writeVCFHeaderJob.toString());
                    graph.addVertex(writeVCFHeaderJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, FastQCCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName);
                    File fastqcR1Output = new File(outputDirectory, r1FastqRootName + ".fastqc.zip");
                    builder.addArgument(FastQCCLI.INPUT, r1FastqFile.getAbsolutePath())
                            .addArgument(FastQCCLI.OUTPUT, fastqcR1Output.getAbsolutePath())
                            .addArgument(FastQCCLI.IGNORE, IgnoreLevelType.ERROR.toString());
                    CondorJob fastQCR1Job = builder.build();
                    logger.info(fastQCR1Job.toString());
                    graph.addVertex(fastQCR1Job);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, FastQCCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName);
                    File fastqcR2Output = new File(outputDirectory, r2FastqRootName + ".fastqc.zip");
                    builder.addArgument(FastQCCLI.INPUT, r2FastqFile.getAbsolutePath())
                            .addArgument(FastQCCLI.OUTPUT, fastqcR2Output.getAbsolutePath())
                            .addArgument(FastQCCLI.IGNORE, IgnoreLevelType.ERROR.toString());
                    CondorJob fastQCR2Job = builder.build();
                    logger.info(fastQCR2Job.toString());
                    graph.addVertex(fastQCR2Job);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, BWAMEMCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName).numberOfProcessors(4);
                    File bwaMemOutFile = new File(outputDirectory, rootFileName + ".mem.sam");
                    builder.addArgument(BWAMEMCLI.FASTADB, referenceSequence).addArgument(BWAMEMCLI.FASTQ1, r1FastqFile.getAbsolutePath())
                            .addArgument(BWAMEMCLI.FASTQ2, r2FastqFile.getAbsolutePath()).addArgument(BWAMEMCLI.THREADS, "4")
                            .addArgument(BWAMEMCLI.VERBOSITY, "1").addArgument(BWAMEMCLI.MARKSHORTERSPLITHITS)
                            .addArgument(BWAMEMCLI.OUTFILE, bwaMemOutFile.getAbsolutePath());
                    CondorJob bwaMemJob = builder.build();
                    logger.info(bwaMemJob.toString());
                    graph.addVertex(bwaMemJob);

                    // new job
                    builder = SequencingWorkflowJobFactory
                            .createJob(++count, PicardAddOrReplaceReadGroupsCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
                    File fixRGOutput = new File(outputDirectory, bwaMemOutFile.getName().replace(".sam", ".fixed-rg.bam"));
                    builder.addArgument(PicardAddOrReplaceReadGroupsCLI.INPUT, bwaMemOutFile.getAbsolutePath())
                            .addArgument(PicardAddOrReplaceReadGroupsCLI.OUTPUT, fixRGOutput.getAbsolutePath())
                            .addArgument(PicardAddOrReplaceReadGroupsCLI.SORTORDER, PicardSortOrderType.COORDINATE.toString().toLowerCase())
                            .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPID,
                                    String.format("%s-%s_L%03d", flowcell.getName(), sample.getBarcode(), sample.getLaneIndex()))
                            .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPLIBRARY, participantId)
                            .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPPLATFORM, readGroupPlatform)
                            .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPPLATFORMUNIT, readGroupPlatformUnit)
                            .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPSAMPLENAME, participantId)
                            .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPCENTERNAME, "UNC");
                    CondorJob picardAddOrReplaceReadGroupsJob = builder.build();
                    logger.info(picardAddOrReplaceReadGroupsJob.toString());
                    graph.addVertex(picardAddOrReplaceReadGroupsJob);
                    graph.addEdge(bwaMemJob, picardAddOrReplaceReadGroupsJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName);
                    File picardAddOrReplaceReadGroupsIndexOut = new File(outputDirectory, fixRGOutput.getName().replace(".bam", ".bai"));
                    builder.addArgument(SAMToolsIndexCLI.INPUT, fixRGOutput.getAbsolutePath()).addArgument(SAMToolsIndexCLI.OUTPUT,
                            picardAddOrReplaceReadGroupsIndexOut.getAbsolutePath());
                    CondorJob samtoolsIndexJob = builder.build();
                    logger.info(samtoolsIndexJob.toString());
                    graph.addVertex(samtoolsIndexJob);
                    graph.addEdge(picardAddOrReplaceReadGroupsJob, samtoolsIndexJob);

                    // new job
                    builder = SequencingWorkflowJobFactory
                            .createJob(++count, PicardMarkDuplicatesCLI.class, attempt.getId(), sample.getId()).siteName(siteName);
                    File picardMarkDuplicatesMetricsFile = new File(outputDirectory,
                            fixRGOutput.getName().replace(".bam", ".deduped.metrics"));
                    File picardMarkDuplicatesOutput = new File(outputDirectory, fixRGOutput.getName().replace(".bam", ".deduped.bam"));
                    builder.addArgument(PicardMarkDuplicatesCLI.INPUT, fixRGOutput.getAbsolutePath())
                            .addArgument(PicardMarkDuplicatesCLI.METRICSFILE, picardMarkDuplicatesMetricsFile.getAbsolutePath())
                            .addArgument(PicardMarkDuplicatesCLI.OUTPUT, picardMarkDuplicatesOutput.getAbsolutePath());
                    CondorJob picardMarkDuplicatesJob = builder.build();
                    logger.info(picardMarkDuplicatesJob.toString());
                    graph.addVertex(picardMarkDuplicatesJob);
                    graph.addEdge(samtoolsIndexJob, picardMarkDuplicatesJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName);
                    File picardMarkDuplicatesIndexOut = new File(outputDirectory,
                            picardMarkDuplicatesOutput.getName().replace(".bam", ".bai"));
                    builder.addArgument(SAMToolsIndexCLI.INPUT, picardMarkDuplicatesOutput.getAbsolutePath())
                            .addArgument(SAMToolsIndexCLI.OUTPUT, picardMarkDuplicatesIndexOut.getAbsolutePath());
                    samtoolsIndexJob = builder.build();
                    logger.info(samtoolsIndexJob.toString());
                    graph.addVertex(samtoolsIndexJob);
                    graph.addEdge(picardMarkDuplicatesJob, samtoolsIndexJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsViewCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName);
                    File samtoolsViewOutput = new File(outputDirectory, picardMarkDuplicatesOutput.getName().replace(".bam", ".view.bam"));
                    builder.addArgument(SAMToolsViewCLI.INPUT, picardMarkDuplicatesOutput.getAbsolutePath())
                            .addArgument(SAMToolsViewCLI.OUTPUT, samtoolsViewOutput.getAbsolutePath())
                            .addArgument(SAMToolsViewCLI.OUTPUTALIGNMENTSWITHBITSPRESENTINFLAG, "0x100")
                            .addArgument(SAMToolsViewCLI.BAMFORMAT);
                    CondorJob samtoolsViewJob = builder.build();
                    logger.info(samtoolsViewJob.toString());
                    graph.addVertex(samtoolsViewJob);
                    graph.addEdge(samtoolsIndexJob, samtoolsViewJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName);
                    File samtoolsViewIndexOutput = new File(outputDirectory, samtoolsViewOutput.getName().replace(".bam", ".bai"));
                    builder.addArgument(SAMToolsIndexCLI.INPUT, samtoolsViewOutput.getAbsolutePath()).addArgument(SAMToolsIndexCLI.OUTPUT,
                            samtoolsViewIndexOutput.getAbsolutePath());
                    samtoolsIndexJob = builder.build();
                    logger.info(samtoolsIndexJob.toString());
                    graph.addVertex(samtoolsIndexJob);
                    graph.addEdge(samtoolsViewJob, samtoolsIndexJob);

                    // new job
                    builder = SequencingWorkflowJobFactory
                            .createJob(++count, GATKRealignerTargetCreatorCLI.class, attempt.getId(), sample.getId()).siteName(siteName)
                            .numberOfProcessors(2);
                    File realignTargetCreatorOut = new File(outputDirectory,
                            picardMarkDuplicatesOutput.getName().replace(".bam", ".targets.intervals"));
                    builder.addArgument(GATKRealignerTargetCreatorCLI.REFERENCESEQUENCE, referenceSequence)
                            .addArgument(GATKRealignerTargetCreatorCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                            .addArgument(GATKRealignerTargetCreatorCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                            .addArgument(GATKRealignerTargetCreatorCLI.KNOWN, knownVCF)
                            .addArgument(GATKRealignerTargetCreatorCLI.INPUTFILE, samtoolsViewOutput.getAbsolutePath())
                            .addArgument(GATKRealignerTargetCreatorCLI.OUT, realignTargetCreatorOut.getAbsolutePath());
                    CondorJob gatkRealignTargetCreatorJob = builder.build();
                    logger.info(gatkRealignTargetCreatorJob.toString());
                    graph.addVertex(gatkRealignTargetCreatorJob);
                    graph.addEdge(samtoolsIndexJob, gatkRealignTargetCreatorJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, GATKIndelRealignerCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName).numberOfProcessors(2);
                    File indelRealignerOut = new File(outputDirectory,
                            picardMarkDuplicatesOutput.getName().replace(".bam", ".realign.bam"));
                    builder.addArgument(GATKIndelRealignerCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                            .addArgument(GATKIndelRealignerCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString().toLowerCase())
                            .addArgument(GATKIndelRealignerCLI.REFERENCESEQUENCE, referenceSequence)
                            .addArgument(GATKIndelRealignerCLI.KNOWNALLELES, knownVCF)
                            .addArgument(GATKIndelRealignerCLI.INPUT, picardMarkDuplicatesOutput.getAbsolutePath())
                            .addArgument(GATKIndelRealignerCLI.TARGETINTERVALS, realignTargetCreatorOut.getAbsolutePath())
                            .addArgument(GATKIndelRealignerCLI.OUT, indelRealignerOut.getAbsolutePath());
                    CondorJob gatkIndelRealignerJob = builder.build();
                    logger.info(gatkIndelRealignerJob.toString());
                    graph.addVertex(gatkIndelRealignerJob);
                    graph.addEdge(gatkRealignTargetCreatorJob, gatkIndelRealignerJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, PicardFixMateCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName);
                    File picardFixMateOutput = new File(outputDirectory, indelRealignerOut.getName().replace(".bam", ".fixmate.bam"));
                    builder.addArgument(PicardFixMateCLI.SORTORDER, PicardSortOrderType.COORDINATE.toString().toLowerCase())
                            .addArgument(PicardFixMateCLI.INPUT, indelRealignerOut.getAbsolutePath())
                            .addArgument(PicardFixMateCLI.OUTPUT, picardFixMateOutput.getAbsolutePath());
                    CondorJob picardFixMateJob = builder.build();
                    logger.info(picardFixMateJob.toString());
                    graph.addVertex(picardFixMateJob);
                    graph.addEdge(gatkIndelRealignerJob, picardFixMateJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName);
                    File picardFixMateIndexOut = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam", ".bai"));
                    builder.addArgument(SAMToolsIndexCLI.INPUT, picardFixMateOutput.getAbsolutePath()).addArgument(SAMToolsIndexCLI.OUTPUT,
                            picardFixMateIndexOut.getAbsolutePath());
                    samtoolsIndexJob = builder.build();
                    logger.info(samtoolsIndexJob.toString());
                    graph.addVertex(samtoolsIndexJob);
                    graph.addEdge(picardFixMateJob, samtoolsIndexJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, GATKCountCovariatesCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName).numberOfProcessors(4);
                    File gatkCountCovariatesRecalFile = new File(outputDirectory,
                            picardFixMateOutput.getName().replace(".bam", ".bam.cov"));
                    builder.addArgument(GATKCountCovariatesCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                            .addArgument(GATKCountCovariatesCLI.REFERENCESEQUENCE, referenceSequence)
                            .addArgument(GATKCountCovariatesCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                            .addArgument(GATKCountCovariatesCLI.KNOWNSITES, knownVCF).addArgument(GATKCountCovariatesCLI.NUMTHREADS, "4")
                            .addArgument(GATKCountCovariatesCLI.STANDARDCOVS)
                            .addArgument(GATKCountCovariatesCLI.INPUTFILE, picardFixMateOutput.getAbsolutePath())
                            .addArgument(GATKCountCovariatesCLI.RECALFILE, gatkCountCovariatesRecalFile.getAbsolutePath());
                    CondorJob gatkCountCovariatesJob = builder.build();
                    logger.info(gatkCountCovariatesJob.toString());
                    graph.addVertex(gatkCountCovariatesJob);
                    graph.addEdge(samtoolsIndexJob, gatkCountCovariatesJob);

                    // new job
                    builder = SequencingWorkflowJobFactory
                            .createJob(++count, GATKTableRecalibrationCLI.class, attempt.getId(), sample.getId()).siteName(siteName)
                            .numberOfProcessors(2);
                    File gatkTableRecalibrationOut = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.bam"));
                    builder.addArgument(GATKTableRecalibrationCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                            .addArgument(GATKTableRecalibrationCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                            .addArgument(GATKTableRecalibrationCLI.REFERENCESEQUENCE, referenceSequence)
                            .addArgument(GATKTableRecalibrationCLI.INPUTFILE, picardFixMateOutput.getAbsolutePath())
                            .addArgument(GATKTableRecalibrationCLI.RECALFILE, gatkCountCovariatesRecalFile.getAbsolutePath())
                            .addArgument(GATKTableRecalibrationCLI.OUT, gatkTableRecalibrationOut.getAbsolutePath());
                    CondorJob gatkTableRecalibrationJob = builder.build();
                    logger.info(gatkTableRecalibrationJob.toString());
                    graph.addVertex(gatkTableRecalibrationJob);
                    graph.addEdge(gatkCountCovariatesJob, gatkTableRecalibrationJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName);
                    File gatkTableRecalibrationIndexOut = new File(outputDirectory,
                            gatkTableRecalibrationOut.getName().replace(".bam", ".bai"));
                    builder.addArgument(SAMToolsIndexCLI.INPUT, gatkTableRecalibrationOut.getAbsolutePath())
                            .addArgument(SAMToolsIndexCLI.OUTPUT, gatkTableRecalibrationIndexOut.getAbsolutePath());
                    samtoolsIndexJob = builder.build();
                    logger.info(samtoolsIndexJob.toString());
                    graph.addVertex(samtoolsIndexJob);
                    graph.addEdge(gatkTableRecalibrationJob, samtoolsIndexJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsFlagstatCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName);
                    File samtoolsFlagstatOut = new File(outputDirectory,
                            gatkTableRecalibrationOut.getName().replace(".bam", ".samtools.flagstat"));
                    builder.addArgument(SAMToolsFlagstatCLI.INPUT, gatkTableRecalibrationOut.getAbsolutePath())
                            .addArgument(SAMToolsFlagstatCLI.OUTPUT, samtoolsFlagstatOut.getAbsolutePath());
                    CondorJob samtoolsFlagstatJob = builder.build();
                    logger.info(samtoolsFlagstatJob.toString());
                    graph.addVertex(samtoolsFlagstatJob);
                    graph.addEdge(samtoolsIndexJob, samtoolsFlagstatJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, GATKFlagStatCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName).numberOfProcessors(2);
                    File gatkFlagstatOut = new File(outputDirectory, gatkTableRecalibrationOut.getName().replace(".bam", ".gatk.flagstat"));
                    builder.addArgument(GATKFlagStatCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                            .addArgument(GATKFlagStatCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                            .addArgument(GATKFlagStatCLI.INPUTFILE, gatkTableRecalibrationOut.getAbsolutePath())
                            .addArgument(GATKFlagStatCLI.REFERENCESEQUENCE, referenceSequence)
                            .addArgument(GATKFlagStatCLI.INTERVALS, flagstatIntervalList)
                            .addArgument(GATKFlagStatCLI.OUT, gatkFlagstatOut.getAbsolutePath());
                    CondorJob gatkFlagstatJob = builder.build();
                    logger.info(gatkFlagstatJob.toString());
                    graph.addVertex(gatkFlagstatJob);
                    graph.addEdge(samtoolsIndexJob, gatkFlagstatJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, GATKDepthOfCoverageCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName).initialDirectory(outputDirectory.getAbsolutePath()).numberOfProcessors(2);
                    builder.addArgument(GATKDepthOfCoverageCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                            .addArgument(GATKDepthOfCoverageCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                            .addArgument(GATKDepthOfCoverageCLI.REFERENCESEQUENCE, referenceSequence)
                            .addArgument(GATKDepthOfCoverageCLI.VALIDATIONSTRICTNESS, "LENIENT")
                            .addArgument(GATKDepthOfCoverageCLI.OMITDEPTHOUTPUTATEACHBASE)
                            .addArgument(GATKDepthOfCoverageCLI.INPUTFILE, gatkTableRecalibrationOut.getAbsolutePath())
                            .addArgument(GATKDepthOfCoverageCLI.INTERVALS, depthOfCoverageIntervalList).addArgument(
                                    GATKDepthOfCoverageCLI.OUTPUTPREFIX, gatkTableRecalibrationOut.getName().replace(".bam", ".coverage"));
                    CondorJob gatkDepthOfCoverageJob = builder.build();
                    logger.info(gatkDepthOfCoverageJob.toString());
                    graph.addVertex(gatkDepthOfCoverageJob);
                    graph.addEdge(samtoolsFlagstatJob, gatkDepthOfCoverageJob);
                    graph.addEdge(gatkFlagstatJob, gatkDepthOfCoverageJob);

                    // new job
                    builder = SequencingWorkflowJobFactory
                            .createJob(++count, GATKUnifiedGenotyperCLI.class, attempt.getId(), sample.getId()).siteName(siteName)
                            .numberOfProcessors(4);
                    File gatkUnifiedGenotyperOut = new File(outputDirectory, gatkTableRecalibrationOut.getName().replace(".bam", ".vcf"));
                    File gatkUnifiedGenotyperMetrics = new File(outputDirectory,
                            gatkTableRecalibrationOut.getName().replace(".bam", ".metrics"));
                    builder.addArgument(GATKUnifiedGenotyperCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                            .addArgument(GATKUnifiedGenotyperCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                            .addArgument(GATKUnifiedGenotyperCLI.REFERENCESEQUENCE, referenceSequence)
                            .addArgument(GATKUnifiedGenotyperCLI.DBSNP, knownVCF).addArgument(GATKUnifiedGenotyperCLI.STANDCALLCONF, "30")
                            .addArgument(GATKUnifiedGenotyperCLI.STANDEMITCONF, "0")
                            .addArgument(GATKUnifiedGenotyperCLI.GENOTYPELIKELIHOODSMODEL, "BOTH")
                            .addArgument(GATKUnifiedGenotyperCLI.INPUTFILE, gatkTableRecalibrationOut.getAbsolutePath())
                            .addArgument(GATKUnifiedGenotyperCLI.NUMTHREADS, "4")
                            .addArgument(GATKUnifiedGenotyperCLI.OUT, gatkUnifiedGenotyperOut.getAbsolutePath())
                            .addArgument(GATKUnifiedGenotyperCLI.INTERVALS, unifiedGenotyperIntervalList)
                            .addArgument(GATKUnifiedGenotyperCLI.OUTPUTMODE, "EMIT_ALL_SITES")
                            .addArgument(GATKUnifiedGenotyperCLI.METRICS, gatkUnifiedGenotyperMetrics.getAbsolutePath())
                            .addArgument(GATKUnifiedGenotyperCLI.DOWNSAMPLETOCOVERAGE, "250")
                            .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "AlleleBalance")
                            .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "DepthOfCoverage")
                            .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "HomopolymerRun")
                            .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "MappingQualityZero")
                            .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "QualByDepth")
                            .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "RMSMappingQuality")
                            .addArgument(GATKUnifiedGenotyperCLI.ANNOTATION, "HaplotypeScore");
                    CondorJob gatkUnifiedGenotyperJob = builder.build();
                    logger.info(gatkUnifiedGenotyperJob.toString());
                    graph.addVertex(gatkUnifiedGenotyperJob);
                    graph.addEdge(gatkDepthOfCoverageJob, gatkUnifiedGenotyperJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, FilterVariantCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName).numberOfProcessors(2);
                    File filterVariant1Output = new File(outputDirectory,
                            gatkTableRecalibrationOut.getName().replace(".bam", ".variant.vcf"));
                    builder.addArgument(FilterVariantCLI.INTERVALLIST, icSNPIntervalList).addArgument(FilterVariantCLI.WITHMISSING)
                            .addArgument(FilterVariantCLI.INPUT, gatkUnifiedGenotyperOut.getAbsolutePath())
                            .addArgument(FilterVariantCLI.OUTPUT, filterVariant1Output.getAbsolutePath());
                    CondorJob filterVariant1Job = builder.build();
                    logger.info(filterVariant1Job.toString());
                    graph.addVertex(filterVariant1Job);
                    graph.addEdge(gatkUnifiedGenotyperJob, filterVariant1Job);

                    // new job
                    builder = SequencingWorkflowJobFactory
                            .createJob(++count, GATKVariantRecalibratorCLI.class, attempt.getId(), sample.getId()).siteName(siteName)
                            .numberOfProcessors(2);
                    File gatkVariantRecalibratorRecalFile = new File(outputDirectory,
                            filterVariant1Output.getName().replace(".vcf", ".recal"));
                    File gatkVariantRecalibratorTranchesFile = new File(outputDirectory,
                            filterVariant1Output.getName().replace(".vcf", ".tranches"));
                    File gatkVariantRecalibratorRScriptFile = new File(outputDirectory,
                            filterVariant1Output.getName().replace(".vcf", ".plots.R"));
                    builder.addArgument(GATKVariantRecalibratorCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                            .addArgument(GATKVariantRecalibratorCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                            .addArgument(GATKVariantRecalibratorCLI.REFERENCESEQUENCE, referenceSequence)
                            .addArgument(GATKVariantRecalibratorCLI.MAXGAUSSIANS, "4")
                            .addArgument(GATKVariantRecalibratorCLI.INPUT, filterVariant1Output.getAbsolutePath())
                            .addArgument(GATKVariantRecalibratorCLI.MODE, "SNP")
                            .addArgument(GATKVariantRecalibratorCLI.RECALFILE, gatkVariantRecalibratorRecalFile.getAbsolutePath())
                            .addArgument(GATKVariantRecalibratorCLI.TRANCHESFILE, gatkVariantRecalibratorTranchesFile.getAbsolutePath())
                            .addArgument(GATKVariantRecalibratorCLI.RESOURCE,
                                    ":hapmap,known=false,training=true,truth=true,prior=15.0^$NCGENES_RESOURCES_DIRECTORY/gatk/bundle/1.2/b37/hapmap_3.3.b37.sites.renci.shortid.vcf")
                            .addArgument(GATKVariantRecalibratorCLI.RESOURCE,
                                    ":omni,known=false,training=true,truth=false,prior=12.0^$NCGENES_RESOURCES_DIRECTORY/gatk/bundle/1.2/b37/1000G_omni2.5.b37.sites.renci.shortid.vcf")
                            .addArgument(GATKVariantRecalibratorCLI.RESOURCE,
                                    ":dbsnp,known=true,training=false,truth=false,prior=8.0^$NCGENES_RESOURCES_DIRECTORY/gatk/bundle/1.2/b37/dbsnp_132.b37.renci.shortid.vcf")
                            .addArgument(GATKVariantRecalibratorCLI.USEANNOTATION, "QD")
                            .addArgument(GATKVariantRecalibratorCLI.USEANNOTATION, "HaplotypeScore")
                            .addArgument(GATKVariantRecalibratorCLI.USEANNOTATION, "MQRankSum")
                            .addArgument(GATKVariantRecalibratorCLI.USEANNOTATION, "ReadPosRankSum")
                            .addArgument(GATKVariantRecalibratorCLI.USEANNOTATION, "MQ")
                            .addArgument(GATKVariantRecalibratorCLI.USEANNOTATION, "FS")
                            .addArgument(GATKVariantRecalibratorCLI.RSCRIPTFILE, gatkVariantRecalibratorRScriptFile.getAbsolutePath());
                    if (StringUtils.isNotEmpty(percentBadVariants)) {
                        builder.addArgument(GATKVariantRecalibratorCLI.PERCENTBADVARIANTS, percentBadVariants);
                    }
                    CondorJob gatkVariantRecalibratorJob = builder.build();
                    logger.info(gatkVariantRecalibratorJob.toString());
                    graph.addVertex(gatkVariantRecalibratorJob);
                    graph.addEdge(filterVariant1Job, gatkVariantRecalibratorJob);

                    // new job
                    builder = SequencingWorkflowJobFactory
                            .createJob(++count, GATKApplyRecalibrationCLI.class, attempt.getId(), sample.getId()).siteName(siteName)
                            .numberOfProcessors(2);
                    File gatkApplyRecalibrationOut = new File(outputDirectory,
                            filterVariant1Output.getName().replace(".vcf", ".recalibrated.filtered.vcf"));
                    builder.addArgument(GATKApplyRecalibrationCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                            .addArgument(GATKApplyRecalibrationCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                            .addArgument(GATKApplyRecalibrationCLI.REFERENCESEQUENCE, referenceSequence)
                            .addArgument(GATKApplyRecalibrationCLI.INPUT, filterVariant1Output.getAbsolutePath())
                            .addArgument(GATKApplyRecalibrationCLI.RECALFILE, gatkVariantRecalibratorRecalFile.getAbsolutePath())
                            .addArgument(GATKApplyRecalibrationCLI.TRANCHESFILE, gatkVariantRecalibratorTranchesFile.getAbsolutePath())
                            .addArgument(GATKApplyRecalibrationCLI.OUT, gatkApplyRecalibrationOut.getAbsolutePath())
                            .addArgument(GATKApplyRecalibrationCLI.TSFILTERLEVEL, "99.0");
                    CondorJob gatkApplyRecalibrationJob = builder.build();
                    logger.info(gatkApplyRecalibrationJob.toString());
                    graph.addVertex(gatkApplyRecalibrationJob);
                    graph.addEdge(gatkVariantRecalibratorJob, gatkApplyRecalibrationJob);

                    // new job
                    builder = SequencingWorkflowJobFactory.createJob(++count, FilterVariantCLI.class, attempt.getId(), sample.getId())
                            .siteName(siteName).numberOfProcessors(2);
                    File filterVariant2Output = new File(outputDirectory, filterVariant1Output.getName().replace(".vcf", ".ic_snps.vcf"));
                    builder.addArgument(FilterVariantCLI.INTERVALLIST, icSNPIntervalList)
                            .addArgument(FilterVariantCLI.INPUT, gatkApplyRecalibrationOut.getAbsolutePath())
                            .addArgument(FilterVariantCLI.OUTPUT, filterVariant2Output.getAbsolutePath());
                    CondorJob filterVariant2Job = builder.build();
                    logger.info(filterVariant2Job.toString());
                    graph.addVertex(filterVariant2Job);
                    graph.addEdge(gatkApplyRecalibrationJob, filterVariant2Job);

                } catch (Exception e) {
                    throw new WorkflowException(e);
                }

            }

        }

        return graph;
    }

    @Override
    public void postRun() throws WorkflowException {
        logger.info("ENTERING postRun()");

        ExecutorService es = Executors.newSingleThreadExecutor();

        try {

            RegisterToIRODSRunnable registerNCGenesToIRODSRunnable = new RegisterToIRODSRunnable(
                    getWorkflowBeanService().getMaPSeqDAOBeanService(), getWorkflowRunAttempt());
            es.submit(registerNCGenesToIRODSRunnable);

            SaveFlagstatAttributesRunnable saveFlagstatAttributeRunnable = new SaveFlagstatAttributesRunnable(
                    getWorkflowBeanService().getMaPSeqDAOBeanService(), getWorkflowRunAttempt());
            es.submit(saveFlagstatAttributeRunnable);

            SaveMarkDuplicatesAttributesRunnable saveMarkDuplicatesAttributesRunnable = new SaveMarkDuplicatesAttributesRunnable(
                    getWorkflowBeanService().getMaPSeqDAOBeanService(), getWorkflowRunAttempt());
            es.submit(saveMarkDuplicatesAttributesRunnable);

            SaveDepthOfCoverageAttributesRunnable saveDepthOfCoverageAttributesRunnable = new SaveDepthOfCoverageAttributesRunnable(
                    getWorkflowBeanService().getMaPSeqDAOBeanService(), getWorkflowRunAttempt());
            es.submit(saveDepthOfCoverageAttributesRunnable);

            es.shutdown();
            es.awaitTermination(1L, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
}
