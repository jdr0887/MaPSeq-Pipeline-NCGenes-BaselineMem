package edu.unc.mapseq.workflow.ncgenes;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.lang.StringUtils;
import org.jgrapht.DirectedGraph;
import org.jgrapht.ext.VertexNameProvider;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.junit.Test;
import org.renci.jlrm.condor.CondorJob;
import org.renci.jlrm.condor.CondorJobBuilder;
import org.renci.jlrm.condor.CondorJobEdge;
import org.renci.jlrm.condor.ext.CondorDOTExporter;

import edu.unc.mapseq.dao.model.Flowcell;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.dao.model.Study;
import edu.unc.mapseq.dao.model.WorkflowSystemType;
import edu.unc.mapseq.module.core.ZipCLI;
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
import edu.unc.mapseq.module.sequencing.picard.PicardSortSAMCLI;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsFlagstatCLI;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsIndexCLI;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsViewCLI;
import edu.unc.mapseq.workflow.WorkflowException;
import edu.unc.mapseq.workflow.core.exporter.CLIScriptExporter;
import edu.unc.mapseq.workflow.sequencing.SequencingWorkflowJobFactory;
import edu.unc.mapseq.workflow.sequencing.SequencingWorkflowUtil;

public class NCGenesWorkflowTest {

    @Test
    public void createDot() {

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(CondorJobEdge.class);

        int count = 0;

        // new job
        CondorJob writeVCFHeaderJob = new CondorJobBuilder().name(String.format("%s_%d", WriteVCFHeaderCLI.class.getSimpleName(), ++count))
                .build();
        graph.addVertex(writeVCFHeaderJob);

        // new job
        CondorJob fastQCR1Job = new CondorJobBuilder().name(String.format("%s_%d", FastQCCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(fastQCR1Job);

        // new job
        CondorJob fastQCR2Job = new CondorJobBuilder().name(String.format("%s_%d", FastQCCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(fastQCR2Job);

        // new job
        CondorJob bwaMemJob = new CondorJobBuilder()
                .name(String.format("%s_%d", BWAMEMCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(bwaMemJob);

        // new job
        CondorJob picardAddOrReplaceReadGroupsJob = new CondorJobBuilder()
                .name(String.format("%s_%d", PicardAddOrReplaceReadGroupsCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(picardAddOrReplaceReadGroupsJob);
        graph.addEdge(bwaMemJob, picardAddOrReplaceReadGroupsJob);

        // new job
        CondorJob samtoolsIndexJob = new CondorJobBuilder().name(String.format("%s_%d", SAMToolsIndexCLI.class.getSimpleName(), ++count))
                .build();
        graph.addVertex(samtoolsIndexJob);
        graph.addEdge(picardAddOrReplaceReadGroupsJob, samtoolsIndexJob);

        // new job
        CondorJob picardMarkDuplicatesJob = new CondorJobBuilder()
                .name(String.format("%s_%d", PicardMarkDuplicatesCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(picardMarkDuplicatesJob);
        graph.addEdge(samtoolsIndexJob, picardMarkDuplicatesJob);

        // new job
        CondorJob picardMarkDuplicatesIndexJob = new CondorJobBuilder()
                .name(String.format("%s_%d", SAMToolsIndexCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(picardMarkDuplicatesIndexJob);
        graph.addEdge(picardMarkDuplicatesJob, picardMarkDuplicatesIndexJob);

        // new job
        CondorJob gatkRealignTargetCreatorJob = new CondorJobBuilder()
                .name(String.format("%s_%d", GATKRealignerTargetCreatorCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(gatkRealignTargetCreatorJob);
        graph.addEdge(picardMarkDuplicatesIndexJob, gatkRealignTargetCreatorJob);

        // new job
        CondorJob gatkIndelRealignerJob = new CondorJobBuilder()
                .name(String.format("%s_%d", GATKIndelRealignerCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(gatkIndelRealignerJob);
        graph.addEdge(gatkRealignTargetCreatorJob, gatkIndelRealignerJob);

        // new job
        CondorJob picardFixMateJob = new CondorJobBuilder().name(String.format("%s_%d", PicardFixMateCLI.class.getSimpleName(), ++count))
                .build();
        graph.addVertex(picardFixMateJob);
        graph.addEdge(gatkIndelRealignerJob, picardFixMateJob);

        // new job
        CondorJob picardFixMateIndexJob = new CondorJobBuilder()
                .name(String.format("%s_%d", SAMToolsIndexCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(picardFixMateIndexJob);
        graph.addEdge(picardFixMateJob, picardFixMateIndexJob);

        // new job
        CondorJob gatkCountCovariatesJob = new CondorJobBuilder()
                .name(String.format("%s_%d", GATKCountCovariatesCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(gatkCountCovariatesJob);
        graph.addEdge(picardFixMateIndexJob, gatkCountCovariatesJob);

        // new job
        CondorJob gatkTableRecalibrationJob = new CondorJobBuilder()
                .name(String.format("%s_%d", GATKTableRecalibrationCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(gatkTableRecalibrationJob);
        graph.addEdge(gatkCountCovariatesJob, gatkTableRecalibrationJob);

        // new job
        CondorJob gatkTableRecalibrationIndexJob = new CondorJobBuilder()
                .name(String.format("%s_%d", SAMToolsIndexCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(gatkTableRecalibrationIndexJob);
        graph.addEdge(gatkTableRecalibrationJob, gatkTableRecalibrationIndexJob);

        // new job
        CondorJob samtoolsFlagstatJob = new CondorJobBuilder()
                .name(String.format("%s_%d", SAMToolsFlagstatCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(samtoolsFlagstatJob);
        graph.addEdge(gatkTableRecalibrationIndexJob, samtoolsFlagstatJob);

        // new job
        CondorJob gatkDepthOfCoverageJob = new CondorJobBuilder()
                .name(String.format("%s_%d", GATKDepthOfCoverageCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(gatkDepthOfCoverageJob);
        graph.addEdge(samtoolsFlagstatJob, gatkDepthOfCoverageJob);

        // new job
        CondorJob gatkGeneDepthOfCoverageJob = new CondorJobBuilder()
                .name(String.format("%s_%d", GATKDepthOfCoverageCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(gatkGeneDepthOfCoverageJob);
        graph.addEdge(samtoolsFlagstatJob, gatkGeneDepthOfCoverageJob);

        // new job
        CondorJob samtoolsViewJob = new CondorJobBuilder().name(String.format("%s_%d", SAMToolsViewCLI.class.getSimpleName(), ++count))
                .build();
        graph.addVertex(samtoolsViewJob);
        graph.addEdge(samtoolsFlagstatJob, samtoolsViewJob);

        // new job
        CondorJob picardSortSAMJob = new CondorJobBuilder().name(String.format("%s_%d", PicardSortSAMCLI.class.getSimpleName(), ++count))
                .build();
        graph.addVertex(picardSortSAMJob);
        graph.addEdge(samtoolsViewJob, picardSortSAMJob);

        // new job
        CondorJob picardSortSAMIndexJob = new CondorJobBuilder()
                .name(String.format("%s_%d", SAMToolsIndexCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(picardSortSAMIndexJob);
        graph.addEdge(picardSortSAMJob, picardSortSAMIndexJob);

        // new job
        CondorJob zipJob = new CondorJobBuilder().name(String.format("%s_%d", ZipCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(zipJob);
        graph.addEdge(picardSortSAMIndexJob, zipJob);

        // new job
        CondorJob gatkUnifiedGenotyperJob = new CondorJobBuilder()
                .name(String.format("%s_%d", GATKUnifiedGenotyperCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(gatkUnifiedGenotyperJob);
        graph.addEdge(gatkDepthOfCoverageJob, gatkUnifiedGenotyperJob);
        graph.addEdge(gatkGeneDepthOfCoverageJob, gatkUnifiedGenotyperJob);
        graph.addEdge(zipJob, gatkUnifiedGenotyperJob);

        // new job
        CondorJob filterVariant1Job = new CondorJobBuilder().name(String.format("%s_%d", FilterVariantCLI.class.getSimpleName(), ++count))
                .build();
        graph.addVertex(filterVariant1Job);
        graph.addEdge(gatkUnifiedGenotyperJob, filterVariant1Job);

        // new job
        CondorJob gatkVariantRecalibratorJob = new CondorJobBuilder()
                .name(String.format("%s_%d", GATKVariantRecalibratorCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(gatkVariantRecalibratorJob);
        graph.addEdge(filterVariant1Job, gatkVariantRecalibratorJob);

        // new job
        CondorJob gatkApplyRecalibrationJob = new CondorJobBuilder()
                .name(String.format("%s_%d", GATKApplyRecalibrationCLI.class.getSimpleName(), ++count)).build();
        graph.addVertex(gatkApplyRecalibrationJob);
        graph.addEdge(gatkVariantRecalibratorJob, gatkApplyRecalibrationJob);

        // new job
        CondorJob filterVariant2Job = new CondorJobBuilder().name(String.format("%s_%d", FilterVariantCLI.class.getSimpleName(), ++count))
                .build();
        graph.addVertex(filterVariant2Job);
        graph.addEdge(gatkApplyRecalibrationJob, filterVariant2Job);

        // new job
        CondorJob filterVariant3Job = new CondorJobBuilder().name(String.format("%s_%d", FilterVariantCLI.class.getSimpleName(), ++count))
                .build();
        graph.addVertex(filterVariant3Job);
        graph.addEdge(gatkApplyRecalibrationJob, filterVariant3Job);

        VertexNameProvider<CondorJob> vnpId = new VertexNameProvider<CondorJob>() {
            @Override
            public String getVertexName(CondorJob job) {
                return job.getName();
            }
        };

        VertexNameProvider<CondorJob> vnpLabel = new VertexNameProvider<CondorJob>() {
            @Override
            public String getVertexName(CondorJob job) {
                return job.getName();
            }
        };

        CondorDOTExporter<CondorJob, CondorJobEdge> dotExporter = new CondorDOTExporter<CondorJob, CondorJobEdge>(vnpId, vnpLabel, null,
                null, null, null);
        File srcSiteResourcesImagesDir = new File("../src/site/resources/images");
        if (!srcSiteResourcesImagesDir.exists()) {
            srcSiteResourcesImagesDir.mkdirs();
        }
        File dotFile = new File(srcSiteResourcesImagesDir, "workflow.dag.dot");
        try {
            FileWriter fw = new FileWriter(dotFile);
            dotExporter.export(fw, graph);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    @Test
    public void createCLI() throws WorkflowException {

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(CondorJobEdge.class);

        File sampleDirectory = new File("/tmp/ncgenes/140205_UNC17-D00216_0141_BH8F37ADXX/L001_ATCACG");
        sampleDirectory.mkdirs();

        File r1FastqFile = new File("/tmp/ncgenes/140205_UNC17-D00216_0141_BH8F37ADXX/L001_ATCACG/CASAVA",
                "140205_UNC17-D00216_0141_BH8F37ADXX_ATCACG_L001_R1.fastq.gz");
        String r1FastqRootName = SequencingWorkflowUtil.getRootFastqName(r1FastqFile.getName());

        File r2FastqFile = new File("/tmp/ncgenes/140205_UNC17-D00216_0141_BH8F37ADXX/L001_ATCACG/CASAVA",
                "140205_UNC17-D00216_0141_BH8F37ADXX_ATCACG_L001_R2.fastq.gz");
        String r2FastqRootName = SequencingWorkflowUtil.getRootFastqName(r2FastqFile.getName());

        String fastqLaneRootName = StringUtils.removeEnd(r2FastqRootName, "_R2");

        String siteName = "Hatteras";
        // start site specific jobs

        Study study = new Study("NCGenes");

        Flowcell flowcell = new Flowcell("140205_UNC17-D00216_0141_BH8F37ADXX");
        Sample sample = new Sample("NCG_00216-LT_1");
        sample.setLaneIndex(1);
        sample.setBarcode("ATCACG");
        sample.setStudy(study);

        String outputDir = System.getenv("MAPSEQ_OUTPUT_DIRECTORY");
        File systemDirectory = new File(outputDir, WorkflowSystemType.PRODUCTION.getValue());
        File studyDirectory = new File(systemDirectory, sample.getStudy().getName());
        File analysisDirectory = new File(studyDirectory, "analysis");
        File flowcellDirectory = new File(analysisDirectory, sample.getFlowcell().getName());
        File sampleOutputDir = new File(flowcellDirectory, String.format("L%03d_%s", sample.getLaneIndex(), sample.getBarcode()));
        File workflowDirectory = new File(sampleOutputDir, "NCGenesBaselineMem");
        if (!workflowDirectory.exists()) {
            workflowDirectory.mkdirs();
        }

        String participantId = "NCG_00216";

        String referenceSequence = "$NCGENES_REFERENCES_DIRECTORY/BUILD.37.1/bwa061sam0118/BUILD.37.1.sorted.shortid.fa";
        String knownVCF = "$NCGENES_RESOURCES_DIRECTORY/gatk/bundle/1.2/b37/dbsnp_132.b37.renci.shortid.vcf";
        String sselProbe = "4";
        String variantHeader = "##targetInterval=/proj/renci/sequence_analysis/resources/intervals/agilent_v4_capture_region_pm_50.txt";
        String flagstatIntervalList = "$NCGENES_RESOURCES_DIRECTORY/intervals/agilent_v4_capture_region_pm_75.shortid.interval_list";
        String depthOfCoverageIntervalList = "$NCGENES_RESOURCES_DIRECTORY/intervals/agilent_v4_capture_region_pm_75.shortid.interval_list";
        String unifiedGenotyperIntervalList = "$NCGENES_RESOURCES_DIRECTORY/intervals/agilent_v4_capture_region_pm_75.shortid.interval_list";
        String icSNPIntervalList = "$NCGENES_RESOURCES_DIRECTORY/ic_snp_v2.list";
        String readGroupPlatform = "ILLUMINA";
        String readGroupPlatformUnit = "Illumina HiSeq 2000";
        String percentBadVariants = "0.05";

        int count = 0;

        // new job
        CondorJobBuilder builder = SequencingWorkflowJobFactory.createJob(++count, WriteVCFHeaderCLI.class, null);
        String flowcellProper = flowcell.getName().substring(flowcell.getName().length() - 9, flowcell.getName().length());
        File writeVCFHeaderOut = new File(workflowDirectory, String.format("%s.vcf.hdr", fastqLaneRootName));
        builder.addArgument(WriteVCFHeaderCLI.VALIDATE, Boolean.FALSE).addArgument(WriteVCFHeaderCLI.DRYRUN)
                .addArgument(WriteVCFHeaderCLI.BARCODE, sample.getBarcode()).addArgument(WriteVCFHeaderCLI.RUN, flowcell.getName())
                .addArgument(WriteVCFHeaderCLI.PARTICIPANTID, participantId)
                .addArgument(WriteVCFHeaderCLI.STUDYNAME, sample.getStudy().getName())
                .addArgument(WriteVCFHeaderCLI.LANE, sample.getLaneIndex().toString())
                .addArgument(WriteVCFHeaderCLI.LABNAME, "Jonathan_Berg").addArgument(WriteVCFHeaderCLI.FLOWCELL, flowcellProper)
                .addArgument(WriteVCFHeaderCLI.OUTPUT, writeVCFHeaderOut.getAbsolutePath());
        CondorJob writeVCFHeaderJob = builder.build();
        graph.addVertex(writeVCFHeaderJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, FastQCCLI.class, null).siteName(siteName);
        File fastqcR1Output = new File(workflowDirectory, r1FastqRootName + ".fastqc.zip");
        builder.addArgument(FastQCCLI.VALIDATE, Boolean.FALSE).addArgument(FastQCCLI.DRYRUN)
                .addArgument(FastQCCLI.INPUT, r1FastqFile.getAbsolutePath()).addArgument(FastQCCLI.OUTPUT, fastqcR1Output.getAbsolutePath())
                .addArgument(FastQCCLI.IGNORE, IgnoreLevelType.ERROR.toString());

        CondorJob fastQCR1Job = builder.build();
        graph.addVertex(fastQCR1Job);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, FastQCCLI.class, null).siteName(siteName);
        File fastqcR2Output = new File(workflowDirectory, r2FastqRootName + ".fastqc.zip");
        builder.addArgument(FastQCCLI.VALIDATE, Boolean.FALSE).addArgument(FastQCCLI.DRYRUN)
                .addArgument(FastQCCLI.INPUT, r2FastqFile.getAbsolutePath()).addArgument(FastQCCLI.OUTPUT, fastqcR2Output.getAbsolutePath())
                .addArgument(FastQCCLI.IGNORE, IgnoreLevelType.ERROR.toString());
        CondorJob fastQCR2Job = builder.build();
        graph.addVertex(fastQCR2Job);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, BWAMEMCLI.class, null)
                .siteName(siteName).numberOfProcessors(4);
        File bwaMemOutFile = new File(workflowDirectory, fastqLaneRootName + ".mem.sam");
        builder.addArgument(BWAMEMCLI.FASTADB, referenceSequence)
                .addArgument(BWAMEMCLI.FASTQ1, r1FastqFile.getAbsolutePath())
                .addArgument(BWAMEMCLI.FASTQ2, r2FastqFile.getAbsolutePath())
                .addArgument(BWAMEMCLI.THREADS, "4")
                .addArgument(BWAMEMCLI.VERBOSITY, "1")
                .addArgument(BWAMEMCLI.MARKSHORTERSPLITHITS)
                .addArgument(BWAMEMCLI.OUTFILE, bwaMemOutFile.getAbsolutePath());
        CondorJob bwaMemJob = builder.build();
        graph.addVertex(bwaMemJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, PicardAddOrReplaceReadGroupsCLI.class, null).siteName(siteName);
        File fixRGOutput = new File(workflowDirectory, bwaMemOutFile.getName().replace(".sam", ".fixed-rg.bam"));
        builder.addArgument(PicardAddOrReplaceReadGroupsCLI.VALIDATE, Boolean.FALSE).addArgument(PicardAddOrReplaceReadGroupsCLI.DRYRUN)
                .addArgument(PicardAddOrReplaceReadGroupsCLI.INPUT, bwaMemOutFile.getAbsolutePath())
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
        graph.addVertex(picardAddOrReplaceReadGroupsJob);
        graph.addEdge(bwaMemJob, picardAddOrReplaceReadGroupsJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, null).siteName(siteName);
        File picardAddOrReplaceReadGroupsIndexOut = new File(workflowDirectory, fixRGOutput.getName().replace(".bam", ".bai"));
        builder.addArgument(SAMToolsIndexCLI.VALIDATE, Boolean.FALSE).addArgument(SAMToolsIndexCLI.DRYRUN)
                .addArgument(SAMToolsIndexCLI.INPUT, fixRGOutput.getAbsolutePath())
                .addArgument(SAMToolsIndexCLI.OUTPUT, picardAddOrReplaceReadGroupsIndexOut.getAbsolutePath());
        CondorJob samtoolsIndexJob = builder.build();
        graph.addVertex(samtoolsIndexJob);
        graph.addEdge(picardAddOrReplaceReadGroupsJob, samtoolsIndexJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, PicardMarkDuplicatesCLI.class, null).siteName(siteName);
        File picardMarkDuplicatesMetricsFile = new File(workflowDirectory, fixRGOutput.getName().replace(".bam", ".deduped.metrics"));
        File picardMarkDuplicatesOutput = new File(workflowDirectory, fixRGOutput.getName().replace(".bam", ".deduped.bam"));
        builder.addArgument(PicardMarkDuplicatesCLI.VALIDATE, Boolean.FALSE).addArgument(PicardMarkDuplicatesCLI.DRYRUN)
                .addArgument(PicardMarkDuplicatesCLI.INPUT, fixRGOutput.getAbsolutePath())
                .addArgument(PicardMarkDuplicatesCLI.METRICSFILE, picardMarkDuplicatesMetricsFile.getAbsolutePath())
                .addArgument(PicardMarkDuplicatesCLI.OUTPUT, picardMarkDuplicatesOutput.getAbsolutePath());
        CondorJob picardMarkDuplicatesJob = builder.build();
        graph.addVertex(picardMarkDuplicatesJob);
        graph.addEdge(samtoolsIndexJob, picardMarkDuplicatesJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, null).siteName(siteName);
        File picardMarkDuplicatesIndexOut = new File(workflowDirectory, picardMarkDuplicatesOutput.getName().replace(".bam", ".bai"));
        builder.addArgument(SAMToolsIndexCLI.VALIDATE, Boolean.FALSE).addArgument(SAMToolsIndexCLI.DRYRUN)
                .addArgument(SAMToolsIndexCLI.INPUT, picardMarkDuplicatesOutput.getAbsolutePath())
                .addArgument(SAMToolsIndexCLI.OUTPUT, picardMarkDuplicatesIndexOut.getAbsolutePath());
        samtoolsIndexJob = builder.build();
        graph.addVertex(samtoolsIndexJob);
        graph.addEdge(picardMarkDuplicatesJob, samtoolsIndexJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, GATKRealignerTargetCreatorCLI.class, null).siteName(siteName)
                .numberOfProcessors(2);
        File realignTargetCreatorOut = new File(workflowDirectory,
                picardMarkDuplicatesOutput.getName().replace(".bam", ".targets.intervals"));
        builder.addArgument(GATKRealignerTargetCreatorCLI.VALIDATE, Boolean.FALSE).addArgument(GATKRealignerTargetCreatorCLI.DRYRUN)
                .addArgument(GATKRealignerTargetCreatorCLI.REFERENCESEQUENCE, referenceSequence)
                .addArgument(GATKRealignerTargetCreatorCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                .addArgument(GATKRealignerTargetCreatorCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                .addArgument(GATKRealignerTargetCreatorCLI.KNOWN, knownVCF)
                .addArgument(GATKRealignerTargetCreatorCLI.INPUTFILE, picardMarkDuplicatesOutput.getAbsolutePath())
                .addArgument(GATKRealignerTargetCreatorCLI.OUT, realignTargetCreatorOut.getAbsolutePath());
        CondorJob gatkRealignTargetCreatorJob = builder.build();
        graph.addVertex(gatkRealignTargetCreatorJob);
        graph.addEdge(samtoolsIndexJob, gatkRealignTargetCreatorJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, GATKIndelRealignerCLI.class, null).siteName(siteName)
                .numberOfProcessors(2);
        File indelRealignerOut = new File(workflowDirectory, picardMarkDuplicatesOutput.getName().replace(".bam", ".realign.bam"));
        builder.addArgument(GATKIndelRealignerCLI.VALIDATE, Boolean.FALSE).addArgument(GATKIndelRealignerCLI.DRYRUN)
                .addArgument(GATKIndelRealignerCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                .addArgument(GATKIndelRealignerCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString().toLowerCase())
                .addArgument(GATKIndelRealignerCLI.REFERENCESEQUENCE, referenceSequence)
                .addArgument(GATKIndelRealignerCLI.KNOWNALLELES, knownVCF)
                .addArgument(GATKIndelRealignerCLI.INPUT, picardMarkDuplicatesOutput.getAbsolutePath())
                .addArgument(GATKIndelRealignerCLI.TARGETINTERVALS, realignTargetCreatorOut.getAbsolutePath())
                .addArgument(GATKIndelRealignerCLI.OUT, indelRealignerOut.getAbsolutePath());
        CondorJob gatkIndelRealignerJob = builder.build();
        graph.addVertex(gatkIndelRealignerJob);
        graph.addEdge(gatkRealignTargetCreatorJob, gatkIndelRealignerJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, PicardFixMateCLI.class, null).siteName(siteName);
        File picardFixMateOutput = new File(workflowDirectory, indelRealignerOut.getName().replace(".bam", ".fixmate.bam"));
        builder.addArgument(PicardFixMateCLI.VALIDATE, Boolean.FALSE).addArgument(PicardFixMateCLI.DRYRUN)
                .addArgument(PicardFixMateCLI.SORTORDER, PicardSortOrderType.COORDINATE.toString().toLowerCase())
                .addArgument(PicardFixMateCLI.INPUT, indelRealignerOut.getAbsolutePath())
                .addArgument(PicardFixMateCLI.OUTPUT, picardFixMateOutput.getAbsolutePath());
        CondorJob picardFixMateJob = builder.build();
        graph.addVertex(picardFixMateJob);
        graph.addEdge(gatkIndelRealignerJob, picardFixMateJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, null).siteName(siteName);
        File picardFixMateIndexOut = new File(workflowDirectory, picardFixMateOutput.getName().replace(".bam", ".bai"));
        builder.addArgument(SAMToolsIndexCLI.VALIDATE, Boolean.FALSE).addArgument(SAMToolsIndexCLI.DRYRUN)
                .addArgument(SAMToolsIndexCLI.INPUT, picardFixMateOutput.getAbsolutePath())
                .addArgument(SAMToolsIndexCLI.OUTPUT, picardFixMateIndexOut.getAbsolutePath());
        samtoolsIndexJob = builder.build();
        graph.addVertex(samtoolsIndexJob);
        graph.addEdge(picardFixMateJob, samtoolsIndexJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, GATKCountCovariatesCLI.class, null).siteName(siteName)
                .numberOfProcessors(4);
        File gatkCountCovariatesRecalFile = new File(workflowDirectory, picardFixMateOutput.getName().replace(".bam", ".bam.cov"));
        builder.addArgument(GATKCountCovariatesCLI.VALIDATE, Boolean.FALSE).addArgument(GATKCountCovariatesCLI.DRYRUN)
                .addArgument(GATKCountCovariatesCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                .addArgument(GATKCountCovariatesCLI.REFERENCESEQUENCE, referenceSequence)
                .addArgument(GATKCountCovariatesCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                .addArgument(GATKCountCovariatesCLI.KNOWNSITES, knownVCF).addArgument(GATKCountCovariatesCLI.NUMTHREADS, "4")
                .addArgument(GATKCountCovariatesCLI.STANDARDCOVS)
                .addArgument(GATKCountCovariatesCLI.INPUTFILE, picardFixMateOutput.getAbsolutePath())
                .addArgument(GATKCountCovariatesCLI.RECALFILE, gatkCountCovariatesRecalFile.getAbsolutePath());
        CondorJob gatkCountCovariatesJob = builder.build();
        graph.addVertex(gatkCountCovariatesJob);
        graph.addEdge(samtoolsIndexJob, gatkCountCovariatesJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, GATKTableRecalibrationCLI.class, null).siteName(siteName)
                .numberOfProcessors(2);
        File gatkTableRecalibrationOut = new File(workflowDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.bam"));
        builder.addArgument(GATKTableRecalibrationCLI.VALIDATE, Boolean.FALSE).addArgument(GATKTableRecalibrationCLI.DRYRUN)
                .addArgument(GATKTableRecalibrationCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                .addArgument(GATKTableRecalibrationCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                .addArgument(GATKTableRecalibrationCLI.REFERENCESEQUENCE, referenceSequence)
                .addArgument(GATKTableRecalibrationCLI.INPUTFILE, picardFixMateOutput.getAbsolutePath())
                .addArgument(GATKTableRecalibrationCLI.RECALFILE, gatkCountCovariatesRecalFile.getAbsolutePath())
                .addArgument(GATKTableRecalibrationCLI.OUT, gatkTableRecalibrationOut.getAbsolutePath());
        CondorJob gatkTableRecalibrationJob = builder.build();
        graph.addVertex(gatkTableRecalibrationJob);
        graph.addEdge(gatkCountCovariatesJob, gatkTableRecalibrationJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class, null).siteName(siteName);
        File gatkTableRecalibrationIndexOut = new File(workflowDirectory, gatkTableRecalibrationOut.getName().replace(".bam", ".bai"));
        builder.addArgument(SAMToolsIndexCLI.VALIDATE, Boolean.FALSE).addArgument(SAMToolsIndexCLI.DRYRUN)
                .addArgument(SAMToolsIndexCLI.INPUT, gatkTableRecalibrationOut.getAbsolutePath())
                .addArgument(SAMToolsIndexCLI.OUTPUT, gatkTableRecalibrationIndexOut.getAbsolutePath());
        samtoolsIndexJob = builder.build();
        graph.addVertex(samtoolsIndexJob);
        graph.addEdge(gatkTableRecalibrationJob, samtoolsIndexJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, SAMToolsFlagstatCLI.class, null).siteName(siteName);
        File samtoolsFlagstatOut = new File(workflowDirectory, gatkTableRecalibrationOut.getName().replace(".bam", ".samtools.flagstat"));
        builder.addArgument(SAMToolsFlagstatCLI.VALIDATE, Boolean.FALSE).addArgument(SAMToolsFlagstatCLI.DRYRUN)
                .addArgument(SAMToolsFlagstatCLI.INPUT, gatkTableRecalibrationOut.getAbsolutePath())
                .addArgument(SAMToolsFlagstatCLI.OUTPUT, samtoolsFlagstatOut.getAbsolutePath());
        CondorJob samtoolsFlagstatJob = builder.build();
        graph.addVertex(samtoolsFlagstatJob);
        graph.addEdge(samtoolsIndexJob, samtoolsFlagstatJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, GATKFlagStatCLI.class, null).siteName(siteName).numberOfProcessors(2);
        File gatkFlagstatOut = new File(workflowDirectory, gatkTableRecalibrationOut.getName().replace(".bam", ".gatk.flagstat"));
        builder.addArgument(GATKFlagStatCLI.VALIDATE, Boolean.FALSE).addArgument(GATKFlagStatCLI.DRYRUN)
                .addArgument(GATKFlagStatCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                .addArgument(GATKFlagStatCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                .addArgument(GATKFlagStatCLI.INPUTFILE, gatkTableRecalibrationOut.getAbsolutePath())
                .addArgument(GATKFlagStatCLI.REFERENCESEQUENCE, referenceSequence)
                .addArgument(GATKFlagStatCLI.INTERVALS, flagstatIntervalList)
                .addArgument(GATKFlagStatCLI.OUT, gatkFlagstatOut.getAbsolutePath());
        CondorJob gatkFlagstatJob = builder.build();
        graph.addVertex(gatkFlagstatJob);
        graph.addEdge(samtoolsIndexJob, gatkFlagstatJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, GATKDepthOfCoverageCLI.class, null).siteName(siteName)
                .initialDirectory(workflowDirectory.getAbsolutePath()).numberOfProcessors(2);
        builder.addArgument(GATKDepthOfCoverageCLI.VALIDATE, Boolean.FALSE).addArgument(GATKDepthOfCoverageCLI.DRYRUN)
                .addArgument(GATKDepthOfCoverageCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                .addArgument(GATKDepthOfCoverageCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                .addArgument(GATKDepthOfCoverageCLI.REFERENCESEQUENCE, referenceSequence)
                .addArgument(GATKDepthOfCoverageCLI.VALIDATIONSTRICTNESS, "LENIENT")
                .addArgument(GATKDepthOfCoverageCLI.OMITDEPTHOUTPUTATEACHBASE)
                .addArgument(GATKDepthOfCoverageCLI.INPUTFILE, gatkTableRecalibrationOut.getAbsolutePath())
                .addArgument(GATKDepthOfCoverageCLI.INTERVALS, depthOfCoverageIntervalList)
                .addArgument(GATKDepthOfCoverageCLI.OUTPUTPREFIX, gatkTableRecalibrationOut.getName().replace(".bam", ".coverage"));
        CondorJob gatkDepthOfCoverageJob = builder.build();
        graph.addVertex(gatkDepthOfCoverageJob);
        graph.addEdge(samtoolsFlagstatJob, gatkDepthOfCoverageJob);
        graph.addEdge(gatkFlagstatJob, gatkDepthOfCoverageJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, GATKUnifiedGenotyperCLI.class, null).siteName(siteName)
                .numberOfProcessors(4);
        File gatkUnifiedGenotyperOut = new File(workflowDirectory, gatkTableRecalibrationOut.getName().replace(".bam", ".vcf"));
        File gatkUnifiedGenotyperMetrics = new File(workflowDirectory, gatkTableRecalibrationOut.getName().replace(".bam", ".metrics"));
        builder.addArgument(GATKUnifiedGenotyperCLI.VALIDATE, Boolean.FALSE).addArgument(GATKUnifiedGenotyperCLI.DRYRUN)
                .addArgument(GATKUnifiedGenotyperCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
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
        graph.addVertex(gatkUnifiedGenotyperJob);
        graph.addEdge(gatkDepthOfCoverageJob, gatkUnifiedGenotyperJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, FilterVariantCLI.class, null).siteName(siteName).numberOfProcessors(2);
        File filterVariant1Output = new File(workflowDirectory, gatkTableRecalibrationOut.getName().replace(".bam", ".variant.vcf"));
        builder.addArgument(FilterVariantCLI.VALIDATE, Boolean.FALSE).addArgument(FilterVariantCLI.DRYRUN)
                .addArgument(FilterVariantCLI.INTERVALLIST, icSNPIntervalList).addArgument(FilterVariantCLI.WITHMISSING)
                .addArgument(FilterVariantCLI.INPUT, gatkUnifiedGenotyperOut.getAbsolutePath())
                .addArgument(FilterVariantCLI.OUTPUT, filterVariant1Output.getAbsolutePath());
        CondorJob filterVariant1Job = builder.build();
        graph.addVertex(filterVariant1Job);
        graph.addEdge(gatkUnifiedGenotyperJob, filterVariant1Job);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, GATKVariantRecalibratorCLI.class, null).siteName(siteName)
                .numberOfProcessors(2);
        File gatkVariantRecalibratorRecalFile = new File(workflowDirectory, filterVariant1Output.getName().replace(".vcf", ".recal"));
        File gatkVariantRecalibratorTranchesFile = new File(workflowDirectory, filterVariant1Output.getName().replace(".vcf", ".tranches"));
        File gatkVariantRecalibratorRScriptFile = new File(workflowDirectory, filterVariant1Output.getName().replace(".vcf", ".plots.R"));
        builder.addArgument(GATKVariantRecalibratorCLI.VALIDATE, Boolean.FALSE).addArgument(GATKVariantRecalibratorCLI.DRYRUN)
                .addArgument(GATKVariantRecalibratorCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
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
                .addArgument(GATKVariantRecalibratorCLI.USEANNOTATION, "MQ").addArgument(GATKVariantRecalibratorCLI.USEANNOTATION, "FS")
                .addArgument(GATKVariantRecalibratorCLI.RSCRIPTFILE, gatkVariantRecalibratorRScriptFile.getAbsolutePath());
        if (StringUtils.isNotEmpty(percentBadVariants)) {
            builder.addArgument(GATKVariantRecalibratorCLI.PERCENTBADVARIANTS, percentBadVariants);
        }
        CondorJob gatkVariantRecalibratorJob = builder.build();
        graph.addVertex(gatkVariantRecalibratorJob);
        graph.addEdge(filterVariant1Job, gatkVariantRecalibratorJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, GATKApplyRecalibrationCLI.class, null).siteName(siteName)
                .numberOfProcessors(2);
        File gatkApplyRecalibrationOut = new File(workflowDirectory,
                filterVariant1Output.getName().replace(".vcf", ".recalibrated.filtered.vcf"));
        builder.addArgument(GATKApplyRecalibrationCLI.VALIDATE, Boolean.FALSE).addArgument(GATKApplyRecalibrationCLI.DRYRUN)
                .addArgument(GATKApplyRecalibrationCLI.PHONEHOME, GATKPhoneHomeType.NO_ET.toString())
                .addArgument(GATKApplyRecalibrationCLI.DOWNSAMPLINGTYPE, GATKDownsamplingType.NONE.toString())
                .addArgument(GATKApplyRecalibrationCLI.REFERENCESEQUENCE, referenceSequence)
                .addArgument(GATKApplyRecalibrationCLI.INPUT, filterVariant1Output.getAbsolutePath())
                .addArgument(GATKApplyRecalibrationCLI.RECALFILE, gatkVariantRecalibratorRecalFile.getAbsolutePath())
                .addArgument(GATKApplyRecalibrationCLI.TRANCHESFILE, gatkVariantRecalibratorTranchesFile.getAbsolutePath())
                .addArgument(GATKApplyRecalibrationCLI.OUT, gatkApplyRecalibrationOut.getAbsolutePath())
                .addArgument(GATKApplyRecalibrationCLI.TSFILTERLEVEL, "99.0");
        CondorJob gatkApplyRecalibrationJob = builder.build();
        graph.addVertex(gatkApplyRecalibrationJob);
        graph.addEdge(gatkVariantRecalibratorJob, gatkApplyRecalibrationJob);

        // new job
        builder = SequencingWorkflowJobFactory.createJob(++count, FilterVariantCLI.class, null).siteName(siteName).numberOfProcessors(2);
        File filterVariant2Output = new File(workflowDirectory, filterVariant1Output.getName().replace(".vcf", ".ic_snps.vcf"));
        builder.addArgument(FilterVariantCLI.VALIDATE, Boolean.FALSE).addArgument(FilterVariantCLI.DRYRUN)
                .addArgument(FilterVariantCLI.INTERVALLIST, icSNPIntervalList)
                .addArgument(FilterVariantCLI.INPUT, gatkApplyRecalibrationOut.getAbsolutePath())
                .addArgument(FilterVariantCLI.OUTPUT, filterVariant2Output.getAbsolutePath());
        CondorJob filterVariant2Job = builder.build();
        graph.addVertex(filterVariant2Job);
        graph.addEdge(gatkApplyRecalibrationJob, filterVariant2Job);

        CLIScriptExporter exporter = new CLIScriptExporter();
        File tmpDir = new File("/tmp/ncgenes");
        if (!tmpDir.exists()) {
            tmpDir.mkdirs();
        }
        exporter.export("test", tmpDir, graph);
    }

}
