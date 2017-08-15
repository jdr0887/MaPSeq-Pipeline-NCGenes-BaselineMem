package edu.unc.mapseq.commons.ncgenes.baselinemem;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.osgi.framework.Bundle;
import org.osgi.framework.BundleContext;
import org.osgi.framework.FrameworkUtil;
import org.renci.common.exec.BashExecutor;
import org.renci.common.exec.CommandInput;
import org.renci.common.exec.CommandOutput;
import org.renci.common.exec.Executor;
import org.renci.common.exec.ExecutorException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.dao.MaPSeqDAOBeanService;
import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.Job;
import edu.unc.mapseq.dao.model.MimeType;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.dao.model.Workflow;
import edu.unc.mapseq.dao.model.WorkflowRun;
import edu.unc.mapseq.dao.model.WorkflowRunAttempt;
import edu.unc.mapseq.module.sequencing.WriteVCFHeader;
import edu.unc.mapseq.module.sequencing.fastqc.FastQC;
import edu.unc.mapseq.module.sequencing.filter.FilterVariant;
import edu.unc.mapseq.module.sequencing.gatk.GATKApplyRecalibration;
import edu.unc.mapseq.module.sequencing.gatk.GATKDepthOfCoverage;
import edu.unc.mapseq.module.sequencing.gatk.GATKFlagStat;
import edu.unc.mapseq.module.sequencing.gatk.GATKTableRecalibration;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsFlagstat;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsIndex;
import edu.unc.mapseq.workflow.WorkflowException;
import edu.unc.mapseq.workflow.sequencing.IRODSBean;
import edu.unc.mapseq.workflow.sequencing.SequencingWorkflowUtil;

public class RegisterToIRODSRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(RegisterToIRODSRunnable.class);

    private MaPSeqDAOBeanService mapseqDAOBeanService;

    private WorkflowRunAttempt workflowRunAttempt;

    public RegisterToIRODSRunnable(MaPSeqDAOBeanService mapseqDAOBeanService, WorkflowRunAttempt workflowRunAttempt) {
        super();
        this.mapseqDAOBeanService = mapseqDAOBeanService;
        this.workflowRunAttempt = workflowRunAttempt;
    }

    @Override
    public void run() {
        logger.info("ENTERING run()");

        final WorkflowRun workflowRun = workflowRunAttempt.getWorkflowRun();
        final Workflow workflow = workflowRun.getWorkflow();

        BundleContext bundleContext = FrameworkUtil.getBundle(getClass()).getBundleContext();
        Bundle bundle = bundleContext.getBundle();
        String version = bundle.getVersion().toString();

        try {

            Set<Sample> sampleSet = SequencingWorkflowUtil.getAggregatedSamples(mapseqDAOBeanService, workflowRunAttempt);

            if (CollectionUtils.isEmpty(sampleSet)) {
                logger.warn("No Samples found");
                return;
            }

            for (Sample sample : sampleSet) {

                File outputDirectory = SequencingWorkflowUtil.createOutputDirectory(sample, workflow);
                File tmpDir = new File(outputDirectory, "tmp");
                if (!tmpDir.exists()) {
                    tmpDir.mkdirs();
                }

                // assumption: a dash is used as a delimiter between a participantId and the external code
                String participantId;

                if (!sample.getName().startsWith("HRC")) {
                    int idx = sample.getName().lastIndexOf("-");
                    participantId = idx != -1 ? sample.getName().substring(0, idx) : sample.getName();
                } else {
                    participantId = sample.getName();
                }

                String irodsDirectory = String.format("/MedGenZone/%s/sequencing/ncgenes/analysis/%s/L%03d_%s/%s",
                        workflow.getSystem().getValue(), sample.getFlowcell().getName(), sample.getLaneIndex(), sample.getBarcode(),
                        workflow.getName());

                CommandOutput commandOutput = null;

                List<CommandInput> commandInputList = new LinkedList<CommandInput>();

                CommandInput commandInput = new CommandInput();
                commandInput.setExitImmediately(Boolean.FALSE);
                StringBuilder sb = new StringBuilder();
                sb.append(String.format("$IRODS_HOME/imkdir -p %s%n", irodsDirectory));
                sb.append(String.format("$IRODS_HOME/imeta add -C %s Project NCGENES%n", irodsDirectory));
                sb.append(String.format("$IRODS_HOME/imeta add -C %s ParticipantId %s NCGENES%n", irodsDirectory, participantId));
                commandInput.setCommand(sb.toString());
                commandInput.setWorkDir(tmpDir);
                commandInputList.add(commandInput);

                List<IRODSBean> files2RegisterToIRODS = new ArrayList<IRODSBean>();

                String rootFileName = String.format("%s_%s_L%03d", sample.getFlowcell().getName(), sample.getBarcode(),
                        sample.getLaneIndex());

                List<ImmutablePair<String, String>> attributeList = Arrays.asList(
                        new ImmutablePair<String, String>("ParticipantId", participantId),
                        new ImmutablePair<String, String>("MaPSeqWorkflowVersion", version),
                        new ImmutablePair<String, String>("MaPSeqWorkflowName", workflow.getName()),
                        new ImmutablePair<String, String>("MaPSeqStudyName", sample.getStudy().getName()),
                        new ImmutablePair<String, String>("MaPSeqSampleId", sample.getId().toString()),
                        new ImmutablePair<String, String>("MaPSeqSystem", workflow.getSystem().getValue()),
                        new ImmutablePair<String, String>("MaPSeqFlowcellId", sample.getFlowcell().getId().toString()));

                List<ImmutablePair<String, String>> attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", WriteVCFHeader.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
                File file = new File(outputDirectory, String.format("%s.vcf.hdr", rootFileName));
                Job job = SequencingWorkflowUtil.findJob(mapseqDAOBeanService, workflowRunAttempt.getId(), WriteVCFHeader.class.getName(),
                        file);
                if (job != null) {
                    attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobId", job.getId().toString()));
                    if (StringUtils.isNotEmpty(job.getCommandLine())) {
                        attributeListWithJob.add(
                                new ImmutablePair<String, String>("MaPSeqJobCommandLine", String.format("'%s'", job.getCommandLine())));
                    }
                } else {
                    logger.warn(String.format("Couldn't find job for: %d, %s", workflowRunAttempt.getId(), WriteVCFHeader.class.getName()));
                }
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", FastQC.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.APPLICATION_ZIP.toString()));
                file = new File(outputDirectory, String.format("%s_R1.fastqc.zip", rootFileName));
                job = SequencingWorkflowUtil.findJob(mapseqDAOBeanService, workflowRunAttempt.getId(), FastQC.class.getName(), file);
                if (job != null) {
                    attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobId", job.getId().toString()));
                    if (StringUtils.isNotEmpty(job.getCommandLine())) {
                        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobCommandLine", job.getCommandLine()));
                    }
                } else {
                    logger.warn(String.format("Couldn't find job for: %d, %s", workflowRunAttempt.getId(), FastQC.class.getName()));
                }
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", FastQC.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.APPLICATION_ZIP.toString()));
                file = new File(outputDirectory, String.format("%s_R2.fastqc.zip", rootFileName));
                job = SequencingWorkflowUtil.findJob(mapseqDAOBeanService, workflowRunAttempt.getId(), FastQC.class.getName(), file);
                if (job != null) {
                    attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobId", job.getId().toString()));
                    if (StringUtils.isNotEmpty(job.getCommandLine())) {
                        attributeListWithJob.add(
                                new ImmutablePair<String, String>("MaPSeqJobCommandLine", String.format("'%s'", job.getCommandLine())));
                    }
                } else {
                    logger.warn(String.format("Couldn't find job for: %d, %s", workflowRunAttempt.getId(), FastQC.class.getName()));
                }
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                File bwaMemOutFile = new File(outputDirectory, String.format("%s.mem.sam", rootFileName));
                File fixRGOutput = new File(outputDirectory, bwaMemOutFile.getName().replace(".sam", ".fixed-rg.bam"));
                File picardMarkDuplicatesOutput = new File(outputDirectory, fixRGOutput.getName().replace(".bam", ".deduped.bam"));
                File indelRealignerOut = new File(outputDirectory, picardMarkDuplicatesOutput.getName().replace(".bam", ".realign.bam"));
                File picardFixMateOutput = new File(outputDirectory, indelRealignerOut.getName().replace(".bam", ".fixmate.bam"));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKTableRecalibration.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.APPLICATION_BAM.toString()));
                file = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.bam"));
                job = SequencingWorkflowUtil.findJob(mapseqDAOBeanService, workflowRunAttempt.getId(),
                        GATKTableRecalibration.class.getName(), file);
                if (job != null) {
                    attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobId", job.getId().toString()));
                    if (StringUtils.isNotEmpty(job.getCommandLine())) {
                        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobCommandLine", job.getCommandLine()));
                    }
                } else {
                    logger.warn(String.format("Couldn't find job for: %d, %s", workflowRunAttempt.getId(),
                            GATKTableRecalibration.class.getName()));
                }
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", SAMToolsIndex.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.APPLICATION_BAM_INDEX.toString()));
                file = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.bai"));
                job = SequencingWorkflowUtil.findJob(mapseqDAOBeanService, workflowRunAttempt.getId(), SAMToolsIndex.class.getName(), file);
                if (job != null) {
                    attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobId", job.getId().toString()));
                    if (StringUtils.isNotEmpty(job.getCommandLine())) {
                        attributeListWithJob.add(
                                new ImmutablePair<String, String>("MaPSeqJobCommandLine", String.format("'%s'", job.getCommandLine())));
                    }
                } else {
                    logger.warn(String.format("Couldn't find job for: %d, %s", workflowRunAttempt.getId(), SAMToolsIndex.class.getName()));
                }
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKDepthOfCoverage.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
                file = new File(outputDirectory,
                        picardFixMateOutput.getName().replace(".bam", ".recal.coverage.sample_cumulative_coverage_counts"));
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKDepthOfCoverage.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
                file = new File(outputDirectory,
                        picardFixMateOutput.getName().replace(".bam", ".recal.coverage.sample_cumulative_coverage_proportions"));
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKDepthOfCoverage.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
                file = new File(outputDirectory,
                        picardFixMateOutput.getName().replace(".bam", ".recal.coverage.sample_interval_statistics"));
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKDepthOfCoverage.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
                file = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.coverage.sample_interval_summary"));
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKDepthOfCoverage.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
                file = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.coverage.sample_statistics"));
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKDepthOfCoverage.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_PLAIN.toString()));
                file = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.coverage.sample_summary"));
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", SAMToolsFlagstat.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_STAT_SUMMARY.toString()));
                file = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.samtools.flagstat"));
                job = SequencingWorkflowUtil.findJob(mapseqDAOBeanService, workflowRunAttempt.getId(), SAMToolsFlagstat.class.getName(),
                        file);
                if (job != null) {
                    attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobId", job.getId().toString()));
                    if (StringUtils.isNotEmpty(job.getCommandLine())) {
                        attributeListWithJob.add(
                                new ImmutablePair<String, String>("MaPSeqJobCommandLine", String.format("'%s'", job.getCommandLine())));
                    }
                } else {
                    logger.warn(
                            String.format("Couldn't find job for: %d, %s", workflowRunAttempt.getId(), SAMToolsFlagstat.class.getName()));
                }
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKFlagStat.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_STAT_SUMMARY.toString()));
                file = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.gatk.flagstat"));
                job = SequencingWorkflowUtil.findJob(mapseqDAOBeanService, workflowRunAttempt.getId(), GATKFlagStat.class.getName(), file);
                if (job != null) {
                    attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobId", job.getId().toString()));
                    if (StringUtils.isNotEmpty(job.getCommandLine())) {
                        attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobCommandLine", job.getCommandLine()));
                    }
                } else {
                    logger.warn(String.format("Couldn't find job for: %d, %s", workflowRunAttempt.getId(), GATKFlagStat.class.getName()));
                }
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", FilterVariant.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_VCF.toString()));
                file = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.variant.vcf"));
                job = SequencingWorkflowUtil.findJob(mapseqDAOBeanService, workflowRunAttempt.getId(), FilterVariant.class.getName(), file);
                if (job != null) {
                    attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobId", job.getId().toString()));
                    if (StringUtils.isNotEmpty(job.getCommandLine())) {
                        attributeListWithJob.add(
                                new ImmutablePair<String, String>("MaPSeqJobCommandLine", String.format("'%s'", job.getCommandLine())));
                    }
                } else {
                    logger.warn(String.format("Couldn't find job for: %d, %s", workflowRunAttempt.getId(), FilterVariant.class.getName()));
                }
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", GATKApplyRecalibration.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_VCF.toString()));
                file = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.variant.recalibrated.filtered.vcf"));
                job = SequencingWorkflowUtil.findJob(mapseqDAOBeanService, workflowRunAttempt.getId(),
                        GATKApplyRecalibration.class.getName(), file);
                if (job != null) {
                    attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobId", job.getId().toString()));
                    if (StringUtils.isNotEmpty(job.getCommandLine())) {
                        attributeListWithJob.add(
                                new ImmutablePair<String, String>("MaPSeqJobCommandLine", String.format("'%s'", job.getCommandLine())));
                    }
                } else {
                    logger.warn(String.format("Couldn't find job for: %d, %s", workflowRunAttempt.getId(),
                            GATKApplyRecalibration.class.getName()));
                }
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                attributeListWithJob = new ArrayList<>(attributeList);
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobName", FilterVariant.class.getSimpleName()));
                attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqMimeType", MimeType.TEXT_VCF.toString()));
                file = new File(outputDirectory, picardFixMateOutput.getName().replace(".bam", ".recal.variant.ic_snps.vcf"));
                job = SequencingWorkflowUtil.findJob(mapseqDAOBeanService, workflowRunAttempt.getId(), FilterVariant.class.getName(), file);
                if (job != null) {
                    attributeListWithJob.add(new ImmutablePair<String, String>("MaPSeqJobId", job.getId().toString()));
                    if (StringUtils.isNotEmpty(job.getCommandLine())) {
                        attributeListWithJob.add(
                                new ImmutablePair<String, String>("MaPSeqJobCommandLine", String.format("'%s'", job.getCommandLine())));
                    }
                } else {
                    logger.warn(String.format("Couldn't find job for: %d, %s", workflowRunAttempt.getId(), FilterVariant.class.getName()));
                }
                files2RegisterToIRODS.add(new IRODSBean(file, attributeListWithJob));

                for (IRODSBean bean : files2RegisterToIRODS) {

                    commandInput = new CommandInput();
                    commandInput.setExitImmediately(Boolean.FALSE);

                    File f = bean.getFile();
                    if (!f.exists()) {
                        logger.warn("file to register doesn't exist: {}", f.getAbsolutePath());
                        continue;
                    }

                    StringBuilder registerCommandSB = new StringBuilder();
                    String registrationCommand = String.format("$IRODS_HOME/ireg -f %s %s/%s", bean.getFile().getAbsolutePath(),
                            irodsDirectory, bean.getFile().getName());
                    String deRegistrationCommand = String.format("$IRODS_HOME/irm -U %s/%s", irodsDirectory, bean.getFile().getName());
                    registerCommandSB.append(registrationCommand).append("\n");
                    registerCommandSB
                            .append(String.format("if [ $? != 0 ]; then %s; %s; fi%n", deRegistrationCommand, registrationCommand));
                    commandInput.setCommand(registerCommandSB.toString());
                    commandInput.setWorkDir(tmpDir);
                    commandInputList.add(commandInput);

                    commandInput = new CommandInput();
                    commandInput.setExitImmediately(Boolean.FALSE);
                    sb = new StringBuilder();
                    for (ImmutablePair<String, String> attribute : bean.getAttributes()) {
                        sb.append(String.format("$IRODS_HOME/imeta add -d %s/%s %s %s NCGenes%n", irodsDirectory, bean.getFile().getName(),
                                attribute.getLeft(), attribute.getRight()));
                    }
                    commandInput.setCommand(sb.toString());
                    commandInput.setWorkDir(tmpDir);
                    commandInputList.add(commandInput);

                }

                File mapseqrc = new File(System.getProperty("user.home"), ".mapseqrc");
                Executor executor = BashExecutor.getInstance();

                for (CommandInput ci : commandInputList) {
                    try {
                        logger.debug("ci.getCommand(): {}", ci.getCommand());
                        commandOutput = executor.execute(ci, mapseqrc);
                        if (commandOutput.getExitCode() != 0) {
                            logger.info("commandOutput.getExitCode(): {}", commandOutput.getExitCode());
                            logger.warn("command failed: {}", ci.getCommand());
                        }
                        logger.debug("commandOutput.getStdout(): {}", commandOutput.getStdout());
                    } catch (ExecutorException e) {
                        if (commandOutput != null) {
                            logger.warn("commandOutput.getStderr(): {}", commandOutput.getStderr());
                        }
                    }
                }

                logger.info("FINISHED PROCESSING: {}", sample.toString());

            }

        } catch (MaPSeqDAOException | WorkflowException e) {
            e.printStackTrace();
        }

    }

}
