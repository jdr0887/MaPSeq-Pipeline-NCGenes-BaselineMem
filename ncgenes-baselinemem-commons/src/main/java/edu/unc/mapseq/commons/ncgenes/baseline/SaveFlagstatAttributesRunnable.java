package edu.unc.mapseq.commons.ncgenes.baseline;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.dao.MaPSeqDAOBeanService;
import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.Attribute;
import edu.unc.mapseq.dao.model.MimeType;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.dao.model.Workflow;
import edu.unc.mapseq.dao.model.WorkflowRun;
import edu.unc.mapseq.dao.model.WorkflowRunAttempt;
import edu.unc.mapseq.module.sequencing.samtools.SAMToolsFlagstat;
import edu.unc.mapseq.workflow.WorkflowException;
import edu.unc.mapseq.workflow.core.WorkflowUtil;
import edu.unc.mapseq.workflow.sequencing.SequencingWorkflowUtil;

public class SaveFlagstatAttributesRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(SaveFlagstatAttributesRunnable.class);

    private MaPSeqDAOBeanService mapseqDAOBeanService;

    private WorkflowRunAttempt workflowRunAttempt;

    public SaveFlagstatAttributesRunnable() {
        super();
    }

    public SaveFlagstatAttributesRunnable(MaPSeqDAOBeanService mapseqDAOBeanService, WorkflowRunAttempt workflowRunAttempt) {
        super();
        this.mapseqDAOBeanService = mapseqDAOBeanService;
        this.workflowRunAttempt = workflowRunAttempt;
    }

    @Override
    public void run() {
        logger.info("ENTERING run()");

        final WorkflowRun workflowRun = workflowRunAttempt.getWorkflowRun();
        final Workflow workflow = workflowRun.getWorkflow();

        try {
            Set<Sample> sampleSet = SequencingWorkflowUtil.getAggregatedSamples(mapseqDAOBeanService, workflowRunAttempt);

            if (CollectionUtils.isEmpty(sampleSet)) {
                logger.warn("No Samples found");
                return;
            }

            for (Sample sample : sampleSet) {

                logger.info(sample.toString());

                File outputDirectory = SequencingWorkflowUtil.createOutputDirectory(sample, workflow);

                Set<Attribute> attributeSet = sample.getAttributes();

                File flagstatFile = WorkflowUtil.findFileByJobAndMimeTypeAndWorkflowId(this.mapseqDAOBeanService, sample.getFileDatas(),
                        SAMToolsFlagstat.class, MimeType.TEXT_STAT_SUMMARY, workflow.getId());

                if (flagstatFile == null) {
                    logger.error("flagstat file to process was not found...checking FS");
                    if (outputDirectory.exists()) {
                        for (File file : outputDirectory.listFiles()) {
                            if (file.getName().endsWith("samtools.flagstat")) {
                                flagstatFile = file;
                                break;
                            }
                        }
                    }
                }

                if (flagstatFile == null) {
                    logger.error("flagstat file to process was still not found");
                    continue;
                }

                Set<String> attributeNameSet = new HashSet<String>();

                for (Attribute attribute : attributeSet) {
                    attributeNameSet.add(attribute.getName());
                }

                Set<String> synchSet = Collections.synchronizedSet(attributeNameSet);

                List<String> lines = FileUtils.readLines(flagstatFile);
                if (CollectionUtils.isNotEmpty(lines)) {

                    for (String line : lines) {

                        if (line.contains("in total")) {
                            String value = line.substring(0, line.indexOf(" ")).trim();
                            if (synchSet.contains("SAMToolsFlagstat.totalPassedReads")) {
                                for (Attribute attribute : attributeSet) {
                                    if (attribute.getName().equals("SAMToolsFlagstat.totalPassedReads")) {
                                        attribute.setValue(value);
                                        try {
                                            mapseqDAOBeanService.getAttributeDAO().save(attribute);
                                        } catch (MaPSeqDAOException e) {
                                            logger.error("MaPSeqDAOException", e);
                                        }
                                        break;
                                    }
                                }
                            } else {
                                attributeSet.add(new Attribute("SAMToolsFlagstat.totalPassedReads", value));
                            }
                        }

                        if (line.contains("mapped (")) {
                            Pattern pattern = Pattern.compile("^.+\\((.+)\\)");
                            Matcher matcher = pattern.matcher(line);
                            if (matcher.matches()) {
                                String value = matcher.group(1);
                                value = value.substring(0, value.indexOf("%")).trim();
                                if (StringUtils.isNotEmpty(value)) {
                                    if (synchSet.contains("SAMToolsFlagstat.aligned")) {
                                        for (Attribute attribute : attributeSet) {
                                            if (attribute.getName().equals("SAMToolsFlagstat.aligned")) {
                                                attribute.setValue(value);
                                                break;
                                            }
                                        }
                                    } else {
                                        attributeSet.add(new Attribute("SAMToolsFlagstat.aligned", value));
                                    }
                                }
                            }
                        }

                        if (line.contains("properly paired (")) {
                            Pattern pattern = Pattern.compile("^.+\\((.+)\\)");
                            Matcher matcher = pattern.matcher(line);
                            if (matcher.matches()) {
                                String value = matcher.group(1);
                                value = value.substring(0, value.indexOf("%"));
                                if (StringUtils.isNotEmpty(value)) {
                                    if (synchSet.contains("SAMToolsFlagstat.paired")) {
                                        for (Attribute attribute : attributeSet) {
                                            if (attribute.getName().equals("SAMToolsFlagstat.paired")) {
                                                attribute.setValue(value);
                                                break;
                                            }
                                        }
                                    } else {
                                        attributeSet.add(new Attribute("SAMToolsFlagstat.paired", value));
                                    }
                                }
                            }
                        }
                    }

                }

                sample.setAttributes(attributeSet);
                mapseqDAOBeanService.getSampleDAO().save(sample);

            }
        } catch (MaPSeqDAOException | WorkflowException | IOException e) {
            e.printStackTrace();
        }

    }

    public MaPSeqDAOBeanService getMapseqDAOBeanService() {
        return mapseqDAOBeanService;
    }

    public void setMapseqDAOBeanService(MaPSeqDAOBeanService mapseqDAOBeanService) {
        this.mapseqDAOBeanService = mapseqDAOBeanService;
    }

    public WorkflowRunAttempt getWorkflowRunAttempt() {
        return workflowRunAttempt;
    }

    public void setWorkflowRunAttempt(WorkflowRunAttempt workflowRunAttempt) {
        this.workflowRunAttempt = workflowRunAttempt;
    }

}
