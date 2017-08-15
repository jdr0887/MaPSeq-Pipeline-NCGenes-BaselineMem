package edu.unc.mapseq.commons.ncgenes.baseline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.dao.MaPSeqDAOBeanService;
import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.Attribute;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.dao.model.Workflow;
import edu.unc.mapseq.dao.model.WorkflowRun;
import edu.unc.mapseq.dao.model.WorkflowRunAttempt;
import edu.unc.mapseq.workflow.WorkflowException;
import edu.unc.mapseq.workflow.sequencing.SequencingWorkflowUtil;

public class SaveDepthOfCoverageAttributesRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(SaveDepthOfCoverageAttributesRunnable.class);

    private MaPSeqDAOBeanService mapseqDAOBeanService;

    private WorkflowRunAttempt workflowRunAttempt;

    public SaveDepthOfCoverageAttributesRunnable() {
        super();
    }

    public SaveDepthOfCoverageAttributesRunnable(MaPSeqDAOBeanService mapseqDAOBeanService, WorkflowRunAttempt workflowRunAttempt) {
        super();
        this.mapseqDAOBeanService = mapseqDAOBeanService;
        this.workflowRunAttempt = workflowRunAttempt;
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        final WorkflowRun workflowRun = workflowRunAttempt.getWorkflowRun();
        final Workflow workflow = workflowRun.getWorkflow();

        try {
            Set<Sample> sampleSet = SequencingWorkflowUtil.getAggregatedSamples(mapseqDAOBeanService, workflowRunAttempt);

            if (CollectionUtils.isEmpty(sampleSet)) {
                logger.warn("No Samples found");
                return;
            }

            for (Sample sample : sampleSet) {

                File outputDirectory = SequencingWorkflowUtil.createOutputDirectory(sample, workflow);

                if (!outputDirectory.exists()) {
                    continue;
                }

                Set<Attribute> attributeSet = sample.getAttributes();

                Set<String> attributeNameSet = new HashSet<String>();

                for (Attribute attribute : attributeSet) {
                    attributeNameSet.add(attribute.getName());
                }

                Set<String> synchSet = Collections.synchronizedSet(attributeNameSet);

                List<File> files = Arrays.asList(outputDirectory.listFiles());

                if (files == null || (files != null && files.isEmpty())) {
                    logger.warn("no files found");
                    continue;
                }

                File sampleSummaryFile = null;
                for (File f : files) {
                    if (f.getName().endsWith(".coverage.sample_summary")) {
                        sampleSummaryFile = f;
                        break;
                    }
                }

                if (sampleSummaryFile != null && sampleSummaryFile.exists()) {
                    List<String> lines = FileUtils.readLines(sampleSummaryFile);
                    if (CollectionUtils.isNotEmpty(lines)) {
                        for (String line : lines) {
                            if (line.contains("Total")) {
                                String[] split = StringUtils.split(line);

                                if (synchSet.contains("GATKDepthOfCoverage.totalCoverage")) {
                                    for (Attribute attribute : attributeSet) {
                                        if (attribute.getName().equals("GATKDepthOfCoverage.totalCoverage")) {
                                            attribute.setValue(split[1]);
                                            break;
                                        }
                                    }
                                } else {
                                    attributeSet.add(new Attribute("GATKDepthOfCoverage.totalCoverage", split[1]));
                                }

                                if (synchSet.contains("GATKDepthOfCoverage.mean")) {
                                    for (Attribute attribute : attributeSet) {
                                        if (attribute.getName().equals("GATKDepthOfCoverage.mean")) {
                                            attribute.setValue(split[1]);
                                            break;
                                        }
                                    }
                                } else {
                                    attributeSet.add(new Attribute("GATKDepthOfCoverage.mean", split[2]));
                                }
                            }
                        }
                    }

                    File sampleIntervalSummaryFile = null;
                    for (File f : files) {
                        if (f.getName().endsWith(".coverage.sample_interval_summary")) {
                            sampleIntervalSummaryFile = f;
                            break;
                        }
                    }

                    if (sampleIntervalSummaryFile != null && sampleIntervalSummaryFile.exists()) {

                        long totalCoverageCount = 0;

                        try (FileReader fr = new FileReader(sampleIntervalSummaryFile); BufferedReader br = new BufferedReader(fr)) {
                            String line;
                            br.readLine();
                            while ((line = br.readLine()) != null) {
                                totalCoverageCount += Long.valueOf(StringUtils.split(line)[1].trim());
                            }
                        } catch (IOException e) {
                            e.printStackTrace();
                        }

                        if (synchSet.contains("GATKDepthOfCoverage.totalCoverageCount")) {
                            for (Attribute attribute : attributeSet) {
                                if (attribute.getName().equals("GATKDepthOfCoverage.totalCoverageCount")) {
                                    attribute.setValue(totalCoverageCount + "");
                                    break;
                                }
                            }
                        } else {
                            attributeSet.add(new Attribute("GATKDepthOfCoverage.totalCoverageCount", totalCoverageCount + ""));
                        }

                        Long totalPassedReads = null;
                        for (Attribute attribute : attributeSet) {
                            if ("SAMToolsFlagstat.totalPassedReads".equals(attribute.getName())) {
                                totalPassedReads = Long.valueOf(attribute.getValue());
                                break;
                            }
                        }

                        if (totalPassedReads != null) {
                            if (synchSet.contains("numberOnTarget")) {
                                for (Attribute attribute : attributeSet) {
                                    if (attribute.getName().equals("numberOnTarget") && totalPassedReads != null) {
                                        attribute.setValue((double) totalCoverageCount / (totalPassedReads * 100) + "");
                                        break;
                                    }
                                }
                            } else {
                                attributeSet
                                        .add(new Attribute("numberOnTarget", (double) totalCoverageCount / (totalPassedReads * 100) + ""));
                            }
                        }

                    }

                }

                sample.setAttributes(attributeSet);
                mapseqDAOBeanService.getSampleDAO().save(sample);

            }
        } catch (MaPSeqDAOException | IOException | NumberFormatException | WorkflowException e) {
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
