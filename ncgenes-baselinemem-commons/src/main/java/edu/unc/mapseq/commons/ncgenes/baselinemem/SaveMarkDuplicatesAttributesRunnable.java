package edu.unc.mapseq.commons.ncgenes.baseline;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.FileFilterUtils;
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

public class SaveMarkDuplicatesAttributesRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(SaveMarkDuplicatesAttributesRunnable.class);

    private MaPSeqDAOBeanService mapseqDAOBeanService;

    private WorkflowRunAttempt workflowRunAttempt;

    public SaveMarkDuplicatesAttributesRunnable() {
        super();
    }

    public SaveMarkDuplicatesAttributesRunnable(MaPSeqDAOBeanService mapseqDAOBeanService, WorkflowRunAttempt workflowRunAttempt) {
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

                File outputDirectory = SequencingWorkflowUtil.createOutputDirectory(sample, workflow);

                Set<Attribute> attributeSet = sample.getAttributes();

                Set<String> attributeNameSet = new HashSet<String>();

                for (Attribute attribute : attributeSet) {
                    attributeNameSet.add(attribute.getName());
                }

                Set<String> synchSet = Collections.synchronizedSet(attributeNameSet);

                Collection<File> fileList = FileUtils.listFiles(outputDirectory, FileFilterUtils.suffixFileFilter(".deduped.metrics"),
                        null);

                if (CollectionUtils.isNotEmpty(fileList)) {
                    File picardMarkDuplicatesMetricsFile = fileList.iterator().next();
                    List<String> lines = FileUtils.readLines(picardMarkDuplicatesMetricsFile);
                    if (CollectionUtils.isNotEmpty(lines)) {
                        Iterator<String> lineIter = lines.iterator();
                        while (lineIter.hasNext()) {
                            String line = lineIter.next();
                            if (line.startsWith("LIBRARY")) {
                                String nextLine = lineIter.next();
                                String[] split = StringUtils.split(nextLine);

                                String readPairDuplicates = split[6];
                                if (synchSet.contains("PicardMarkDuplicates.readPairDuplicates")) {
                                    for (Attribute attribute : attributeSet) {
                                        if (attribute.getName().equals("PicardMarkDuplicates.readPairDuplicates")) {
                                            attribute.setValue(readPairDuplicates);
                                            break;
                                        }
                                    }
                                } else {
                                    attributeSet.add(new Attribute("PicardMarkDuplicates.readPairDuplicates", readPairDuplicates));
                                }

                                String percentDuplication = split[7];
                                if (synchSet.contains("PicardMarkDuplicates.percentDuplication")) {
                                    for (Attribute attribute : attributeSet) {
                                        if (attribute.getName().equals("PicardMarkDuplicates.percentDuplication")) {
                                            attribute.setValue(percentDuplication);
                                            break;
                                        }
                                    }
                                } else {
                                    attributeSet.add(new Attribute("PicardMarkDuplicates.percentDuplication", percentDuplication));
                                }

                                break;
                            }

                        }
                    }

                    sample.setAttributes(attributeSet);
                    mapseqDAOBeanService.getSampleDAO().save(sample);
                }
            }
        } catch (WorkflowException | IOException | MaPSeqDAOException e) {
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
