package edu.unc.mapseq.ws.ncgenes.baselinemem.impl;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.StringUtils;
import org.renci.vcf.VCFParser;
import org.renci.vcf.VCFResult;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.SampleDAO;
import edu.unc.mapseq.dao.WorkflowDAO;
import edu.unc.mapseq.dao.model.FileData;
import edu.unc.mapseq.dao.model.MimeType;
import edu.unc.mapseq.dao.model.Sample;
import edu.unc.mapseq.dao.model.Workflow;
import edu.unc.mapseq.workflow.sequencing.SequencingWorkflowUtil;
import edu.unc.mapseq.ws.ncgenes.baselinemem.NCGenesBaselineMemService;
import edu.unc.mapseq.ws.ncgenes.baselinemem.QualityControlInfo;

public class NCGenesBaselineMemServiceImpl implements NCGenesBaselineMemService {

    private static final Logger logger = LoggerFactory.getLogger(NCGenesBaselineMemServiceImpl.class);

    private SampleDAO sampleDAO;

    private WorkflowDAO workflowDAO;

    @Override
    public QualityControlInfo lookupQuantificationResults(Long sampleId) {
        logger.debug("ENTERING lookupQuantificationResults(Long)");
        if (sampleId == null) {
            logger.warn("sampleId is null");
            return null;
        }

        Sample sample = null;
        try {
            sample = sampleDAO.findById(sampleId);
            logger.info(sample.toString());
        } catch (MaPSeqDAOException e) {
            logger.error("Failed to find Sample", e);
        }

        if (sample == null) {
            return null;
        }

        logger.debug(sample.toString());

        Set<FileData> sampleFileDataSet = sample.getFileDatas();

        QualityControlInfo ret = new QualityControlInfo();

        try {

            if (sampleFileDataSet != null) {

                File flagstatFile = fileToFind(sample, MimeType.TEXT_STAT_SUMMARY, "samtools.flagstat");

                if (flagstatFile == null) {
                    logger.error("flagstat file to process was still not found");
                    return ret;
                }

                logger.info("flagstat file is: {}", flagstatFile.getAbsolutePath());
                if (flagstatFile.exists()) {
                    List<String> lines = null;
                    try {
                        lines = FileUtils.readLines(flagstatFile);
                    } catch (IOException e1) {
                        e1.printStackTrace();
                    }

                    if (lines != null) {
                        for (String line : lines) {

                            if (line.contains("in total")) {
                                String value = line.substring(0, line.indexOf(" ")).trim();
                                try {
                                    ret.setPassedReads(Integer.valueOf(value));
                                } catch (Exception e) {
                                    logger.error("problem getting passedReads, value: {}", value);
                                }
                            }

                            if (line.contains("mapped (")) {
                                Pattern pattern = Pattern.compile("^.+\\((.+)\\)");
                                Matcher matcher = pattern.matcher(line);
                                if (matcher.matches()) {
                                    String value = matcher.group(1);
                                    value = value.substring(0, value.indexOf("%")).trim();
                                    if (StringUtils.isNotEmpty(value)) {
                                        try {
                                            ret.setAligned(Float.valueOf(value));
                                        } catch (Exception e) {
                                            logger.error("problem getting mapped, value: {}", value);
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
                                        try {
                                            ret.setPaired(Float.valueOf(value));
                                        } catch (Exception e) {
                                            logger.error("problem getting paired, value: {}", value);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                File vcfFile = fileToFind(sample, MimeType.TEXT_VCF, "ic_snps.vcf");

                if (vcfFile == null) {
                    logger.error("vcf file to process was still not found");
                    logger.info(sample.toString());
                    return ret;
                }

                logger.info("vcfFile file is: {}", vcfFile.getAbsolutePath());
                if (vcfFile.exists()) {
                    VCFParser parser = VCFParser.getInstance();
                    VCFResult results = parser.parse(vcfFile);
                    ret.setIcSNPResultList(results);
                }

                File depthOfCoverageSummaryFile = fileToFind(sample, MimeType.TEXT_DEPTH_OF_COVERAGE_SUMMARY, "coverage.sample_summary");

                if (depthOfCoverageSummaryFile == null) {
                    logger.error("depthOfCoverageSummaryFile to process was still not found");
                    logger.info(sample.toString());
                    return ret;
                }

                logger.info("depthOfCoverageSummaryFile file is: {}", depthOfCoverageSummaryFile.getAbsolutePath());

                if (depthOfCoverageSummaryFile.exists()) {
                    List<String> lines = FileUtils.readLines(depthOfCoverageSummaryFile);
                    for (String line : lines) {
                        if (line.contains("Total")) {
                            String[] split = line.split("\t");
                            ret.setTotalCoverage(Long.valueOf(split[1]));
                            ret.setMean(Double.valueOf(split[2]));
                        }
                    }
                }

            }
        } catch (Exception e) {
            logger.error(e.getMessage(), e);
        }

        return ret;
    }

    private File fileToFind(Sample sample, MimeType mimeType, String suffix) throws MaPSeqDAOException {
        Set<FileData> sampleFileDataSet = sample.getFileDatas();

        Workflow workflow = null;
        List<Workflow> workflowList = workflowDAO.findByName("NCGenesBaselineMem");
        if (CollectionUtils.isNotEmpty(workflowList)) {
            workflow = workflowList.get(0);
        }

        if (CollectionUtils.isNotEmpty(sampleFileDataSet)) {
            for (FileData fileData : sampleFileDataSet) {
                if (fileData.getName().endsWith(suffix) && mimeType.equals(fileData.getMimeType())) {
                    return fileData.toFile();
                }
            }
        }

        File outputDirectory = SequencingWorkflowUtil.createOutputDirectory(sample, workflow);
        if (outputDirectory.exists()) {
            for (File file : outputDirectory.listFiles()) {
                if (file.getName().endsWith(suffix)) {
                    return file;
                }
            }
        }

        return null;
    }

    @Override
    public VCFResult lookupIdentityInfoFromVCF(Long sampleId) {
        logger.debug("ENTERING lookupIdentityInfoFromVCF(Long)");
        if (sampleId == null) {
            logger.warn("sampleId is null");
            return null;
        }

        Sample sample = null;
        try {
            sample = sampleDAO.findById(sampleId);
            logger.info(sample.toString());
        } catch (MaPSeqDAOException e) {
            logger.error("Failed to find Sample", e);
        }

        if (sample == null) {
            return null;
        }

        Set<FileData> sampleFileDataSet = sample.getFileDatas();

        VCFParser parser = VCFParser.getInstance();
        VCFResult ret = null;
        if (sampleFileDataSet != null) {
            for (FileData fileData : sampleFileDataSet) {
                if (MimeType.TEXT_VCF.equals(fileData.getMimeType()) && fileData.getName().endsWith(".ic_snps.vcf")) {
                    File icSNPVCFFile = new File(fileData.getPath(), fileData.getName());
                    logger.info("icSNPVCFFile file is: {}", icSNPVCFFile.getAbsolutePath());
                    if (icSNPVCFFile.exists()) {
                        ret = parser.parse(icSNPVCFFile);
                    }
                }
            }
        }
        return ret;
    }

    public SampleDAO getSampleDAO() {
        return sampleDAO;
    }

    public void setSampleDAO(SampleDAO sampleDAO) {
        this.sampleDAO = sampleDAO;
    }

    public WorkflowDAO getWorkflowDAO() {
        return workflowDAO;
    }

    public void setWorkflowDAO(WorkflowDAO workflowDAO) {
        this.workflowDAO = workflowDAO;
    }

}
