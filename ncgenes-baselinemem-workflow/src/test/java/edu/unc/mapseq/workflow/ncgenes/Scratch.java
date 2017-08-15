package edu.unc.mapseq.workflow.ncgenes;

import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

public class Scratch {

    @Test
    public void testSampleNameExternalCodeParsing() {
        String value = "MCF-7-SW2";
        int idx = value.lastIndexOf("-");
        assertTrue(value.substring(0, idx).equals("MCF-7"));

        value = "PEROU_LAB_REF-SW2";
        idx = value.lastIndexOf("-");
        assertTrue(value.substring(0, idx).equals("PEROU_LAB_REF"));

        value = "NCG_00007-LT_1";
        idx = value.lastIndexOf("-");
        assertTrue(value.substring(0, idx).equals("NCG_00007"));

        value = "NCG_00007_LT_1";
        idx = value.lastIndexOf("-");
        assertTrue(idx == -1);

    }

    @Test
    public void scratch() {
        try {
            List<String> lines = FileUtils.readLines(new File("/tmp/asdf.txt"));
            double total = 0;
            for (String line : lines) {
                total += Double.valueOf(line.trim());
            }
            System.out.println(total / 1000);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    @Test
    public void splitTest() {
        try {
            List<String> lines = FileUtils.readLines(new File("/tmp/asdf.txt"));
            assertTrue(Integer.valueOf(lines.iterator().next().split("\\t")[1]) == 1);
            assertTrue(Integer.valueOf(StringUtils.split(lines.iterator().next())[2]) == 2);
        } catch (NumberFormatException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void calculateNumberOnTarget() {

        File docFile = new File(
                "/home/jdr0887/tmp/120517_UNC15-SN850_0212_AD0VD1ACXX_ACAGTG_L004.fixed-rg.deduped.realign.fixmate.recal.coverage.sample_interval_summary");
        try {

            BufferedReader br = new BufferedReader(new FileReader(docFile));
            String line;
            long totalCoverageCount = 0;
            br.readLine();
            while ((line = br.readLine()) != null) {
                totalCoverageCount += Long.valueOf(StringUtils.split(line)[1].trim());
            }
            br.close();
            System.out.println((double) totalCoverageCount / (52983590 * 100));
        } catch (IOException e1) {
            e1.printStackTrace();
        }
    }

}
