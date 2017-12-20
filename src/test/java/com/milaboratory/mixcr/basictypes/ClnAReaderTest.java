package com.milaboratory.mixcr.basictypes;

import cc.redberry.pipe.CUtils;
import cc.redberry.pipe.OutputPort;
import com.milaboratory.mixcr.assembler.AlignmentsMappingMerger;
import com.milaboratory.mixcr.util.RunMiXCR;
import com.milaboratory.util.TempFileManager;
import io.repseq.core.VDJCLibraryRegistry;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertEquals;

public class ClnAReaderTest {

    @Test
    public void test4() throws Exception {
        RunMiXCR.RunMiXCRAnalysis params = new RunMiXCR.RunMiXCRAnalysis(
                RunMiXCR.class.getResource("/sequences/test_R1.fastq").getFile(),
                RunMiXCR.class.getResource("/sequences/test_R2.fastq").getFile());

//        RunMiXCR.RunMiXCRAnalysis params = new RunMiXCR.RunMiXCRAnalysis(
//                "/Users/poslavskysv/Projects/milab/temp/clean2/synth_R1.fastq",
//                "/Users/poslavskysv/Projects/milab/temp/clean2/synth_R2.fastq");


        params.cloneAssemblerParameters.setAddReadsCountOnClustering(true);
        RunMiXCR.AlignResult align = RunMiXCR.align(params);
        RunMiXCR.AssembleResult assemble = RunMiXCR.assemble(align, false);

        AlignmentsMappingMerger merged = new AlignmentsMappingMerger(align.resultReader(), assemble.cloneAssembler.getAssembledReadsPort());

        File file = TempFileManager.getTempFile();
        ClnAWriter writer = new ClnAWriter(file, align.parameters.alignerParameters);
        writer.writeClones(assemble.cloneSet);
        writer.sortAlignments(merged, align.alignments.size());
        writer.writeAlignmentsAndIndex();

        writer.close();

        ClnAReader reader = new ClnAReader(file.toPath(), VDJCLibraryRegistry.createDefaultRegistry());

        assertEquals(assemble.cloneSet.size(), reader.numberOfClones());

        for (ClnAReader.CloneAlignments c : CUtils.it(reader.clonesAndAlignments())) {
            assertEquals("" + c.cloneId, c.clone.count, count(c.alignments()));
            CUtils.it(c.alignments()).forEach(a -> assertEquals(c.cloneId, a.getCloneIndex()));
        }
    }

    static int count(OutputPort port) {
        int c = 0;
        while (port.take() != null)
            ++c;
        return c;
    }
}