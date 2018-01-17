package com.milaboratory.mixcr.assembler.fullseq;

import cc.redberry.pipe.CUtils;
import com.milaboratory.core.io.sequence.PairedRead;
import com.milaboratory.core.io.sequence.SingleReadImpl;
import com.milaboratory.core.sequence.NSequenceWithQuality;
import com.milaboratory.core.sequence.NucleotideSequence;
import com.milaboratory.mixcr.util.RunMiXCR;
import com.milaboratory.mixcr.vdjaligners.VDJCParametersPresets;
import org.junit.Test;

import java.util.Arrays;
import java.util.stream.Collectors;

/**
 *
 */
public class FullSeqAggregatorTest {

    static final class MasterSequence {
        final int vPart, cdr3Part, jPart;
        final NucleotideSequence masterSequence;

        MasterSequence(NucleotideSequence vPart, NucleotideSequence cdr3Part, NucleotideSequence jPart) {
            this.vPart = vPart.size();
            this.cdr3Part = cdr3Part.size();
            this.jPart = jPart.size();
            this.masterSequence = vPart.concatenate(cdr3Part).concatenate(jPart);
        }

        MasterSequence(String vPart, String cdr3Part, String jPart) {
            this(new NucleotideSequence(vPart), new NucleotideSequence(cdr3Part), new NucleotideSequence(jPart));
        }

        NucleotideSequence getRange(int vPadd, int jPadd) {
            return masterSequence.getRange(vPart + vPadd, vPart + cdr3Part + jPadd);
        }

        NucleotideSequence getRangeFromV(int vPadd, int len) {
            return masterSequence.getRange(vPart + vPadd, vPart + vPadd + len);
        }

        NucleotideSequence getRangeFromJ(int jPadd, int len) {
            return masterSequence.getRange(vPart + cdr3Part + jPadd, vPart + cdr3Part + jPadd + len);
        }
    }

    static final MasterSequence masterSeq1 = new MasterSequence(
            "CTGAAGAAAACCAGCCCTGCAGCTCTGGGAGAGGAGCCCCAGCCCTGGGATTCCCAGCTGTTTCTGCTTGCTGATCAGGACTGCACACAGAGAACTCACCATGGAGTTTGGGCTGAGCTGGGTTTTCCTTGTTGCTATTTTAAAAGGTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCCGGGGGAGGCTTAGTTCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTACTGGATGCACTGGGTCCGCCAAGCTCCAGGGAAGGGGCTGGTGTGGGTCTCACGTATTAATAGTGATGGGAGTAGCACAAGCTACGCGGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTGTATCTGCAAATGAACAGTCTGAGAGCCGAGGACACGGCTGTGTATTTTAC",
            "TGTGCAAGAGGGCCCCAAGAAAATAGTGGTTATTACTACGGGTTTGACTACTGG",
            "GGCCAGGGAACCCTGGTCCGTCTCCTCAG");

    @Test
    public void test1() throws Exception {

        int len = 140;
        PairedRead read1 = new PairedRead(
                new SingleReadImpl(0, new NSequenceWithQuality(masterSeq1.getRangeFromV(-20, len)), "R1"),
                new SingleReadImpl(0, new NSequenceWithQuality(masterSeq1.getRangeFromV(-200, len).getReverseComplement()), "R2"));

        PairedRead read2 = new PairedRead(
                new SingleReadImpl(1, new NSequenceWithQuality(masterSeq1.getRangeFromV(-30, len)), "R1"),
                new SingleReadImpl(1, new NSequenceWithQuality(masterSeq1.getRangeFromV(-150, len).getReverseComplement()), "R2"));

        RunMiXCR.RunMiXCRAnalysis params = new RunMiXCR.RunMiXCRAnalysis(read1, read2);

        params.alignerParameters = VDJCParametersPresets.getByName("rna-seq");
        params.alignerParameters.setSaveOriginalReads(true);

        RunMiXCR.AlignResult align = RunMiXCR.align(params);
        RunMiXCR.AssembleResult assemble = RunMiXCR.assemble(align);

        FullSeqAggregator agg = new FullSeqAggregator(assemble.cloneSet.get(0));
        PointSequence[] sss = agg.convertToPointSequences(align.alignments.get(0));
        System.out.println(Arrays.toString(Arrays.stream(sss).filter(p -> p.sequence.size() != 1).toArray()));

        FullSeqAggregator.PreparedData prep = agg.initialRound(() -> CUtils.asOutputPort(align.alignments));
        System.out.println(prep);
        for (int[] ints : CUtils.it(prep.createPort())) {
            System.out.println(Arrays.stream(ints).mapToObj(s -> String.format("%08X", s)).collect(Collectors.joining(", ")));
        }
//        for (Clone clone : assemble.cloneSet.getClones()) {
//            System.out.println(clone);
//        }
    }
}