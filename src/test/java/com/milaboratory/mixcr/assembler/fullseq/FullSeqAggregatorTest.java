package com.milaboratory.mixcr.assembler.fullseq;

import cc.redberry.pipe.CUtils;
import com.milaboratory.core.io.sequence.PairedRead;
import com.milaboratory.core.io.sequence.SingleReadImpl;
import com.milaboratory.core.sequence.NSequenceWithQuality;
import com.milaboratory.core.sequence.NucleotideSequence;
import com.milaboratory.mixcr.util.RunMiXCR;
import com.milaboratory.mixcr.vdjaligners.VDJCParametersPresets;
import gnu.trove.set.hash.TIntHashSet;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.stream.StreamSupport;

/**
 *
 */
public class FullSeqAggregatorTest {

    static final class MasterSequence {
        final int vPart, cdr3Part, jPart, cPart;
        final NucleotideSequence masterSequence;

        MasterSequence(NucleotideSequence vPart, NucleotideSequence cdr3Part, NucleotideSequence jPart, NucleotideSequence cPart) {
            this.vPart = vPart.size();
            this.cdr3Part = cdr3Part.size();
            this.jPart = jPart.size();
            this.cPart = cPart.size();
            this.masterSequence = vPart.concatenate(cdr3Part).concatenate(jPart).concatenate(cPart);
        }

        MasterSequence(String vPart, String cdr3Part, String jPart, String cPart) {
            this(new NucleotideSequence(vPart), new NucleotideSequence(cdr3Part),
                    new NucleotideSequence(jPart), new NucleotideSequence(cPart));
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
            "GGCCAGGGAACCCTGGTCCGTCTCCTCAG",
            "CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAACCGGTGACGGTGTCGTGGAACTCAGGCGCTCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCAGCGTGGTGACCGTGCCCTCCAGCAACTTCGGCACCCAGACCTACACCTGCAACGTAGATCACAAGCCCAGCAACACCAAGGTGGACAAGACAGTTGGTGAGAGGCCAGCTCAGGGAGGGAGGGTGTCTGCTGGAAGCCAGGCTCAGCCCTCCTGCCTGGACGCACCCCGGCTGTGCAGCCCCAGCCCAGGGCAGCAAGGCAGGCCCCATCTGTCTCCTCACCCGGAGGCCTCTGCCCGCCCCACTCATGCTCAGGGAGAGGGTCTTCTGGCTTTTTCCACCAGGCTCCAGGCAGGCACAGGCTGGGTGCCCCTACCCCAGGCCCTTCACACACAGGGGCAGGTGCTTGGCTCAGACCTGCCAAAAGCCATATCCGG");

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

        // [-200, -60]  [-20, 120]
        //      [-150, 110]
        //
        // [-200, -150], [110, 120] = 60
        // [-60, -20] = 40
        params.alignerParameters = VDJCParametersPresets.getByName("rna-seq");
        params.alignerParameters.setSaveOriginalReads(true);

        RunMiXCR.AlignResult align = RunMiXCR.align(params);
//        for (VDJCAlignments al : align.alignments) {
//            for (int i = 0; i < al.numberOfTargets(); i++)
//                System.out.println(VDJCAlignmentsFormatter.getTargetAsMultiAlignment(al, i));
//            System.out.println();
//        }

        RunMiXCR.AssembleResult assemble = RunMiXCR.assemble(align);

        FullSeqAggregator agg = new FullSeqAggregator(assemble.cloneSet.get(0));

        PointSequence[] r2s = agg.convertToPointSequences(align.alignments.get(1));
        TIntHashSet p2 = new TIntHashSet(Arrays.stream(r2s).mapToInt(s -> s.point).toArray());
        Assert.assertEquals(260 - masterSeq1.cdr3Part, p2.size());

        PointSequence[] r1s = agg.convertToPointSequences(align.alignments.get(0));
        TIntHashSet p1 = new TIntHashSet(Arrays.stream(r1s).mapToInt(s -> s.point).toArray());
        Assert.assertEquals(280 - masterSeq1.cdr3Part, p1.size());

        FullSeqAggregator.PreparedData prep = agg.initialRound(() -> CUtils.asOutputPort(align.alignments));

        long uniq1 = StreamSupport.stream(CUtils.it(prep.createPort()).spliterator(), false)
                .mapToInt(l -> l[0])
                .filter(c -> c == 0xFFFFFFFF).count();
        long uniq2 = StreamSupport.stream(CUtils.it(prep.createPort()).spliterator(), false)
                .mapToInt(l -> l[1])
                .filter(c -> c == 0xFFFFFFFF).count();

        Assert.assertEquals(40, uniq1);
        Assert.assertEquals(60, uniq2);


//        int i = 0;
//        for (int[] ints : CUtils.it(prep.createPort())) {
//            System.out.println(Arrays.stream(ints).mapToObj(s -> String.format("%08X", s)).collect(Collectors.joining(", ")) + "    ->" + prep.points[i]);
//            ++i;
//        }
    }
}