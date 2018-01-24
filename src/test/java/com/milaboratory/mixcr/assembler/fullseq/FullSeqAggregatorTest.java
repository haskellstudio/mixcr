package com.milaboratory.mixcr.assembler.fullseq;

import cc.redberry.pipe.CUtils;
import com.milaboratory.core.io.sequence.PairedRead;
import com.milaboratory.core.io.sequence.SequenceRead;
import com.milaboratory.core.io.sequence.SequenceReaderCloseable;
import com.milaboratory.core.io.sequence.SingleReadImpl;
import com.milaboratory.core.sequence.NSequenceWithQuality;
import com.milaboratory.core.sequence.NucleotideSequence;
import com.milaboratory.mixcr.basictypes.Clone;
import com.milaboratory.mixcr.basictypes.CloneSet;
import com.milaboratory.mixcr.basictypes.VDJCAlignments;
import com.milaboratory.mixcr.basictypes.VDJCAlignmentsFormatter;
import com.milaboratory.mixcr.cli.ActionExportClonesPretty;
import com.milaboratory.mixcr.util.RunMiXCR;
import com.milaboratory.mixcr.vdjaligners.VDJCParametersPresets;
import gnu.trove.set.hash.TIntHashSet;
import io.repseq.core.GeneFeature;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
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

        NucleotideSequence getRangeFromCDR3Begin(int vPadd, int len) {
            return masterSequence.getRange(vPart + vPadd, vPart + vPadd + len);
        }

        NucleotideSequence getRangeFromCDR3End(int jPadd, int len) {
            return masterSequence.getRange(vPart + cdr3Part + jPadd, vPart + cdr3Part + jPadd + len);
        }
    }

    static final MasterSequence masterSeq1WT = new MasterSequence(
            "CTGAAGAAAACCAGCCCTGCAGCTCTGGGAGAGGAGCCCCAGCCCTGGGATTCCCAGCTGTTTCTGCTTGCTGATCAGGACTGCACACAGAGAACTCACC" +
                    "ATGGAGTTTGGGCTGAGCTGGGTTTTCCTTGTTGCTATTTTAAAAGGTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCCGGGGGAGGCTTAGTTCAGCC" +
                    "TGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTACTGGATGCACTGGGTCCGCCAAGCTCCAGGGAAGGGGCTGGTGT" +
                    "GGGTCTCACGTATTAATAGTGATGGGAGTAGCACAAGCTACGCGGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTG" +
                    "TATCTGCAAATGAACAGTCTGAGAGCCGAGGACACGGCTGTGTATTAC",
            "TGTGCAAGAGGGCCCCAAGAAAATAGTGGTTATTACTACGGGTTTGACTACTGG",
            "GGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
            "CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCC" +
                    "CGAACCGGTGACGGTGTCGTGGAACTCAGGCGCTCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCA" +
                    "GCGTGGTGACCGTGCCCTCCAGCAACTTCGGCACCCAGACCTACACCTGCAACGTAGATCACAAGCCCAGCAACACCAAGGTGGACAAGACAGTTGGT" +
                    "GAGAGGCCAGCTCAGGGAGGGAGGGTGTCTGCTGGAAGCCAGGCTCAGCCCTCCTGCCTGGACGCACCCCGGCTGTGCAGCCCCAGCCCAGGGCAGCA" +
                    "AGGCAGGCCCCATCTGTCTCCTCACCCGGAGGCCTCTGCCCGCCCCACTCATGCTCAGGGAGAGGGTCTTCTGGCTTTTTCCACCAGGCTCCAGGCAG" +
                    "GCACAGGCTGGGTGCCCCTACCCCAGGCCCTTCACACACAGGGGCAGGTGCTTGGCTCAGACCTGCCAAAAGCCATATCCGG");

    static final MasterSequence masterSeq1VIns1 = new MasterSequence(
            "CTGAAGAAAACCAGCCCTGCAGCTCTGGGAGAGGAGCCCCAGCCCTGGGATTCCCAGCTGTTTCTGCTTGCTGATCAGGACTGCACACAGAGAACTCACC" +
                    "ATGGAGTTTGGGCTGAGCTGGGTTTTCCTTGTTGCTATTTTAAAAGGTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCCGGGGGAGGCTTAGTTCAGCC" +
                    "TGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTACTGGATGCACTGGGTCCGCCAAGCTCCAGGGAAGGGGCTGGTGT" +
                    "GGGTCTCACGTATTAATAGTGATGGGAGTAGCACAAGCTtattACGCGGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTG" +
                    "TATCTGCAAATGAACAGTCTGAGAGCCGAGGACACGGCTGTGTATTAC",
            "TGTGCAAGAGGGCCCCAAGAAAATAGTGGTTATTACTACGGGTTTGACTACTGG",
            "GGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
            "CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCC" +
                    "CGAACCGGTGACGGTGTCGTGGAACTCAGGCGCTCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCA" +
                    "GCGTGGTGACCGTGCCCTCCAGCAACTTCGGCACCCAGACCTACACCTGCAACGTAGATCACAAGCCCAGCAACACCAAGGTGGACAAGACAGTTGGT" +
                    "GAGAGGCCAGCTCAGGGAGGGAGGGTGTCTGCTGGAAGCCAGGCTCAGCCCTCCTGCCTGGACGCACCCCGGCTGTGCAGCCCCAGCCCAGGGCAGCA" +
                    "AGGCAGGCCCCATCTGTCTCCTCACCCGGAGGCCTCTGCCCGCCCCACTCATGCTCAGGGAGAGGGTCTTCTGGCTTTTTCCACCAGGCTCCAGGCAG" +
                    "GCACAGGCTGGGTGCCCCTACCCCAGGCCCTTCACACACAGGGGCAGGTGCTTGGCTCAGACCTGCCAAAAGCCATATCCGG");

    static final MasterSequence masterSeq1VSub1 = new MasterSequence(
            "CTGAAGAAAACCAGCCCTGCAGCTCTGGGAGAGGAGCCCCAGCCCTGGGATTCCCAGCTGTTTCTGCTTGCTGATCAGGACTGCACACAGAGAACTCACC" +
                    "ATGGAGTTTGGGCTGAGCTGGGTTTTCCTTGTTGCTATTTTAAAAGGTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCCGGGGGAGGCTTAGTTCAGCC" +
                    "TGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTACTGGATGCACTGGGTCCGCCAAGCTCCAGGGAAGGGGCTGGTGT" +
                    "GGGTCTCACGTATTAATAGTGATGGGAGTAGCACAAGCTACGCGGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTG" +
                    "TATCTGCAAATGAACAGTCTcAGAGCCGAGGACACGGCTGTGTATTAC",
            "TGTGCAAGAGGGCCCCAAGAAAATAGTGGTTATTACTACGGGTTTGACTACTGG",
            "GGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
            "CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCC" +
                    "CGAACCGGTGACGGTGTCGTGGAACTCAGGCGCTCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCA" +
                    "GCGTGGTGACCGTGCCCTCCAGCAACTTCGGCACCCAGACCTACACCTGCAACGTAGATCACAAGCCCAGCAACACCAAGGTGGACAAGACAGTTGGT" +
                    "GAGAGGCCAGCTCAGGGAGGGAGGGTGTCTGCTGGAAGCCAGGCTCAGCCCTCCTGCCTGGACGCACCCCGGCTGTGCAGCCCCAGCCCAGGGCAGCA" +
                    "AGGCAGGCCCCATCTGTCTCCTCACCCGGAGGCCTCTGCCCGCCCCACTCATGCTCAGGGAGAGGGTCTTCTGGCTTTTTCCACCAGGCTCCAGGCAG" +
                    "GCACAGGCTGGGTGCCCCTACCCCAGGCCCTTCACACACAGGGGCAGGTGCTTGGCTCAGACCTGCCAAAAGCCATATCCGG");

    static final MasterSequence masterSeq1VIns1Sub1 = new MasterSequence(
            "CTGAAGAAAACCAGCCCTGCAGCTCTGGGAGAGGAGCCCCAGCCCTGGGATTCCCAGCTGTTTCTGCTTGCTGATCAGGACTGCACACAGAGAACTCACC" +
                    "ATGGAGTTTGGGCTGAGCTGGGTTTTCCTTGTTGCTATTTTAAAAGGTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCCGGGGGAGGCTTAGTTCAGCC" +
                    "TGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTACTGGATGCACTGGGTCCGCCAAGCTCCAGGGAAGGGGCTGGTGT" +
                    "GGGTCTCACGTATTAATAGTGATGGGAGTAGCACAAGCTtattACGCGGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTG" +
                    "TATCTGCAAATGAACAGTCTcAGAGCCGAGGACACGGCTGTGTATTAC",
            "TGTGCAAGAGGGCCCCAAGAAAATAGTGGTTATTACTACGGGTTTGACTACTGG",
            "GGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
            "CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCC" +
                    "CGAACCGGTGACGGTGTCGTGGAACTCAGGCGCTCTGACCAGCGGCGTGCACACCTTCCCGGCTGTCCTACAGTCCTCAGGACTCTACTCCCTCAGCA" +
                    "GCGTGGTGACCGTGCCCTCCAGCAACTTCGGCACCCAGACCTACACCTGCAACGTAGATCACAAGCCCAGCAACACCAAGGTGGACAAGACAGTTGGT" +
                    "GAGAGGCCAGCTCAGGGAGGGAGGGTGTCTGCTGGAAGCCAGGCTCAGCCCTCCTGCCTGGACGCACCCCGGCTGTGCAGCCCCAGCCCAGGGCAGCA" +
                    "AGGCAGGCCCCATCTGTCTCCTCACCCGGAGGCCTCTGCCCGCCCCACTCATGCTCAGGGAGAGGGTCTTCTGGCTTTTTCCACCAGGCTCCAGGCAG" +
                    "GCACAGGCTGGGTGCCCCTACCCCAGGCCCTTCACACACAGGGGCAGGTGCTTGGCTCAGACCTGCCAAAAGCCATATCCGG");

    @Test
    public void testRandom1() throws Exception {
        CloneFraction[] clones = {
                // new CloneFraction(16, masterSeq1WT),
                // new CloneFraction(32, masterSeq1VSub1),
                new CloneFraction(64, masterSeq1VIns1),
                // new CloneFraction(128, masterSeq1VIns1Sub1),
        };

        Well19937c rand = new Well19937c();
        RandomDataGenerator rdg = new RandomDataGenerator(rand);

        List<PairedRead> readsOrig = new ArrayList<>();

        int readLength = 100;

        int id = -1;

        for (CloneFraction clone : clones) {
            for (int i = 0; i < clone.count; i++) {
                // Left read with CDR3
                ++id;
                readsOrig.add(new PairedRead(
                        new SingleReadImpl(id, new NSequenceWithQuality(clone.seq.getRangeFromCDR3Begin(-rand.nextInt(readLength - clone.seq.cdr3Part), readLength)), "R1_" + id),
                        new SingleReadImpl(id, new NSequenceWithQuality(clone.seq.getRangeFromCDR3End(rdg.nextInt(-clone.seq.cdr3Part / 2, clone.seq.jPart),
                                readLength).getReverseComplement()), "R2_" + id)));

                ++id;
                readsOrig.add(new PairedRead(
                        new SingleReadImpl(id, new NSequenceWithQuality(clone.seq.getRangeFromCDR3Begin(rdg.nextInt(-clone.seq.vPart, clone.seq.cdr3Part / 2 - readLength),
                                readLength)), "R1_" + id),
                        new SingleReadImpl(id, new NSequenceWithQuality(clone.seq.getRangeFromCDR3Begin(-rand.nextInt(readLength - clone.seq.cdr3Part),
                                readLength)).getReverseComplement(), "R2_" + id)));
            }
        }

        int[] perm = rdg.nextPermutation(readsOrig.size(), readsOrig.size());
        List<PairedRead> reads = new ArrayList<>();
        for (int i = 0; i < readsOrig.size(); i++)
            reads.add(readsOrig.get(perm[i]));

        RunMiXCR.RunMiXCRAnalysis params = new RunMiXCR.RunMiXCRAnalysis(
                new SequenceReaderCloseable<SequenceRead>() {
                    int counter = 0;

                    @Override
                    public void close() {
                    }

                    @Override
                    public long getNumberOfReads() {
                        return counter;
                    }

                    @Override
                    public synchronized SequenceRead take() {
                        if (counter == reads.size())
                            return null;
                        return reads.get(counter++);
                    }
                }, true);

        params.alignerParameters = VDJCParametersPresets.getByName("rna-seq");
        params.alignerParameters.setSaveOriginalReads(true);

        RunMiXCR.AlignResult align = RunMiXCR.align(params);

        // TODO exception for translation
        for (VDJCAlignments al : align.alignments) {
            for (int i = 0; i < al.numberOfTargets(); i++) {
                System.out.println(VDJCAlignmentsFormatter.getTargetAsMultiAlignment(al, i));
                System.out.println();
            }
            System.out.println();
            System.out.println(" ================================================ ");
            System.out.println();
        }


        RunMiXCR.AssembleResult assemble = RunMiXCR.assemble(align);

        FullSeqAggregator agg = new FullSeqAggregator(assemble.cloneSet.get(0), align.parameters.alignerParameters);

        FullSeqAggregator.PreparedData prep = agg.initialRound(() -> CUtils.asOutputPort(
                align.alignments.stream().filter(a -> a.getFeature(GeneFeature.CDR3) != null).collect(Collectors.toList())
        ));

        for (Clone clone : new CloneSet(Arrays.asList(agg.process(prep))).getClones())
            ActionExportClonesPretty.outputCompact(System.out, clone);
    }

    public static class CloneFraction {
        final int count;
        final MasterSequence seq;

        public CloneFraction(int count, MasterSequence seq) {
            this.count = count;
            this.seq = seq;
        }
    }

    @Test
    public void test1() throws Exception {
        int len = 140;
        PairedRead read1 = new PairedRead(
                new SingleReadImpl(0, new NSequenceWithQuality(masterSeq1WT.getRangeFromCDR3Begin(-20, len)), "R1"),
                new SingleReadImpl(0, new NSequenceWithQuality(masterSeq1WT.getRangeFromCDR3Begin(-200, len).getReverseComplement()), "R2"));

        PairedRead read2 = new PairedRead(
                new SingleReadImpl(1, new NSequenceWithQuality(masterSeq1WT.getRangeFromCDR3Begin(-30, len)), "R1"),
                new SingleReadImpl(1, new NSequenceWithQuality(masterSeq1WT.getRangeFromCDR3Begin(-150, len).getReverseComplement()), "R2"));

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

        FullSeqAggregator agg = new FullSeqAggregator(assemble.cloneSet.get(0), align.parameters.alignerParameters);

        PointSequence[] r2s = agg.convertToPointSequences(align.alignments.get(1));
        TIntHashSet p2 = new TIntHashSet(Arrays.stream(r2s).mapToInt(s -> s.point).toArray());
        Assert.assertEquals(260 - masterSeq1WT.cdr3Part, p2.size());

        PointSequence[] r1s = agg.convertToPointSequences(align.alignments.get(0));
        TIntHashSet p1 = new TIntHashSet(Arrays.stream(r1s).mapToInt(s -> s.point).toArray());
        Assert.assertEquals(280 - masterSeq1WT.cdr3Part, p1.size());

        FullSeqAggregator.PreparedData prep = agg.initialRound(() -> CUtils.asOutputPort(align.alignments));

        long uniq1 = StreamSupport.stream(CUtils.it(prep.createPort()).spliterator(), false)
                .mapToInt(l -> l[0])
                .filter(c -> c == 0xFFFFFFFF).count();
        long uniq2 = StreamSupport.stream(CUtils.it(prep.createPort()).spliterator(), false)
                .mapToInt(l -> l[1])
                .filter(c -> c == 0xFFFFFFFF).count();

        Assert.assertEquals(40, uniq1);
        Assert.assertEquals(60, uniq2);

        for (Clone clone : new CloneSet(Arrays.asList(agg.process(prep))).getClones()) {
            ActionExportClonesPretty.outputCompact(System.out, clone);
            System.out.println();
            System.out.println(" ================================================ ");
            System.out.println();
        }
    }
}