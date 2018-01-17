package com.milaboratory.mixcr.assembler.fullseq;

import cc.redberry.pipe.CUtils;
import cc.redberry.pipe.OutputPort;
import com.milaboratory.core.Range;
import com.milaboratory.core.alignment.Alignment;
import com.milaboratory.core.sequence.NSequenceWithQuality;
import com.milaboratory.core.sequence.NucleotideSequence;
import com.milaboratory.mixcr.basictypes.Clone;
import com.milaboratory.mixcr.basictypes.VDJCAlignments;
import com.milaboratory.mixcr.basictypes.VDJCHit;
import com.milaboratory.mixcr.basictypes.VDJCPartitionedSequence;
import gnu.trove.impl.Constants;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import io.repseq.core.GeneFeature;
import io.repseq.core.ReferencePoint;
import io.repseq.core.VDJCGene;
import io.repseq.gen.VDJCGenes;

import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static io.repseq.core.GeneType.Joining;
import static io.repseq.core.GeneType.Variable;

/**
 */
final class FullSeqAggregator {
    final Clone clone;
    final GeneFeature assemblingFeature;
    final VDJCGenes genes;
    final boolean hasV, hasJ;
    final int nLeftDummies, lengthV, assemblingFeatureLength, jOffset;

    final double relativeSumQualityThreshold = 0.1;
    final long absoluteSumQualityThreshold = 120;

    FullSeqAggregator(Clone clone) {
        this.clone = clone;
        GeneFeature[] assemblingFeatures = clone.getAssemblingFeatures();
        if (assemblingFeatures.length != 1)
            throw new IllegalArgumentException();

        if (assemblingFeatures[0].isComposite())
            throw new IllegalArgumentException();

        this.assemblingFeature = assemblingFeatures[0];
        this.genes = clone.getBestHitGenes();

        ReferencePoint
                start = assemblingFeature.getFirstPoint(),
                end = assemblingFeature.getLastPoint();

        this.hasV = start.getGeneType() == Variable;
        this.hasJ = end.getGeneType() == Joining;


        //  nLeftDummies     assemblingFeatureLength
        //  ------|--------------|--------------|------------------------>
        //        ↓              ↓              ↓
        //  0000000vvvvvvvvvvvvvvCDR3CDR3CDR3CDR3jjjjjjjjjjjjjjjjCCCCCCCCC

        this.nLeftDummies = 1024; // fixme

        if (hasV) {
            VDJCHit vHit = clone.getBestHit(Variable);
            GeneFeature vFeature = vHit.getAlignedFeature();
            VDJCGene gene = vHit.getGene();
            this.lengthV =
                    gene.getFeature(vFeature).size()
                            - gene.getFeature(GeneFeature.intersection(assemblingFeature, vFeature)).size();
        } else
            this.lengthV = 0;

        this.assemblingFeatureLength = clone.getFeature(assemblingFeature).size();

        if (hasJ) {
            VDJCHit jHit = clone.getBestHit(Joining);
            VDJCGene gene = jHit.getGene();
            this.jOffset = gene.getPartitioning().getPosition(assemblingFeature.getLastPoint());
        } else
            this.jOffset = 0;
    }

    final TObjectIntHashMap<NucleotideSequence> possibleSequences
            = new TObjectIntHashMap<>(Constants.DEFAULT_CAPACITY, Constants.DEFAULT_LOAD_FACTOR, -1);

    PreparedData initialRound(Supplier<OutputPort<VDJCAlignments>> alignments) {
        for (byte letter = 0; letter < NucleotideSequence.ALPHABET.basicSize(); letter++)
            possibleSequences.put(new NucleotideSequence(new byte[]{letter}), letter);
        possibleSequences.put(NucleotideSequence.EMPTY, NucleotideSequence.ALPHABET.basicSize());

        TIntIntHashMap coverage = new TIntIntHashMap();
        TIntObjectHashMap<TIntObjectHashMap<VariantAggregator>> variants = new TIntObjectHashMap<>();

        int nAlignments = 0;
        for (VDJCAlignments al : CUtils.it(alignments.get())) {
            ++nAlignments;
            for (PointSequence point : convertToPointSequences(al)) {

                int seqIndex;
                if (point.sequence.size() == 0)
                    seqIndex = NucleotideSequence.ALPHABET.basicSize();
                else if (point.sequence.size() == 1)
                    seqIndex = point.sequence.getSequence().codeAt(0);
                else {
                    seqIndex = possibleSequences.putIfAbsent(point.sequence.getSequence(), possibleSequences.size());
                    if (seqIndex == -1)
                        seqIndex = possibleSequences.size() - 1;
                }

                coverage.adjustOrPutValue(point.point, 1, 1);

                TIntObjectHashMap<VariantAggregator> map = variants.get(point.point);
                if (map == null)
                    variants.put(point.point, map = new TIntObjectHashMap<>());

                VariantAggregator var = map.get(seqIndex);
                if (var == null)
                    map.put(point.point, var = new VariantAggregator());

                var.count += 1;
                var.sumQuality += point.sequence.getQuality().minValue();
            }
        }

        /*TIntObjectIterator<TIntObjectHashMap<VariantAggregator>> iterator = variants.iterator();
        while (iterator.hasNext()) {
            iterator.advance();

            int point = iterator.key();
            TIntObjectHashMap<VariantAggregator> vars = iterator.value();
            long totalQuality = vars.valueCollection().stream().mapToLong(c -> c.sumQuality).sum();
            if (vars.size() > 1) {
                TIntObjectIterator<VariantAggregator> it = vars.iterator();
                while (it.hasNext()) {
                    it.advance();
                    if (it.value().sumQuality < absoluteSumQualityThreshold
                            && 1.0 * it.value().sumQuality / totalQuality < relativeSumQualityThreshold)
                        it.remove();
                }
            }
        }*/

        long[] forSort = new long[coverage.size()];
        TIntIntIterator iterator = coverage.iterator();
        int i = 0;
        while (iterator.hasNext()) {
            iterator.advance();
            forSort[i++] = -((((long) iterator.value()) << 32) | iterator.key());
        }

        Arrays.sort(forSort);
        int[] pointsArray = Arrays.stream(forSort).mapToInt(l -> (int) (-l)).toArray();
        TIntIntHashMap revIndex = new TIntIntHashMap();
        for (int j = 0; j < pointsArray.length; j++)
            revIndex.put(pointsArray[j], j);

        int[] coverageArray = Arrays.stream(forSort).mapToInt(l -> (int) ((-l) >> 32)).toArray();

        int[][] packedData = new int[pointsArray.length][nAlignments];
        for (int[] aPackedData : packedData)
            Arrays.fill(aPackedData, -1);

        i = 0;
        for (VDJCAlignments al : CUtils.it(alignments.get())) {
            for (PointSequence point : convertToPointSequences(al)) {
                int pointIndex = revIndex.get(point.point);
                packedData[pointIndex][i] =
                        (possibleSequences.get(point.sequence.getSequence()) << 8)
                                | point.sequence.getQuality().minValue();
            }
            i++;
        }

        System.out.println(possibleSequences);
        return new PreparedData(pointsArray, coverageArray) {
            @Override
            OutputPort<int[]> createPort() {
                return CUtils.asOutputPort(Arrays.asList(packedData));
            }
        };
    }

    static abstract class PreparedData {
        final int[] points;
        final int[] coverage;

        PreparedData(int[] points, int[] coverage) {
            this.points = points;
            this.coverage = coverage;
        }

        // array[readId] = (variantId << 8) | minQuality
        abstract OutputPort<int[]> createPort();

        void destroy() {}

        @Override
        public String toString() {
            return "PreparedData{" +
                    "\npoints  =" + Arrays.toString(points) +
                    "\ncoverage=" + Arrays.toString(coverage) +
                    '}';
        }
    }


    static final class VariantAggregator {
        long sumQuality = 0;
        int count = 0;
    }

    PointSequence[] convertToPointSequences(VDJCAlignments alignments) {
        return IntStream.range(0, alignments.getNumberOfReads())
                .mapToObj(i -> convertToPointSequences(alignments, i))
                .flatMap(Collection::stream)
                .collect(Collectors.groupingBy(s -> s.point))
                .values().stream()
                .map(l -> l.stream().max(Comparator.comparingInt(a -> a.sequence.getQuality().meanValue())).get())
                .toArray(PointSequence[]::new);
    }

    void convertToPointSequences(List<PointSequence> points,
                                 Alignment<NucleotideSequence> alignments,
                                 NSequenceWithQuality seq2,
                                 Range seq2Range,
                                 int offset) {

        Range
                alSeq2Range = alignments.getSequence2Range(),
                alSeq2RangeIntersection = alSeq2Range.intersection(seq2Range),
                alSeq1RangeIntersection = alignments.convertToSeq1Range(alSeq2RangeIntersection);

        assert alSeq2RangeIntersection != null;
        assert alSeq1RangeIntersection != null;

        int shift;

        // left
        shift = offset + alignments.getSequence1Range().getFrom() - alignments.getSequence2Range().getFrom();
        for (int i = seq2Range.getFrom(); i < alSeq2RangeIntersection.getFrom(); ++i) {
            NSequenceWithQuality seq = seq2.getRange(i, i + 1);
            points.add(new PointSequence(i + shift, seq));
        }

        // central
        for (int i = alSeq1RangeIntersection.getFrom(); i < alSeq1RangeIntersection.getTo(); ++i) {
            NSequenceWithQuality seq = seq2.getRange(
                    Alignment.aabs(alignments.convertToSeq2Position(i)),
                    Alignment.aabs(alignments.convertToSeq2Position(i + 1)));
            points.add(new PointSequence(i + offset, seq));
        }

        // right
        shift = offset + alignments.getSequence1Range().getTo() - alignments.getSequence2Range().getTo();
        for (int i = alSeq2RangeIntersection.getTo(); i < seq2Range.getTo(); ++i) {
            NSequenceWithQuality seq = seq2.getRange(i, i + 1);
            points.add(new PointSequence(i + shift, seq));
        }
    }

    void convertToPointSequences(List<PointSequence> points,
                                 NSequenceWithQuality seq2,
                                 Range seq2Range,
                                 int offset) {

        for (int i = seq2Range.getFrom(); i < seq2Range.getTo(); ++i) {
            NSequenceWithQuality seq = seq2.getRange(i, i + 1);
            points.add(new PointSequence(i + offset, seq));
        }
    }

    List<PointSequence> convertToPointSequences(VDJCAlignments alignments, int iTarget) {

        //
        //  nLeftDummies     assemblingFeatureLength
        //  ------|--------------|--------------|------------------------>
        //        ↓              ↓              ↓
        //  0000000vvvvvvvvvvvvvvCDR3CDR3CDR3CDR3jjjjjjjjjjjjjjjjCCCCCCCCC

        VDJCPartitionedSequence target = alignments.getPartitionedTarget(iTarget);
        NSequenceWithQuality targetSeq = alignments.getTarget(iTarget);

        Optional<Alignment<NucleotideSequence>> vHitOpt = Arrays.stream(alignments.getHits(Variable))
                .filter(c -> Objects.equals(genes.v, c.getGene()))
                .map(c -> c.getAlignment(iTarget))
                .filter(Objects::nonNull)
                .findAny();

        Optional<Alignment<NucleotideSequence>> jHitOpt = Arrays.stream(alignments.getHits(Joining))
                .filter(c -> Objects.equals(genes.j, c.getGene()))
                .map(c -> c.getAlignment(iTarget))
                .filter(Objects::nonNull)
                .findAny();

        List<PointSequence> points = new ArrayList<>();
        if (target.getPartitioning().isAvailable(assemblingFeature)) {
            int leftStop = target.getPartitioning().getPosition(assemblingFeature.getFirstPoint());
            if (hasV) {
                assert vHitOpt.isPresent();
                convertToPointSequences(points,
                        vHitOpt.get(),
                        targetSeq,
                        new Range(0, leftStop),
                        nLeftDummies);
            } else
                convertToPointSequences(points, targetSeq, new Range(0, leftStop), nLeftDummies - leftStop);

            int rightStart = target.getPartitioning().getPosition(assemblingFeature.getLastPoint());
            if (hasJ) {
                assert jHitOpt.isPresent();
                convertToPointSequences(points,
                        jHitOpt.get(),
                        targetSeq,
                        new Range(rightStart, targetSeq.size()),
                        nLeftDummies + lengthV + assemblingFeatureLength);
            } else
                convertToPointSequences(points, targetSeq, new Range(rightStart, targetSeq.size()), nLeftDummies + lengthV + assemblingFeatureLength - rightStart);
        } else if (hasV && vHitOpt.isPresent())
            convertToPointSequences(points,
                    vHitOpt.get(),
                    targetSeq,
                    new Range(0, vHitOpt.get().getSequence2Range().getTo()),
                    nLeftDummies);
        else if (hasJ && jHitOpt.isPresent())
            convertToPointSequences(points,
                    jHitOpt.get(),
                    targetSeq,
                    new Range(jHitOpt.get().getSequence2Range().getFrom(), targetSeq.size()),
                    nLeftDummies + lengthV + assemblingFeatureLength);

        return points;
    }
}
