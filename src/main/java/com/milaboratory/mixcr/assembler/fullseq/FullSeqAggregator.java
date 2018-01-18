package com.milaboratory.mixcr.assembler.fullseq;

import cc.redberry.pipe.CUtils;
import cc.redberry.pipe.OutputPort;
import com.milaboratory.core.Range;
import com.milaboratory.core.alignment.Alignment;
import com.milaboratory.core.sequence.NSequenceWithQuality;
import com.milaboratory.core.sequence.NSequenceWithQualityBuilder;
import com.milaboratory.core.sequence.NucleotideSequence;
import com.milaboratory.core.sequence.SequenceQuality;
import com.milaboratory.mixcr.assembler.CloneFactoryParameters;
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

    final double minimalQualityShare = 0.1;
    final long minimalSumQuality = 20;
    /**
     * Even if other criteria are not met, but sum quality reaches this threshold, variant is marked as significant.
     */
    final long decisiveSumQualityThreshold = 120;

    final CloneFactoryParameters cloneFactoryParameters;

    final TObjectIntHashMap<NucleotideSequence> sequenceToVariantId
            = new TObjectIntHashMap<>(Constants.DEFAULT_CAPACITY, Constants.DEFAULT_LOAD_FACTOR, -1);
    final TIntObjectHashMap<NucleotideSequence> variantIdToSequence = new TIntObjectHashMap<>();

    FullSeqAggregator(Clone clone, CloneFactoryParameters cloneFactoryParameters) {
        this.clone = clone;
        this.cloneFactoryParameters = cloneFactoryParameters;
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
                    gene.getPartitioning().getLength(vFeature)
                            - gene.getFeature(GeneFeature.intersection(assemblingFeature, vFeature)).size();
        } else
            this.lengthV = 0;

        this.assemblingFeatureLength = clone.getFeature(assemblingFeature).size();

        if (hasJ) {
            VDJCHit jHit = clone.getBestHit(Joining);
            VDJCGene gene = jHit.getGene();
            this.jOffset = gene.getPartitioning().getRelativePosition(jHit.getAlignedFeature(), assemblingFeature.getLastPoint());
        } else
            this.jOffset = 0;
    }

    Clone[] process(PreparedData data) {
        OutputPort<int[]> port = data.createPort();
        List<VariantBranch> branches = new ArrayList<>();
        BitSet allReads = new BitSet();
        allReads.set(0, data.nReads);
        branches.add(new VariantBranch(clone.getCount(), allReads));
        for (int i = 0; i < data.points.length; ++i) {
            int[] variantInfos = port.take();
            List<VariantBranch> newBranches = new ArrayList<>();
            for (VariantBranch branch : branches) {
                List<Variant> variants = callVariants(variantInfos, branch.reads);
                int sumSignificant = 0;
                for (Variant variant : variants)
                    sumSignificant += variant.nSignificant;
                for (Variant variant : variants)
                    newBranches.add(branch.addVariant(variant, sumSignificant));
            }
            branches = newBranches;
        }


        return null;
    }

    AssembledSequences assembleSequences(int[] points, VariantBranch branch) {
        long[] positionedStates = new long[points.length];
        for (int i = 0; i < points.length; i++)
            positionedStates[i] = ((long) points[i]) << 32 | branch.pointStates[i];
        Arrays.sort(positionedStates);

        List<NSequenceWithQuality> sequences = new ArrayList<>();
        List<Range> ranges = new ArrayList<>();
        NSequenceWithQualityBuilder sequenceBuilder = new NSequenceWithQualityBuilder();
        int blockStartPosition = extractPosition(positionedStates[0]);
        for (int i = 0; i < positionedStates.length; ++i) {
            int currentPosition = extractPosition(positionedStates[i]);
            int nextPosition = i == positionedStates.length - 1
                    ? Integer.MAX_VALUE
                    : extractPosition(positionedStates[i + 1]);

            assert currentPosition != nextPosition;

            sequenceBuilder.append(
                    new NSequenceWithQuality(
                            variantIdToSequence.get((int) (positionedStates[i] >>> 8)),
                            (byte) positionedStates[i]));

            if (currentPosition != nextPosition - 1) {
                sequences.add(sequenceBuilder.createAndDestroy());
                ranges.add(new Range(blockStartPosition, currentPosition + 1));
                sequenceBuilder = new NSequenceWithQualityBuilder();
                blockStartPosition = nextPosition;
            }
        }

        int assemblingFeatureTargetId = -1;
        for (int i = 0; i < ranges.size(); i++) {
            if (ranges.get(i).getTo() == nLeftDummies + lengthV) {
                if (i < ranges.size() - 1
                        && ranges.get(i + 1).getFrom() == nLeftDummies + lengthV + assemblingFeatureLength) {
                    // seq[i]-AssemblingFeature-seq[i+1]
                    ranges.set(i, new Range(ranges.get(i).getFrom(), ranges.get(i + 1).getTo()));
                    ranges.remove(i + 1);
                    sequences.set(i, sequences.get(i).concatenate(clone.getTarget(0)).concatenate(sequences.get(i + 1)));
                    sequences.remove(i + 1);
                } else {
                    // seq[i]-AssemblingFeature
                    ranges.set(i, new Range(ranges.get(i).getFrom(), nLeftDummies + lengthV + assemblingFeatureLength));
                    sequences.set(i, sequences.get(i).concatenate(clone.getTarget(0)));
                }
                assemblingFeatureTargetId = i;
                break;
            }

            if (ranges.get(i).getFrom() == nLeftDummies + lengthV + assemblingFeatureLength) {
                // AssemblingFeature-seq[i]
                ranges.set(i, new Range(nLeftDummies + lengthV, ranges.get(i).getTo()));
                sequences.set(i, clone.getTarget(0).concatenate(sequences.get(i)));
                assemblingFeatureTargetId = i;
                break;
            }

            if (ranges.get(i).getFrom() > nLeftDummies + lengthV + assemblingFeatureLength) {
                // seq[i-1]    AssemblingFeature    seq[i]
                ranges.add(i, new Range(nLeftDummies + lengthV, nLeftDummies + lengthV + assemblingFeatureLength));
                sequences.add(i, clone.getTarget(0));
                assemblingFeatureTargetId = i;
                break;
            }
        }

        if (assemblingFeatureTargetId == -1) {
            // seq[last]   AssemblingFeature
            ranges.add(new Range(nLeftDummies + lengthV, nLeftDummies + lengthV + assemblingFeatureLength));
            sequences.add(clone.getTarget(0));
            assemblingFeatureTargetId = ranges.size() - 1;
        }

        return new AssembledSequences(
                assemblingFeatureTargetId,
                ranges.toArray(new Range[ranges.size()]),
                sequences.toArray(new NSequenceWithQuality[sequences.size()]));
    }

    private static int extractPosition(long positionedState) {
        return (int) (positionedState >>> 32);
    }

    private final class AssembledSequences {
        final int assemblingFeatureTargetId;
        final Range[] ranges;
        final NSequenceWithQuality[] sequences;

        public AssembledSequences(int assemblingFeatureTargetId, Range[] ranges, NSequenceWithQuality[] sequences) {
            this.assemblingFeatureTargetId = assemblingFeatureTargetId;
            this.ranges = ranges;
            this.sequences = sequences;
        }
    }

    private static int ABSENT_PACKED_VARIANT_INFO = -1;

    private byte getQuality(int packedVariantInfo) {
        return (byte) (0xFF & packedVariantInfo);
    }

    private int getVariantId(int packedVariantInfo) {
        return packedVariantInfo >>> 8;
    }

    private int packPointVariantInfo(int variantId, byte quality) {
        return variantId << 8 | quality;
    }

    private static class VariantBranch {
        final double count;
        final int[] pointStates;
        final BitSet reads;

        public VariantBranch(double count, BitSet reads) {
            this(count, new int[0], reads);
        }

        public VariantBranch(double count, int[] pointStates, BitSet reads) {
            this.count = count;
            this.pointStates = pointStates;
            this.reads = reads;
        }

        public VariantBranch addVariant(Variant variant, int sumSignificant) {
            int[] newStates = Arrays.copyOf(pointStates, pointStates.length + 1);
            newStates[newStates.length - 1] = variant.variantInfo;
            return new VariantBranch(count * variant.nSignificant / sumSignificant, newStates, variant.reads);
        }
    }

    private static class Variant {
        final int variantInfo;
        final BitSet reads;
        final int nSignificant;

        public Variant(int variantInfo, BitSet reads, int nSignificant) {
            this.variantInfo = variantInfo;
            this.reads = reads;
            this.nSignificant = nSignificant;
        }
    }

    private List<Variant> callVariants(int[] pointVariantInfos, BitSet targetReads) {
        // Pre-calculating number of present variants
        int count = 0;
        for (int readId = targetReads.nextSetBit(0);
             readId >= 0;
             readId = targetReads.nextSetBit(readId + 1)) {
            if (pointVariantInfos[readId] != ABSENT_PACKED_VARIANT_INFO)
                ++count;
        }

        // List of readIds of reads that either:
        //   - don't cover this point
        //   - has insignificant variant in this position
        BitSet unassignedVariants = new BitSet();

        long totalSumQuality = 0;

        // Sorting to GroupBy variantId
        long[] targets = new long[count];
        int i = 0;
        for (int readId = targetReads.nextSetBit(0);
             readId >= 0;
             readId = targetReads.nextSetBit(readId + 1)) {
            if (pointVariantInfos[readId] != ABSENT_PACKED_VARIANT_INFO) {
                targets[i++] = ((long) pointVariantInfos[readId]) << 32 | readId;
                totalSumQuality += 0xFF & pointVariantInfos[readId];
            } else
                unassignedVariants.set(readId);
        }
        Arrays.sort(targets);

        // Collecting measures for each variant
        int blockBegin = 0;
        int currentVariant = (int) (targets[blockBegin] >>> 40);
        int currentIndex = 0;
        long variantSumQuality = 0;

        // Will be used if no significant variant is found
        int bestVariant = -1;
        int bestVariantSumQuality = -1;

        ArrayList<Variant> variants = new ArrayList<>();

        do {
            if (currentIndex == count
                    || currentVariant != (int) (targets[currentIndex] >>> 40)) {

                // Checking significance conditions
                if ((variantSumQuality >= minimalSumQuality
                        && variantSumQuality >= minimalQualityShare * totalSumQuality)
                        || variantSumQuality >= decisiveSumQualityThreshold) {
                    // Variant is significant
                    BitSet reads = new BitSet();
                    for (int j = currentIndex - 1; j >= blockBegin; --j)
                        reads.set((int) targets[j]);
                    variants.add(new Variant((currentVariant << 8) | (int) Math.min(SequenceQuality.MAX_QUALITY_VALUE, variantSumQuality),
                            reads, currentIndex - blockBegin));
                } else {
                    // Variant is not significant
                    for (int j = currentIndex - 1; j >= blockBegin; --j)
                        unassignedVariants.set((int) targets[j]);
                    if (variantSumQuality > bestVariantSumQuality) {
                        bestVariant = currentVariant;
                        // totalSumQuality is definitely less than Long because variantSumQuality < decisiveSumQualityThreshold
                        bestVariantSumQuality = (int) totalSumQuality;
                    }
                }

                if (currentIndex != count) {
                    variantSumQuality = (0xFF & (targets[currentIndex] >>> 32));
                    blockBegin = currentIndex;
                    currentVariant = (int) (targets[blockBegin] >>> 40);
                }
            } else
                variantSumQuality += (0xFF & (targets[currentIndex] >>> 32));
        } while (++currentIndex <= count);

        if (variants.isEmpty()) {
            assert bestVariant != -1;
            BitSet reads = new BitSet();
            for (int j = 0; j < targets.length; j++)
                reads.set((int) targets[j]);
            // nSignificant = 1 (will not be practically used, only one variant, don't care)
            return Collections.singletonList(
                    new Variant(
                            bestVariant << 8 | Math.min(SequenceQuality.MAX_QUALITY_VALUE, bestVariantSumQuality),
                            reads, 1));
        } else {
            for (Variant variant : variants)
                variant.reads.or(unassignedVariants);
            return variants;
        }
    }

    PreparedData initialRound(Supplier<OutputPort<VDJCAlignments>> alignments) {
        if (!sequenceToVariantId.isEmpty())
            throw new IllegalStateException();

        for (byte letter = 0; letter < NucleotideSequence.ALPHABET.basicSize(); letter++) {
            NucleotideSequence seq = new NucleotideSequence(new byte[]{letter});
            sequenceToVariantId.put(seq, letter);
            variantIdToSequence.put(letter, seq);
        }
        sequenceToVariantId.put(NucleotideSequence.EMPTY, NucleotideSequence.ALPHABET.basicSize());
        variantIdToSequence.put(NucleotideSequence.ALPHABET.basicSize(), NucleotideSequence.EMPTY);

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
                    seqIndex = sequenceToVariantId.putIfAbsent(point.sequence.getSequence(), sequenceToVariantId.size());
                    if (seqIndex == -1) {
                        seqIndex = sequenceToVariantId.size() - 1;
                        variantIdToSequence.put(seqIndex, point.sequence.getSequence());
                    }
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
                        (sequenceToVariantId.get(point.sequence.getSequence()) << 8)
                                | point.sequence.getQuality().minValue();
            }
            i++;
        }

        return new PreparedData(nAlignments, pointsArray, coverageArray) {
            @Override
            OutputPort<int[]> createPort() {
                return CUtils.asOutputPort(Arrays.asList(packedData));
            }
        };
    }

    static abstract class PreparedData {
        final int nReads;
        final int[] points;
        final int[] coverage;

        public PreparedData(int nReads, int[] points, int[] coverage) {
            this.nReads = nReads;
            this.points = points;
            this.coverage = coverage;
        }

        // array[readId] = (variantId << 8) | minQuality
        abstract OutputPort<int[]> createPort();

        void destroy() {
        }

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
        return IntStream.range(0, alignments.numberOfTargets())
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
                        nLeftDummies + lengthV + assemblingFeatureLength - jOffset);
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
                    nLeftDummies + lengthV + assemblingFeatureLength - jOffset);

        return points;
    }
}
