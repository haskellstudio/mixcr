package com.milaboratory.mixcr.assembler.fullseq;

import com.milaboratory.core.sequence.NSequenceWithQuality;

/**
 *
 */
final class PointSequence {
    final int point;
    final NSequenceWithQuality sequence;

    PointSequence(int point, NSequenceWithQuality sequence) {
        this.point = point;
        this.sequence = sequence;
    }

    @Override
    public String toString() {
        return String.format("%s->(%s, %s)", point, sequence.getSequence(), sequence.getQuality().meanValue());
    }
}
