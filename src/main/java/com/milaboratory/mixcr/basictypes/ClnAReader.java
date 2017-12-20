/*
 * Copyright (c) 2014-2017, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */
package com.milaboratory.mixcr.basictypes;

import cc.redberry.pipe.OutputPort;
import com.milaboratory.primitivio.PipeDataInputReader;
import com.milaboratory.primitivio.PrimitivI;
import io.repseq.core.GeneFeature;
import io.repseq.core.VDJCGene;
import io.repseq.core.VDJCLibraryRegistry;

import java.io.DataInput;
import java.io.DataInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.List;

public final class ClnAReader implements AutoCloseable {
    final int chunk = 262144;
    final FileChannel channel;
    final long firstClonePosition;
    final long[] index;
    final long[] counts;
    final VDJCLibraryRegistry libraryRegistry;
    final GeneFeature[] assemblingFeatures;
    final CloneSetIO.GT2GFAdapter alignedFeatures;
    final List<VDJCGene> genes;

    public ClnAReader(Path path, VDJCLibraryRegistry libraryRegistry) throws IOException {
        this.channel = FileChannel.open(path, StandardOpenOption.READ);
        this.libraryRegistry = libraryRegistry;

        PrimitivI input;
        // reading index
        ByteBuffer buf = ByteBuffer.allocate(ClnAWriter.MAGIC_LENGTH + 4);
        channel.read(buf, 0L);

        byte[] magicBytes = new byte[ClnAWriter.MAGIC_LENGTH];
        buf.get(magicBytes);
        String magicString = new String(magicBytes, StandardCharsets.UTF_8);

        if (!magicString.equals(ClnAWriter.MAGIC))
            throw new IllegalArgumentException("Wrong file type. Magic = " + magicString +
                    ", expected = " + ClnAWriter.MAGIC);

        int numberOfClones = buf.getInt();

        buf.reset();
        buf.limit(16);
        long fSize = channel.size();
        channel.read(buf, fSize - 16);

        this.firstClonePosition = buf.getLong();
        long indexBegin = buf.getLong();

        input = new PrimitivI(new InputDataStream(indexBegin, fSize - 8));

        this.index = new long[numberOfClones + 2];
        this.counts = new long[numberOfClones + 2];
        long previousValue = 0;
        for (int i = 0; i < numberOfClones + 2; i++) {
            previousValue = index[i] = previousValue + input.readVarInt();
            counts[i] = input.readVarLong();
        }

        // reading gene features

        input = new PrimitivI(new InputDataStream(ClnAWriter.MAGIC_LENGTH + 4, firstClonePosition));
        this.assemblingFeatures = input.readObject(GeneFeature[].class);
        this.alignedFeatures = new CloneSetIO.GT2GFAdapter(IO.readGF2GTMap(input));
        this.genes = IOUtil.readGeneReferences(input, libraryRegistry);
    }

    public int numberOfClones() {
        return index.length - 2;
    }

    public CloneSet readCloneSet() throws IOException {
        PrimitivI input = new PrimitivI(new InputDataStream(firstClonePosition, index[0]));
        IOUtil.putGeneReferences(input, genes, alignedFeatures);

        int count = index.length - 2;
        List<Clone> clones = new ArrayList<>(count);
        for (int i = 0; i < count; i++)
            clones.add(input.readObject(Clone.class));

        return new CloneSet(clones, genes, alignedFeatures.map, assemblingFeatures);
    }

    public OutputPort<Clone> readClones() throws IOException {
        PrimitivI input = new PrimitivI(new InputDataStream(firstClonePosition, index[0]));
        IOUtil.putGeneReferences(input, genes, alignedFeatures);

        return new PipeDataInputReader<>(Clone.class, input, numberOfClones());
    }

    public OutputPort<VDJCAlignments> alignmentsPort(int clone) throws IOException {
        PrimitivI input = new PrimitivI(new InputDataStream(index[clone + 1], index[clone + 2]));
        IOUtil.putGeneReferences(input, genes, alignedFeatures);
        return new PipeDataInputReader<>(VDJCAlignments.class, input, counts[clone + 1]);
    }

    public OutputPort<CloneAlignments> clonesAndAlignments() throws IOException {
        PrimitivI input = new PrimitivI(new InputDataStream(firstClonePosition, index[0]));
        IOUtil.putGeneReferences(input, genes, alignedFeatures);

        PipeDataInputReader<Clone> clones = new PipeDataInputReader<>(Clone.class, input, numberOfClones());

        return new OutputPort<CloneAlignments>() {
            private int cloneIndex = 0;

            @Override
            public CloneAlignments take() {
                Clone clone = clones.take();
                if (clone == null)
                    return null;
                CloneAlignments result = new CloneAlignments(clone, cloneIndex);
                ++cloneIndex;
                return result;
            }
        };
    }

    public final class CloneAlignments {
        public final Clone clone;
        final int cloneId;

        private CloneAlignments(Clone clone, int cloneId) {
            this.clone = clone;
            this.cloneId = cloneId;
        }

        public OutputPort<VDJCAlignments> alignments() throws IOException {
            return alignmentsPort(cloneId);
        }
    }

    @Override
    public void close() throws Exception {
        channel.close();
    }

    private class InputDataStream implements DataInput {
        private final long to;
        private final ByteBuffer buffer;
        private long lastPosition;

        public InputDataStream(long from, long to) throws IOException {
            this.to = to;
            this.buffer = ByteBuffer.allocate(chunk);
            this.lastPosition = from;
            fillBuffer(true);
        }

        void fillBuffer(boolean first) throws IOException {
            int size = (int) Math.min(chunk, to - lastPosition);
            buffer.position(chunk - size);
            int read = channel.read(buffer, first
                    ? lastPosition
                    : lastPosition - buffer.remaining());
            if (read != size)
                throw new IOException("Wrong block positions.");
            if (read != 0)
                throw new IOException("Wrong block size.");
        }

        void ensureBuffer(int requiredSize) throws IOException {
            if (buffer.remaining() < requiredSize)
                fillBuffer(false);
        }

        @Override
        public void readFully(byte[] b) throws IOException {
            ensureBuffer(b.length);
            buffer.get(b);
        }

        @Override
        public void readFully(byte[] b, int off, int len) throws IOException {
            ensureBuffer(len);
            buffer.get(b, off, len);
        }

        @Override
        public int skipBytes(int n) throws IOException {
            ensureBuffer(n);
            buffer.position(buffer.position() + n);
            return n;
        }

        @Override
        public boolean readBoolean() throws IOException {
            byte b = buffer.get();
            if (b == 1)
                return true;
            else if (b == 0)
                return false;
            else
                throw new IOException("Illegal file format, can't deserialize boolean.");
        }

        @Override
        public byte readByte() throws IOException {
            ensureBuffer(1);
            return buffer.get();
        }

        @Override
        public int readUnsignedByte() throws IOException {
            ensureBuffer(1);
            return 0xFF & buffer.get();
        }

        @Override
        public short readShort() throws IOException {
            ensureBuffer(2);
            return (short) buffer.getChar();
        }

        @Override
        public int readUnsignedShort() throws IOException {
            ensureBuffer(2);
            return 0xFFFF & buffer.getChar();
        }

        @Override
        public char readChar() throws IOException {
            ensureBuffer(2);
            return buffer.getChar();
        }

        @Override
        public int readInt() throws IOException {
            ensureBuffer(4);
            return buffer.getInt();
        }

        @Override
        public long readLong() throws IOException {
            ensureBuffer(8);
            return buffer.getLong();
        }

        @Override
        public float readFloat() throws IOException {
            ensureBuffer(4);
            return buffer.getFloat();
        }

        @Override
        public double readDouble() throws IOException {
            ensureBuffer(8);
            return buffer.getDouble();
        }

        @Override
        public String readLine() throws IOException {
            throw new UnsupportedOperationException();
        }

        @Override
        public String readUTF() throws IOException {
            return DataInputStream.readUTF(this);
        }
    }
}
