/*
 * Copyright (c) 2014-2018, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
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
package com.milaboratory.mixcr.cli;

import cc.redberry.pipe.CUtils;
import cc.redberry.pipe.OutputPort;
import cc.redberry.pipe.util.FlatteningOutputPort;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.IntegerConverter;
import com.milaboratory.cli.*;
import com.milaboratory.mixcr.basictypes.*;
import gnu.trove.map.hash.TIntIntHashMap;
import io.repseq.core.VDJCLibraryRegistry;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

public class ActionSliceClnA implements Action {
    final Params params = new Params();

    @Override
    public void go(ActionHelper helper) throws Exception {
        Collections.sort(params.cloneIds);

        try (ClnAReader reader = new ClnAReader(params.getInputFileName(), VDJCLibraryRegistry.getDefault());
             ClnAWriter writer = new ClnAWriter(params.getOutputFileName())) {

            // Getting full clone set
            CloneSet cloneSet = reader.readCloneSet();

            // old clone id -> new clone id
            TIntIntHashMap idMapping = new TIntIntHashMap();

            long newNumberOfAlignments = 0;

            // Creating new cloneset
            List<Clone> clones = new ArrayList<>();
            int i = 0;
            List<OutputPort<VDJCAlignments>> allAlignmentsList = new ArrayList<>();
            for (Integer cloneId : params.cloneIds) {
                newNumberOfAlignments += reader.numberOfAlignmentsInClone(cloneId);
                Clone clone = cloneSet.get(cloneId);
                idMapping.put(clone.getId(), i);
                clones.add(clone.setId(i).resetParentCloneSet());
                OutputPort<VDJCAlignments> als = reader.readAlignmentsOfClone(cloneId);
                final int ii = i;
                allAlignmentsList.add(() -> {
                    VDJCAlignments al = als.take();
                    if (al == null)
                        return null;
                    return al.updateCloneIndex(ii);
                });
                i++;
            }

            CloneSet newCloneSet = new CloneSet(clones, cloneSet.getUsedGenes(),
                    cloneSet.getAlignedFeatures(), cloneSet.getAlignmentParameters(), cloneSet.getAssemblerParameters());

            OutputPort<VDJCAlignments> allAlignmentsPortRaw = new FlatteningOutputPort<>(CUtils.asOutputPort(allAlignmentsList));
            AtomicLong idGen = new AtomicLong();
            OutputPort<VDJCAlignments> allAlignmentsPort = () -> {
                VDJCAlignments al = allAlignmentsPortRaw.take();
                if (al == null)
                    return null;
                return al.setAlignmentsIndex(idGen.getAndIncrement());
            };

            writer.writeClones(newCloneSet);

            writer.sortAlignments(allAlignmentsPort, newNumberOfAlignments);

            writer.writeAlignmentsAndIndex();
        }
    }

    @Override
    public String command() {
        return "sliceClnA";
    }

    @Override
    public ActionParameters params() {
        return params;
    }

    @Parameters(commandDescription = "Slice ClnA file.")
    @HiddenAction
    public static final class Params extends ActionParametersWithOutput {
        @Parameter(description = "[input_file1.vdjca[.gz] [input_file2.vdjca[.gz] ....]] output_file.vdjca[.gz]")
        List<String> parameters;

        @Parameter(description = "List of read ids to export",
                names = {"-i", "--clone-id"},
                converter = IntegerConverter.class)
        List<Integer> cloneIds = new ArrayList<>();

        String getInputFileName() {
            return parameters.get(0);
        }

        String getOutputFileName() {
            return parameters.get(1);
        }

        @Override
        protected List<String> getOutputFiles() {
            return Collections.singletonList(getOutputFileName());
        }

        @Override
        public void validate() {
            if (parameters.size() != 2)
                throw new ParameterException("Wrong number of parameters.");
        }
    }
}
