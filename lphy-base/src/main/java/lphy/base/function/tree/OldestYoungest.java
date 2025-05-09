package lphy.base.function.tree;

import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Arrays;

import static lphy.base.evolution.EvolutionConstants.treeParamName;

/**
 * A function to remove any fossils from a tree that are not the oldest and youngest fossil of each species
 */

public class OldestYoungest extends DeterministicFunction<TimeTree> {

    public OldestYoungest(@ParameterInfo(name = treeParamName, description = "the full tree to remove fossils from") Value<TimeTree> tree) {
        setParam(treeParamName, tree);
    }

    @GeneratorInfo(name = "oldestYoungest",
            category = GeneratorCategory.TREE, examples = {"simSRFBD.lphy"},
            description = "A tree with all fossils removed except for the oldest and youngest of each species, to create a stratigraphic range")
    public Value<TimeTree> apply() {

        Value<TimeTree> tree = getParams().get(treeParamName);

        // do deep copy
        TimeTree prunedTree = new TimeTree(tree.value());
        int[] nextFossilNumber = new int[prunedTree.getMaxLineage()];
        int[] minFossilNumber = new int[prunedTree.getMaxLineage()];
        Arrays.fill(minFossilNumber, -1);



        for (TimeTreeNode node : prunedTree.getNodes()) {
            String id = node.getId();
            String[] idVals;
            int lineage;
            int fossilNo;
            if (id != null) {
                idVals = id.split("_");
                if (idVals.length == 2) {
                    lineage = Integer.parseInt(idVals[0].substring(1));
                    fossilNo = Integer.parseInt(idVals[1]);
                    node.setLineage(lineage);
                    node.setFossilNo(fossilNo);
                    if (nextFossilNumber[lineage-1] < fossilNo) {
                        nextFossilNumber[lineage-1] = fossilNo;
                    }
                    if (minFossilNumber[lineage-1] == -1 || minFossilNumber[lineage-1] > fossilNo) {
                        minFossilNumber[lineage-1] = fossilNo;
                    }
                }
//                else {
//                    System.out.println("id value is not standard\n");
//                }
            }

            }
        List treeNodes = prunedTree.getNodes();
        ArrayList<TimeTreeNode> nodes = new ArrayList<>(treeNodes.size());
        nodes.addAll(treeNodes);
        List<Integer> removeNodes = new ArrayList<>();
        for (Iterator<TimeTreeNode> iterator = nodes.iterator(); iterator.hasNext(); ) {
            TimeTreeNode node = iterator.next();
            if (node.getId() != null) {
                int fossilNo = node.getFossilNo();
                int lineage = node.getLineage();
                if (fossilNo == nextFossilNumber[lineage - 1]) {
                    System.out.println("setting " + node.getId() + " to " + "f" + lineage + "_first");
                    node.setId("f" + lineage + "_first");
                }
                else {
                    if (fossilNo == minFossilNumber[lineage-1]) {
                        node.setId("f"+lineage+"_last");
                        System.out.println("setting " + node.getId() + " to " + "f" + lineage + "_last");
                    } else {
                        TimeTreeNode parent = node.getParent();
                        TimeTreeNode grandParent = parent.getParent();
                        TimeTreeNode sibling = parent.getChild(1);
                        parent.removeChild(node);
                        parent.removeChild(sibling);

                        grandParent.removeChild(parent);
                        grandParent.addChild(sibling);
                        removeNodes.add(nodes.indexOf(node));
                        removeNodes.add(nodes.indexOf(parent));
                    }
                }
            }
        }
        Collections.sort(removeNodes, Collections.reverseOrder());
        for (int i = 0; i < removeNodes.size(); i++) {
            nodes.remove(removeNodes.get(i).intValue());
        }
        nodes.trimToSize();
        prunedTree = new TimeTree(nodes);


        return new Value<>(null, prunedTree, this);
    }
}
