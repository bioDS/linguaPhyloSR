package lphy.base.evolution.birthdeath;

import lphy.base.evolution.Taxon;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.ValueUtils;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import static lphy.base.evolution.EvolutionConstants.treeParamName;
import static lphy.base.evolution.birthdeath.BirthDeathConstants.psiParamName;

/**
 * A java implementation of sim.fossils.poisson in https://github.com/fossilsim/fossilsim/blob/master/R/sim.fossils.R
 */
public class SimFossilsPoisson implements GenerativeDistribution<TimeTree> {

    private Value<TimeTree> tree;
    private Value<Number> psi;

    RandomGenerator random;

    static final boolean generateSampledAncestorsAsLeafNodes = true;

    public SimFossilsPoisson(@ParameterInfo(name = treeParamName, description = "Tree to add simulated fossils to.") Value<TimeTree> tree,
                             @ParameterInfo(name = psiParamName, description = "The fossilization rate per unit time per lineage.") Value<Number> psi) {

        this.tree = tree;
        this.psi = psi;
        this.random = RandomUtils.getRandom();
    }


    @GeneratorInfo(name = "SimFossilsPoisson",
            category = GeneratorCategory.BD_TREE, examples = {"simFossils.lphy"},
            description = "A tree with fossils added to the given tree at rate psi.")
    public RandomVariable<TimeTree> sample() {

        double samplingRate = ValueUtils.doubleValue(psi);

        TimeTree treeCopy = new TimeTree(tree.value());

        simulateFossils(treeCopy, samplingRate);
        relabelExtant(treeCopy);

        return new RandomVariable<>(null, treeCopy, this);
    }

    private void relabelExtant(TimeTree tree) {
        List<TimeTreeNode> extant = tree.getExtantNodes();
        for( TimeTreeNode node : extant ) {
            int lineage = node.getLineage();
//            System.out.println(node);
//            System.out.println(node.getLineage());
//            System.out.println(nextFossilNumber[node.getLineage() - 1]);
            node.setId("f"+lineage+"_0");
        }
    }

    private void simulateFossils(TimeTree tree, double psi) {

//        int nextFossilNumber = 0;
//        List<Integer> nextFossilNumbers = new ArrayList<>();
        int[] nextFossilNumber = new int[tree.getMaxLineage()];



        for (TimeTreeNode node : tree.getNodes()) {

            if (!node.isRoot()) {
                int lineage = node.getLineage();
                double min = node.getAge();
                double max = node.getParent().getAge();
                double expectedFossils = (max - min) * psi;

                PoissonDistribution poissonDistribution = new PoissonDistribution(random,expectedFossils,1e-8, 100);
                int fossils = poissonDistribution.sample();
                if (fossils > 0) {
                    double[] fossilTimes = new double[fossils];
                    for (int i = 0; i < fossils; i++) {
                        fossilTimes[i] = random.nextDouble() * (max - min) + min;
                    }
                    addFossils(fossilTimes, node.getParent(), node, lineage, nextFossilNumber[lineage-1], tree);
                }
                nextFossilNumber[lineage-1] += fossils;
//                nextFossilNumber += fossils;
            }
        }

        tree.setRoot(tree.getRoot(), true);
    }

    private void addFossils(double[] times, TimeTreeNode parent, TimeTreeNode child, int lineage, int nextFossilNumber, TimeTree tree) {
        Arrays.sort(times);

        nextFossilNumber = nextFossilNumber + times.length;

        for (int i = times.length - 1; i >= 0; i--) {
            Taxon fossilTaxon = new Taxon("f"+lineage+"_"+nextFossilNumber, times[i]);

            parent.removeChild(child);

            TimeTreeNode fossilNode;
            if (generateSampledAncestorsAsLeafNodes) {
                TimeTreeNode fossilLeafNode = new TimeTreeNode(fossilTaxon,tree);
                fossilNode = new TimeTreeNode(fossilTaxon.getAge());
                fossilNode.addChild(fossilLeafNode);
            } else {
                fossilNode = new TimeTreeNode(fossilTaxon,tree);
            }

            parent.addChild(fossilNode);
            fossilNode.addChild(child);
            parent = fossilNode;
            nextFossilNumber -= 1;
        }
    }
    
    @Override
    public double logDensity(TimeTree timeTree) {
        throw new UnsupportedOperationException("Not implemented!");
    }

    @Override
    public Map<String, Value> getParams() {
        return new TreeMap<>() {{
            put(treeParamName, tree);
            put(psiParamName, psi);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(treeParamName)) tree = value;
        else if (paramName.equals(psiParamName)) psi = value;
        else throw new IllegalArgumentException("Expected either " + treeParamName + " or " + psiParamName + " as parameter name but got " + paramName );
    }

    public String toString() {
        return getName();
    }
}
