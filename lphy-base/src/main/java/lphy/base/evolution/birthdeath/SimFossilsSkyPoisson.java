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

import java.util.*;

import static lphy.base.evolution.EvolutionConstants.treeParamName;
import static lphy.base.evolution.birthdeath.BirthDeathConstants.psiParamName;
import static lphy.base.evolution.birthdeath.BirthDeathConstants.skyTimesParamName;

/**
 * A java implementation of sim.fossils.poisson in https://github.com/fossilsim/fossilsim/blob/master/R/sim.fossils.R
 */
public class SimFossilsSkyPoisson implements GenerativeDistribution<TimeTree> {

    private Value<TimeTree> tree;
    private Value<Number[]> psi;
    private Value<Number[]> skyTimes;


    private int[] nextFossilNumber;
    private double[] samplingRate;
    private double[] intervals;

    RandomGenerator random;

    static final boolean generateSampledAncestorsAsLeafNodes = true;

    public SimFossilsSkyPoisson(@ParameterInfo(name = treeParamName, description = "Tree to add simulated fossils to.") Value<TimeTree> tree,
                                @ParameterInfo(name = psiParamName, description = "The fossilization rate per unit time per lineage.") Value<Number[]> psi,
                                @ParameterInfo(name = skyTimesParamName, description = "interval times") Value<Number[]> skyTimes) {

        this.tree = tree;
        this.psi = psi;
        this.skyTimes = skyTimes;
        this.random = RandomUtils.getRandom();
    }


    @GeneratorInfo(name = "SimFossilsPoisson",
            category = GeneratorCategory.BD_TREE, examples = {"simFossils.lphy"},
            description = "A tree with fossils added to the given tree at rate psi.")
    public RandomVariable<TimeTree> sample() {

        intervals = ValueUtils.doubleArrayValue(skyTimes);
        samplingRate = ValueUtils.doubleArrayValue(psi);


        TimeTree treeCopy = new TimeTree(tree.value());

        simulateFossils(treeCopy, samplingRate);
        relabelPast(treeCopy);
        relabelExtant(treeCopy);

        return new RandomVariable<>(null, treeCopy, this);
    }

    private void relabelPast(TimeTree tree) {
        for (int i = 0; i < tree.getMaxLineage(); i++){
            List<TimeTreeNode> species_nodes =tree.getLineageNodes(i);
            ArrayList<Double> ages = new ArrayList<>();
            if (!species_nodes.isEmpty()){
                for (TimeTreeNode node : species_nodes){
                    ages.add(node.getAge());
                }
                Collections.sort(ages);
                for (TimeTreeNode node : species_nodes){
                    int fossilNo = ages.indexOf(node.getAge())+1;
                    node.setId("f"+i+"_"+ fossilNo);
                }

            }
        }
    }

    private void relabelExtant(TimeTree tree) {
        List<TimeTreeNode> extant = tree.getExtantNodes();
        for( TimeTreeNode node : extant ) {
            int lineage = node.getLineage();
            node.setId("f"+lineage+"_0");
        }
    }

    private void simulateFossils(TimeTree tree, double[] psi) {
        nextFossilNumber = new int[tree.getMaxLineage()];

        for (TimeTreeNode node : tree.getNodes()) {

            if (!node.isRoot()) {
                int lineage = node.getLineage();
                double min = node.getAge();
                double max = node.getParent().getAge();

                int minInterval = 0;
                int maxInterval = 0;
                System.out.println("max: " + max + " min: " + min);
                for(int i = 0; i < intervals.length; i++) {
                    if (maxInterval == 0 &&  max >= intervals[i]) {
                        maxInterval = i;
                    }
                    if (minInterval == 0 && min >= intervals[i]) {
                        minInterval = i;
                    }
                }
                System.out.println("maxInterval: " + maxInterval +  " minInterval: " + minInterval);
                if (maxInterval == minInterval) {
                    System.out.println("A simulating fossils from " + min + " to " + max);
                    nextFossilNumber[lineage - 1] += getFossils(min, max, psi[maxInterval], node, lineage, tree);
                }
                else {
                    nextFossilNumber[lineage - 1] += getFossils(intervals[maxInterval], max, psi[maxInterval], node, lineage, tree);
                    System.out.println("B simulating fossils from " + intervals[maxInterval] + " to " + max + " psi: " + psi[maxInterval]);
                    nextFossilNumber[lineage - 1] += getFossils(min, intervals[minInterval-1], psi[minInterval], node, lineage, tree);
                    System.out.println("C simulating fossils from " + min + " to " + intervals[minInterval-1] + " psi: " + psi[minInterval]);
                    for (int i = maxInterval + 1; i < minInterval-1; i++) {
                        nextFossilNumber[lineage - 1] += getFossils(intervals[i], intervals[i+1], psi[i], node, lineage, tree);
                        System.out.println("D simulating fossils from " + intervals[i] + " to " + intervals[i+1] + " psi: " + psi[i]);
                    }
                }
            }
        }

        tree.setRoot(tree.getRoot(), true);
    }

    private int getFossils(double min, double max, double psi, TimeTreeNode node, int lineage,TimeTree tree){
        double expectedFossils = (max - min) * psi;

        PoissonDistribution poissonDistribution = new PoissonDistribution(random, expectedFossils, 1e-8, 100);
        int fossils = poissonDistribution.sample();
        if (fossils > 0) {
            double[] fossilTimes = new double[fossils];
            for (int i = 0; i < fossils; i++) {
                fossilTimes[i] = random.nextDouble() * (max - min) + min;
            }
            addFossils(fossilTimes, node.getParent(), node, lineage, tree);
        }
        return fossils;
    }

    private void addFossils(double[] times, TimeTreeNode parent, TimeTreeNode child, int lineage, TimeTree tree) {
        Arrays.sort(times);

        nextFossilNumber[lineage - 1] = nextFossilNumber[lineage - 1] + times.length;

        for (int i = times.length - 1; i >= 0; i--) {
            Taxon fossilTaxon = new Taxon("f"+lineage+"_"+nextFossilNumber[lineage - 1], times[i]);

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
            nextFossilNumber[lineage - 1] -= 1;
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
            put(skyTimesParamName, skyTimes);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(treeParamName)) tree = value;
        else if (paramName.equals(psiParamName)) psi = value;
        else if (paramName.equals(skyTimesParamName)) skyTimes = value;
        else throw new IllegalArgumentException("Expected either " + treeParamName + " or " + psiParamName + " or " + skyTimesParamName + " as parameter name but got " + paramName );
    }

    public String toString() {
        return getName();
    }
}
