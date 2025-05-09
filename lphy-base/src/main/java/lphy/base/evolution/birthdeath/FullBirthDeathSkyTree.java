package lphy.base.evolution.birthdeath;

import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.ValueUtils;
import lphy.core.model.annotation.Citation;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Arrays;
import java.util.TreeMap;
import java.util.Collections;
import static lphy.base.evolution.birthdeath.BirthDeathConstants.*;

/**
 * A Birth-death tree generative distribution
 */
@Citation(value="David G. Kendall. On the Generalized \"Birth-and-Death\" Process, " +
        "The Annals of Mathematical Statistics, Ann. Math. Statist. 19(1), 1-15, March, 1948.",
        year=1948,
        title="On the Generalized \"Birth-and-Death\" Process",
        authors={"Kendall"},
        DOI="https://doi.org/10.1214/aoms/1177730285")
public class FullBirthDeathSkyTree implements GenerativeDistribution<TimeTree> {

    private Value<Number[]> birthRate;
    private Value<Number[]> deathRate;
    private Value<Number[]> skyTimes;
    private Value<Number> rootAge;
    private Value<Number> originAge;

    private double[] lambda;
    private double[] mu;
    private double[] intervals;

    private List<TimeTreeNode> activeNodes;
    private int maxLineage = 1;

    RandomGenerator random;

    private static final int MAX_ATTEMPTS = 1000;

    public FullBirthDeathSkyTree(@ParameterInfo(name = lambdaParamName, description = "per-lineage birth rate.") Value<Number[]> birthRate,
                                 @ParameterInfo(name = muParamName, description = "per-lineage death rate.") Value<Number[]> deathRate,
                                 @ParameterInfo(name = skyTimesParamName, description = "skyline interval times") Value<Number[]> skyTimes,
                                 @ParameterInfo(name = rootAgeParamName, description = "the age of the root of the tree (only one of rootAge and originAge may be specified).", optional=true) Value<Number> rootAge,
                                 @ParameterInfo(name = originAgeParamName, description = "the age of the origin of the tree  (only one of rootAge and originAge may be specified).", optional=true) Value<Number> originAge) {
        this.birthRate = birthRate;
        this.deathRate = deathRate;
        this.rootAge = rootAge;
        this.originAge = originAge;
        this.skyTimes = skyTimes;
        this.random = RandomUtils.getRandom();



        if (rootAge != null && originAge != null) throw new IllegalArgumentException("Only one of rootAge and originAge may be specified!");
        if (rootAge == null && originAge == null) throw new IllegalArgumentException("One of rootAge and originAge must be specified!");

        lambda = ValueUtils.doubleArrayValue(birthRate);
        mu = ValueUtils.doubleArrayValue(deathRate);
        intervals = ValueUtils.doubleArrayValue(skyTimes);
        int l = lambda.length;
        int m = mu.length;
        int i = intervals.length;
        double[] intervalsSorted = new double[i];
        System.arraycopy( intervals, 0, intervalsSorted, 0, intervals.length );
        java.util.Arrays.sort(intervalsSorted);
        intervalsSorted = reverse(intervalsSorted);

        if (l != m || l != i || m != i) {
            throw new IllegalArgumentException("Birth, death and skyline intervals must have the same length");
        }
        if (!(Arrays.equals(intervalsSorted, intervals))){
            throw new IllegalArgumentException("Interval times must be ordered from largest to smallest");
        }

        System.out.println("intervals " + Arrays.toString(intervals));


//        for (int a = 0; a < i; a++){
//            intervals[a] = ValueUtils.doubleValue(originAge) - intervals[a];
//        }


        activeNodes = new ArrayList<>();
    }
    // reverse code from https://www.geeksforgeeks.org/reverse-an-array-in-java/
    private double[] reverse(double[] array){
        for (int j = 0; j < array.length / 2; j++) {
            double t = array[j];
            array[j] = array[array.length - 1 - j];
            array[array.length - 1 - j] = t;
        }
        return array;
    }


    @GeneratorInfo(name = "FullBirthDeath",
            category = GeneratorCategory.BD_TREE, examples = {"simpleFullBirthDeath.lphy"},
            description = "A birth-death tree with both extant and extinct species.<br>" +
            "Conditioned on age of root or origin.")
    public RandomVariable<TimeTree> sample() {

        boolean success = false;
        TimeTree tree = new TimeTree();
        TimeTreeNode root = null;
        int currentInt = 0;


        int attempts = 0;

        while (!success && attempts < MAX_ATTEMPTS) {
            currentInt = 0;

            activeNodes.clear();

            root = new TimeTreeNode((String)null, tree);
            if (rootAge != null) {
                root.setAge(ValueUtils.doubleValue(rootAge));
            } else {
                root.setAge(ValueUtils.doubleValue(originAge));
            }
            root.setLineage(1);


            double time = root.getAge();

            if (rootAge != null) {
                activeNodes.add(root);
                doBirth(activeNodes, time, tree);
            } else {
                TimeTreeNode origin = root;
                root = new TimeTreeNode((String)null, tree);
                root.setAge(origin.getAge());
                root.setLineage(1);
                origin.addChild(root);
                activeNodes.add(root);
                root = root.getParent();
            }

            while (time > 0.0 && activeNodes.size() > 0) {


                while(intervals[currentInt] >= time){
                    currentInt++;
                }
//                System.out.println("time: " + time + " lower interval bound: " + intervals[currentInt]);
//                System.out.println("birth: " + lambda[currentInt] + " death: " + mu[currentInt]);

                int k = activeNodes.size();

                double totalRate = (lambda[currentInt] + mu[currentInt]) * (double) k;

                // random exponential variate
                double x = -Math.log(random.nextDouble()) / totalRate;
                time -= x;

                if(time <= intervals[currentInt]){
                    time = intervals[currentInt];
//                    System.out.println("next interval");
                    continue;
                }

                if (time < 0) break;

                double U = random.nextDouble();
                if (U < lambda[currentInt] / (lambda[currentInt] + mu[currentInt])) {
                    doBirth(activeNodes, time, tree);
                } else {
                    doDeath(activeNodes, time);
                }
            }

            int number = 0;
            for (TimeTreeNode node : activeNodes) {
                node.setAge(0.0);
                node.setId(number+"");
                number += 1;
            }

            success = activeNodes.size() > 0;
            attempts += 1;
        }

        if (!success) {
            throw new RuntimeException("Failed to simulated FullBirthSkyDeathTree after " + MAX_ATTEMPTS + " attempts.");
        }

        tree.setRoot(root, true);

        return new RandomVariable<>(null, tree, this);
    }

    private void doBirth(List<TimeTreeNode> activeNodes, double age, TimeTree tree) {
        TimeTreeNode parent = activeNodes.remove(random.nextInt(activeNodes.size()));
        parent.setAge(age);
        TimeTreeNode child1 = new TimeTreeNode((String)null, tree);
        TimeTreeNode child2 = new TimeTreeNode((String)null, tree);
        child1.setAge(age);
        child2.setAge(age);
        parent.addChild(child1);
        parent.addChild(child2);
        child1.setLineage(parent.getLineage());
        maxLineage = maxLineage + 1;
        child2.setLineage(maxLineage);
        activeNodes.add(child1);
        activeNodes.add(child2);
    }

    private void doDeath(List<TimeTreeNode> activeNodes, double age) {
        TimeTreeNode deadNode = activeNodes.remove(random.nextInt(activeNodes.size()));
        deadNode.setAge(age);
    }


    @Override
    public double logDensity(TimeTree timeTree) {

        throw new UnsupportedOperationException("Not implemented!");
    }

    @Override
    public Map<String, Value> getParams() {
        return new TreeMap<>() {{
            put(lambdaParamName, birthRate);
            put(muParamName, deathRate);
            put(skyTimesParamName, skyTimes);
            if (rootAge != null) put(rootAgeParamName, rootAge);
            if (originAge != null) put(originAgeParamName, originAge);
        }};
    }



    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(lambdaParamName)) birthRate = value;
        else if (paramName.equals(muParamName)) deathRate = value;
        else if (paramName.equals(skyTimesParamName)) skyTimes = value;
        else if (paramName.equals(rootAgeParamName)) rootAge = value;
        else if (paramName.equals(originAgeParamName)) originAge = value;
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);
    }

    public String toString() {
        return getName();
    }

}
