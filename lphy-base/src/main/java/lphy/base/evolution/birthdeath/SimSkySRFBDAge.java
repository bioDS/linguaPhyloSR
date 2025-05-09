package lphy.base.evolution.birthdeath;

import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.base.function.tree.OldestYoungest;
import lphy.base.function.tree.PruneTree;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.ValueUtils;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import static lphy.base.evolution.birthdeath.BirthDeathConstants.*;

/**
 * A Birth-death tree generative distribution
 */
public class SimSkySRFBDAge implements GenerativeDistribution<TimeTree> {

    private Value<Number[]> birthRate;
    private Value<Number[]> deathRate;
    private Value<Number[]> psiVal;
    private Value<Number[]> skyTimes;
    private Value<Double> fracVal;
    private Value<Number> originAge;
    private Value<Number> maxSamples;
    private Value<Number> minExtant;

    RandomGenerator random;

    private static final int MAX_ATTEMPTS = 1000;

    public SimSkySRFBDAge(@ParameterInfo(name = lambdaParamName, description = "per-lineage birth rate.") Value<Number[]> birthRate,
                          @ParameterInfo(name = muParamName, description = "per-lineage death rate.") Value<Number[]> deathRate,
                          @ParameterInfo(name = skyTimesParamName, description = "skyline intervals") Value<Number[]> skyTimes,
                          @ParameterInfo(name = fracParamName, description = "fraction of extant taxa sampled.") Value<Double> fracVal,
                          @ParameterInfo(name = psiParamName, description = "per-lineage sampling-through-time rate.") Value<Number[]> psiVal,
                          @ParameterInfo(name = originAgeParamName, description = "the age of the origin.") Value<Number> originAge,
                          @ParameterInfo(name = maxSamplesParamName, description = "maximum number of samples (extant and fossils) to allow in the tree") Value<Number> maxSamples,
                          @ParameterInfo(name = minExtantParamName, description = "minimum number of extant samples permitted in the tree") Value<Number> minExtant) {


        this.birthRate = birthRate;
        this.deathRate = deathRate;
        this.skyTimes = skyTimes;
        this.fracVal = fracVal;
        this.psiVal = psiVal;
        this.originAge = originAge;
        this.maxSamples = maxSamples;
        this.minExtant = minExtant;

        random = RandomUtils.getRandom();
    }

    @GeneratorInfo(name = "SimSkySRFBDAge",
            category = GeneratorCategory.BD_TREE, examples = {"simSRFBDAge.lphy"},
            description = "A skyline stratigraphic range tree of extant species and those sampled through time (in stratigraphic ranges), which is conceptually embedded in a full species tree produced by a speciation-extinction (birth-death) branching process.<br>" +
                    "Conditioned on origin age.")
    public RandomVariable<TimeTree> sample() {
        System.out.println("beginning of sample");
        int nonNullLeafCount = 0;
        TimeTree sampleTree = null;

        int attempts = 0;

        while (nonNullLeafCount < 1 && attempts < MAX_ATTEMPTS) {
            System.out.println("attempts: " + attempts);
            FullBirthDeathSkyTree fullBirthDeathSkyTree = new FullBirthDeathSkyTree(birthRate, deathRate, skyTimes, null, originAge);
            System.out.println("bd tree\n");
            RandomVariable<TimeTree> fullTree = fullBirthDeathSkyTree.sample();
            System.out.println("full tree\n");

            System.out.println("no fossils:");
            System.out.println(fullTree.value().toNewick(true));

            SimFossilsSkyPoisson simFossilsSkyPoisson = new SimFossilsSkyPoisson(fullTree, psiVal, skyTimes);

            Value<TimeTree> fullTreeWithFossils = simFossilsSkyPoisson.sample();

            sampleTree = new TimeTree(fullTreeWithFossils.value());


            List<TimeTreeNode> leafNodes = new ArrayList<>();

            for (TimeTreeNode node : sampleTree.getNodes()) {
                if (node.isLeaf() && node.getAge() == 0.0) {
                    leafNodes.add(node);
                }
            }

            int toNull = (int)Math.round(leafNodes.size()* (1.0-fracVal.value()));
            List<TimeTreeNode> nullList = new ArrayList<>();
            for (int i =0; i < toNull; i++) {
                nullList.add(leafNodes.remove(random.nextInt(leafNodes.size())));
            }
            for (TimeTreeNode node : nullList) {
                node.setId(null);
            }

            List<TimeTreeNode> extant = sampleTree.getObservedExtantNodes();
            List<TimeTreeNode> fossils = sampleTree.getFossilNodes();
            if (ValueUtils.doubleValue(maxSamples) != -1 && ValueUtils.doubleValue(minExtant) != -1) {
                if (fossils.size() + extant.size() > ValueUtils.doubleValue(maxSamples) || extant.size() < ValueUtils.doubleValue(minExtant)) {
                    attempts += 1;
                    continue;
                }
            }

            nonNullLeafCount = leafNodes.size();
            attempts += 1;
        }

        if (attempts == MAX_ATTEMPTS) throw new RuntimeException("Failed to simulate SimSRFBDAge after " + MAX_ATTEMPTS + " attempts.");

        System.out.println("original tree:");
        System.out.println(sampleTree.toNewick(true));

        OldestYoungest oldestYoungest = new OldestYoungest(new Value<>(null, sampleTree));
        TimeTree tree = oldestYoungest.apply().value();
        System.out.println("oldest and youngest:");
        System.out.println(tree.toNewick(true));
        PruneTree pt = new PruneTree(new Value<>(null, tree));
        tree = pt.apply().value();
        System.out.println("observed:");
        System.out.println(tree.toNewick(true));
        return new RandomVariable<>(null, tree, this);

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
            put(fracParamName, fracVal);
            put(psiParamName, psiVal);
            put(originAgeParamName, originAge);
            put(maxSamplesParamName, maxSamples);
            put(minExtantParamName, minExtant);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        switch (paramName) {
            case lambdaParamName:
                birthRate = value;
                break;
            case muParamName:
                deathRate = value;
                break;
            case skyTimesParamName:
                skyTimes = value;
                break;
            case fracParamName:
                fracVal = value;
                break;
            case psiParamName:
                psiVal = value;
                break;
            case originAgeParamName:
                originAge = value;
                break;
            case maxSamplesParamName:
                maxSamples = value;
                break;
            case minExtantParamName:
                minExtant = value;
                break;
            default:
                throw new RuntimeException("Unexpected parameter " + paramName);
        }
    }

    public Value<Number[]> getBirthRate() {
        return birthRate;
    }

    public Value<Number[]> getDeathRate() {
        return deathRate;
    }

    public Value<Double> getRho() {
        return fracVal;
    }

    public Value<Number[]> getPsi() {
        return psiVal;
    }

    public Value<Number[]> getIntervals() {
        return skyTimes;
    }
}