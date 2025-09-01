import javax.swing.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Processor {

    SequenceGenerator seqGen = new SequenceGenerator();

    ArrayList<AminoAcid> aminoOG = seqGen.createNewFolding();
    ArrayList<AminoAcid> aminoCopy = new ArrayList<>();
    ArrayList<AminoAcid> aminoChampion = new ArrayList<>();

    ArrayList<Double> fitnessValues = new ArrayList<>();
    ArrayList<ArrayList<AminoAcid>> sequenceCopies = new ArrayList<>();

    int championHBonds = 0;
    int championOverlaps = 0;

    double bestGenerationFitness =0;
    double generationFitness =0;
    double bestFitness =0;
    int iterations = 200;

    double mutationRate = 0.2;

    PrintWriter csvWriter;

    public void initialize()
    {

        try {
            // Open the file for writing and initialize the PrintWriter
            csvWriter = new PrintWriter(new FileWriter("amino_acid_folding.csv"));

            // Write the header to the CSV file
            csvWriter.println("Generation,Average Fitness,Best Generational Fitness,Best Overall Fitness,Mutation Rate,H-Bonds,Overlaps");

            int generations = 100;
            int tournamentSize = 5;

            aminoSequenceFolding();

            display(0);

            generationFitness =0;
            bestFitness =0;

            for (int i = 0; i < generations; i++) {

                //selectSequenceByRouletteWheel();
                selectSequenceByTournament(tournamentSize,i,generations);
                display(i);

                generationFitness = 0;
                bestGenerationFitness = 0;
            }


            visualizeAminoSequence(aminoChampion);

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (csvWriter != null) {
                csvWriter.close(); // Close the writer after the data writing is done
            }
        }

    }

    public void display(int generation)
    {

        System.out.println("Generation: " + generation +
                ", Average fitness: " + generationFitness/iterations +
                ", Best generational fitness: " + bestGenerationFitness +
                ", Best overall fitness: " + bestFitness +
                ", Mutation rate: " + mutationRate +
                ", Number of H-Bonds: " + championHBonds +
                ", Number of Overlaps: " + championOverlaps);

        csvWriter.println(generation + "," +
                (generationFitness / iterations) + "," +
                bestGenerationFitness + "," +
                bestFitness + "," +
                mutationRate + "," +
                championHBonds + "," +
                championOverlaps);

    }

    public void visualizeAminoSequence(ArrayList<AminoAcid> aminoFolding)
    {
        Playground playground = new Playground(aminoFolding);

        JFrame frame = new JFrame("Amino Acid Sequence Visualization");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(500, 500);
        frame.add(playground);
        frame.setVisible(true);
    }

    public double calculateFitness(ArrayList<AminoAcid> aminoSequence)
    {
        int finalHydrogenBonds = 0;
        int finalOverlaps = 0;
        double fitness = 0;
        for(AminoAcid a : aminoSequence)
        {
            if(a.hydrophobic)
            {
                finalHydrogenBonds += checkHydrogenBonds(aminoSequence, a);
            }

            finalOverlaps += countOverlaps(aminoSequence, a);
        }

        fitness = (finalHydrogenBonds * 3) / (finalOverlaps +1);

        if(fitness > bestFitness)
        {
            bestFitness = fitness;
            championOverlaps = finalOverlaps;
            championHBonds = finalHydrogenBonds;
            aminoChampion = deepCopyAminoSequence(aminoSequence);
        }

        //generationFitness += fitness;

        return fitness;

    }

    public int checkHydrogenBonds(ArrayList<AminoAcid> aminoSequence, AminoAcid aminoAcid)
    {
        int hydrogenBonds = 0;
        int currentX = aminoAcid.xCoordinate;
        int currentY = aminoAcid.yCoordinate;

        int northBonds = checkNorth(aminoSequence, aminoAcid, currentX, currentY);
        int eastBonds = checkEast(aminoSequence, aminoAcid, currentX, currentY);
        int southBonds = checkSouth(aminoSequence, aminoAcid, currentX, currentY);
        int westBonds = checkWest(aminoSequence, aminoAcid, currentX, currentY);

        hydrogenBonds = northBonds + eastBonds + southBonds + westBonds;

        return hydrogenBonds;
    }

    public int checkNorth(ArrayList<AminoAcid> aminoSequence, AminoAcid aminoAcid, int hostX, int hostY)
    {
        int bonds=0;
        int northX = hostX;
        int northY = hostY -1;

        bonds = countHydrogenBonds(aminoSequence, aminoAcid, northX, northY);

        return bonds;
    }

    public int checkEast(ArrayList<AminoAcid> aminoSequence, AminoAcid aminoAcid, int hostX, int hostY)
    {
        int bonds=0;
        int eastX = hostX +1;
        int eastY = hostY;

        bonds = countHydrogenBonds(aminoSequence, aminoAcid, eastX, eastY);

        return bonds;
    }

    public int checkSouth(ArrayList<AminoAcid> aminoSequence, AminoAcid aminoAcid, int hostX, int hostY)
    {
        int bonds=0;
        int southX = hostX;
        int southY = hostY +1;

        bonds = countHydrogenBonds(aminoSequence, aminoAcid, southX, southY);

        return bonds;
    }

    public int checkWest(ArrayList<AminoAcid> aminoSequence, AminoAcid aminoAcid, int hostX, int hostY)
    {
        int bonds=0;
        int westX = hostX -1;
        int westY = hostY;

        bonds = countHydrogenBonds(aminoSequence, aminoAcid, westX, westY);

        return bonds;
    }

    public int countHydrogenBonds(ArrayList<AminoAcid> aminoSequence, AminoAcid aminoAcid, int friendX, int friendY)
    {
        int bonds=0;
        for(AminoAcid a : aminoSequence)
        {
            boolean neighbour = a.xCoordinate == friendX && a.yCoordinate == friendY;
            boolean linked = a.previousAA == aminoAcid || a.subsequentAA == aminoAcid;
            boolean counted = a.hydrogenBonds.contains(aminoAcid);
            boolean hydrophobic = a.hydrophobic;

            if(neighbour && !linked && !counted && hydrophobic)
            {
                bonds++;
                storeLinkedAminoAcids(aminoAcid, a);
            }
        }

        return bonds;
    }
    public void storeLinkedAminoAcids(AminoAcid aminoAcidOne, AminoAcid aminoAcidTwo)
    {
        aminoAcidOne.hydrogenBonds.add(aminoAcidTwo);
        aminoAcidTwo.hydrogenBonds.add(aminoAcidOne);
    }

    public int countOverlaps(ArrayList<AminoAcid> aminoSequence, AminoAcid aminoAcid)
    {
        int myX = aminoAcid.xCoordinate;
        int myY = aminoAcid.yCoordinate;

        int overlaps = 0;

        for(AminoAcid a : aminoSequence)
        {
            boolean overlap = a.xCoordinate == myX && a.yCoordinate == myY;
            boolean counted = a.overlaps.contains(aminoAcid);
            boolean me = a == aminoAcid;

            if(overlap && !counted && !me)
            {
                overlaps++;
                storeOverlappedAminoAcids(aminoAcid, a);
            }
        }

        return overlaps;
    }
    public void storeOverlappedAminoAcids(AminoAcid aminoAcidOne, AminoAcid aminoAcidTwo)
    {
        aminoAcidOne.overlaps.add(aminoAcidTwo);
        aminoAcidTwo.overlaps.add(aminoAcidOne);
    }

    public void aminoSequenceFolding()
    {

        for(int i=0; i<iterations; i++)
        {

            //create a deep copy, so the subsequent foldings have still the same origin
            aminoCopy= deepCopyAminoSequence(aminoOG);

            //fold the sequence
            foldSequence();

            //fitness value of the folding gets calculated and the aminoChampion picked
            double fitness = calculateFitness(aminoCopy);

            //fitness value as well as the folding get stored in their respective array
            fitnessValues.add(fitness);
            sequenceCopies.add(deepCopyAminoSequence(aminoCopy));

            generationFitness += fitness;

            //sets the best generational fitness
            collectData(fitness);
        }
    }

    public void selectSequenceByTournament(int tournamentSize, int generation, int maxGenerations) {
        //Random random = new Random();
        ArrayList<Double> selectedFitnessValues = new ArrayList<>();
        ArrayList<ArrayList<AminoAcid>> selectedSequenceCopies = new ArrayList<>();

        double initialMutationRate = 0.1;
        double finalMutationRate = 0.3;

        mutationRate = getMutationRate(generation, maxGenerations, initialMutationRate, finalMutationRate);

        // Perform tournament selection for each candidate in the next generation
        for (int i = 0; i < iterations; i++) {
            // Conduct a tournament to pick one candidate
            ArrayList<AminoAcid> winner = tournament(tournamentSize);

            // Perform mutation on the winner
            ArrayList<AminoAcid> mutant = mutation(winner);

            // Calculate fitness of the mutant
            double newFitness = calculateFitness(mutant);
            collectData(newFitness);
            generationFitness +=newFitness;

            // Store the new candidate in the next generation
            selectedFitnessValues.add(newFitness);
            selectedSequenceCopies.add(deepCopyAminoSequence(mutant));
        }

        // Replace the current generation with the new generation
        fitnessValues = selectedFitnessValues;
        sequenceCopies = selectedSequenceCopies;

        doCrossover(selectedFitnessValues, selectedSequenceCopies);
    }

    public ArrayList<AminoAcid> tournament(int tournamentSize) {
        Random random = new Random();
        ArrayList<AminoAcid> bestCandidate = null;
        double bestFitness = Double.NEGATIVE_INFINITY;

        // Randomly select `tournamentSize` candidates and choose the best one
        for (int i = 0; i < tournamentSize; i++) {
            // Randomly select a candidate from the population
            int randomIndex = random.nextInt(sequenceCopies.size());
            ArrayList<AminoAcid> candidate = sequenceCopies.get(randomIndex);

            // Evaluate the fitness of the candidate
            double candidateFitness = calculateFitness(candidate);

            // Track the best candidate found in the tournament
            if (candidateFitness > bestFitness) {
                bestFitness = candidateFitness;
                bestCandidate = candidate;
            }
        }

        // Return the best candidate from the tournament
        return deepCopyAminoSequence(bestCandidate);
    }

    public double getMutationRate(int generation, int maxGenerations, double initialRate, double finalRate) {
        // Linearly decay the mutation rate from initialRate to finalRate
        return initialRate - ((initialRate - finalRate) * generation / maxGenerations);
    }

    public void selectSequenceByRouletteWheel()
    {
        Random random = new Random();

        ArrayList<Double> selectedFitnessValues = new ArrayList<>();
        ArrayList<ArrayList<AminoAcid>> selectedSequenceCopies = new ArrayList<>();

        boolean crossover = false;

        ArrayList<AminoAcid> crossoverParent = new ArrayList<>();

        for(double fitness : fitnessValues)
        {
            generationFitness += fitness;
        }

        mutationRate = 0.1;

        // Select candidates for the next generation using roulette wheel selection
        for (int i = 0; i < iterations; i++) {
            // Generate a random number between 0 and totalFitness
            double r = random.nextDouble() * generationFitness;

            // Perform roulette wheel selection based on the random number
            double cumulativeFitness = 0;
            for (int j = 0; j < fitnessValues.size(); j++) {
                cumulativeFitness += fitnessValues.get(j);
                if (cumulativeFitness >= r) {

                    //ArrayList<AminoAcid> folding = sequenceCopies.get(j);
/*
                    if(!crossover && getRandomBoolean())
                    {
                        crossoverParent = sequenceCopies.get(j);
                        crossover = true;
                        break;
                    }

                    if(crossover)
                    {
                        performCrossover(sequenceCopies.get(j), crossoverParent, selectedFitnessValues, selectedSequenceCopies);
                        crossover = false;
                        crossoverParent.clear();
                        break;
                    }
 */
                    // Folding gets mutated
                    ArrayList<AminoAcid> mutant = mutation(sequenceCopies.get(j));

                    double newFitness = calculateFitness(mutant);
                    collectData(newFitness);

                    selectedFitnessValues.add(newFitness);
                    selectedSequenceCopies.add(deepCopyAminoSequence(sequenceCopies.get(j)));
                }
            }
        }

        fitnessValues = selectedFitnessValues;
        sequenceCopies = selectedSequenceCopies;

    }

    public void doCrossover(ArrayList<Double> selectedFitnessValues, ArrayList<ArrayList<AminoAcid>> selectedSequenceCopies)
    {
        for(int i=0; i<selectedSequenceCopies.size()-1; i++)
        {
            if(getRandomBoolean())
            {
                performCrossover(selectedSequenceCopies.get(i),selectedSequenceCopies.get(i+1),selectedFitnessValues,selectedSequenceCopies);
                selectedSequenceCopies.remove(i+1);
                selectedSequenceCopies.remove(i);
            }
        }
    }

    public void performCrossover(ArrayList<AminoAcid> parentOne, ArrayList<AminoAcid> parentTwo, ArrayList<Double> selectedFitnessValues, ArrayList<ArrayList<AminoAcid>> selectedSequenceCopies)
    {

        ArrayList<AminoAcid> childOne = new ArrayList<>();
        ArrayList<AminoAcid> childTwo = new ArrayList<>();

        //pick 0 - 19 (lets say 4)
        int length = getRandomPosition(parentOne.size()-1);
        length += 1;
        //Exchange genetic material
        //(4 elements copied)
        for(int i=0; i<length; i++)
        {
            if(i==0)
            {
                AminoAcid aminoAcidOne = new AminoAcid(
                        parentOne.get(i).hydrophobic,
                        parentOne.get(i).enumeration,
                        parentOne.get(i).xCoordinate,
                        parentOne.get(i).yCoordinate,
                        null);
                childOne.add(aminoAcidOne);

                AminoAcid aminoAcidTwo = new AminoAcid(
                        parentTwo.get(i).hydrophobic,
                        parentTwo.get(i).enumeration,
                        parentTwo.get(i).xCoordinate,
                        parentTwo.get(i).yCoordinate,
                        null);
                childTwo.add(aminoAcidTwo);
            }
            else
            {
                AminoAcid aminoAcidOne = new AminoAcid(
                        parentOne.get(i).hydrophobic,
                        parentOne.get(i).enumeration,
                        parentOne.get(i).xCoordinate,
                        parentOne.get(i).yCoordinate,
                        childOne.get(i-1));
                childOne.add(aminoAcidOne);

                AminoAcid aminoAcidTwo = new AminoAcid(
                        parentTwo.get(i).hydrophobic,
                        parentTwo.get(i).enumeration,
                        parentTwo.get(i).xCoordinate,
                        parentTwo.get(i).yCoordinate,
                        childTwo.get(i-1)
                );
                childTwo.add(aminoAcidTwo);
            }
        }

        int offsetXChildOne = childOne.get(length-1).xCoordinate+1 - parentOne.get(length).xCoordinate;
        int offsetYChildOne = childOne.get(length-1).yCoordinate - parentOne.get(length).yCoordinate;

        int offsetXChildTwo = childTwo.get(length-1).xCoordinate+1 - parentTwo.get(length).xCoordinate;
        int offsetYChildTwo = childTwo.get(length-1).yCoordinate - parentTwo.get(length).yCoordinate;

        //Complete child with parent genes
        for(int i= childOne.size(); i<parentOne.size(); i++)
        {
            AminoAcid aminoAcidOne = new AminoAcid(
                    parentOne.get(i).hydrophobic,
                    parentOne.get(i).enumeration,
                    parentOne.get(i).xCoordinate+offsetXChildOne,
                    parentOne.get(i).yCoordinate+offsetYChildOne,
                    childOne.get(i-1));
            childOne.add(aminoAcidOne);

            AminoAcid aminoAcidTwo = new AminoAcid(
                    parentTwo.get(i).hydrophobic,
                    parentTwo.get(i).enumeration,
                    parentTwo.get(i).xCoordinate+offsetXChildTwo,
                    parentTwo.get(i).yCoordinate+offsetYChildTwo,
                    childTwo.get(i-1)
            );
            childTwo.add(aminoAcidTwo);
        }

        for(int i=0; i<childOne.size()-1; i++)
        {
            childOne.get(i).subsequentAA = childOne.get(i+1);
            childTwo.get(i).subsequentAA = childTwo.get(i+1);
        }

        ArrayList<AminoAcid> mutant = mutation(childOne);

        double newFitness = calculateFitness(mutant);
        collectData(newFitness);

        selectedFitnessValues.add(newFitness);
        selectedSequenceCopies.add(deepCopyAminoSequence(mutant));

        mutant = mutation(childTwo);

        newFitness = calculateFitness(mutant);
        collectData(newFitness);

        selectedFitnessValues.add(newFitness);
        selectedSequenceCopies.add(deepCopyAminoSequence(mutant));
    }

    public ArrayList<AminoAcid> mutation(ArrayList<AminoAcid> aminoFolding)
    {
        Random random = new Random();

        for (AminoAcid aminoAcid : aminoFolding)
        {
            if (random.nextDouble() < mutationRate)
            {
                aminoAcid.hydrophobic = !aminoAcid.hydrophobic;
            }
        }

        return aminoFolding;
    }

    public int getRandomPosition(int number)
    {
        Random random = new Random();

        return random.nextInt(number);
    }


    public boolean getRandomBoolean()
    {
        Random random = new Random();

        boolean chance;

        chance = 0 == random.nextInt(100);
        return chance;
    }

    public void foldSequence()
    {
        for(AminoAcid a : aminoCopy)
        {
            if(a.previousAA != null)
            {
                switch (rotateDirection())
                {
                    case 0: break;
                    case 1: rotateClockwise(a);
                        break;
                    case 2: rotateCounterclockwise(a);
                        break;
                }
            }
        }
    }

    public int rotateDirection()
    {
        return new Random().nextInt(3);
    }

    public void rotateCounterclockwise(AminoAcid aminoAcid)
    {
        // Get the pivot point (previous amino acid)
        int pivotX = aminoAcid.previousAA.xCoordinate;
        int pivotY = aminoAcid.previousAA.yCoordinate;

        while (aminoAcid != null) {

            // Calculate the relative position to the pivot
            int deltaX = aminoAcid.xCoordinate - pivotX;
            int deltaY = aminoAcid.yCoordinate - pivotY;

            // Apply counterclockwise rotation
            aminoAcid.xCoordinate = pivotX - deltaY;
            aminoAcid.yCoordinate = pivotY + deltaX;

            // Move to the next amino acid in the chain
            aminoAcid = aminoAcid.subsequentAA;
        }
    }

    public void rotateClockwise(AminoAcid aminoAcid)
    {

        // Get the pivot point (previous amino acid)
        int pivotX = aminoAcid.previousAA.xCoordinate;
        int pivotY = aminoAcid.previousAA.yCoordinate;

        while (aminoAcid != null) {

            // Calculate the relative position to the pivot
            int deltaX = aminoAcid.xCoordinate - pivotX;
            int deltaY = aminoAcid.yCoordinate - pivotY;

            // Apply clockwise 90-degree rotation
            aminoAcid.xCoordinate = pivotX + deltaY;
            aminoAcid.yCoordinate = pivotY - deltaX;

            // Move to the next amino acid in the chain
            aminoAcid = aminoAcid.subsequentAA;
        }
    }

    public void collectData(double localFitness)
    {
        if(localFitness >= bestGenerationFitness)
        {
            bestGenerationFitness = localFitness;
        }
    }

    public ArrayList<AminoAcid> deepCopyAminoSequence(ArrayList<AminoAcid> originalSequence) {
        ArrayList<AminoAcid> copy = new ArrayList<>();
        HashMap<AminoAcid, AminoAcid> aminoMap = new HashMap<>(); // Map to link original and copy for setting previous and subsequent relationships

        // First pass: create new AminoAcid objects with the same properties, but without linking the chain
        for (AminoAcid aa : originalSequence) {
            AminoAcid copyAA = new AminoAcid(aa.hydrophobic, aa.enumeration, aa.xCoordinate, aa.yCoordinate, null);
            copy.add(copyAA);
            aminoMap.put(aa, copyAA);
        }

        // Second pass: set the correct previousAA and subsequentAA relationships in the copied sequence
        for (int i = 0; i < originalSequence.size(); i++) {
            AminoAcid originalAA = originalSequence.get(i);
            AminoAcid copyAA = copy.get(i);

            if (originalAA.previousAA != null) {
                copyAA.previousAA = aminoMap.get(originalAA.previousAA);
            }
            if (originalAA.subsequentAA != null) {
                copyAA.setSubsequentAA(aminoMap.get(originalAA.subsequentAA));
            }
        }
        return copy;
    }
}
