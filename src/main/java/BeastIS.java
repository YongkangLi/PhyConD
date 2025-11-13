import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.yongkangl.phylogenetics.io.BeastTreeLineParser;
import com.yongkangl.phylogenetics.io.TreeNode;
import com.yongkangl.phylogenetics.model.HiddenNodesSMC;
import com.yongkangl.phylogenetics.sampler.HiddenNodes;
import com.yongkangl.phylogenetics.sampler.TreeNodesIS;
import com.yongkangl.phylogenetics.substitution.*;
import org.apache.commons.cli.*;

public class BeastIS {
    public static void main(String[] args) {
        Options options = new Options();
        options.addOption("s", "steps", true, "Number of steps");
        options.addOption("p", "particles", true, "Number of particles");
        options.addOption("m", "mutationSteps", true, "Number of mutation steps");
        options.addOption("l", "line", true, "Line of input");
        options.addOption("f", "file", true, "File path");

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("Error parsing command line: " + e.getMessage());
            System.exit(1);
        }

        int nSteps = 8;
        int nParticles = 16;
        int mutationSteps = 8;
        int line = 10000;
        String filePath = "treeLines.txt";

        if (cmd.hasOption("steps")) {
            try {
                nSteps = Integer.parseInt(cmd.getOptionValue("steps"));
            } catch (NumberFormatException e) {
                System.err.println("Invalid number for steps");
                System.exit(1);
            }
        }
        if (cmd.hasOption("particles")) {
            try {
                nParticles = Integer.parseInt(cmd.getOptionValue("particles"));
            } catch (NumberFormatException e) {
                System.err.println("Invalid number for particles");
                System.exit(1);
            }
        }
        if (cmd.hasOption("mutationSteps")) {
            try {
                mutationSteps = Integer.parseInt(cmd.getOptionValue("mutationSteps"));
            } catch (NumberFormatException e) {
                System.err.println("Invalid number for mutation steps");
                System.exit(1);
            }
        }
        if (cmd.hasOption("line")) {
            try {
                line = Integer.parseInt(cmd.getOptionValue("line"));
            } catch (NumberFormatException e) {
                System.err.println("Invalid number for line");
                System.exit(1);
            }
        }
        if (cmd.hasOption("file")) {
            filePath = cmd.getOptionValue("file");
        }

        String inputLine = "";
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            for (int i = 0; i < line; i++) {
                inputLine = reader.readLine();
                if (inputLine == null)
                    break;
            }
        } catch (IOException e) {
            System.err.println("Error reading file: " + e.getMessage());
            System.exit(1);
        }


//        double[][] Q = {{-1.13964572, 0.33807326, 0.54065275, 0.26091972}, {0.20687479, -0.74294866, 0.19805718, 0.33801669}, {0.51107285, 0.27280825, -0.992562, 0.2086809}, {0.18542334, 0.44322683, 0.22394699, -0.85259716}};
        double[][] scaledQ = {{-1.299225557498793, 0.3854122453465731, 0.6163580944860869, 0.297455217666133}, {0.2358426020428258, -0.8469806551036532, 0.22579029897487418, 0.38534775408595323}, {0.5826362447422561, 0.3487971468891999, -1.169335005707481, 0.23790161407602503}, {0.2113873951730493, 0.4475211646353517, 0.17885063667912704, -0.8377591964875281}};
        IndependentSiteModel ism = new IndependentSiteModel(scaledQ);
        SiteModel target = new ARMADiLLOSiteModel(1.14002582);

        BeastTreeLineParser treeLineParser = new BeastTreeLineParser(inputLine);
        TreeNode root = treeLineParser.parseBeast(true, true);
        root.setUpMessages(ism, target.getContextLength());

        HiddenNodesSMC SMC = new HiddenNodesSMC(root, ism, target, nSteps, nParticles, mutationSteps);
        SMC.run();
        System.out.println(SMC.getLogWeight());
//        System.out.println();
//        System.out.println(SMC.sequences());

//        root.calculateExactly(target, false, null);
//        System.out.println("DONE!");

//        HiddenNodes hiddenNodes = new HiddenNodes(1, root, 1, new SiteModel[]{ism, target});
//        for (int i = 0; i < 10000; i++) {
//            hiddenNodes.mutate();
//            if (i % 200 == 0) {
//                System.out.println(hiddenNodes);
//            }
//        }
    }
}