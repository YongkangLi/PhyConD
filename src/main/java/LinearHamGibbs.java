import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import com.yongkangl.phylogenetics.io.LinearHamTreeLineParser;
import com.yongkangl.phylogenetics.io.TreeNode;
import com.yongkangl.phylogenetics.model.HiddenNodesSMC;
import com.yongkangl.phylogenetics.sampler.HiddenNodes;
import com.yongkangl.phylogenetics.sampler.TreeNodesIS;
import com.yongkangl.phylogenetics.substitution.*;
import org.apache.commons.cli.*;

public class LinearHamGibbs {
    public static void main(String[] args) {
        Options options = new Options();
        options.addOption("l", "line", true, "Line of input");
        options.addOption("f", "file", true, "file path");

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("Error parsing command line: " + e.getMessage());
            System.exit(1);
        }

        int line = 12;
        String filePath = "treeLines.txt";

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

        double[][] scaledQ = {{-1.299225557498793, 0.3854122453465731, 0.6163580944860869, 0.297455217666133}, {0.2358426020428258, -0.8469806551036532, 0.22579029897487418, 0.38534775408595323}, {0.5826362447422561, 0.3487971468891999, -1.169335005707481, 0.23790161407602503}, {0.2113873951730493, 0.4475211646353517, 0.17885063667912704, -0.8377591964875281}};
        IndependentSiteModel ism = new IndependentSiteModel(scaledQ);
        SiteModel target = new ARMADiLLOSiteModel(1.14002582);

        LinearHamTreeLineParser treeLineParser = new LinearHamTreeLineParser(inputLine);
        TreeNode root = treeLineParser.parseLinearHam();

        HiddenNodes hiddenNodes = new HiddenNodes(1, root, 1, new SiteModel[]{ism, target});
        for (int i = 0; i < 20000; i++) {
            hiddenNodes.mutate();
            if (i % 100 == 0) {
                System.out.println(hiddenNodes);
            }
        }
    }
}