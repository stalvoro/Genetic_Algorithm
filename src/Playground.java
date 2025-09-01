import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.HashMap;

public class Playground extends JPanel {

    private ArrayList<AminoAcid> aminoSequence;

    public Playground(ArrayList<AminoAcid> aminoSequence) {
        this.aminoSequence = aminoSequence;
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);

        // Calculate the bounding box for the amino acid coordinates
        int minX = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;
        int minY = Integer.MAX_VALUE;
        int maxY = Integer.MIN_VALUE;

        for (AminoAcid aa : aminoSequence) {
            if (aa.xCoordinate < minX) minX = aa.xCoordinate;
            if (aa.xCoordinate > maxX) maxX = aa.xCoordinate;
            if (aa.yCoordinate < minY) minY = aa.yCoordinate;
            if (aa.yCoordinate > maxY) maxY = aa.yCoordinate;
        }

        // Calculate the width and height of the bounding box
        int boxWidth = (maxX - minX + 1) * 50;
        int boxHeight = (maxY - minY + 1) * 50;

        // Calculate the offsets to center the drawing in the panel
        int panelWidth = getWidth();
        int panelHeight = getHeight();
        int xOffset = (panelWidth - boxWidth) / 2 - minX * 50;
        int yOffset = (panelHeight - boxHeight) / 2 - minY * 50;

        // A map to store overlapping nodes by their coordinates
        HashMap<Point, ArrayList<AminoAcid>> positionMap = new HashMap<>();

        // Populate the map with amino acids grouped by their coordinates
        for (AminoAcid aa : aminoSequence) {
            Point p = new Point(aa.xCoordinate, aa.yCoordinate);
            positionMap.putIfAbsent(p, new ArrayList<>());
            positionMap.get(p).add(aa);
        }

        // First, draw the connecting lines for each amino acid
        for (AminoAcid aa : aminoSequence) {
            // Get current coordinates
            int x = aa.xCoordinate * 50 + xOffset;
            int y = aa.yCoordinate * 50 + yOffset;

            // Draw a line to the previous amino acid
            if (aa.previousAA != null) {
                int prevX = aa.previousAA.xCoordinate * 50 + xOffset;
                int prevY = aa.previousAA.yCoordinate * 50 + yOffset;
                g.setColor(Color.BLACK);
                g.drawLine(x, y, prevX, prevY);
            }
        }

        // Then, draw the nodes and stack the enumeration numbers for overlapping nodes
        for (Point p : positionMap.keySet()) {
            ArrayList<AminoAcid> aminoAcidsAtSamePosition = positionMap.get(p);

            // Scale and shift coordinates for drawing
            int x = p.x * 50 + xOffset -20;
            int y = p.y * 50 + yOffset -20;

            // Draw the circle (node) once
            AminoAcid firstAA = aminoAcidsAtSamePosition.get(0);
            if (firstAA.hydrophobic) {
                g.setColor(Color.BLACK);
            } else {
                g.setColor(Color.WHITE);
            }
            g.fillOval(x, y, 40, 40); // Draw a filled circle

            // Draw the border of the circle
            g.setColor(Color.BLACK);
            g.drawOval(x, y, 40, 40);

            // Stack the enumeration numbers
            g.setColor(Color.RED);
            int numYOffset = 15; // Vertical offset between numbers
            for (int i = 0; i < aminoAcidsAtSamePosition.size(); i++) {
                int numY = y + numYOffset + (i * 12); // Adjust the Y position for each number
                g.drawString(String.valueOf(aminoAcidsAtSamePosition.get(i).enumeration), x + 15, numY);
            }
        }
    }
}
