import java.util.ArrayList;

public class AminoAcid {

    int enumeration;

    boolean hydrophobic;

    int xCoordinate;
    int yCoordinate;

    AminoAcid previousAA;
    AminoAcid subsequentAA;

    ArrayList<AminoAcid> hydrogenBonds = new ArrayList<>();
    ArrayList<AminoAcid> overlaps = new ArrayList<>();

    AminoAcid(boolean phob, int number, int x, int y, AminoAcid prev)
    {
        hydrophobic = phob;
        enumeration = number;
        xCoordinate = x;
        yCoordinate = y;
        previousAA = prev;
    }

    public void setSubsequentAA(AminoAcid nextAA){
        subsequentAA = nextAA;
    }

}
