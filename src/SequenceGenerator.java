import java.util.ArrayList;

public class SequenceGenerator {

    Examples examples;

    public void createFirstSequence() {

        ArrayList<AminoAcid> AminoSequence = new ArrayList<>();

        AminoAcid aa0 = new AminoAcid(true,0,0,0,null);
        AminoSequence.add(aa0);

        AminoAcid aa1 = new AminoAcid(false,1,0,-1, aa0);
        aa0.setSubsequentAA(aa1);
        AminoSequence.add(aa1);

        AminoAcid aa2 = new AminoAcid(true,2,-1,-1, aa1);
        aa1.setSubsequentAA(aa2);
        AminoSequence.add(aa2);

        AminoAcid aa3 = new AminoAcid(true,3, -1, 0, aa2);
        aa2.setSubsequentAA(aa3);
        AminoSequence.add(aa3);

        AminoAcid aa4 = new AminoAcid(false, 4, -2, 0, aa3);
        aa3.setSubsequentAA(aa4);
        AminoSequence.add(aa4);

        AminoAcid aa5 = new AminoAcid(false,5,-2,-1, aa4);
        aa4.setSubsequentAA(aa5);
        AminoSequence.add(aa5);

        AminoAcid aa6 = new AminoAcid(false,6,-2,-2, aa5);
        aa5.setSubsequentAA(aa6);
        AminoSequence.add(aa6);

        AminoAcid aa7 = new AminoAcid(true,7,-1,-2, aa6);
        aa6.setSubsequentAA(aa7);
        AminoSequence.add(aa7);

        //ass.storeAminoSequence(AminoSequence);

    }

    public void createSecondSequence() {

        ArrayList<AminoAcid> AminoSequence = new ArrayList<>();

        AminoAcid aa0 = new AminoAcid(true,0,0,0,null);
        AminoSequence.add(aa0);

        AminoAcid aa1 = new AminoAcid(false,1,0,-1, aa0);
        aa0.setSubsequentAA(aa1);
        AminoSequence.add(aa1);

        AminoAcid aa2 = new AminoAcid(true,2,-1,-1, aa1);
        aa1.setSubsequentAA(aa2);
        AminoSequence.add(aa2);

        AminoAcid aa3 = new AminoAcid(true,3,-1,0, aa2);
        aa2.setSubsequentAA(aa3);
        AminoSequence.add(aa3);

        AminoAcid aa4 = new AminoAcid(false,4,-2,0, aa3);
        aa3.setSubsequentAA(aa4);
        AminoSequence.add(aa4);

        AminoAcid aa5 = new AminoAcid(false,5,-2,1, aa4);
        aa4.setSubsequentAA(aa5);
        AminoSequence.add(aa5);

        AminoAcid aa6 = new AminoAcid(false,6,-1,1, aa5);
        aa5.setSubsequentAA(aa6);
        AminoSequence.add(aa6);

        AminoAcid aa7 = new AminoAcid(true,7,-1,0, aa6);
        aa6.setSubsequentAA(aa7);
        AminoSequence.add(aa7);

        //ass.storeAminoSequence(AminoSequence);
    }

    public ArrayList<AminoAcid> createRandomSequence(int length)
    {
        ArrayList<AminoAcid> AminoSequence = new ArrayList<>();


        AminoAcid aa0 = new AminoAcid(getRandomBoolean(),0,0,0,null);
        AminoSequence.add(aa0);

        AminoAcid prevAa = aa0;

        for(int i=1; i<length; i++)
        {
            AminoAcid aa = new AminoAcid(getRandomBoolean(),i,i,0,prevAa);
            AminoSequence.add(aa);
            prevAa = aa;

            if(!AminoSequence.isEmpty())
            {
                AminoSequence.get(i-1).setSubsequentAA(aa);
            }
        }
        aa0.setSubsequentAA(AminoSequence.get(1));
        return AminoSequence;

    }

    public boolean getRandomBoolean()
    {
        return Math.random() > 0.5;
    }

    public ArrayList<AminoAcid> createNewFolding()
    {
        String SEQ64 = "1111111111110101001100110010011001100100110011001010111111111111";

        ArrayList<AminoAcid> AminoSequence = new ArrayList<>();

        AminoAcid aa0 = new AminoAcid(getInteger(SEQ64.charAt(0)),0,0,0,null);
        AminoSequence.add(aa0);

        AminoAcid prevAa = aa0;
        for(int i=1; i<SEQ64.length(); i++)
        {
            AminoAcid aa = new AminoAcid(getInteger(SEQ64.charAt(i)),i,i,0,prevAa);
            AminoSequence.add(aa);
            prevAa = aa;

            if(!AminoSequence.isEmpty())
                {
                AminoSequence.get(i-1).setSubsequentAA(aa);
                }
        }
        aa0.setSubsequentAA(AminoSequence.get(1));

        return AminoSequence;
    }

    public boolean getInteger(char string)
    {
        int value = Integer.parseInt(String.valueOf(string));
        return value > 0;
    }


}
