package Models.FusionModelEtienne.GenomeInterface;

public class LongGenome implements Genome{
    long data=0;
    @Override
    public void Set(int i, boolean val) {
        if(val==true){
            data=data|(1l<<i);
        }
        if(val==false){
            data=data&~(1l<<i);
        }
    }

    @Override
    public boolean Get(int i) {
        return ((data>>i)&1)!=0;
    }
}
