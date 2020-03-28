package Models.FusionModelEtienne.FusionModel3DHeterogeneityContactInhibition;

import Framework.Extensions.SphericalAgent3D;
import Framework.GridsAndAgents.AgentGrid3D;
import Framework.Gui.Window3DOpenGL;
//import Framework.Tools.BitGenome;
import Framework.Gui.Window3DOpenGL;
import Framework.Rand;
import Framework.Util;

import java.io.File;
import java.io.FileWriter;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;

import static Framework.Util.*;

class Dish extends AgentGrid3D<Cell> {
    final static int BLACK=RGB(0,0,0),RED=RGB(1,0,0),GREEN=RGB(0,1,0),YELLOW=RGB(1,1,0),BLUE=RGB(0,0,1),WHITE=RGB(1,1,1),CYTOPLASM=RGB256(191,156,147);
    //GLOBAL CONSTANTS
    double DIVISION_PROB=0.67/24;
    double DEATH_PROB=0.12/24;
    //double FUSION_PROB=0.00001;
    //double FUSION_PROB=0.005;
//    double FUSION_PROB=0.00035;
    double FUSION_PROB=0;
//    double FUSION_PROB=0.0015;
//    double FUSION_PROB=0.02;
    //double FUSION_PROB=0.1;

    //BitSet bs = new BitSet();
    //long bits=1+1<<1;
    double CELL_RAD=0.3;
    double MAX_RAD=Math.sqrt(2)*CELL_RAD;
    double FRICTION=0.5;
    double STEADY_STATE_FORCE=0;
    double MAX_STEADY_STATE_LOOPS=2;
    double DIV_RADIUS=CELL_RAD*(2.0/3.0);
    double FORCE_EXPONENT=2;
    double FORCE_SCALER=0.7;
    int GENE_NUMBER=20;
    //double mutationProb=0.00001/24;
    double mutationProb=0.005;
//    int[] WILDTYPE={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int [] WILDTYPE=new int[GENE_NUMBER];
    int MAXSCORE;
    int[] TypeRepresentants;

    //double MAX_VEL=1000000000;

    int fusionCt=0;
    //INTERNAL VARIABLES
    Rand rn=new Rand();
    ArrayList<Cell> cellScratch=new ArrayList<>();
    double[] divCoordScratch=new double[3];

    public Dish(int sideLen,int startingPop,double startingRadius){
        super(sideLen,sideLen,sideLen,Cell.class,true,true,true);
        double[] startCoords=new double[3];
        for (int i=0;i<GENE_NUMBER;i++) {
            WILDTYPE[i]=0;
        }
        MAXSCORE=0;
        for (int i=0;i<GENE_NUMBER;i++) {
            MAXSCORE=MAXSCORE+(int)Math.pow(2,i);
        }
        TypeRepresentants=new int[MAXSCORE];
        for (int i=0;i<MAXSCORE;i++) {
            TypeRepresentants[i]=0;
        }
        for (int i = 0; i < startingPop; i++) {
            rn.RandomPointInCircle(startingRadius, startCoords);
            Cell c=NewAgentPT(startCoords[0]+xDim/2.0,startCoords[1]+yDim/2.0,startCoords[2]+zDim/2.0);
            if(i%2==0) {
                //for (int j=0;j<GENE_NUMBER;j++) {
                //        WILDTYPE[j]=0;
                //    }
                c.Init(RED,WILDTYPE);

                //c.Mutate();
                //for (int j=0;j<GENE_NUMBER;j++){
                //    System.out.println("WILDTYPE="+j+WILDTYPE[j]);
                //}
            }
            else {
                //for (int j=0;j<GENE_NUMBER;j++) {
                //    WILDTYPE[j]=0;
                //}
                c.Init(RED,WILDTYPE);
            }
        }
    }


    //int[] Blendind(int[] Genotype1,int[] Genotype2){
    //    int[] BlendedGenome=new int[GENE_NUMBER];
    //    for (int i=0;i<GENE_NUMBER;i++){
    //        if(rn.nextDouble()<0.5) {
    //            BlendedGenome[i]=Genotype1[i];
    //            System.out.println("vrai");
    //        }
    //        else{
    //            BlendedGenome[i]=Genotype2[i];
    //            System.out.println("faux");
    //        }
    //    }
    //    return BlendedGenome;
    //}


    void Blendind(int[] Genotype1,int[] Genotype2){
//        for (int i=0;i<GENE_NUMBER;i++){
//            System.out.println("entries="+i+Genotype1[i]+Genotype2[i]);
//        }
        //int[] ScratchGenome=new int[GENE_NUMBER];
        int tempBit=0;
        for (int i=0;i<GENE_NUMBER;i++){
            if(rn.Double()<0.5) {
                tempBit = Genotype1[i];
                //System.out.println("vrai");
                Genotype1[i] = Genotype2[i];
                Genotype2[i] = tempBit;
            }
            else{
                //System.out.println("faux");
            }
        }
//        for (int i=0;i<GENE_NUMBER;i++){
//            System.out.println("results="+i+Genotype1[i]+Genotype2[i]);
//        }
    }


    void Fusion(Cell c1,Cell c2){
        if (c1!=null && c1!=null) {
            int[] ScratchGenome1 = c1.Genotype;
            int[] ScratchGenome2 = c2.Genotype;

            Blendind(ScratchGenome1, ScratchGenome2);

            double X1 = c1.Xpt();
            double X2 = c2.Xpt();

            double Y1 = c1.Ypt();
            double Y2 = c2.Ypt();

            double Z1 = c1.Zpt();
            double Z2 = c2.Zpt();

            c1.Dispose();
            if (c2!=c1){
                c2.Dispose();
                Cell H2 = NewAgentPT(X2, Y2,Z2);
                H2.Init(RED, ScratchGenome2);


            }

            Cell H1 = NewAgentPT(X1, Y1, Z1);

            H1.Init(RED, ScratchGenome1);
        }
    }

    int SteadyStateMovement(){
        int loopCt=0;
        while(loopCt<MAX_STEADY_STATE_LOOPS) {
            double maxForce=0;
            for (Cell c : this) {
                maxForce=Math.max(maxForce,c.Observe());
            }
            for (Cell c : this) {
                c.Act();
            }
            loopCt++;
            if(maxForce<STEADY_STATE_FORCE){
                //System.out.println(loopCt+","+maxForce);
                break;
            }
        }
        return loopCt;
    }
    void Step(){
        SteadyStateMovement();
        //for (Cell c:this) {
        //    c.Observe();
        //}
        //for (Cell c:this){
        //    c.Act();
        //}
        for (Cell c:this) {
            c.Step();
        }
        for (Cell c:this){
            fusionCt+=c.Fuse()?1:0;
        }
        //System.out.println(fusionCt);
        IncTick();
    }


    int ConvertBinary(int[] BinaryVector){
        int total=0;
        for (int i=0;i<BinaryVector.length;i++){
            total=total+BinaryVector[i]*(int)Math.pow(2,i);
        }
        return total;
    }


    int Score(int[] BinaryVector){
        int total=0;
        for (int i=0;i<BinaryVector.length;i++){
            if (BinaryVector[i]>0){
                total=total+BinaryVector[i];
            }
        }
        return total;
    }


    double ComputeTotal(){
        int TotalCell=0;

        for (Cell c:this){
            TotalCell=TotalCell+1;
        }

        System.out.println("checkPop="+TotalCell);
        return TotalCell;
    }

    double ComputeShannon(){
        double Shanon=0;
        int GenoInt=0;
        int TotalCell=0;
        double sumProp=0;

        for (Cell c:this){
            TotalCell=TotalCell+1;
            GenoInt=ConvertBinary(c.Genotype);
            TypeRepresentants[GenoInt]=TypeRepresentants[GenoInt]+1;
        }
        for (int i=0;i<MAXSCORE;i++){
            if(TypeRepresentants[i]>0){
                Shanon=Shanon-(TypeRepresentants[i]*1.0/TotalCell)*Math.log(TypeRepresentants[i]*1.0/TotalCell);
                System.out.println("checkProp="+i+(float)TypeRepresentants[i]/TotalCell);
                //System.out.println("PUTAIN DE LOG="+Math.log(TypeRepresentants[i]*1.0/TotalCell));
                sumProp=sumProp+(float)TypeRepresentants[i]/TotalCell;
            }
            TypeRepresentants[i]=0;
        }
        System.out.println("checkPop="+TotalCell+","+sumProp);
        return Shanon;
    }

    double ComputeRichness(){
        int Richness=0;
        int GenoInt=0;
        int TotalCell=0;

        for (Cell c:this){
            GenoInt=ConvertBinary(c.Genotype);
            TypeRepresentants[GenoInt]=TypeRepresentants[GenoInt]+1;
            TotalCell=TotalCell+1;
        }
        for (int i=0;i<MAXSCORE;i++){
            if(TypeRepresentants[i]>0){
                Richness=Richness+1;
            }
            TypeRepresentants[i]=0;
        }

        System.out.println("checkPop="+TotalCell);
        return Richness;
    }


    double ComputeMaxScore() {
        int MaxScore = 0;
        int Score = 0;

        for (Cell c : this) {
            Score = Score(c.Genotype);
            if (Score(c.Genotype) > MaxScore) {
                MaxScore = Score;
            }
        }
        return MaxScore;
    }

}

class Cell extends SphericalAgent3D<Cell,Dish> {
    int color;
    boolean hybrid;
    int[] Genotype;
    //double forceSum;//used with contact inhibition calculation
    double DIV_BIAS;
    double INHIB_WEIGHT;


    void Init(int InitialColor,int[] GenotypeIni){
        radius=G().CELL_RAD;
        DIV_BIAS=0.0279;
        INHIB_WEIGHT=0.45;
        xVel=0;
        yVel=0;
        zVel=0;
        //color=InitialColor;
        hybrid=false;
        Genotype=new int[G().GENE_NUMBER];
        for (int i=0;i<G().GENE_NUMBER;i++){
            Genotype[i]=GenotypeIni[i];
        }
        //for (int i=0;i<G().GENE_NUMBER;i++){
        //    Genotype[i]=0;
        //}
        SetCellColor();
    }
    void Init(int InitialColor,boolean IsHybrid,int[] GenotypeIni){
        xVel=0;
        yVel=0;
        zVel=0;

        //Genotype=new int[G().GENE_NUMBER];
        Genotype=GenotypeIni;
        //for (int i=0;i<G().GENE_NUMBER;i++){
        //    Genotype[i]=0;
        //}
        hybrid=IsHybrid;
        if(hybrid==false){
            radius=G().CELL_RAD;
        }
        else{
            radius=Math.sqrt(2)*G().CELL_RAD;
        }
        //color=InitialColor;
        SetCellColor();
    }

    //void SetCellColor(int newColor){
    //    color=newColor;
    //}

    void SetCellColor(){
        double binary=0;
        binary=G().ConvertBinary(this.Genotype);
//        for (int i=0;i<G().GENE_NUMBER;i++){
//            sum=sum+Genotype[i];
//        }
        color=LongRainbowMap(1-1.0*(binary)/(Math.pow(2,G().GENE_NUMBER)-1));
        //if (sum==0){
        //    color=G().GREEN;
        //}
        //else if (sum==1.0){
        //    color=G().RED;
        //}
        //else if (sum==2.0){
        //    color=G().BLUE;
        //}
    }


    public boolean CanDivide(double div_bias,double inhib_weight){
        return G().rn.Double()<Math.tanh(div_bias-this.Observe()*inhib_weight);
    }

    void Mutate(){
//        for (int i=0;i<G().GENE_NUMBER;i++){
//            System.out.println("entries="+i+Genotype[i]);
//        }
        int indexMut=G().rn.Int(G().GENE_NUMBER);
        if (Genotype[indexMut]==0){
            Genotype[indexMut]=1;
        }
//        else{
//            Genotype[indexMut]=0;
//        }
//        for (int i=0;i<G().GENE_NUMBER;i++){
//            System.out.println("sorties="+i+Genotype[i]);
//        }
        SetCellColor();
    }

    double OverlapToForce(double overlap){
        if(overlap<0){
            return 0;
        }
        return Math.pow(G().FORCE_SCALER*overlap,G().FORCE_EXPONENT);
        //return G().FORCE_SCALER*overlap;
    }
    boolean Fuse(){
        if(hybrid){return false;}
        //listing all cells in the area
        G().cellScratch.clear();
        G().AgentsInRad(G().cellScratch,Xpt(),Ypt(),Zpt(),G().CELL_RAD*2);
        int neighborCt=0;
        //getting valid fusion neighbors
        for (int i=0;i<G().cellScratch.size();i++) {
            Cell c=G().cellScratch.get(i);
            if(!c.hybrid&&c!=this&&Util.DistSquared(Xpt(),Ypt(),Zpt(),c.Xpt(),c.Ypt(),c.Zpt())<G().CELL_RAD*2){
                G().cellScratch.set(neighborCt,c);
                neighborCt++;
            }
        }
        //fusing
        if(neighborCt>0&&G().rn.Double()<Util.ProbScale(G().FUSION_PROB,neighborCt)){
            G().Fusion(this,G().cellScratch.get(G().rn.Int(neighborCt)));
            return true;
        }
        return false;
    }


    double Observe(){

        return SumForces(radius+G().MAX_RAD,G().cellScratch,this::OverlapToForce);
        //return forceSum;
        /*
        if(ret>G().MAX_VEL){
            xVel*=G().MAX_VEL/ret;
            yVel*=G().MAX_VEL/ret;
        }
        return ret;
        */
    }
    void Act(){
        ForceMove();
        ApplyFriction(G().FRICTION);
    }
    void Step(){

        if(CanDivide(this.DIV_BIAS,this.INHIB_WEIGHT)){
            Cell child=Divide(G().DIV_RADIUS,G().divCoordScratch,G().rn);
            if(G().rn.Double()<G().mutationProb){
                child.Init(this.color,this.Genotype);
                child.Mutate();
                //child.MutColorChange(this,G().rn,0.1);
            }
            else{
                child.Init(this.color,this.Genotype);
            }
        }
        if(G().rn.Double()<G().DEATH_PROB){
            Dispose();
            return;
        }

    }
}

public class FusionModelContactInhibition {
    static int SIDE_LEN = 40;
    static int STARTING_POP = 1;
    static double STARTING_RADIUS = 2;
    static int TIMESTEPS = 8760;//(365 Days)
    static float[] circleCoords = Util.GenCirclePoints(1, 10);

    public static void main(String[] args) {
        //TickTimer trt=new TickRateTimer();
        for (int j = 19; j < 20; j++) {
            Window3DOpenGL vis = new Window3DOpenGL("Cell Fusion Visualization", 1000, 1000, SIDE_LEN, SIDE_LEN,SIDE_LEN);
            Dish d = new Dish(SIDE_LEN, STARTING_POP, STARTING_RADIUS);
            double Shanon = 0;
            double Richness=0;
            double MaxMut=0;
            double Total=0;
            //double Shanon = 0;
            //int[] vectorExample={1,0,0,0,1,0};
            //int[] vectorExample2={0,1,1,1,1,0};
            //for (int i=0;i<6;i++){
            //    System.out.println("entries="+i+vectorExample[i]+vectorExample2[i]);
            //}
            //d.Blendind(vectorExample,vectorExample2);
            //for (int i=0;i<6;i++){
            //    System.out.println("results="+i+vectorExample[i]+vectorExample2[i]);
            //}

            //d.SetCellsColor("red");


//        for (int i = 0; i < TIMESTEPS; i++) {
//            System.out.println("Time=" + i);
//            vis.TickPause(0);
//            d.Step();
//            DrawCells(vis,d,"/Users/baratcEA/work/Moffitt/Phoenix/Fusion/ResultsSpatial/Heterogeneity2D/Mutation/",i);
//            Shanon=0;
//            Shanon=d.ComputeShannon();
//
//            //System.out.println("Shanon index=" + Shanon);
//        }
//        System.out.println("Shanon index=" + Shanon);

            try {
                String fileName="Shanon3D"+j+".txt";
                File file = new File(fileName);
                String fileName2="Richness3D"+j+".txt";
                File file2 = new File(fileName2);
//            File file = new File("Shanon.txt");
                PrintWriter printWriter = new PrintWriter(file);
                PrintWriter printWriter2 = new PrintWriter(file2);

                for (int i = 0; i < TIMESTEPS; i++) {
                    System.out.println("Time=" + i);
                    vis.TickPause(0);
                    d.Step();

                    DrawCells(vis, d, "./Results/Heterogeneity/3DContactInhibition/3DMut/", i);

                    Total=d.ComputeTotal();
                    //Shanon = d.ComputeShannon();
//                    Richness=d.ComputeRichness();
                    MaxMut=d.ComputeMaxScore();
                    System.out.println("MaxMut=" + MaxMut);
//                    System.out.println("Richness=" + Richness);
//                    System.out.println("Shanon=" + Shanon);

                    String line = "";
                    line = "" + i + "," + Shanon + "\n";
                    String line2 = "";
                    line2 = "" + i + "," + Richness + "\n";
                    if (i % 100 == 0) {
//                        System.out.println("Shanon index=" + Shanon);
//                        printWriter.print(line);
//                        System.out.println("Richness=" + Richness);
//                        printWriter2.print(line2);
                    }
                }

                printWriter.close();
                printWriter2.close();
            }// end try block
            catch (Exception e) {
                System.out.println(e.getClass());
            }
        }
    }


    static void DrawCells(Window3DOpenGL vis, Dish d, String path, int i){
        vis.Clear(Dish.WHITE);
        for (Cell c:d) {
            //color "cytoplasm"
            vis.CelSphere(c.Xpt(),c.Ypt(),c.Zpt(),c.radius, c.color);
        }
//        for (Cell c:d) {
//            //color "nucleus"
//            vis.CelSphere(c.Xpt(),c.Ypt(),c.Zpt(),c.radius/2.0,c.color);
//        }
        vis.Show();
        if (i%10==0) {
            vis.ToPNG((path.concat(Integer.toString(i))).concat(".png"));
        }
    }
}
