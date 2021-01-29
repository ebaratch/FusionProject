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

//Tissue Class
class Tissue extends AgentGrid3D<Cell> {
    //Tissue class parameters
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
    int [] WILDTYPE=new int[GENE_NUMBER];
    int MAXSCORE;
    int[] TypeRepresentants;



    int fusionCt=0;
    //INTERNAL VARIABLES
    Rand rn=new Rand();
    ArrayList<Cell> cellScratch=new ArrayList<>();
    double[] divCoordScratch=new double[3];

    public Tissue(int sideLen,int startingPop,double startingRadius){
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

                c.Init(RED,WILDTYPE);

            }
            else {

                c.Init(RED,WILDTYPE);
            }
        }
    }


    //Genomes recombination function
    void Blendind(int[] Genotype1,int[] Genotype2){

        int tempBit=0;
        for (int i=0;i<GENE_NUMBER;i++){
            if(rn.Double()<0.5) {
                tempBit = Genotype1[i];
                Genotype1[i] = Genotype2[i];
                Genotype2[i] = tempBit;
            }
        }
    }

    //Fusion function
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
                break;
            }
        }
        return loopCt;
    }
    void Step(){
        SteadyStateMovement();

        for (Cell c:this) {
            c.Step();
        }
        for (Cell c:this){
            fusionCt+=c.Fuse()?1:0;
        }

        IncTick();
    }


    //Function Converting a bit array into vector into the numerical value in base 2
    int ConvertBinary(int[] BinaryVector){
        int total=0;
        for (int i=0;i<BinaryVector.length;i++){
            total=total+BinaryVector[i]*(int)Math.pow(2,i);
        }
        return total;
    }


    //Function computing the amount of mutations in a genome vector
    int Score(int[] BinaryVector){
        int total=0;
        for (int i=0;i<BinaryVector.length;i++){
            if (BinaryVector[i]>0){
                total=total+BinaryVector[i];
            }
        }
        return total;
    }


    //Function computing Total number of cancer cells
    double ComputeTotal(){
        int TotalCell=0;

        for (Cell c:this){
            TotalCell=TotalCell+1;
        }

        System.out.println("checkPop="+TotalCell);
        return TotalCell;
    }


    //Function computing Shanon Diversity Index of the population
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
                sumProp=sumProp+(float)TypeRepresentants[i]/TotalCell;
            }
            TypeRepresentants[i]=0;
        }
        System.out.println("checkPop="+TotalCell+","+sumProp);
        return Shanon;
    }


    //Function computing population richness
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


    //Function computing Maximal amount of mutation in a cell
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


//Cell Class
class Cell extends SphericalAgent3D<Cell,Tissue> {
    //Cell attributes
    int color;
    boolean hybrid;
    int[] Genotype;
    double DIV_BIAS;
    double INHIB_WEIGHT;


    //Cell Initializing function
    void Init(int InitialColor,int[] GenotypeIni){
        radius=G().CELL_RAD;
        DIV_BIAS=0.0279;
        INHIB_WEIGHT=0.45;
        xVel=0;
        yVel=0;
        zVel=0;

        hybrid=false;
        Genotype=new int[G().GENE_NUMBER];
        for (int i=0;i<G().GENE_NUMBER;i++){
            Genotype[i]=GenotypeIni[i];
        }

        SetCellColor();
    }
    void Init(int InitialColor,boolean IsHybrid,int[] GenotypeIni){
        xVel=0;
        yVel=0;
        zVel=0;

        Genotype=GenotypeIni;

        hybrid=IsHybrid;
        if(hybrid==false){
            radius=G().CELL_RAD;
        }
        else{
            radius=Math.sqrt(2)*G().CELL_RAD;
        }

        SetCellColor();
    }


    //Cell color setting function
    void SetCellColor(){
        double binary=0;
        binary=G().ConvertBinary(this.Genotype);

        color=LongRainbowMap(1-1.0*(binary)/(Math.pow(2,G().GENE_NUMBER)-1));

    }


    //Returns whether cell divides or not
    public boolean CanDivide(double div_bias,double inhib_weight){
        return G().rn.Double()<Math.tanh(div_bias-this.Observe()*inhib_weight);
    }


    //Function changing a random 0 into a 1 in the genome vector
    void Mutate(){

        int indexMut=G().rn.Int(G().GENE_NUMBER);
        if (Genotype[indexMut]==0){
            Genotype[indexMut]=1;
        }

        SetCellColor();
    }


    //Function changing a random 0 into a 1 in the genome vector
    double OverlapToForce(double overlap){
        if(overlap<0){
            return 0;
        }
        return Math.pow(G().FORCE_SCALER*overlap,G().FORCE_EXPONENT);
    }


    //Function making cell agent fuse with a random neighbor and then divide
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


    //Returns Repulsion force exerted on the cells by neighbors
    double Observe(){
        return SumForces(radius+G().MAX_RAD,G().cellScratch,this::OverlapToForce);
    }


    //Applies force and friction to cell agent
    void Act(){
        ForceMove();
        ApplyFriction(G().FRICTION);
    }

    //Step function applying division, mutation and death
    void Step(){

        if(CanDivide(this.DIV_BIAS,this.INHIB_WEIGHT)){
            Cell child=Divide(G().DIV_RADIUS,G().divCoordScratch,G().rn);
            if(G().rn.Double()<G().mutationProb){
                child.Init(this.color,this.Genotype);
                child.Mutate();
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


//Model Run class
public class FusionModelContactInhibition {
    //Simulation parameters: Domaine size, starting ppulation, number of time steps
    static int SIDE_LEN = 40;
    static int STARTING_POP = 1;
    static double STARTING_RADIUS = 2;
    static int TIMESTEPS = 8760;//(365 Days)
    static float[] circleCoords = Util.GenCirclePoints(1, 10);


    //Model Running function
    public static void main(String[] args) {
        for (int j = 19; j < 20; j++) {
            Window3DOpenGL vis = new Window3DOpenGL("Cell Fusion Visualization", 1000, 1000, SIDE_LEN, SIDE_LEN,SIDE_LEN);
            Tissue t = new Tissue(SIDE_LEN, STARTING_POP, STARTING_RADIUS);
            double Shanon = 0;
            double Richness=0;
            double MaxMut=0;
            double Total=0;

            try {
                String fileName="Shanon3D"+j+".txt";
                File file = new File(fileName);
                String fileName2="Richness3D"+j+".txt";
                File file2 = new File(fileName2);

                //Results files
                PrintWriter printWriter = new PrintWriter(file);
                PrintWriter printWriter2 = new PrintWriter(file2);

                //Time Loop
                for (int i = 0; i < TIMESTEPS; i++) {
                    System.out.println("Time=" + i);
                    vis.TickPause(0);
                    t.Step();

                    DrawCells(vis, t, "./Results/Heterogeneity/3DContactInhibition/3DMut/", i);

                    Total=t.ComputeTotal();
                    MaxMut=t.ComputeMaxScore();
                    System.out.println("MaxMut=" + MaxMut);

                    String line = "";
                    line = "" + i + "," + Shanon + "\n";
                    String line2 = "";
                    line2 = "" + i + "," + Richness + "\n";
                    if (i % 100 == 0) {

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


    //Cells Drawing Function
    static void DrawCells(Window3DOpenGL vis, Tissue t, String path, int i){
        vis.Clear(Tissue.WHITE);
        for (Cell c:t) {
            //color "cytoplasm"
            vis.CelSphere(c.Xpt(),c.Ypt(),c.Zpt(),c.radius, c.color);
        }

        vis.Show();
        if (i%10==0) {
            vis.ToPNG((path.concat(Integer.toString(i))).concat(".png"));
        }
    }
}
