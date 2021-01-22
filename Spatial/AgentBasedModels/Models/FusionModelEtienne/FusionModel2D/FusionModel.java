package Models.FusionModelEtienne.FusionModel2D;

import Framework.Extensions.SphericalAgent2D;
import Framework.GridsAndAgents.AgentGrid2D;
import Framework.Gui.Window2DOpenGL;
import Framework.Rand;
import Framework.Util;

import java.awt.*;
import java.util.ArrayList;
import java.util.Random;

import static Framework.Util.*;

class Dish extends AgentGrid2D<Cell> {
    final static int BLACK=RGB(0,0,0),RED=RGB(1,0,0),GREEN=RGB(0,1,0),YELLOW=RGB(1,1,0),BLUE=RGB(0,0,1),WHITE=RGB(1,1,1),GREY=RGB(0.5,0.5,0.5);
    //GLOBAL CONSTANTS
    double DEATH_PROB=0.12/24;
    double DIVISION_PROB=DEATH_PROB+0.67/24;
    double FUSION_PROB=0.00035;
    //double FUSION_PROB=0.00011;
    //Spaceless fusion rates [0.0015, 0.0015, 0.0004,0.00045,0.00015,0.00015]
    //Fusion rates [0.00035, 0.00035, 0.000085, 0.00009, 0.00004, 0.00004, 0.0001, 0.00043, 0.00004, 0.00013, 0.0002, 0.00021, 0.00006, 0.000025, 0.00003]
    //Death rates [0.12/24, 0.075/24, 0.06/24, 0.06/24, 0.4/24, 0.4/24, 0.5/24, 2.5/24, 0.75/24, 1.2/24, 0.21/24, 0.2/24, 0.22/24, 0.17/24, 0.2/24]
    double CELL_RAD=0.3;
    double MAX_RAD=Math.sqrt(2)*CELL_RAD;
    double FRICTION=0.8;
    double STEADY_STATE_FORCE=0;
    double MAX_STEADY_STATE_LOOPS=100;
    double DIV_RADIUS=CELL_RAD*(2.0/3.0);
    double FORCE_EXPONENT=2;
    double FORCE_SCALER=0.7;
    int fusionCt=0;
    //INTERNAL VARIABLES
    Rand rn=new Rand();
    ArrayList<Cell> cellScratch=new ArrayList<>();
    double[] divCoordScratch=new double[2];

    public Dish(int sideLen,int startingRed,int startingGreen,double startingRadius){
        super(sideLen,sideLen,Cell.class);
        double[] startCoords=new double[2];
        for (int i = 0; i < startingRed; i++) {
            rn.RandomPointInCircle(startingRadius, startCoords);
            Cell c=NewAgentPT(startCoords[0]+xDim/2.0,startCoords[1]+yDim/2.0);
            c.Init(RED);
        }
        for (int i = 0; i < startingGreen; i++) {
            rn.RandomPointInCircle(startingRadius, startCoords);
            Cell c=NewAgentPT(startCoords[0]+xDim/2.0,startCoords[1]+yDim/2.0);
            c.Init(GREEN);
        }
    }
    void Fusion(Cell c1,Cell c2){
        if (!c1.dead && !c2.dead && !c1.hybrid && !c2.hybrid) {
            Cell H = NewAgentPT((c1.Xpt() + c2.Xpt()) / 2, (c1.Ypt() + c2.Ypt()) / 2);
            if (c1.color == GREEN && c2.color == GREEN) {
                H.Init(GREEN, true);
            } else if ((c1.color == GREEN && c2.color == RED) || (c1.color == RED && c2.color == GREEN)) {
                H.Init(YELLOW, true);
            } else if (c1.color == RED && c2.color == RED) {
                H.Init(RED, true);
            }
            c1.Dispose();
            c2.Dispose();
        }
    }

    void SpacelessFusion(Cell c1,Cell c2){
        if (!c1.dead && !c2.dead && !c1.hybrid && !c2.hybrid) {
            Cell H = NewAgentPT(c1.Xpt(), c1.Ypt());
            if (c1.color == GREEN && c2.color == GREEN) {
                H.Init(GREEN, true);
            } else if ((c1.color == GREEN && c2.color == RED) || (c1.color == RED && c2.color == GREEN)) {
                H.Init(YELLOW, true);
            } else if (c1.color == RED && c2.color == RED) {
                H.Init(RED, true);
            }
            c1.Dispose();
            c2.Dispose();
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

        cellScratch.clear();
        for (Cell c:this) {
            c.Step();
        }
        for (Cell c:this){
            fusionCt+=c.Fuse()?1:0;
//            fusionCt+=c.SpacelessFuse()?1:0;
        }
        IncTick();
    }


}

class Cell extends SphericalAgent2D<Cell,Dish> {
    int color;
    boolean hybrid;
    boolean dead;

    void Init(int InitialColor){
        radius=G().CELL_RAD;
        xVel=0;
        yVel=0;
        color=InitialColor;
        hybrid=false;
        dead=false;
    }
    void Init(int InitialColor,boolean IsHybrid){
        xVel=0;
        yVel=0;
        color=InitialColor;
        hybrid=IsHybrid;
        dead=false;
        if(hybrid==false){
            radius=G().CELL_RAD;
        }
        else{
            radius=Math.sqrt(2)*G().CELL_RAD;
        }
    }

    void SetCellColor(int newColor){
        color=newColor;
    }

    double OverlapToForce(double overlap){
        if(overlap<0){
            return 0;
        }
        return Math.pow(G().FORCE_SCALER*overlap,G().FORCE_EXPONENT);
    }
    boolean Fuse(){
        if(hybrid || dead){return false;}
        //listing all cells in the area
        G().cellScratch.clear();
        G().AgentsInRad(G().cellScratch,Xpt(),Ypt(),G().CELL_RAD*2);
        int neighborCt=0;
        //getting valid fusion neighbors

        for (int i=0;i<G().cellScratch.size();i++) {
            Cell c=G().cellScratch.get(i);
            if(!c.hybrid&&!c.dead&&c!=this&&Util.DistSquared(Xpt(),Ypt(),c.Xpt(),c.Ypt())<G().CELL_RAD*2){
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

    boolean SpacelessFuse(){
        if(hybrid || dead){return false;}

        //fusing
        if(G().rn.Double()<G().FUSION_PROB){
            G().SpacelessFusion(this,G().RandomAgent2(G().rn));
            return true;
        }
        return false;
    }

    boolean Die(){
        this.dead=true;
        this.color=RGB(0.5,0.5,0.5);

        return dead;
    }

    double Observe(){

        return SumForces(radius+G().MAX_RAD,G().cellScratch,this::OverlapToForce);

    }
    void Act(){
        ForceMove();
        ApplyFriction(G().FRICTION);
    }
    void Step() {
        if (dead == false) {
            if (G().rn.Double() < G().DIVISION_PROB && hybrid == false) {
                Cell child = Divide(G().DIV_RADIUS, G().divCoordScratch, G().rn);
                child.Init(this.color);
                Init(this.color);
            }
            if (G().rn.Double() < G().DEATH_PROB && hybrid == false) {
                Die();
            }
        }
    }
}

public class FusionModel {
    static int SIDE_LEN=100;
    static int STARTING_RED=500;
    static int STARTING_GREEN=500;
    static double STARTING_RADIUS=7;
    static int TIMESTEPS=24*3;
    static float[] circleCoords=Util.GenCirclePoints(1,10);
    public static void main(String[] args) {
        int hybridCount = 0;
        int totalCount = 0;
        int totalLiveCount=0;
        int redCount = 0;
        int greenCount = 0;
        int yellowCount = 0;
        int deadCount= 0 ;
        String outputh_path="Results/FusionGreenRed/";
        Window2DOpenGL vis = new Window2DOpenGL("Cell Fusion Visualization", 1000, 1000, SIDE_LEN, SIDE_LEN);
        Dish d = new Dish(SIDE_LEN, STARTING_RED, STARTING_GREEN, STARTING_RADIUS);

        for (int i = 0; i < TIMESTEPS; i++) {
            vis.TickPause(0);
            d.Step();
            DrawCells(vis, d,i,outputh_path);
        }
        for (Cell c : d) {
            totalCount = totalCount + 1;
            totalLiveCount = totalLiveCount + 1;
            if (c.dead){
                deadCount=deadCount+1;
                totalLiveCount=totalLiveCount-1;
            }
            if (c.hybrid == true) {
                hybridCount = hybridCount + 1;
                if (c.color == Dish.YELLOW) {
                    yellowCount = yellowCount + 1;
                }
            } else {
                if (c.color == Dish.RED) {
                    redCount = redCount + 1;
                }
                if (c.color == Dish.GREEN) {
                    greenCount = greenCount + 1;
                }
            }
        }

        System.out.println("Dead cells prct="+(float) 100 * deadCount * 1.0 / totalCount);
        System.out.println("Hybrid cells prct="+(float) 100 * hybridCount * 1.0 / totalLiveCount);
        System.out.println("Yellow cells prct="+(float) 100 * yellowCount * 1.0 / totalLiveCount);
        System.out.println("Red cells prct="+(float) 100 * redCount * 1.0 / totalLiveCount);
        System.out.println("Green cells prct="+(float) 100 * greenCount * 1.0 / totalLiveCount);
        System.out.println("Total cells prct="+(float) 100 * (hybridCount + redCount + greenCount) / totalLiveCount);


        System.out.println("Dead cells="+deadCount);
        System.out.println("Hybrid cells="+hybridCount);
        System.out.println("Yellow cells="+yellowCount);
        System.out.println("Red cells="+redCount);
        System.out.println("Green cells="+greenCount);
        System.out.println("Total live cells="+(hybridCount + redCount + greenCount));
        System.out.println("Total cells="+(hybridCount + redCount + greenCount + deadCount));

    }

    static void DrawCells(Window2DOpenGL vis, Dish d,int i, String path){
        vis.Clear(Dish.WHITE);
        for (Cell c:d) {
            //color "cytoplasm"
            vis.Circle(c.Xpt(),c.Ypt(),c.radius,c.color);
        }
        for (Cell c:d) {
            //color "nucleus"
            vis.Circle(c.Xpt(),c.Ypt(),c.radius/2.0, Dish.BLUE);
        }
        vis.Show();
        vis.ToPNG(path.concat(Integer.toString(i)).concat(".png"));
    }
}
