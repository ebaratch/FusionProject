package LEARN_HERE.Agents;

import Framework.GridsAndAgents.AgentGrid3D;
import Framework.GridsAndAgents.AgentSQ3Dunstackable;
import Framework.Gui.GridWindow;
import Framework.Gui.Window3DOpenGL;
import Framework.Util;

import java.util.Random;

/**
 * Created by Rafael on 10/29/2017.
 */
class MoveAgent extends AgentSQ3Dunstackable<Mover3D> {

}

public class Mover3D extends AgentGrid3D<MoveAgent> {
    public Mover3D(int x, int y, int z) {
        super(x, y,z, MoveAgent.class,false,false,false);//set to true for wraparound
    }

    public static void main(String[] args) {
        int GREEN= Util.RGB(0,0.2,0);
        int WHITE= Util.RGB(1,1,1);

        Mover3D test=new Mover3D(10,10,10);
        Window3DOpenGL win3D=new Window3DOpenGL("3D",500,500,test.xDim,test.yDim,test.zDim);
        GridWindow win2D=new GridWindow("2D",test.xDim,test.yDim,20);
        MoveAgent ourHero=test.NewAgentSQ(5,5,5);
        Random rn=new Random();

        for (int i = 0; i < 10000; i++) {
            win2D.TickPause(10);
            ourHero.MoveSafeSQ(ourHero.Xsq()+(rn.nextInt(3)-1),ourHero.Ysq()+(rn.nextInt(3)-1),ourHero.Zsq()+(rn.nextInt(3)-1));//random movement

            win2D.Clear(GREEN);
            win2D.SetPix(ourHero.Xsq(),ourHero.Ysq(), Util.HeatMapRGB(ourHero.Zsq(),test.zDim,0));//draw 2d

            win3D.Clear(GREEN);
            win3D.Circle(ourHero.Xsq(),ourHero.Ysq(),ourHero.Zsq(),0.5,WHITE);//draw 3d
            win3D.Show();
            if(win3D.CheckClosed()){//quit if close button is clicked
                break;
            }
        }
        win3D.Dispose();//destroy guis
        win2D.Dispose();
    }
}
