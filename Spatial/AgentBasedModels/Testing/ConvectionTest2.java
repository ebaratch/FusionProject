package Testing;

import Framework.GridsAndAgents.PDEGrid2D;
import Framework.Gui.GuiGrid;
import Framework.Gui.GuiLabel;
import Framework.Gui.GuiWindow;
import Framework.Gui.TickTimer;
import Framework.Util;

public class ConvectionTest2 {
    static int sideLen=1000;
    static int steps=100000;
    public static void Setup(PDEGrid2D g, double[]diffRates){
        for (int x = 1; x < 11; x++) {
            for (int y = 1; y < 11; y++) {
                g.Set(x,y,1);
            }
        }
        for (int i = 0; i < diffRates.length; i++) {
            if(g.ItoX(i)>g.xDim/2) {
                diffRates[i] = 0.1;
            }
            else{
                diffRates[i]=0.2;
            }
        }
    }
    public static void main(String[] args) {
        PDEGrid2D g=new PDEGrid2D(sideLen,sideLen,false,false);
        GuiWindow win=new GuiWindow("Advection1stOrder Example",true);
        GuiGrid vis=new GuiGrid(sideLen,sideLen,1);
        GuiLabel lbl=new GuiLabel("tick");
        TickTimer trt=new TickTimer();
        double[]diffRates=new double[g.length];
        win.AddCol(0, lbl);
        win.AddCol(0, vis);
        win.RunGui();
        Setup(g,diffRates);
        for (int i = 0; i < steps; i++) {
            trt.TickPause(0);
            //g.Advection1stOrder(0.01,0);
            g.Diffusion(0.1,1);
            //g.MultiThread(4,(x,y,j)->{
            //    Util.Diffusion2(x,y,g.GetField(),g.GetSwapField(),g.xDim,g.yDim,0.1,true,1,false,false);
            //});
            g.Diffusion(0.1,1);
            vis.DrawGridDiff(g,(val)->{
                return Util.HeatMapRBG(Util.ScaleMinToMax(val, (double) 0, (double) 1));
            });
            g.IncTick();
            lbl.SetText("tick "+g.GetTick());
        }
    }
}
