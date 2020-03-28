package LEARN_HERE.Diffusbiles;

import Framework.GridsAndAgents.PDEGrid3D;
import Framework.Gui.GridWindow;
import Framework.Util;

/**
 * Created by Rafael on 10/28/2017.
 */
public class DiffusionAdvection3D {
    public static void main(String[] args) {
        GridWindow win = new GridWindow("advection",20,20,10,true);
        PDEGrid3D grid=new PDEGrid3D(20,20,20,true,true,true);//last booleans are for wraparound
        grid.Set(grid.xDim/2,grid.yDim/2,grid.zDim/2,1);
        while(true){
            win.TickPause(100);
            grid.Advection(0.1,0,0.1);//advection
            grid.Diffusion(0.13);//diffusion
            win.DrawGridDiffXY(grid, (val)->(Util.HeatMapBGR(val*1000)));
            //win.DrawGridDiffXZ(grid, (val)->(Util.HeatMapBGR(val*1000)));//uncomment to view from different angles
            //win.DrawGridDiffYZ(grid, (val)->(Util.HeatMapBGR(val*1000)));
        }
    }
}
