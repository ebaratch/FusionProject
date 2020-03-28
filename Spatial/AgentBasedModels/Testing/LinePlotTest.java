package Testing;

import Framework.Gui.GridWindow;

import static Framework.Util.RGB;

/**
 * Created by Rafael on 10/28/2017.
 */
public class LinePlotTest {
    public static void main(String[] args) {
        GridWindow win = new GridWindow(200,200,1);
        win.PlotLine(new double[]{10,8,24},new double[]{30,46,40}, RGB(1,0,0), 0,1,0.5);
    }
}
