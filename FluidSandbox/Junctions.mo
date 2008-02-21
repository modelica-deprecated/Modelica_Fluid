within FluidSandbox;
package Junctions "Junction and splitter models" 
  extends Icons.VariantLibrary;
  
  model IdealJunction 
    "Volume with inlet and outlet ports (flow reversal is allowed)" 
    
    extends Interfaces.PartialComponent;
    
    extends FluidInterface.PartialIdealJunction;
    
    annotation (Icon(
        Rectangle(extent=[-100,20; 100,-20], style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[-20,100; 20,18], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=1,
            fillColor=8,
            rgbfillColor={192,192,192})),
        Rectangle(extent=[-100,16; 100,-16], style(
            color=69,
            gradient=2,
            fillColor=69)),
        Rectangle(extent=[-16,100; 16,14], style(
            color=69,
            rgbcolor={0,127,255},
            gradient=1,
            fillColor=69,
            rgbfillColor={0,127,255})),
        Rectangle(extent=[-6,94; 6,4], style(
            color=69,
            rgbcolor={0,127,255},
            fillColor=69,
            rgbfillColor={0,127,255}))));
    
  end IdealJunction;
end Junctions;
