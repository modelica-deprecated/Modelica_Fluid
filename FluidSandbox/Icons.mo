within FluidSandbox;
package Icons "Library of resuable icons" 
  extends Modelica.Icons.Library;
  
  annotation (Documentation(info="<html>
 
</html>"));
  partial class BaseClassLibrary 
    "Same icons as in Modelica_Fluid but a CLASS, which avoids GUI confusion" 
    
    annotation (Coordsys(
        extent=[-100,-100; 100,100],
        grid=[1,1],
        component=[20,20]), Icon(
        Rectangle(extent=[-100,-100; 80,50], style(
            color=9,
            rgbcolor={175,175,175},
            fillColor=31,
            rgbfillColor={235,255,253},
            fillPattern=1)),
        Polygon(points=[-100,50; -80,70; 100,70; 80,50; -100,50], style(
            color=9,
            rgbcolor={175,175,175},
            fillColor=31,
            rgbfillColor={235,255,253},
            fillPattern=1)),
        Polygon(points=[100,70; 100,-80; 80,-100; 80,50; 100,70], style(
            color=9,
            rgbcolor={175,175,175},
            fillColor=31,
            rgbfillColor={235,255,253},
            fillPattern=1)),
        Text(
          extent=[-94,15; 73,-16],
          style(color=3),
          string="Library of"),
        Text(
          extent=[-120,122; 120,73],
          string="%name",
          style(color=1)),
        Text(
          extent=[-92,-44; 73,-72],
          style(color=3),
          string="Base classes")));
  end BaseClassLibrary;
  
  partial class GenericPackage 
    "Icon for generic packages (not specific to any modelling approach)" 
    annotation (Coordsys(
        extent=[-100,-100; 100,100],
        grid=[1,1],
        component=[20,20]), Icon(
        Rectangle(extent=[-100,-100; 80,50], style(
            color=76,
            rgbcolor={159,159,255},
            fillColor=30,
            rgbfillColor={235,235,235},
            fillPattern=1)),
        Polygon(points=[-100,50; -80,70; 100,70; 80,50; -100,50], style(
            color=76,
            rgbcolor={159,159,255},
            fillColor=30,
            rgbfillColor={235,235,235},
            fillPattern=1)),
        Polygon(points=[100,70; 100,-80; 80,-100; 80,50; 100,70], style(
            color=76,
            rgbcolor={159,159,255},
            fillColor=30,
            rgbfillColor={235,235,235},
            fillPattern=1))));
  equation 
    
  end GenericPackage;
  
  partial class GenericVariantPackage 
    "Icon for generic variant packages (not specific to any modelling approach)" 
    annotation (Coordsys(
        extent=[-100,-100; 100,100],
        grid=[1,1],
        component=[20,20]), Icon(
        Rectangle(extent=[-40,-40; 100,100], style(
            color=9,
            rgbcolor={175,175,175},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Rectangle(extent=[-70,-70; 70,70], style(
            color=9,
            rgbcolor={175,175,175},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Rectangle(extent=[-100,-100; 40,40], style(
            color=9,
            rgbcolor={175,175,175},
            fillColor=7,
            rgbfillColor={255,255,255}))));
  equation 
    
  end GenericVariantPackage;
  
  partial class ModellingApproach "Icon for a basic modelling approach" 
    annotation (Coordsys(
        extent=[-100,-100; 100,100],
        grid=[1,1],
        component=[20,20]), Icon(
        Rectangle(extent=[-40,-40; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=68,
            rgbfillColor={190,223,255})),
        Rectangle(extent=[-70,-70; 70,70], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=68,
            rgbfillColor={190,223,255})),
        Rectangle(extent=[-100,-100; 40,40], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=68,
            rgbfillColor={190,223,255}))));
  equation 
    
  end ModellingApproach;
  
  partial class FluidInterface "Icon for a fluid interface implementation" 
    annotation (Coordsys(
        extent=[-100,-100; 100,100],
        grid=[1,1],
        component=[20,20]), Icon(
        Rectangle(extent=[-40,-40; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=30,
            rgbfillColor={255,255,210})),
        Rectangle(extent=[-70,-70; 70,70], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=30,
            rgbfillColor={255,255,210})),
        Rectangle(extent=[-100,-100; 40,40], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=30,
            rgbfillColor={255,255,210}))));
  equation 
    
  end FluidInterface;
  
  partial class FluidDiscretization 
    "Icon for a fluid discretization implementation" 
    annotation (Coordsys(
        extent=[-100,-100; 100,100],
        grid=[1,1],
        component=[20,20]), Icon(
        Rectangle(extent=[-40,-40; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=52,
            rgbfillColor={213,255,170})),
        Rectangle(extent=[-70,-70; 70,70], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=52,
            rgbfillColor={213,255,170})),
        Rectangle(extent=[-100,-100; 40,40], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=52,
            rgbfillColor={213,255,170}))));
  equation 
    
  end FluidDiscretization;
  
  partial class VariantLibrary 
    "Icon for a library that contains several variants of one component" 
    
    annotation (Coordsys(
        extent=[-100,-100; 100,100],
        grid=[1,1],
        component=[20,20]), Icon(
        Rectangle(extent=[-40,-40; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Rectangle(extent=[-70,-70; 70,70], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Rectangle(extent=[-100,-100; 40,40], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Text(
          extent=[-125,158; 115,103],
          string="%name",
          style(color=1))));
  end VariantLibrary;
  
  partial class Examples 
    annotation (Icon(
        Rectangle(extent=[-90,90; 90,-90], style(
            color=58,
            rgbcolor={0,127,0},
            fillColor=59,
            rgbfillColor={118,200,118})),
        Polygon(points=[-28,42; 56,0; -28,-42; -28,42], style(
            pattern=0,
            fillColor=0,
            rgbfillColor={0,0,0}))));
  end Examples;
  
  partial class Example 
    annotation (Icon(
        Ellipse(extent=[-100,100; 100,-100], style(
            color=59,
            rgbcolor={150,255,150},
            gradient=3,
            fillColor=59,
            rgbfillColor={150,255,150})),
        Polygon(points=[-28,42; 56,0; -28,-42; -28,42], style(
            pattern=0,
            fillColor=0,
            rgbfillColor={0,0,0}))));
  end Example;
end Icons;
