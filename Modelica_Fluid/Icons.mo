package Icons 
  extends Modelica.Icons.Library;
  partial class VariantLibrary 
    "Icon for a library that contains several variants of one component" 
    
    annotation (Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[1, 1],
        component=[20, 20]), Icon(
        Rectangle(extent=[-40, -40; 100, 100], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Rectangle(extent=[-70, -70; 70, 70], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Rectangle(extent=[-100, -100; 40, 40], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Text(
          extent=[-125, 158; 115, 103],
          string="%name",
          style(color=1))));
  end VariantLibrary;
end Icons;
