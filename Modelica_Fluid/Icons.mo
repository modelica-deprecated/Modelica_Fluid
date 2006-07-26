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
  
  partial package BaseClassLibrary "Icon for library" 
    annotation (Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[1, 1],
        component=[20, 20]), Icon(
        Rectangle(extent=[-100, -100; 80, 50], style(
            color=9,
            rgbcolor={175,175,175},
            fillColor=31,
            rgbfillColor={235,255,253},
            fillPattern=1)),
        Polygon(points=[-100, 50; -80, 70; 100, 70; 80, 50; -100, 50], style(
            color=9,
            rgbcolor={175,175,175},
            fillColor=31,
            rgbfillColor={235,255,253},
            fillPattern=1)),
        Polygon(points=[100, 70; 100, -80; 80, -100; 80, 50; 100, 70], style(
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
          extent=[-120, 122; 120, 73],
          string="%name",
          style(color=1)),
        Text(
          extent=[-92,-44; 73,-72],
          style(color=3),
          string="Base classes")));
  end BaseClassLibrary;
end Icons;
