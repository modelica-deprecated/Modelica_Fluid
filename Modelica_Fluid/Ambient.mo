model Ambient "Ambient field component" 
  import Modelica.SIunits.Conversions;
  
  parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure 
    default_p_ambient =                                                                  101325 
    "Default ambient pressure" 
      annotation(Dialog(group="Defaults"));
  parameter Modelica.Media.Interfaces.PartialMedium.Temperature 
    default_T_ambient=
      Conversions.from_degC(20) "Default ambient temperature" 
      annotation(Dialog(group="Defaults"));
  parameter SI.Acceleration g=Modelica.Constants.g_n 
    "Constant gravity acceleration"                                                  annotation(Dialog(group="Environment"));
  
  annotation (
    preferedView="info",
    defaultComponentName="ambient",
    defaultComponentPrefixes="inner",
    missingInnerMessage="Your model is using an outer \"ambient\" component. An inner \"ambient\" component is not defined. For simulation drag Modelica_Fluid.Components.Ambient into your model and specify ambient conditions.",
    Icon(
      Rectangle(extent=[-100,100; 100,-100], style(
          color=3,
          rgbcolor={0,0,255},
          fillColor=7,
          rgbfillColor={255,255,255})),
      Text(
        extent=[-160,160; 160,110],
        style(color=3, rgbcolor={0,0,255}),
        string="%name"),
      Line(points=[-86,-30; 82,-30], style(color=0, rgbcolor={0,0,0})),
      Line(points=[-82,-68; -52,-30], style(color=0, rgbcolor={0,0,0})),
      Line(points=[-48,-68; -18,-30], style(color=0, rgbcolor={0,0,0})),
      Line(points=[-14,-68; 16,-30], style(color=0, rgbcolor={0,0,0})),
      Line(points=[22,-68; 52,-30], style(color=0, rgbcolor={0,0,0})),
      Line(points=[74,84; 74,14], style(color=0, rgbcolor={0,0,0})),
      Polygon(points=[60,14; 88,14; 74,-18; 60,14], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Text(
        extent=[16,20; 60,-18],
        style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1),
        string="g"),
        Text(
          extent=[-90,82; 74,50],
          style(color=0, rgbcolor={0,0,0}),
          string="defaults")),
    Diagram,
    Documentation(info="<HTML>
<p>
This models defines <b>default ambient conditions</b> (such as ambient pressure) 
for all components that are on the same 
or on a lower level as this component, as well as the constant 
<b>gravity acceleration</b> (on earth). Dragging this component in a model results
in the following declaration:
</p>
<pre>
   <b>inner</b> Modelica_Fluid.Components.Ambient ambient;
</pre>
<p>
The parameters of this instance can be 
then accessed via a corresponding outer declaration:
</p>
<pre>
   <b>outer</b> Modelica_Fluid.Components.Ambient ambient;
</pre>
 
</HTML>
"));
  
    Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape dummyShape 
    "Just temporarily, to force Dymola to open an animation window (only then animation setup is available for diagram animation)"
      annotation (extent=[-60,20; -40,40]);
end Ambient;
