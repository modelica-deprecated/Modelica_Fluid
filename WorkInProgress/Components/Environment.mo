model Environment "Environment conditions and global options" 
  import SI = Modelica.SIunits;
  import Modelica_Fluid.WorkInProgress.Types;
  
  replaceable package Medium = PackageMedium extends 
    Modelica.Media.Interfaces.PartialMedium "Ambient medium model" 
     annotation (choicesAllMatching=true);
  parameter Medium.AbsolutePressure p = Medium.reference_p "Ambient pressure";
  parameter Boolean use_T = true "Use T if true, otherwise h" annotation(Dialog(Evaluate=true));
  parameter Medium.Temperature T = if use_T then 293.15 else 
            Medium.T_phX(p,h,X[1:Medium.nXi]) "Ambient temperature" annotation(Dialog(enable=use_T));
  parameter Medium.SpecificEnthalpy h=
            if use_T then Medium.h_pTX(p, T, X[1:Medium.nXi]) else 1e4 
    "Ambient specific enthalpy" annotation(Dialog(enable=not use_T));
  parameter Medium.MassFraction X[Medium.nX] = Medium.reference_X 
    "Ambient mass fractions m_i/m"  annotation (Dialog(enable=Medium.nXi > 0));
  parameter SI.Acceleration g = Modelica.Constants.g_n "Gravity constant";
  parameter Types.Init.Temp initOption=Modelica_Fluid.WorkInProgress.Types.Init.SteadyState 
    "Default initialization";
  parameter Types.Flow.Temp flowOption=Modelica_Fluid.WorkInProgress.Types.Flow.Bidirectional 
    "Default flow type";
  Medium.BaseProperties medium(p(start=p), T(start=T), h(start=h),
                               Xi(start=X[1:Medium.nXi])) "Ambient medium";
  annotation (
    preferedView="info",
    defaultComponentName="environment",
    defaultComponentPrefixes="inner",
    missingInnerMessage="An \"environment\" component is not defined. A default 
environment component with initOptions = SteadyState will be used. If this is not desired, 
drag Modelica_Fluid.Components.Environment into the top level of your model.",
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
      Line(points=[40,84; 40,14], style(color=0, rgbcolor={0,0,0})),
      Polygon(points=[26,14; 54,14; 40,-18; 26,14], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Text(
        extent=[46,86; 94,42],
        style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1),
        string="g"),
      Text(
        extent=[-86,84; 20,38],
        style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1),
        string="p,T,X"),
      Text(
        extent=[-142,-112; 154,-142],
        string="p=%p",
        style(color=0, rgbcolor={0,0,0}))),
    Diagram,
    Documentation(info="<HTML>
<p>
This models defines ambient environment conditions and global
options for all components that are on the same or on a lower level
as this component. Dragging this component in a model results
in the following declaration:
</p>
<pre>
   <b>inner</b> Modelica_Fluid.WorkInProgress.Components.Environment environment;
</pre>
<p>
The parameters and variables of this component can be 
then accessed via a corresponding outer declaration:
</p>
<pre>
   <b>outer</b> Modelica_Fluid.WorkInProgress.Components.Environment environment;
</pre>
<p>
Note, the parameters \"InitOption\" and \"FlowOption\" are used as 
default setting by the Modelica_Fluid components but can
be individually redefined in every Modelica_Fluid component.
</p>
</HTML>
"));
equation 
  medium.p = p;
  if use_T then
    medium.T = T;
  else
    medium.h = h;
  end if;
  medium.Xi = X[1:Medium.nXi];
end Environment;
