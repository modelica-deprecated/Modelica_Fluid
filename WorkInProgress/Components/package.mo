package Components 
model GlobalOptions "Global options" 
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
end GlobalOptions;


model Environment "Environment conditions and global options (should be improved)" 
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

model PressureLoss "Generic pressure loss component" 
  extends Modelica_Fluid.WorkInProgress.Interfaces.PressureLossWithoutIcon;
  annotation (
    defaultComponentName="orifice",
    Diagram,
    Icon(
      Text(
        extent=[-120, 130; 116, 64],
        string="%name",
        style(gradient=2, fillColor=69)),
      Line(points=[-60, -50; -60, 50; 60, -50; 60, 50], style(color=0,
            thickness=2)),
      Line(points=[-60, 0; -100, 0], style(color=69)),
      Line(points=[60, 0; 100, 0], style(color=69))),
    Documentation(info="<html>
<p>
Generic pressure loss component. The loss factors are defined
with record \"PressureLossFactors\".
</p>
</html>"));
end PressureLoss;

model WallFriction 
  "Pressure loss due to friction in a straight pipe with walls of nonuniform roughness (commercial pipes)" 
  import SI = Modelica.SIunits;
  extends Modelica_Fluid.WorkInProgress.Interfaces.PressureLossWithoutIcon(
     final lossFactors = Modelica_Fluid.WorkInProgress.Utilities.PressureLossFactors.wallFriction(length, diameter, roughness));
  parameter SI.Length length "Length of pipe";
  parameter SI.Diameter diameter "Inner diameter of pipe";
  parameter SI.Length roughness(min=1e-10) 
    "Absolute roughness of pipe (> 0 required, details see info layer)";
  annotation (defaultComponentName="pipe",
    Documentation(info="<html>
<p>
The absolute roughness <font face=\"Symbol\">d</font> has usually to
be estimated. In <i>[Idelchik 1994, pp. 105-109,
Table 2-5; Miller 1990, p. 190, Table 8-1]</i> many examples are given.
As a short summary:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><b>Smooth pipes</b></td>
      <td>Drawn brass, coper, aluminium, glass, etc.</td>
      <td><font face=\"Symbol\">d</font> = 0.0025 mm</td>
  </tr>
  <tr><td rowspan=\"3\"><b>Steel pipes</b></td>
      <td>New smooth pipes</td>
      <td><font face=\"Symbol\">d</font> = 0.025 mm</td>
  </tr>
  <tr><td>Mortar lined, average finish</td>
      <td><font face=\"Symbol\">d</font> = 0.1 mm</td>
  </tr>
  <tr><td>Heavy rust</td>
      <td><font face=\"Symbol\">d</font> = 1 mm</td>
  </tr>
  <tr><td rowspan=\"3\"><b>Concrete pipes</b></td>
      <td>Steel forms, first class workmanship</td>
      <td><font face=\"Symbol\">d</font> = 0.025 mm</td>
  </tr>
  <tr><td>Steel forms, average workmanship</td>
      <td><font face=\"Symbol\">d</font> = 0.1 mm</td>
  </tr>
  <tr><td>Block linings</td>
      <td><font face=\"Symbol\">d</font> = 1 mm</td>
  </tr>
</table>
</html>"),
    Icon(
      Text(
        extent=[-150,140; 150,80],
        string="%name",
        style(gradient=2, fillColor=69)),
      Rectangle(extent=[-100,60; 100,-60],   style(
          color=0,
          gradient=2,
          fillColor=8)),
      Rectangle(extent=[-100,34; 100,-36],   style(
          color=69,
          gradient=2,
          fillColor=69))),
    Diagram(
      Rectangle(extent=[-100,50; 100,-50], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=7,
          rgbfillColor={255,255,255})),
      Line(points=[-50,-50; -50,50], style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1)),
      Text(
        extent=[-40,26; 28,14],
        style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1),
        string="diameter"),
      Line(points=[-100,60; 100,60], style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1)),
      Text(
        extent=[-30,74; 38,62],
        style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1),
        string="length")));
end WallFriction;

model SuddenExpansion "Pressure drop in pipe due to suddenly expanding area" 
  import SI = Modelica.SIunits;
  extends Modelica_Fluid.WorkInProgress.Interfaces.PressureLossWithoutIcon(
     final lossFactors = Modelica_Fluid.WorkInProgress.Utilities.PressureLossFactors.suddenExpansion(D_a, D_b));
  parameter SI.Diameter D_a "Inner diameter of pipe at port_a";
  parameter SI.Diameter D_b "Inner diameter of pipe at port_b";
  
  annotation (
    defaultComponentName="suddenExpansion",
    Diagram(
      Line(points=[0,40; -100,40; -100,-40; 0,-40; 0,-100; 100,-100; 100,100; 0,
            100; 0,40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1)),
      Rectangle(extent=[-100,40; 0,-40], style(
          color=7,
          rgbcolor={255,255,255},
          fillColor=7,
          rgbfillColor={255,255,255})),
      Rectangle(extent=[0,100; 100,-100], style(
          color=7,
          rgbcolor={255,255,255},
          fillColor=7,
          rgbfillColor={255,255,255})),
      Line(points=[0,40; -100,40; -100,-40; 0,-40; 0,-100; 100,-100; 100,100; 0,
            100; 0,40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1)),
      Line(points=[-60,-40; -60,40], style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1)),
      Text(
        extent=[-50,16; -26,-10],
        style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1),
        string="D_a"),
      Line(points=[34,-100; 34,100], style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1)),
      Text(
        extent=[54,16; 78,-10],
        style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1),
        string="D_b")),
    Icon(
      Text(
        extent=[-116,154; 120,88],
        string="%name",
        style(gradient=2, fillColor=69)),
      Rectangle(extent=[-100,80; 100,-80],   style(
          color=0,
          gradient=2,
          fillColor=8)),
      Rectangle(extent=[-100,20; 0,-20],     style(
          color=69,
          gradient=2,
          fillColor=69)),
      Rectangle(extent=[0,60; 100,-60],      style(
          color=69,
          gradient=2,
          fillColor=69))));
end SuddenExpansion;

model SharpEdgedOrifice 
  "Pressure loss due to sharp edged orifice (for both flow directions)" 
  import SI = Modelica.SIunits;
  import NonSI = Modelica.SIunits.Conversions.NonSIunits;
  extends Modelica_Fluid.WorkInProgress.Interfaces.PressureLossWithoutIcon(
     final lossFactors = Modelica_Fluid.WorkInProgress.Utilities.PressureLossFactors.sharpEdgedOrifice(D_pipe, D_min, L, alpha));
  parameter SI.Diameter D_pipe 
    "Inner diameter of pipe (= same at port_a and port_b)";
  parameter SI.Diameter D_min "Smallest diameter of orifice";
  parameter SI.Diameter L "Length of orifice";
  parameter NonSI.Angle_deg alpha "Angle of orifice";
  annotation (defaultComponentName="orifice",
    Documentation(info="<html>
</html>"),
    Icon(
      Text(
        extent=[-150,140; 150,80],
        string="%name",
        style(gradient=2, fillColor=69)),
      Rectangle(extent=[-100,80; 100,-80],   style(
          color=0,
          gradient=2,
          fillColor=8)),
      Rectangle(extent=[-100,60; 100,-60],   style(
          color=69,
          gradient=2,
          fillColor=69)),
      Polygon(points=[-24,60; -24,12; 36,50; 36,60; -24,60], style(
          color=0,
          rgbcolor={0,0,0},
          gradient=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=8)),
      Polygon(points=[-22,-10; -22,-60; 38,-60; 38,-50; -22,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=8))),
    Diagram(       Rectangle(extent=[-100,60; 100,-60], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=7,
          rgbfillColor={255,255,255})),
      Polygon(points=[-30,60; -30,12; 30,50; 30,60; -30,60], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=8)),
      Polygon(points=[-30,-10; -30,-60; 30,-60; 30,-50; -30,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=8)),
      Line(points=[-82,-60; -82,60], style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1)),
      Text(
        extent=[-78,16; -44,-8],
        style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1),
        string="D_pipe"),
      Line(points=[-30,-10; -30,12], style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1)),
      Text(
        extent=[-24,14; 8,-10],
        style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1),
        string="D_min"),
      Text(
        extent=[-20,84; 18,70],
        style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1),
        string="L"),
      Line(points=[30,68; -30,68],   style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=1)),
      Line(points=[16,40; 32,18; 36,-2; 34,-20; 20,-42], style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=8)),
      Text(
        extent=[38,8; 92,-6],
        style(
          color=3,
          rgbcolor={0,0,255},
          arrow=3,
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=8),
        string="alpha")));
end SharpEdgedOrifice;

model IsolatedPipe 
  "Model of an isolated pipe consisting of n pipe segments/FiniteVolumes" 
  import SI = Modelica.SIunits;
  
  replaceable package Medium = PackageMedium extends 
    Modelica_Media.Interfaces.PartialMedium "Medium in the component"  annotation (
      choicesAllMatching =                                       true);
  
  extends Modelica_Fluid.Interfaces.PartialMenuInitialization;
  
  parameter Integer nVolumes(min=1)=1 "Number of pipe segments/finite volumes";
  
  parameter SI.Length L "Length of pipe";
  parameter SI.AbsolutePressure dp_nominal(min=1.e-10) = 1 
    "|frictionType = ConstantLaminar or ConstantTurbulent| Nominal pressure drop";
  
  parameter SI.MassFlowRate m_flow_nominal = 1E-3 
    "Nominal mass flow rate at nominal pressure drop";
  
  parameter SI.Area A_a;
  parameter SI.Area A_b=A_a;
  
  parameter SI.Length Z_a=0;
  parameter SI.Length Z_b=Z_a;
  
  parameter Boolean dynamicMomentumBalance=false 
    "If false, der(m_flow) is neglected in momentum balance" 
                                                 annotation(Evaluate=true,
      Dialog(tab="Level of Detail"));
  parameter Boolean includeKineticTerm=false 
    "If false, d*v^2 is neglected in momentum balance" 
                                             annotation(Evaluate=true,
      Dialog(tab="Level of Detail"));
  parameter Boolean includeViscosity=false 
    "If false, artifical viscosity is neglected" 
                                          annotation(Evaluate=true, Dialog(tab=
          "Level of Detail"));
  parameter Real viscosityFactor1=0 annotation(Dialog(enable=includeViscosity,tab="Level of Detail"));
  parameter Real viscosityFactor2=1 annotation(Dialog(enable=includeViscosity,tab="Level of Detail"));
  
  Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare model Medium = Medium) 
              annotation (extent=[-120, -10; -100, 10]);
  Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare model Medium = Medium) 
              annotation (extent=[120, -10; 100, 10]);
  Modelica_Fluid.Utilities.PipeSegment pipeSegment[nVolumes](
      redeclare package Medium = Medium,
      each initType = initType,
      each init_p = init_p,
      each p_start = p_start,
      each d_start = d_start,
      each init_T = init_T,
      each T_start = T_start,
      each h_start = h_start,
      each X_start = X_start,
      each L=L/nVolumes,
      each dp_nominal=dp_nominal/nVolumes,
      each A_a=A_a "has to be corrected: linear distribution of A",
      each Z_a=Z_a "has to be corrected: linear distribution of Z",
      each m_flow_nominal=m_flow_nominal,
      each dynamicMomentumBalance=dynamicMomentumBalance,
      each includeKineticTerm=includeKineticTerm,
      each includeViscosity=includeViscosity,
      each viscosityFactor1=viscosityFactor1,
      each viscosityFactor2=viscosityFactor2);
  
annotation (Icon(
    Rectangle(extent=[-100, 60; 100, -60], style(color=0, fillColor=8)),
    Rectangle(extent=[-100, 34; 100, -36], style(
        color=69,
        gradient=2,
        fillColor=69)),
    Text(
      extent=[-150, 125; 150, 65],
      string="%name",
      style(gradient=2, fillColor=69)),
      Ellipse(extent=[-58,14; -28,-14], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
      Ellipse(extent=[22,14; 52,-14], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0}))));
equation 
  connect(port_a, pipeSegment[1].port_a);
  connect(port_b, pipeSegment[nVolumes].port_b);
  for i in 1:nVolumes - 1 loop
    connect(pipeSegment[i].port_b, pipeSegment[i + 1].port_a);
  end for;
end IsolatedPipe;

end Components;
