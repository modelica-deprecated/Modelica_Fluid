package Components 
model GlobalOptions "Global options" 
    import SI = Modelica.SIunits;
    import Modelica_Fluid.WorkInProgress.Types;
    
  parameter SI.Acceleration g = Modelica.Constants.g_n "Gravity constant";
  parameter Types.Init.Temp initOption=Modelica_Fluid.WorkInProgress.Types.Init.SteadyState 
      "Default initialization";
  parameter Types.Flow.Temp flowOption=Modelica_Fluid.WorkInProgress.Types.Flow.Bidirectional 
      "Default flow type";
  annotation (
    preferedView="info",
    defaultComponentName="globalOptions",
    defaultComponentPrefixes="inner",
    missingInnerMessage="An \"globalOptions\" component is not defined. A default 
globalOptions component with initOptions = SteadyState will be used. If this is not desired, 
drag Modelica_Fluid.Components.GlobalOptions into the top level of your model.",
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
          extent=[-98,76; 30,46],
          string="options",
          style(color=0, rgbcolor={0,0,0}))),
    Diagram,
    Documentation(info="<HTML>
<p>
This models defines global
options for all components that are on the same or on a lower level
as this component. Dragging this component in a model results
in the following declaration:
</p>
<pre>
   <b>inner</b> Modelica_Fluid.WorkInProgress.Components.GlobalOptions globalOptions;
</pre>
<p>
The parameters and variables of this component can be 
then accessed via a corresponding outer declaration:
</p>
<pre>
   <b>outer</b> Modelica_Fluid.WorkInProgress.Components.GlobalOptions globalOptions;
</pre>
<p>
Note, the parameters \"InitOption\" and \"FlowOption\" are used as 
default setting by the Modelica_Fluid components but can
be individually redefined in every Modelica_Fluid component.
</p>
</HTML>
"));
    
end GlobalOptions;
  
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
This model computes the pressure loss of a pipe
segment (orifice, bending etc.) with a minimum amount of data
provided via parameter <b>lossFactors</b>.
If available, data should be provided for <b>both flow directions</b>,
i.e., flow from port_a to port_b and from port_b to port_a, 
as well as for the <b>laminar</b> and the <b>turbulent</b> region.
It is also an option to provide the loss factor <b>only</b> for the
<b>turbulent</b> region for a flow from port_a to port_b.
</p>
<p>
The following equations are used:
</p>
<pre>   &Delta;p = 0.5*&zeta;*&rho;*v*|v|
      = 0.5*&zeta;/A^2 * (1/&rho;) * m_flow*|m_flow|
        Re = |v|*D*&rho;/&eta;
</pre>
<table border=1 cellspacing=0 cellpadding=2>
<tr><td><b>flow type</b></td>
    <td><b>&zeta;</b> = </td>
    <td><b>flow region</b></td></tr>
<tr><td>turbulent</td>
    <td><b>zeta1</b> = const.</td>
    <td>Re &ge;  Re_turbulent, v &ge; 0</td></tr>
<tr><td></td>
    <td><b>zeta2</b> = const.</td>
    <td>Re &ge; Re_turbulent, v &lt; 0</td></tr>
<tr><td>laminar</td>
    <td><b>c0</b>/Re</td>
    <td>both flow directions, Re small; c0 = const.</td></tr>
</table>
<p>
where
</p>
<ul>
<li> &Delta;p is the pressure drop: &Delta;p = port_a.p - port_b.p</li>
<li> v is the mean velocity.</li>
<li> &rho; is the density.</li>
<li> &zeta; is the loss factor that depends on the geometry of
     the pipe. In the turbulent flow regime, it is assumed that
     &zeta; is constant and is given by \"zeta1\" and
     \"zeta2\" depending on the flow direction.<br>
     When the Reynolds number Re is below \"Re_turbulent\", the
     flow is laminar for small flow velocities. For higher 
     velocities there is a transition region from 
     laminar to turbulent flow. The loss factor for
     laminar flow at small velocities is defined by the often occuring
     approximation c0/Re. If c0 is different for the two
     flow directions, the mean value has to be used 
     (c0 = (c0_ab + c0_ba)/2).<li>
<li> The equation \"&Delta;p = 0.5*&zeta;*&rho;*v*|v|\" is either with
     respect to port_a or to port_b, depending on the definition
     of the particular loss factor &zeta; (in some references loss
     factors are defined with respect to port_a, in other references
     with respect to port_b).</li>
 
<li> Re = |v|*D_Re*&rho;/&eta; = |m_flow|*D_Re/(A_Re*&eta;) 
     is the Reynolds number at the smallest cross
     section area. This is often at port_a or at port_b, but can
     also be between the two ports. In the record, the diameter
     D_Re of this smallest cross section area has to be provided, as
     well, as Re_turbulent, the absolute value of the 
     Reynolds number at which
     the turbulent flow starts. If Re_turbulent is different for
     the two flow directions, use the smaller value as Re_turbulent.</li>
<li> D is the diameter of the pipe. If the pipe has not a 
     circular cross section, D = 4*A/P, where A is the cross section
     area and P is the wetted perimeter.</li>
<li> A is the cross section area with A = &pi;(D/2)^2.
<li> &eta; is the dynamic viscosity.</li>
</ul>
<p>
The laminar and the transition region is usually of
not much technical interest because the operating point is
mostly in the turbulent regime. For simplification and for
numercial reasons, this whole region is described by two
polynomials of third order, one polynomial for m_flow &ge; 0 
and one for m_flow &lt; 0. The polynomials start at 
Re = |m_flow|*4/(&pi;*D_Re*&eta;), where D_Re is the
smallest diameter between port_a and port_b.
The common derivative
of the two polynomials at Re = 0 is
computed from the equation \"c0/Re\". Note, the pressure drop
equation above in the laminar region is always defined
with respect to the smallest diameter D_Re.
</p>
<p>
If no data for c0 is available, the derivative at Re = 0 is computed in such
a way, that the second derivatives of the two polynomials
are identical at Re = 0. The polynomials are constructed, such that
they smoothly touch the characteristic curves in the turbulent
regions. The whole characteristic is therefore <b>continuous</b>
and has a <b>finite</b>, <b>continuous first derivative everywhere</b>.
In some cases, the constructed polynomials would \"vibrate\". This is 
avoided by reducing the derivative at Re=0 in such a way that
the polynomials are guaranteed to be monotonically increasing.
The used sufficient criteria for monotonicity follows from:
</p>

<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <b>Monotone piecewise cubic interpolation</b>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
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
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching = true);
    
  extends Modelica_Fluid.Interfaces.PartialInitializationParameters;
    
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
  Modelica_Fluid.WorkInProgress.Utilities.PipeSegment pipeSegment[nVolumes](
      redeclare package Medium = Medium,
      each initOption = initOption,
      each p_start = p_start,
      each use_T_start = use_T_start,
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
  
model Pipe 
    
  extends Modelica_Fluid.WorkInProgress.Interfaces.Flow1D(
    Qs_flow=heat.Q_flow,
    ms_flow=zeros(n),
    msXi_flow=zeros(n, Medium.nXi),
    dp=friction.dp);
  parameter SI.Area area_h = P_inner*length "Heat transfer area" annotation(Dialog(tab="Advanced", group="Heat transfer and friction loss"));
  parameter SI.Diameter d_h = 4*A_inner/P_inner "Hydraulic diameter" annotation(Dialog(tab="Advanced", group="Heat transfer and friction loss"));
  parameter Boolean use_wall=false 
      "= true, use wall component between fluid and thermalPort" 
                                                                annotation(Dialog(tab="General", group="Wall - optional"),Evaluate=true);
  parameter SI.Diameter d_outer "Outer diameter of circular pipe" annotation(Dialog(tab="General", group="Wall - optional", enable=(crossSectionType==1 and use_wall)));
  parameter SI.Length h_outer "Outer height of rectangular pipe"  annotation(Dialog(tab="General", group="Wall - optional", enable=(crossSectionType==2 and use_wall)));
  parameter SI.Length w_outer "Outer width of rectangular pipe"  annotation(Dialog(tab="General", group="Wall - optional", enable=(crossSectionType==2 and use_wall)));
  parameter SI.Length A_outer = if crossSectionType == 1 then Modelica.Constants.pi*d_outer*d_outer/4 else if crossSectionType == 2 then h_outer*w_outer else 1 
      "Outer cross section area" 
                               annotation(Dialog(tab="General", group="Wall - optional", enable=(use_wall and crossSectionType==3)));
  inner Medium.ThermodynamicState[n] state = medium.state;
    
  replaceable 
      Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PipeHT_constAlpha
      heat(
    redeclare final package Medium = Medium,
    final n=n,
    final d_h=d_h,
    final A_h=area_h,
    T=medium.T) extends 
      Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PartialPipeHeatTransfer(
    redeclare final package Medium = Medium,
    final n=n,
    final d_h=d_h,
    final A_h=area_h,
    T=medium.T) "Convective heat transfer" 
                annotation (Dialog(tab="Advanced", group="Heat transfer and friction loss"),choicesAllMatching, extent=[-20,-20;
        20,20]);
  replaceable 
      Modelica_Fluid.WorkInProgress.Utilities.PipeFriction.PipeFriction_RoughSurface
      friction(
    redeclare final package Medium = Medium,
    final d_h = d_h,
    final length = length,
    final n=n,
    final np=np,
    final lumped_dp=lumped_dp,
    m_flow=if lumped_dp then ones(np)*m_flow[2] else m_flow,
    d_1 = if lumped_dp then ones(np)*d_port_a else {if i>1 then medium[i-1].d else d_port_a for i in 1:np},
    d_2 = if lumped_dp then ones(np)*d_port_b else {if i<np then medium[i].d else d_port_b for i in 1:np}) extends 
      Modelica_Fluid.WorkInProgress.Utilities.PipeFriction.PartialPipeFriction(
    redeclare final package Medium = Medium,
    final d_h = d_h,
    final length = length,
    final n=n,
    final np=np,
    final lumped_dp=lumped_dp,
    m_flow=if lumped_dp then ones(np)*m_flow[2] else m_flow,
    d_1 = if lumped_dp then ones(np)*d_port_a else {if i>1 then medium[i-1].d else d_port_a for i in 1:np},
    d_2 = if lumped_dp then ones(np)*d_port_b else {if i<np then medium[i].d else d_port_b for i in 1:np}) 
      "Pressure drop due to friction" 
                                    annotation(Dialog(tab="Advanced", group="Heat transfer and friction loss"), choicesAllMatching, extent=[-58,-20;
        -18,20]);
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
      "Thermal port" 
    annotation (extent=[-20,60; 20,80]);
  replaceable model Wall = 
        Modelica_Fluid.WorkInProgress.Components.Wall_constProps                    extends 
      Modelica_Fluid.WorkInProgress.Interfaces.PartialPipeWall "Wall model" 
                                                                   annotation(choicesAllMatching, Dialog(enable=use_wall, tab="General", group="Wall - optional"));
  Wall wall(final n=n, final a_inner=A_inner, final a_outer=A_outer, final 
        length=length, initOption=initOption) if use_wall 
                           annotation (extent=[10,20; 50,60]);
  annotation (Icon(Rectangle(extent=[-100,60; 100,40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)), Rectangle(extent=[-100,-40; 100,-60], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)),
      Text(
        extent=[-100,-60; 100,-100],
        string="%name",
        style(color=3, rgbcolor={0,0,255}))),
                            Diagram, 
      Documentation(info="<html>
<p>
From Katrins email, Nov. 28, 2005:
</p>

<p>
extends Interfaces.1DFlow. Pressure drop and heat transfer are added in terms of replaceable components. The main problem here is to make all required variables and parameters available to the respective component (medium state, pipe geometry, Medium functions, empirical parameters). Only those shared by all future replaceable models (the simple one parameter model and the highly sophisticated (fictitious) two phase Nusselt correlation) can be set by modifiers (which is not straightforward in Dymola at the moment if a contsraining clause is used).  Those not required by all models as i.e. viscosity and conductivitiy must be computed inside the component from medium properties made available via inner and outer. I always try to avoid this as it it as bit like free climbing, but in this case I see no better solution.
</p>

<p>
Martin, I have not tested your latest pressure drop implementation with this model, but will do so as soon as possible. However, it is used in a completely different way, that means as an array of components, not as a  base class, in order to be able to handle distributed flow. I will check if another implementation would be more practical.
</p>

<p>
The pipe model contains a Boolean flag useWall which determines if a wall component is added. Unfortunately the icon does not represent the difference. In this way a heat exchanger can be created using two instances of the pipe model, one with a wall and one without. If interested in transients it could also make sense to include a wall in an insulated pipe. 
</p>

</html>"));
equation 
if use_wall then
  connect(wall.thermalPort_a, thermalPort) annotation (points=[30,50; 30,60; 0,
          60; 0,70],
      style(
      color=10,
      rgbcolor={95,95,95},
      fillColor=0,
      rgbfillColor={0,0,0},
      fillPattern=7));
  connect(wall.thermalPort_b, heat.thermalPort) 
                                              annotation (points=[30,30; 30,22;
          0,22; 0,14],
      style(
      color=10,
      rgbcolor={95,95,95},
      fillColor=0,
      rgbfillColor={0,0,0},
      fillPattern=7));
else
  connect(thermalPort, heat.thermalPort) 
                                       annotation (points=[0,70; 0,14],
      style(
      color=10,
      rgbcolor={95,95,95},
      fillColor=0,
      rgbfillColor={0,0,0},
      fillPattern=7));
end if;
end Pipe;
  
model HeatExchanger "Double pipe heat exchanger with neglectible outer wall" 
  //General
  parameter Integer n(min=1) "Spatial segmentation";
  replaceable package Medium_1 = Modelica.Media.Water.StandardWater extends 
      Modelica.Media.Interfaces.PartialMedium "Inner pipe medium" 
                                                    annotation(choicesAllMatching);
  replaceable package Medium_2 = Modelica.Media.Water.StandardWater extends 
      Modelica.Media.Interfaces.PartialMedium "Outer pipe medium" 
                                                    annotation(choicesAllMatching);
  parameter SI.Length di_1(min=0) "Inner diameter of inner pipe"     annotation(Dialog(tab="General", group="Dimensions"));
  parameter SI.Length da_1(min=0) "Outer diameter of inner pipe"     annotation(Dialog(tab="General", group="Dimensions"));
  parameter SI.Length di_2(min=0) "Inner diameter of outer pipe"     annotation(Dialog(tab="General", group="Dimensions"));
  parameter SI.Length length(min=0) "Length of both pipes" annotation(Dialog(tab="General", group="Dimensions"));
  //Wall
  parameter SI.Density d_wall "Density of wall material" annotation(Dialog(tab="General", group="Constant material properties"));
  parameter SI.SpecificHeatCapacity c_wall 
      "Specific heat capacity of wall material"                                      annotation(Dialog(tab="General", group="Constant material properties"));
  parameter SI.Temperature T_start_w "Start value of wall termperature" annotation(Dialog(tab="Initialization", group="Wall"));
  final parameter SI.Mass m_wall=sum(pipe_1.wall.m) "Wall mass";
  //Initialization pipe 1
  parameter Modelica_Fluid.Types.InitTypes.Temp initOption_1 
      "Initialization option" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe"));
  parameter Boolean use_T_start_1=true "Use T_start if true, otherwise h_start"
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe"));
  parameter Medium_1.AbsolutePressure p_start_1=Medium_1.reference_p 
      "Start value of pressure" 
    annotation(Dialog(tab = "Initialization", group = "Inner pipe"));
  parameter Medium_1.Temperature T_start_1=if use_T_start_1 then 293.15 else 
      Medium_1.T_phX(p_start_1, h_start_1, X_start_1) 
      "Start value of temperature" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe", enable = use_T_start_1));
  parameter Medium_1.SpecificEnthalpy h_start_1=if use_T_start_1 then 
      Medium_1.h_pTX(p_start_1, T_start_1, X_start_1[1:Medium_1.nXi]) else 1e4 
      "Start value of specific enthalpy" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe", enable = not use_T_start_1));
  parameter Medium_1.MassFraction X_start_1[Medium_1.nX]=Medium_1.reference_X 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Initialization", group = "Inner pipe", enable=(Medium_1.nXi > 0)));
  parameter Medium_1.MassFlowRate mflow_start_1 "Start value of mass flow rate"
                                    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe"));
  //Initialization pipe 2
  parameter Modelica_Fluid.Types.InitTypes.Temp initOption_2 
      "Initialization option" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe"));
  parameter Boolean use_T_start_2=true "Use T_start if true, otherwise h_start"
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe"));
  parameter Medium_2.AbsolutePressure p_start_2=Medium_2.reference_p 
      "Start value of pressure" 
    annotation(Dialog(tab = "Initialization", group = "Outer pipe"));
  parameter Medium_2.Temperature T_start_2=if use_T_start_2 then 293.15 else 
      Medium_2.T_phX(p_start_2, h_start_2, X_start_2) 
      "Start value of temperature" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe", enable = use_T_start_2));
  parameter Medium_2.SpecificEnthalpy h_start_2=if use_T_start_2 then 
      Medium_2.h_pTX(p_start_2, T_start_2, X_start_2[1:Medium_2.nXi]) else 1e4 
      "Start value of specific enthalpy" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe", enable = not use_T_start_2));
  parameter Medium_2.MassFraction X_start_2[Medium_2.nX]=Medium_2.reference_X 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Initialization", group = "Outer pipe", enable=Medium_2.nXi>0));
  parameter Medium_2.MassFlowRate mflow_start_2 "Start value of mass flow rate"
                                       annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe"));
  //Advanced
  parameter Boolean lumped_dp 
      " = true, lumped pressure drop, reduces number of pressure states to one"
                              annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
  parameter Boolean static "= true, use quasistatic mass and energy balances" 
                           annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
  parameter Boolean kineticTerm 
      " = true, include kinetic term in momentum balance" 
                                annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
  parameter Real K1=1 "Enhancement factor for outer pipe heat transfer area" 
                                                                          annotation(Dialog(tab="Advanced", group="Heat transfer and friction loss"));
  replaceable model PipeFriction = 
      Modelica_Fluid.WorkInProgress.Utilities.PipeFriction.PipeFriction_SimpleLinear
                                                                                        extends 
      Modelica_Fluid.WorkInProgress.Utilities.PipeFriction.PartialPipeFriction 
      "Friction model"                                                                                  annotation(choicesAllMatching, Dialog(tab="Advanced", group="Heat transfer and friction loss"));
  replaceable model HeatTransfer = 
      Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PipeHT_constAlpha
      extends 
      Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PartialPipeHeatTransfer
      "Heat transfer model"                                                                       annotation(choicesAllMatching, Dialog(tab="Advanced", group="Heat transfer and friction loss"));
    
  //Display variables
  SI.HeatFlowRate Q_flow_1 "Total heat flow rate of inner pipe";
  SI.HeatFlowRate Q_flow_2 "Total heat flow rate of outer pipe";
    
  Modelica_Fluid.WorkInProgress.Components.Pipe pipe_1(
    redeclare package Medium = Medium_1,
    n=n,
    H_ab=0,
    static=static,
    lumped_dp=lumped_dp,
    kineticTerm=kineticTerm,
    gravityTerm=false,
    dynamicTerm=false,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.Circular,
    length=length,
    use_wall=true,
    redeclare model Wall = 
        Modelica_Fluid.WorkInProgress.Components.Wall_constProps (
        d_wall=d_wall,
        c_wall=c_wall,
        T_start=T_start_w),
    redeclare HeatTransfer heat,
    redeclare PipeFriction friction,
    initOption=initOption_1,
    use_T_start=use_T_start_1,
    p_start=p_start_1,
    T_start=T_start_1,
    h_start=h_start_1,
    X_start=X_start_1,
    mflow_start=mflow_start_1,
    d_inner=di_1,
    d_outer=da_1)            annotation (extent=[-40,-60; 20,0]);
    
  Modelica_Fluid.WorkInProgress.Components.Pipe pipe_2(
    redeclare package Medium = Medium_2,
    n=n,
    H_ab=0,
    static=static,
    lumped_dp=lumped_dp,
    kineticTerm=kineticTerm,
    gravityTerm=false,
    dynamicTerm=false,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.General,
    length=length,
    redeclare HeatTransfer heat,
    redeclare PipeFriction friction,
    use_T_start=use_T_start_2,
    p_start=p_start_2,
    T_start=T_start_2,
    h_start=h_start_2,
    X_start=X_start_2,
    initOption=initOption_2,
    mflow_start=mflow_start_2,
    P_inner=Modelica.Constants.pi*(da_1 + di_2),
    A_inner=Modelica.Constants.pi/4*(di_2*di_2 - da_1*da_1),
    area_h=Modelica.Constants.pi*da_1*length*K1) 
              annotation (extent=[-40,88; 20,28]);
  annotation (Diagram(Line(points=[-10,36; -10,-8], style(
          color=1,
          rgbcolor={255,0,0},
          thickness=2))), Icon(
      Rectangle(extent=[-100,-26; 100,-30], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)),
      Rectangle(extent=[-100,30; 100,26], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)),
      Rectangle(extent=[-100,60; 100,30], style(
          color=0,
          rgbcolor={0,0,0},
          gradient=2,
          fillColor=70,
          rgbfillColor={0,63,125})),
      Rectangle(extent=[-100,-30; 100,-60], style(
          color=0,
          rgbcolor={0,0,0},
          gradient=2,
          fillColor=70,
          rgbfillColor={0,63,125})),
      Rectangle(extent=[-100,26; 100,-26], style(
          color=69,
          rgbcolor={0,128,255},
          gradient=2,
          fillColor=69,
          rgbfillColor={0,128,255})),
      Text(
        extent=[-100,-60; 100,-100],
        string="%name",
        style(color=3, rgbcolor={0,0,255}))));
  Modelica_Fluid.Interfaces.FluidPort_b port_b1(redeclare package Medium = 
        Medium_1) annotation (extent=[100,-12; 120,8]);
  Modelica_Fluid.Interfaces.FluidPort_a port_a1(redeclare package Medium = 
        Medium_1) annotation (extent=[-120,-12; -100,8]);
  Modelica_Fluid.Interfaces.FluidPort_a port_a2(redeclare package Medium = 
        Medium_2) annotation (extent=[-120,36; -100,56]);
  Modelica_Fluid.Interfaces.FluidPort_b port_b2(redeclare package Medium = 
        Medium_2) annotation (extent=[100,-56; 120,-36]);
equation 
  Q_flow_1=sum(pipe_1.Qs_flow);
  Q_flow_2=sum(pipe_2.Qs_flow);
  connect(pipe_2.thermalPort, pipe_1.thermalPort);
  connect(pipe_2.port_b, port_b2) annotation (points=[20,58; 60,58; 60,-46; 110,
        -46], style(
      color=69,
      rgbcolor={0,127,255},
      thickness=2,
      gradient=2,
      fillColor=42,
      rgbfillColor={213,0,0}));
  connect(pipe_1.port_b, port_b1) annotation (points=[20,-30; 42,-30; 42,-2;
        110,-2], style(
      color=69,
      rgbcolor={0,127,255},
      thickness=2,
      gradient=2,
      fillColor=42,
      rgbfillColor={213,0,0}));
  connect(pipe_1.port_a, port_a1) annotation (points=[-40.6,-30; -75.3,-30;
          -75.3,-2; -110,-2],
                            style(
      color=69,
      rgbcolor={0,127,255},
      thickness=2,
      gradient=2,
      fillColor=42,
      rgbfillColor={213,0,0}));
  connect(pipe_2.port_a, port_a2) annotation (points=[-40.6,58; -76,58; -76,46;
          -110,46],
                  style(
      color=69,
      rgbcolor={0,127,255},
      thickness=2,
      gradient=2,
      fillColor=42,
      rgbfillColor={213,0,0}));
end HeatExchanger;
  
model Wall_constProps 
    "Pipe wall, assuming ideal 1D-conduction and constant material properties" 
  extends Modelica_Fluid.WorkInProgress.Interfaces.PartialPipeWall;
  parameter SI.Density d_wall "Density of wall material";
  parameter SI.SpecificHeatCapacity c_wall 
      "Specific heat capacity of wall material";
  parameter SI.Temperature T_start "Start value for wall temperature";
  parameter SI.Mass[n] m=ones(n)*(a_outer-a_inner)*length*d_wall/n "Wall mass";
  parameter Modelica_Fluid.Types.InitTypes.Temp initOption;
  SI.Temperature[n] T(start=ones(n)*T_start, stateSelect=StateSelect.prefer) 
      "Wall temperature";
initial equation 
  if initOption==2 then //2: full steady state initialization
    der(T)=zeros(n);
  else
    T=ones(n)*T_start;
  end if;
equation 
  for i in 1:n loop
   c_wall*m[i]*der(T[i]) = thermalPort_a[i].Q_flow + thermalPort_b[i].Q_flow;
  end for;
  //assuming ideal heat conduction perpendicular to fluid flow, conduction in remaining two dimensions is neglected
  thermalPort_a.T=T;
  thermalPort_b.T=T;
end Wall_constProps;
  
end Components;
