package Utilities 
model TankAttachment "Equations to attach pipe at tank" 
    import SI = Modelica.SIunits;
     replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
     annotation (choicesAllMatching=true);
    
    Modelica.Fluid.Interfaces.FluidPort_a port(redeclare package Medium = Medium) 
    annotation (extent=[-10,-112; 10,-92],    rotation=90);
   // Real mXi_flow;
    parameter Boolean onlyInFlow = false 
      "= true, if flow only into the tank (e.g. top ports)" 
                                                          annotation(Evaluate=true);
    parameter SI.Diameter pipeDiameter=0.0 "Inner (hydraulic) pipe diameter" 
                                                                         annotation(Dialog(enable=not onlyInFlow));
    parameter SI.Height pipeHeight "Height of pipe";
    parameter Medium.SpecificEnthalpy h_start 
      "Start value of specific enthalpy (used as enthalpy for topPorts if back flow)";
    parameter Medium.MassFraction X_start[Medium.nX] 
      "Start value of mass fractions m_i/m";
    parameter SI.Height level_start "Start value of tank level";
    
    parameter SI.AbsolutePressure p_ambient "Tank surface pressure" 
                                                                annotation(Dialog);
    input SI.Height level "Actual tank level" annotation(Dialog);
    input Medium.SpecificEnthalpy h "Actual specific enthalpy of fluid in tank"
                                                                annotation(Dialog);
    input Medium.Density d "Actual specific density of fluid in tank" 
                                                      annotation(Dialog);
    input Medium.MassFraction Xi[Medium.nXi] 
      "Actual mass fractions of fluid in tank"                  annotation(Dialog);
    parameter Real k_small(min=0) = 1e-5 
      "Small regularization range if tank level is below bottom_height or side_height; k_small = 0 gives ideal switch"
              annotation(Evaluate=true);
    parameter Real s_start = 0;
    
    output Medium.EnthalpyFlowRate H_flow 
      "= port.H_flow (used to transform vector of connectors in vector of Real numbers)";
    output Medium.MassFlowRate m_flow 
      "= port.m_flow (used to transform vector of connectors in vector of Real numbers)";
    output Medium.MassFlowRate mXi_flow[Medium.nXi] 
      "= port.mXi_flow (used to transform vector of connectors in vector of Real numbers)";
    
  annotation (Documentation(info="<html>
<p>
This component contains the equations that attach the pipe
to the tank. The main reason to introduce this component is
that Dymola has currently limitations for connector arrays
when the dimension is zero. Without this utility component
it would not be possible to set, e.g., n_topPorts to zero.
</p>
 
 
</html>"),
         Icon(Rectangle(extent=[-100,0; 100,-100], style(
          color=3,
          rgbcolor={0,0,255},
          fillColor=7,
          rgbfillColor={255,255,255})), Text(
        extent=[-122,48; 132,6],
        style(
          color=3,
          rgbcolor={0,0,255},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1),
        string="%name")));
  protected 
  outer Modelica_Fluid.Ambient ambient "Global default options";
  parameter SI.Area pipeArea = Modelica.Constants.pi*(pipeDiameter/2)^2;
  parameter Medium.MassFlowRate m_flow_nominal = 1 
      "Nominal mass flow rate used for scaling (has only an effect if k_small > 0)";
  parameter Medium.AbsolutePressure p_nominal = p_ambient;
  SI.Length aboveLevel = level - pipeHeight;
  Boolean m_flow_out(start=true,fixed=true) "true= massflow out of tank";
  Real s(start=s_start) 
      "path parameter of parameterized curve description (either m_flow/m_flow_nominal or (port.p-p_ambient)/p_ambient)";
equation 
  H_flow = port.H_flow;
  m_flow = port.m_flow;
  mXi_flow = port.mXi_flow;
    
  if onlyInFlow then
     m_flow_out = false "Dummy value in this if branch";
     port.p = p_ambient;
     /* flow should never out of the port. However, if this occurs in a 
        small time interval (e.g. during initialization), the start values of
        h and X are provided, since otherwise there is a singular
        system 
     */
     port.H_flow = semiLinear(port.m_flow, port.h, h_start);
     port.mXi_flow = semiLinear(port.m_flow, port.Xi, X_start[1:Medium.nXi]);
     assert(port.m_flow > -1e-6, "Mass flows out of tank via topPort. This indicates a wrong model");
     s = 0;
  else
     port.H_flow = semiLinear(port.m_flow, port.h, h);
     port.mXi_flow = semiLinear(port.m_flow, port.Xi, Xi);
      
/* Original equations from Poschlad/Remelhe:
*/
     s = 0;
     m_flow_out = (pre(m_flow_out) and not port.p>p_ambient) or (port.m_flow < -1e-6);
      
     if (aboveLevel > 0) then
       port.p = aboveLevel*ambient.g*d + p_ambient - smooth(2,noEvent(if m_flow < 0 then m_flow^2/(2*d*pipeArea^2) else 0));
     else
       if pre(m_flow_out) then
          m_flow = 0;
       else
          port.p = p_ambient;
       end if;
     end if;
      
/* Martin Otter: The following equations are a declarative form 
   (parameterized curve description) of the above equations and
   should theoretically work better. However, some examples with
   IF97 water fail, whereas the above works. Therefore, not used. 
       Add the following text to OpenTank, once the initialization
   with this solution works for Modelica_Fluid.Examples.Tanks.ThreeOpenTanks:
 
OpenTank:
<p>
The situation when the tank level is below bottom_heights[i] or side_heights[i]
is handeled properly. Details are described
<a href="Modelica:Modelica_Fluid.Utilities.TankAttachment">here</a>
</p> 
 
TankAttachment:
<p>
If a bottom or side connector is above the actual tank level, the
following characteristic is used to compute the mass flow rate port.m_flow
from the connector to the tank and the absolute pressure port.p
in the port:
</p>
 
<img src="../Images/Components/Tank_PipeAboveTankLevel.png">   
 
 
     m_flow_out = s <= 0;
     if aboveLevel >= 0 then
        m_flow = m_flow_nominal*s "equation to compute s, which is a dummy in this branch";
        port.p - p_ambient = aboveLevel*ambient.g*d  -
                             smooth(2,if m_flow_out then s*abs(s)*m_flow_nominal^2/(2*d*pipeArea^2) else k_small*m_flow_nominal*s);
     else
        m_flow = (if m_flow_out then k_small*p_nominal else m_flow_nominal)*s;
        port.p - p_ambient = (if m_flow_out then p_nominal else k_small*m_flow_nominal)*s;
     end if;
*/
      
  end if;
    
  /*
  More precise equations (introduce them later; need to transform
  them from energy balance form 1 to form 2):
 
  Momentum balance:
  (integrated momentum equation for frictionless fluid with density that is
   independent of the level, i.e., the unsteady Bernoulli equation for incompressible fluid)
  v_level = der(level);
  v = -port.m_flow/(rho*A_outlet);
  level*der(v_level) + (v^2 - v_level^2)/2 - g*level + (p - p_ambient)/rho = 0; or
  rho*level*der(v_level) + rho*(v^2 - v_level^2)/2 - rho*g*level + (p - p_ambient) = 0;
 
  Energy balance:
  Potential energy: E_pot = integ(dm*g*s)
                          = g*integ(rho*A*s*ds)
                          = g*rho*A*z^2/2
  Kinetic energy  : E_kin = integ(dm*v^2/2)
                          = integ(rho*A*v^2/2*ds)
                          = rho*A*v^2/2*integ(ds)
                          = rho*A*v^2/2*z
                          = M*v^2/2
  E = U + M*g*z/2 + M*v_level^2/2
  der(E) = port.H_flow + port.m_flow*v^2/2 - p_ambient*area*der(level)
*/
    
end TankAttachment;
  
 record PressureLossFactors 
    "Data structure defining constant loss factors zeta for dp = zeta*rho*v*|v|/2 and functions providing the data for some loss types" 
    import SI = Modelica.SIunits;
  extends Modelica.Icons.Record;
    
  parameter SI.Diameter D_a "Diameter at port_a";
  parameter SI.Diameter D_b = D_a "Diameter at port_b";
  parameter Real zeta1 "Loss factor for flow port_a -> port_b";
  parameter Real zeta2=zeta1 "Loss factor for flow port_b -> port_a";
  parameter SI.ReynoldsNumber Re_turbulent = 1000 
      "Loss factors suited for |Re| >= Re_turbulent";
  parameter SI.Diameter D_Re = min(D_a,D_b) "Diameter used to compute Re";
  parameter Boolean zeta1_at_a = true 
      "dp = zeta1*(if zeta1_at_a then d_a*v_a^2/2 else d_b*v_b^2/2)";
  parameter Boolean zeta2_at_a = false 
      "dp = -zeta2*(if zeta2_at_a then d_a*v_a^2/2 else d_b*v_b^2/2)";
  parameter Boolean zetaLaminarKnown = false 
      "= true, if zeta = c0/Re in laminar region";
  parameter Real c0 = 1 
      "zeta = c0/Re; dp = zeta*d_Re*v_Re^2/2, Re=v_Re*D_Re*d_Re/eta_Re)"                   annotation(Dialog(enable=zetaLaminarKnown));
    
  annotation (preferedView="info", Documentation(info="<html>
<p>
This record defines the pressure loss factors of a pipe
segment (orifice, bending etc.) with a minimum amount of data.
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
    
   encapsulated function wallFriction 
      "Compute pressure loss factors due to friction in a straight pipe with walls of nonuniform roughness (commercial pipes)" 
      import SI = Modelica.SIunits;
      import Modelica_Fluid.WorkInProgress.Utilities.PressureLossFactors;
      import lg = Modelica.Math.log10;
      
     input SI.Length length "Length of pipe";
     input SI.Diameter diameter "Inner diameter of pipe";
     input SI.Length roughness(min=1e-10) 
        "Absolute roughness of pipe (> 0 required, details see info layer)";
     output PressureLossFactors data 
        "Pressure loss factors for both flow directions";
     annotation (Icon(Rectangle(extent=[-100,48; 100,-50], style(
             color=0,
             rgbcolor={0,0,0},
             fillColor=7,
             rgbfillColor={255,255,255}))),
                               Diagram(
         Rectangle(extent=[-100,48; 100,-50], style(
               color=0,
               rgbcolor={0,0,0},
               fillColor=7,
               rgbfillColor={255,255,255},
               fillPattern=1)),
         Line(points=[-60,-50; -60,48], style(
             color=3,
             rgbcolor={0,0,255},
             arrow=3,
             fillColor=3,
             rgbfillColor={0,0,255},
             fillPattern=1)),
         Text(
           extent=[-50,16; 6,-10],
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
           extent=[-34,80; 34,62],
           style(
             color=3,
             rgbcolor={0,0,255},
             arrow=3,
             fillColor=3,
             rgbfillColor={0,0,255},
             fillPattern=1),
           string="length")),
       Documentation(info="<html>
<p>
Friction in straight pipe with walls of nonuniform roughness 
(commercial pipes)
</p>
<p>
The loss factors are given for mass flow rates from 
port_a to port_b as:
</p>
<pre>
  turbulent flow (Idelchik 1994, diagram 2-5, p. 117)
     zeta = (L/D)/(2*lg(3.7 / &Delta;))^2, for Re >= 560/&Delta;
     for Re >= 560/&Delta; the loss factor does not depend on the
     Reynolds number. For Re >= 4000, the flow is turbulent,
     but depends both on &Delta; and slightly on Re.
&nbsp;
  laminar flow (Idelchick 1994, diagram 2-1, p. 110):
     zeta = 64*(L/D)/Re
</pre>
<p>
where
</p>
<ul>
<li> D is the inner pipe diameter</li>
<li> L is the lenght of the pipe</li>
<li> &Delta; = &delta;/D is the relative roughness where &delta; is
     the absolute \"roughness\", i.e., the averaged height of asperities in the pipe.
     (&delta; may change over time due to growth of surface asperities during
      service, see [Idelchick 1994, p. 85, Tables 2-1, 2-2]).</li>
</ul>
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
</html>"));
    protected 
     Real Delta = roughness/diameter "relative roughness";
   algorithm 
     data.D_a          := diameter;
     data.D_b          := diameter;
     data.zeta1        := (length/diameter)/(2*lg(3.7 /Delta))^2;
     data.zeta2        := data.zeta1;
     data.Re_turbulent := 4000 
        ">= 560/Delta flow does not depend on Re, but interpolation is bad";
     data.D_Re         := diameter;
     data.zeta1_at_a   := true;
     data.zeta2_at_a   := false;
     data.zetaLaminarKnown := true;
     data.c0               := 64*(length/diameter);
   end wallFriction;
    
   encapsulated function suddenExpansion 
      "Compute pressure loss factors for sudden expansion or contraction in a pipe (for both flow directions)" 
      import SI = Modelica.SIunits;
      import Modelica_Fluid.WorkInProgress.Utilities.PressureLossFactors;
      
     input SI.Diameter D_a "Inner diameter of pipe at port_a";
     input SI.Diameter D_b "Inner diameter of pipe at port_b";
     output PressureLossFactors data 
        "Pressure loss factors for both flow directions";
     annotation (Icon(
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
             fillPattern=1))), Diagram(
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
       Documentation(info="<html>
<p>
The loss factors are given for mass flow rates from 
port_a to port_b as:
</p>
<pre>
   A_a &lt; A_b (Idelchik 1994, diagram 4-1, p. 208):
      zeta = dp/(d_a*v_a^2/2)
           = (1 - A_a/A_b)^2 for Re_a &ge; 3.3e3 (turbulent flow)
      zeta = 30/Re           for Re_a &lt; 10    (laminar flow)
&nbsp;
   A_a &gt; A_b (Idelchik 1994, diagram 4-9, p. 216 and diagram 4-10, p. 217)
      zeta = dp/(d_b*v_b^2/2)
           = 0.5*(1 - A_b/A_a)^0.75 for Re_b &ge; 1e4 (turbulent flow)
      zeta = 30/Re                  for Re_a &lt; 10  (laminar flow)
</pre>
</html>"));
    protected 
     Real A_rel;
   algorithm 
     data.D_a          := D_a;
     data.D_b          := D_b;
     data.Re_turbulent := 100;
     data.zetaLaminarKnown := true;
     data.c0 := 30;
      
     if D_a <= D_b then
        A_rel :=(D_a/D_b)^2;
        data.zeta1 :=(1 - A_rel)^2;
        data.zeta2 :=0.5*(1 - A_rel)^0.75;
        data.zeta1_at_a :=true;
        data.zeta2_at_a :=true;
        data.D_Re := D_a;
     else
        A_rel :=(D_b/D_a)^2;
        data.zeta1 :=0.5*(1 - A_rel)^0.75;
        data.zeta2 :=(1 - A_rel)^2;
        data.zeta1_at_a :=false;
        data.zeta2_at_a :=false;
        data.D_Re := D_b;
     end if;
   end suddenExpansion;
    
   encapsulated function sharpEdgedOrifice 
      "Compute pressure loss factors for sharp edged orifice (for both flow directions)" 
      import SI = Modelica.SIunits;
      import NonSI = Modelica.SIunits.Conversions.NonSIunits;
      import Modelica_Fluid.WorkInProgress.Utilities.PressureLossFactors;
      
     input SI.Diameter D_pipe 
        "Inner diameter of pipe (= same at port_a and port_b)";
     input SI.Diameter D_min "Smallest diameter of orifice";
     input SI.Diameter L "Length of orifice";
     input NonSI.Angle_deg alpha "Angle of orifice";
     output PressureLossFactors data 
        "Pressure loss factors for both flow directions";
     annotation (Icon(Rectangle(extent=[-100,60; 100,-60], style(
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
               fillPattern=8))),
                               Diagram(
                      Rectangle(extent=[-100,60; 100,-60], style(
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
             string="alpha")),
       Documentation(info="<html>
<p>
Loss factor for mass flow rate from port_a to port_b
(Idelchik 1994, diagram 4-14, p. 221):
</p>
<pre>
   zeta = [(1-A0/A1) + 0.707*(1-A0/A1)^0.375]^2*(A1/A0)^2 
          for Re(A0) >= 1e5,  independent of alpha
</pre>
<p>
Loss factor for mass flow rate from port_b to port_a
(Idelchik 1994, diagram 4-13, p. 220, with A2=A1):
</p>
<pre>
   zeta = k*(1 - A0/A1)^0.75 + (1 - A0/A1)^2 + 2*sqrt(k*(1-A0/A1)^0.375) + (1- A0/A1)
          k  = 0.13 + 0.34*10^(-(3.4*LD+88.4*LD^2.3)) 
               (there is a typing error in the formula in diagram 4-13, the above
                equation corresponds to table (a) in diagram 4-12)
          LD = L/D0
          for Re(A0) >= 1e4, 40 deg &le; alpha &le; 60 deg 
                             for other values of alpha, k is given as table
                             in diagram 3-7 (this is not yet included in the function)
</pre
</html>"));
    protected 
     Real D_rel = D_min/D_pipe;
     Real LD = L/D_min;
     Real k = 0.13 + 0.34*10^(-(3.4*LD+88.4*LD^2.3));
   algorithm 
     data.D_a          := D_pipe;
     data.D_b          := D_pipe;
     data.zeta1        := ((1-D_rel) + 0.707*(1-D_rel)^0.375)^2*(1/D_rel)^2;
     data.zeta2        := k*(1 - D_rel)^0.75 + (1 - D_rel)^2 +
                          2*sqrt(k*(1-D_rel)^0.375) + (1- D_rel);
     data.Re_turbulent := 1e4;
     data.D_Re         := D_min;
     data.zeta1_at_a   := true;
     data.zeta2_at_a   := false;
     data.zetaLaminarKnown := false;
     data.c0               := 0;
   end sharpEdgedOrifice;
 end PressureLossFactors;
  
  function diameter_of_squarePipe 
    "Determine hydraulic diameter of pipe with square cross sectional area" 
    import SI = Modelica.SIunits;
    extends Modelica.Icons.Function;
    input SI.Length width "Inner width of pipe";
    input SI.Length height " Inner height of pipe";
    output SI.Diameter D "Inner (hydraulic) diameter of pipe";
  algorithm 
    D :=4*width*height/(2*(width+height));
  end diameter_of_squarePipe;
  
    model FiniteVolume 
    "One dimensional volume according to the finite volume method with 1 mass, 1 energy and 2 momentum balances" 
    
    import SI = Modelica.SIunits;
    import Modelica.Math;
    
      replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
          annotation (choicesAllMatching = true);
    
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare model Medium = Medium) 
        annotation(extent=[-120, -10; -100, 10]);
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare model Medium =  Medium) 
        annotation(extent=[120, -10; 100, 10]);
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
        annotation(extent=[-10, 60; 10, 80], rotation=-90);
      Medium.BaseProperties medium(preferredMediumStates=true) 
      "Medium properties in the middle of the finite volume";
      SI.Mass M "Total mass in volume";
      SI.Mass[Medium.nXi] MXi "Independent component masses";
      SI.Energy U "Inner energy";
      parameter SI.Length L "Length of volume";
    
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
    
      input SI.Pressure dp_a 
      "Pressure loss due to pipe friction between port_a and middle of pipe";
      input SI.Pressure dp_b 
      "Pressure loss due to pipe friction between middle of pipe and port_b";
    
      Real my "Artifical viscosity";
    
      annotation (
        Documentation(info="<html>
<p>
Model <b>FiniteVolume</b> is a generic finite volume for 1-dim. thermo-fluid flow
in piping networks which has the following properties:
</p>
<ul>
<li> A FiniteVolume model is <b>independent</b> of the <b>medium</b> model. 
     The only requirement is that the medium model has to be
     a subclass of Modelica.Media.Interfaces.<b>PartialMedium</b>.
     As a consequence, the FiniteVolume model can be used
     for incompressible or compressible media, fluids with one 
     and multiple substances as well as for one and multiple phases. 
     The FiniteVolume model depends also not on the independent
     variables of the medium model.</li>
<li> A FiniteVolume model contains a <b>staggered grid</b> consisting
     of one volume from port_a to port_b for which the mass and
     energy balances are formulated and from two volumes from port_a 
     to the middle of the FiniteVolume and from the middle to
     port_b for which momentum
     balances are provided. When a FiniteVolume model is connected
     to another FiniteVolume, the adjacent momentum balance volumes
     are automatically merged together to form one volume. Thus, 
     a connected network of FiniteVolumes consists of a staggered
     grid for the balance equations.</li>
<li> For the intensive properties, such as density, temperature,
     an <b>upwind scheme</b> is used, i.e., depending on the direction
     of the mass flow rate, the upwind value is used in the equations.
     Also zero mass flow rate is handeled appropriately.</li>
<li> In order that the FiniteVolume model can be utilized, equations
     for the pipe friction have to be added via the input variables
     dp_a (pressure loss from port_a to the middle of the FiniteVolume)
     and dp_b (pressure loss from port_b to the middle of the FiniteVolume).</li>
<li> A FiniteVolume component contains <b>one medium</b> model in the 
     <b>middle</b> of the FiniteVolume.</li>
</ul>
</html>",     revisions="<html>
<ul>
<li><i>Aug. 30, 2004</i>
    by Hilding Elmqvist, Dynasim:<br>
    Further improvements + artifical viscosity introduced to remove
    unphysical oscillations in shock waves.</li>
<li><i>May 28, 2004</i>
    by Hilding Elmqvist, Dynasim:<br>
    Implemented.</li>
</ul>
</html>"),
        Diagram,
        Icon(Rectangle(extent=[-100, -60; 100, 60], style(
              color=3,
              rgbcolor={0,0,255},
              gradient=2,
              fillColor=76,
              rgbfillColor={170,170,255})),
          Ellipse(extent=[-16,16; 14,-12],    style(fillColor=0)),
          Rectangle(extent=[-90,46; 92,-46], style(color=0, rgbcolor={0,0,0})),
          Rectangle(extent=[-80,36; -6,-36], style(color=0, rgbcolor={0,0,0})),
          Rectangle(extent=[8,36; 82,-36], style(color=0, rgbcolor={0,0,0}))));
  protected 
      SI.MassFlowRate m_flow_a;
      SI.MassFlowRate m_flow_b;
      SI.MassFlowRate m_flow_middle;
      constant Real pi=Modelica.Constants.pi;
      constant Real g=Modelica.Constants.g_n;
      parameter SI.Area A_m=(A_a + A_b)/2;
      parameter SI.Length dx=L;
    equation 
      //Extensive properties
        M=medium.d*A_m*dx;
        MXi=M*medium.Xi;
        U=M*medium.u;
    
      // Mass balance over the interval a to b
      //der(medium.d)*A_m*dx = port_a.m_flow + port_b.m_flow;
      der(M)=port_a.m_flow + port_b.m_flow;
    
      // Substance mass balances over the interval a to b
      // der(medium.d*medium.X)*A_m*dx = port_a.mXi_flow + port_b.mXi_flow;
      //(der(medium.d)*medium.X + medium.d*der(medium.X))*A_m*dx = port_a.mXi_flow + port_b.mXi_flow;
      der(MXi)= port_a.mXi_flow + port_b.mXi_flow;
    
      // Energy balance over the interval a to b
      // der(medium.d*medium.u)*A_m*dx = port_a.H_flow + port_b.H_flow + m_flow_middle/
      //   medium.d*(port_b.p - port_a.p) + heatPort.Q_flow;
      //(der(medium.d)*medium.u + medium.d*der(medium.u))*A_m*dx = port_a.H_flow + port_b.H_flow - m_flow_middle/
      //  medium.d*(port_a.p - port_b.p - dp_a - dp_b) + heatPort.Q_flow;
      der(U)= port_a.H_flow + port_b.H_flow - m_flow_middle/medium.d*(port_a.p - port_b.p - dp_a - dp_b) + heatPort.Q_flow;
    
      m_flow_middle = (port_a.m_flow - port_b.m_flow)/2 
      "since assumed same density in entire interval a to b";
    
      // Momentum balance over interval a to dx/2
      (if dynamicMomentumBalance then der(m_flow_a)*dx/2 else 0) =
        A_m*(port_a.p - medium.p - dp_a) +
        (if includeKineticTerm then 
          - m_flow_middle^2/(A_m*medium.d) else 0)
        - A_m*medium.d/2*g*(Z_b - Z_a) +
        (if includeViscosity then my*((-m_flow_b-m_flow_a)/dx - 0) else 0);
        /* Removed: port_a.m_flow^2/(A_a*medium_a.d) */
    
      // Momentum balance over interval dx/2 to b
      (if dynamicMomentumBalance then -der(m_flow_b)*dx/2 else 0) =
        A_m*(medium.p - port_b.p - dp_b) +
        (if includeKineticTerm then 
          m_flow_middle^2/(A_m*medium.d) else 0)
        - A_m*medium.d/2*g*(Z_b - Z_a) +
        (if includeViscosity then my*(0 - (-m_flow_b-m_flow_a)/dx) else 0);
        /* Removed: - port_b.m_flow^2/(A_b*medium_b.d) */
    
       if includeViscosity then
         my = viscosityFactor1 + viscosityFactor2*dx*(if m_flow_middle*(-m_flow_b-m_flow_a) < 0 then 
            abs(-m_flow_b-m_flow_a)/(A_m*medium.d) else 0);
       else
         my = 0;
       end if;
    
      // Coupling to environment  
      m_flow_a = port_a.m_flow 
      "Due to problem with non-aliasing and semiLinear";
      m_flow_b = port_b.m_flow;
    
      // Upwind scheme (use properties from upwind port and handle zero flow)  
      port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium.h);
      port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium.h);
      port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium.Xi);
      port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium.Xi);
    
      // Heat port has the medium temperature
      heatPort.T = medium.T;
    
    end FiniteVolume;
  
model PipeSegment 
    "One segment of a pipe with 1 mass, 1 energy, 2 momementum balances and pipe friction" 
    
    import SI = Modelica.SIunits;
    import Modelica_Fluid.Types.Init.*;
  extends FiniteVolume(medium(
             p(start=p_start),
             T(start=T_start), h(start=h_start), Xi(start=X_start[1:Medium.nXi])));
  extends Modelica_Fluid.Interfaces.PartialInitializationParameters;
    
  parameter Boolean linearPressureDrop=true;
  parameter SI.AbsolutePressure dp_nominal(min=1.e-10) = 1 
      "Nominal pressure drop";
  parameter SI.MassFlowRate m_flow_nominal = 1E-3 
      "Nominal mass flow rate at nominal pressure drop";
    
  annotation (Documentation(info="<html>
<p>
Model <b>PipeSegment</b> describes one segment of a pipe.
It consists of the following parts:
</p>
<ul>
<li> One <a href=\"Modelica:Modelica_Fluid.Utilities.FiniteVolume\">FiniteVolume</a>
     model described by 1 mass, 1 energy, and 2 momemtum balances.</li>
<li> Different types of methods to initialize the FiniteVolume.</li>
<li> Different pipe friction models (ConstantLaminar, ConstantTurbulent,
     DetailedFriction) to describe the pressure loss due to the wall friction.</li>
</ul>
</html>"));
initial equation 
  // Initial conditions
  if initType == NoInit then
    // no initial equations
  elseif initType == InitialValues then
    if not Medium.singleState then
      medium.p = p_start;
    end if;
    if use_T_start then
      medium.T = T_start;
    else
      medium.h = h_start;
    end if;
    medium.Xi = X_start[1:Medium.nXi];
  elseif initType == SteadyState then
    if not Medium.singleState then
       der(medium.p) = 0;
    end if;
    der(medium.h) = 0;
    der(medium.Xi) = zeros(Medium.nXi);
  elseif initType == SteadyStateHydraulic then
    if not Medium.singleState then
       der(medium.p) = 0;
    end if;
    if use_T_start then
      medium.T = T_start;
    else
      medium.h = h_start;
    end if;
    medium.Xi = X_start[1:Medium.nXi];
  else
    assert(false, "Unsupported initialization option");
  end if;
equation 
  /*
  LongPipes.Components.PipeFriction friction[pipe.n](
    each from_dp=false, 
    each dp_nominal=500/pipe.n, 
    each roughness=1, 
    each diameter=30, 
    each length=length/pipe.n);
*/
    
  // Simple linear pressure drop in each segment
  dp_a = dp_nominal*(if linearPressureDrop then m_flow_a/m_flow_nominal else 
                        abs(m_flow_a)*m_flow_a/m_flow_nominal^2);
  dp_b = dp_nominal*(if linearPressureDrop then -m_flow_b/m_flow_nominal else 
                        abs(-m_flow_b)*(-m_flow_b)/m_flow_nominal^2);
end PipeSegment;
  
model PipeFriction 
    "Computes different types of pressure losses in pipes due to friction" 
    
    import SI = Modelica.SIunits;
    import FT = Modelica_Fluid.Types.FrictionTypes;
    import CT = Modelica_Fluid.Types.CrossSectionTypes;
    import Modelica.Math;
    
/* This model requires eta and d as input and provides
   an equation m_flow = f1 (dp) or dp = f2(m_flow)
*/
  input SI.DynamicViscosity eta 
      "Dummy or upstream dynamic viscosity for detailed friction model used for pressure loss calculation";
  input SI.Density d 
      "Dummy or upstream density for detailed friction model used for pressure loss calculation";
  SI.Pressure dp(start=0) "Pressure loss due to pipe friction";
  SI.MassFlowRate m_flow(start=0) "Mass flow rate from port_a to port_b";
    
  parameter Modelica_Fluid.Types.FrictionTypes.Temp frictionType=Modelica_Fluid.Types.
      FrictionTypes.ConstantTurbulent 
      "Type of friction to determine pressure loss";
  parameter SI.AbsolutePressure dp_nominal(min=1.e-10)=
    Modelica.SIunits.Conversions.from_bar(1.0) " Nominal pressure drop" 
    annotation (Dialog(enable=frictionType==FT.ConstantLaminar or frictionType==FT.ConstantTurbulent, group=
          "frictionType = ConstantLaminar or ConstantTurbulent"));
    
  parameter SI.MassFlowRate m_flow_nominal(min=1.e-10) = 1 
      " Nominal mass flow rate at nominal pressure drop" 
                                                       annotation (Dialog(
         enable=frictionType==FT.ConstantLaminar or frictionType==FT.ConstantTurbulent, group=
         "frictionType = ConstantLaminar or ConstantTurbulent"));
  parameter SI.Length length=1 " Length of pipe" 
    annotation (Dialog(enable=frictionType==FT.DetailedFriction, group="frictionType = DetailedFriction"));
  parameter SI.Length roughness=0 " Roughness of pipe" 
    annotation (Dialog(enable=frictionType==FT.DetailedFriction, group="frictionType = DetailedFriction"));
  parameter Modelica_Fluid.Types.CrossSectionTypes.Temp crossSectionType=
                     Modelica_Fluid.Types.CrossSectionTypes.Circular 
      " Type of cross section of pipe" 
    annotation (Dialog(enable=frictionType==FT.DetailedFriction, group="frictionType = DetailedFriction"));
  parameter SI.Diameter diameter=0.1 " Inner diameter of pipe" 
    annotation (Dialog(enable=frictionType==FT.DetailedFriction and crossSectionType==CT.Circular, group="frictionType = DetailedFriction"));
  parameter SI.Length width=0.05 " Inner width of pipe" 
    annotation (Dialog(enable=frictionType==FT.DetailedFriction and crossSectionType==CT.Rectangular, group="frictionType = DetailedFriction"));
  parameter SI.Length height=0.02 " Inner height of pipe" 
    annotation (Dialog(enable=frictionType==FT.DetailedFriction and crossSectionType==CT.Rectangular, group="frictionType = DetailedFriction"));
  parameter SI.Area area=0.01 " Cross sectional area of pipe" 
    annotation (Dialog(enable=frictionType==FT.DetailedFriction and crossSectionType==CT.General, group="frictionType = DetailedFriction"));
  parameter SI.Length perimeter=0.1 " Wetted perimeter of cross sectional area"
    annotation (Dialog(enable=frictionType==FT.DetailedFriction and crossSectionType==CT.General, group="frictionType = DetailedFriction"));
  parameter Boolean from_dp=true 
      " = true, use m_flow = f(dp) otherwise use dp = f(m_flow), i.e., inverse equation"
    annotation (Evaluate=true, Dialog(tab="Advanced"));
  parameter SI.Pressure p_small(min=1.e-10) = 1 
      " A small laminar region is introduced around p_small" 
                                                           annotation (Dialog(
        tab="Advanced", group="Only for frictionType = ConstantTurbulent"));
    
  annotation (
Images(Parameters(group="frictionType = ConstantLaminar or ConstantTurbulent", source=""),
       Parameters(group="frictionType = DetailedFriction", source="Images/PipeFriction1_small.png")),
structurallyIncomplete,
preferedView="info",
    Diagram,
    Icon,
    Documentation(info="<html>
<p>
This component models the pressure loss in a short pipe
due to friction under the assumption of quasi steady state flow (i.e., the
mass flow rate varies only slowly). This model is not complete
but may be used in a pipe model to provide an equation to compute
the friction pressure loss from the mass flow rate through
the pipe (see, e.g., <a href=\"Modelica://Modelica_Fluid.Components.ShortPipe\">Modelica_Fluid.Components.ShortPipe</a>).
</p>
<p>
Three loss models can be selected via
parameter <b>frictionType</b>:
</p>
<pre>
   frictionType = <b>ConstantLaminar</b>  :  dp =  k*m_flow
                = <b>ConstantTurbulent</b>:  dp =  k*m_flow^2  if m_flow &gt; 0
                                         = -k*m_flow^2  if m_flow &lt; 0
                = <b>DetailedFriction</b> :  dp = lambda(Re,Delta)*(L*rho/D)*v^2/2
                                         = lambda2(Re,Delta)*L*eta^2/(2*D^3*rho^3)
</pre>
<p>
where dp = \"port_a.p - port_b.p\" is the pressure loss and
m_flow is the mass flow rate from port_a to port_b.
</p>
<h3>ConstantLaminar and ConstantTurbulent</h3>
<p>
The pressure loss factor \"k\" is computed by providing the
mass flow rate \"m_flow_nominal\" and the corresponding
pressure loss \"dp_nominal\" for one flow condition
(usually the desired nominal flow condition). These factors might
be estimated or determined by measurements.
</p>
<p>
For \"ConstantTurbulent\" a small laminar region
is introduced around zero mass flow rate by interpolating
with a cubic polynomial (this technique is copied from the
ThermoFluid library).
</p>
<p>
The first two formulations are useful, if the pipe data is directly
measured and the main operating points are either fully in the
laminar or fully in the turbulent region. It would be better
for \"ConstantTurbulent\" to use the \"real\" laminar region. However,
then more data is required, especially the viscosity and the
diameter of the pipe.
</p>
<h3>DetailedFriction</h3>
<p>
The \"DetailedFriction\" option provides a detailed model
of frictional losses for commercial pipes with
<b>nonuniform roughness</b> (including the smooth pipe
as a special case). For pipes with circular cross section
the pressure loss is computed as:
</p>
<pre>
   dp = lambda*(L/D)*rho*v^2/2
      = lambda2*(L/(2*D^3))*(eta^2/rho)
        (with lambda2 = lambda*Re^2)
</pre>
<p>
where
</p>
<ul>
<li> L is the length of the pipe,</li>
<li> D is the diameter of the pipe,</li>
<li> lambda = lambda(Re,<font face=\"Symbol\">D</font>) is the \"usual\" friction coefficient,</li>
<li> lambda2 = lambda*Re^2 is the friction coefficient used in this model,</li>
<li> Re = v*D*rho/eta is the Reynolds number</li>
<li> <font face=\"Symbol\">D</font> = <font face=\"Symbol\">d</font>/D is the relative roughness where
     \"<font face=\"Symbol\">d</font>\" is
     the absolute \"roughness\", i.e., the averaged height of asperities in the pipe
     (<font face=\"Symbol\">d</font> may change over time due to growth of surface asperities during
      service, see <i>[Idelchick 1994, p. 85, Tables 2-1, 2-2])</i>,</li>
<li> rho is the density,</li>
<li> eta is the dynamic viscosity, and </li>
<li> v is the mean velocity.</li>
</ul>
<p>
The first form is usually given in books but is not suited
for a simulation program since lambda is infinity for zero mass flow rate.
The second form is the one implemented
in this model (lambda2=0 for zero mass flow rate).
The friction coefficient <b>lambda</b> is shown in the next figure:
</p>
<IMG SRC=\"../Images/Components/PipeFriction1.png\" ALT=\"PipeFriction1\">
<p>
More useful for a simulation model is the slightly
differently defined friction coefficient <b>lambda2</b> = lambda*Re^2,
as shown in the next figure:
</p>
<IMG SRC=\"../Images/Components/PipeFriction2.png\" ALT=\"PipeFriction2\">
<p>
<ul>
<li> For <b>Re &le; 2000</b>, the flow is <b>laminar</b> and the exact solution of the
     3-dim. Navier-Stokes equations (momentum and mass balance) is used under the
     assumptions of steady flow, constant pressure gradient and constant
     density and viscosity (= Hagen-Poiseuille flow). </li>
<li> For <b>Re &ge; 4000</b>, the flow is <b>turbulent</b>.
     Depending on the calculation direction (see \"Inverse formulation\"
     below) either of two explicite equations are used. If the pressure drop is assumed
     known (and therefore implicitly also lambda2), then the
     corresponding Reynolds number is computed with the Colebrook-White equation
     <i>[Colebrook 1939; Idelchik 1994, p. 83, eq. (2-9)]</i>.
     These are the <b>red</b> curves in the diagrams above.
     If the mass flow rate is assumed known (and therefore implicitly
     also the Reynolds number), then lambda2 is computed by an approximation of the
     inverse of the Colebrook-White equation <i>[Swamee and Jain 1976;
     Miller 1990, p. 191, eq.(8.4)]</i>.</li>
<li> For <b>2000 &le; Re &le; 4000</b> there is a transition region between laminar
     and turbulent flow. The value of lambda2 depends on more factors as just
     the Reynolds number and the relative roughness, therefore only crude approximations
     are possible in this area.<br>
     The deviation from the laminar region depends on the
     relative roughness. A laminar flow at Re=2000 is only reached for smooth pipes.
     The deviation Reynolds number Re1 is computed according to
     <i>[Samoilenko 1968; Idelchik 1994, p. 81, sect. 2.1.21].</i>
     These are the <b>blue</b> curves in the diagrams above.<br>
     Between Re1=Re1(<font face=\"Symbol\">d</font>/D) and Re2=4000, lambda2 is approximated by a cubic
     polynomial in the \"lg(lambda2) - lg(Re)\" chart (see figure above) such that the
     first derivative is continuous at these two points. In order to avoid
     the solution of non-linear equations, two different cubic polynomials are used
     for the direct and the inverse formulation. This leads to some discrepancies
     in lambda2 (= differences between the red and the blue curves).
     This is acceptable, because the transition region is anyway not
     precisely known since the actual friction coefficient depends on
     additional factors and since the operating points are usually
     not in this region.</li>
</ul>
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
<p>
The equations above are valid for incompressible flow.
They can also be applied for <b>compressible</b> flow up to about <b>Ma = 0.6</b>
(Ma is the Mach number) with a maximum error in lambda of about 3 %.
The effect of gas compressibility in a wide region can be taken into
account by the following formula derived by Voronin
<i>[Voronin 1959; Idelchick 1994, p. 97, sect. 2.1.81]</i>:
</p>
<pre>
  lambda_comp = lambda*(1 + (kappa-1)/2 * Ma^2)^(-0.47)
        kappa = cp/cv // specific heat ratio
</pre>
<p>
An appreciable decrease in the coefficent \"lambda_comp\" is observed
only in a narrow transonic region and also at supersonic flow velocities
by about 15% <i>[Idelchick 1994, p. 97, sect. 2.1.81]</i>.
</p>
<h3>Inverse formulation</h3>
<p>
In the \"Advanced menu\" it is possible via parameter
\"from_dp\" to define in which form the
loss equation is actually evaluated (<b>default</b> is from_dp = <b>true</b>):
</p>
<pre>
   from_dp = <b>true</b>:   m_flow = f1(dp)
           = <b>false</b>:  dp    = f2(m_flow)
</pre>
<p>
\"from_dp\" can be useful to avoid nonlinear systems of equations
in cases where the inverse pressure loss function is needed.
</p>
<p>
At the 34th Modelica meeting in Vienna it was discussed to introduce
a language element for alternatives, such that the tool can
figure out what alternative to use. If this would be available,
parameter from_dp could be removed and the equations would
be written as:
</p>
<pre>
  alternative
    // m_flow = f1(dp);
  or
    // dp = f2(m_flow);
  end alternative;
</pre>
<p>
The tool has then \"somehow\" to select the better alternative.
Further research is needed to develop appropriate symbolic
transformation algorithms.
If you have examples where this is an issue, please provide
them, in order that it is possible to experiment with.
</p>
<h3>References</h3>
<dl><dt>Colebrook F. (1939):</dt>
    <dd><b>Turbulent flow in pipes with particular reference to the transition
         region between the smooth and rough pipe laws</b>.
         J. Inst. Civ. Eng. no. 4, 14-25.</dd>
    <dt>Idelchik I.E. (1994):</dt>
    <dd><a href=\"http://www.begellhouse.com/books/00c0f05b040d2ec0.html\"><b>Handbook
        of Hydraulic Resistance</b></a>. 3rd edition, Begell House, ISBN
        0-8493-9908-4</dd>
    <dt>Miller D. S. (1990):</dt>
    <dd><b>Internal flow systems</b>.
    2nd edition. Cranfield:BHRA(Information Services).</dd>
    <dt>Samoilenko L.A. (1968):</dt>
    <dd><b>Investigation of the Hydraulic Resistance of Pipelines in the
        Zone of Transition from Laminar into Turbulent Motion</b>.
        Thesis (Cand. of Technical Science), Leningrad.</dd>
    <dt>Swamee P.K. and Jain A.K. (1976):</dt>
    <dd><b>Explicit equations for pipe-flow problems</b>.
         Proc. ASCE, J.Hydraul. Div., 102 (HY5), pp. 657-664.</dd>
    <dt>Voronin F.S. (1959):</dt>
    <dd><b>Effect of contraction on the friction coefficient in a
           turbulent gas flow</b>.
           Inzh. Fiz. Zh., vol. 2, no. 11, pp. 81-85.</dd>
</dl>
</html>", revisions="<html>
<h3>Author</h3>
<p>
<a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a><br>
Deutsches Zentrum f&uuml;r Luft und Raumfahrt e.V. (DLR)<br>
Institut f&uuml;r Robotik und Mechatronik<br>
Postfach 1116<br>
D-82230 Wessling<br>
Germany<br>
email: <A HREF=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</A><br>
</p>
</html>"));
  SI.ReynoldsNumber Re 
      "Dummy or Reynolds number of flow, if frictionType = DetailedFriction";
  Real lambda 
      "Dummy or friction coefficient, if frictionType = DetailedFriction";
  Real lambda2 
      "Dummy or non-standard friction coefficient, if frictionType = DetailedFriction (= lambda*Re^2)";
  final parameter Real Delta=roughness/D "Relative roughness";
    
  // Auxiliary variables for ConstantLaminar and ConstantTurbulent
  protected 
  parameter Real k=if frictionType == FT.ConstantLaminar then 
      dp_nominal/m_flow_nominal else (if frictionType == FT.ConstantTurbulent then 
     dp_nominal/m_flow_nominal^2 else length/(2*D*D*D)) 
      "Pressure loss coefficient (dp = k*f(m_flow))";
  parameter Real delta=if from_dp then p_small else sqrt(dp_nominal/k);
  parameter Real C1=if from_dp then 0.5/sqrt(delta) - 3.0*C3*delta^2 else 0.5
      *delta "Coefficient 1 of cubic polynomial in the laminar region";
  parameter Real C3=if from_dp then -0.25/(sqrt(delta)*delta^2) else 0.5/
      delta "Coefficient 3 of cubic polynomial in the laminar region";
    
  // Auxiliary variables for DetailedFriction model
  parameter SI.Diameter D=if crossSectionType == CT.Circular then 
            diameter else (if crossSectionType == CT.Rectangular then 
            4*width*height/(2*(width+height)) else 4*area/
      perimeter) "Diameter of pipe in SI units";
  parameter SI.ReynoldsNumber Re1=(745*exp(if Delta <= 0.0065 then 1 else 
      0.0065/Delta))^(if from_dp then 0.97 else 1) "Re leaving laminar curve";
  parameter SI.ReynoldsNumber Re2=4000 "Re entering turbulent curve";
    
  // point lg(lambda2(Re1)) with derivative at lg(Re1)
  parameter Real x1=if from_dp then Math.log10(64*Re1) else Math.log10(Re1);
  parameter Real y1=if from_dp then Math.log10(Re1) else Math.log10(64*Re1);
  parameter Real yd1=1;
    
  // Point lg(lambda2(Re2)) with derivative at lg(Re2)
  parameter Real aux1=(0.5/Math.log(10))*5.74*0.9;
  parameter Real aux2=Delta/3.7 + 5.74/Re2^0.9;
  parameter Real aux3=Math.log10(aux2);
  parameter Real L2=0.25*(Re2/aux3)^2;
  parameter Real aux4=2.51/sqrt(L2) + 0.27*Delta;
  parameter Real aux5=-2*sqrt(L2)*Math.log10(aux4);
  parameter Real x2=if from_dp then Math.log10(L2) else Math.log10(Re2);
  parameter Real y2=if from_dp then Math.log10(aux5) else Math.log10(L2);
  parameter Real yd2=if from_dp then 0.5 + (2.51/Math.log(10))/(aux5*aux4) else 
            2 + 4*aux1/(aux2*aux3*(Re2)^0.9);
    
  // Constants: Cubic polynomial between lg(Re1) and lg(Re2)
  parameter Real diff_x=x2 - x1;
  parameter Real m=(y2 - y1)/diff_x;
  parameter Real c2=(3*m - 2*yd1 - yd2)/diff_x;
  parameter Real c3=(yd1 + yd2 - 2*m)/(diff_x*diff_x);
  parameter Real lambda2_1=64*Re1;
  constant Real pi=Modelica.Constants.pi;
  Real dx;
  Real aux7;
equation 
  if frictionType <> FT.DetailedFriction then
    // Assign dummy values for auxiliary variables
    Re = 0;
    dx = 0;
    lambda = 0;
    lambda2 = 0;
    aux7 = 0;
  else
    lambda = noEvent(if Re < 64 then 1 else lambda2/(Re*Re));
  end if;
    
  if from_dp then
    // equations in the form m_flow = m_flow(dp)
    if frictionType == FT.ConstantLaminar then
      m_flow = dp/k;
    elseif frictionType == FT.ConstantTurbulent then
      m_flow = noEvent(if dp > delta then sqrt(dp) else (if dp < -delta then -
        sqrt(-dp) else (C1 + C3*dp*dp)*dp))/sqrt(k);
    else
      lambda2 = noEvent(d*abs(dp)/(k*eta*eta));
      if noEvent(lambda2/64 <= Re1) then
        aux7 = 0;
        dx = 0;
        Re = lambda2/64;
      else
        aux7 = -2*sqrt(lambda2)*Math.log10(2.51/sqrt(lambda2) + 0.27*Delta);
        dx = if noEvent(aux7 >= Re2) then 0 else Math.log10(lambda2/lambda2_1);
        Re = if noEvent(aux7 >= Re2) then aux7 else Re1*(lambda2/lambda2_1)^(
          1 + dx*(c2 + dx*c3));
      end if;
      m_flow = noEvent((pi*D/4)*eta*Re*(if dp >= 0 then +1 else -1));
    end if;
  else
    // equations in the form dp = dp(m_flow)
    if frictionType == FT.ConstantLaminar then
      dp = k*m_flow;
    elseif frictionType == FT.ConstantTurbulent then
      dp = k*noEvent(if m_flow > delta then m_flow*m_flow else (if m_flow < -
        delta then -m_flow*m_flow else (C1 + C3*m_flow*m_flow)*m_flow));
    else
      Re = noEvent((4/pi)*abs(m_flow)/(D*eta));
      dx = noEvent(if Re < Re1 or Re > Re2 then 0 else Math.log10(Re/Re1));
      lambda2 = noEvent(if Re <= Re1 then 64*Re else (if Re >= Re2 then 0.25*
        (Re/Math.log10(Delta/3.7 + 5.74/Re^0.9))^2 else 64*Re1*(Re/Re1)^(1 +
        dx*(c2 + dx*c3))));
      aux7 = 0;
      dp = noEvent(k*lambda2*eta*eta/d*(if m_flow >= 0 then 1 else -1));
    end if;
  end if;
end PipeFriction;
end Utilities;
