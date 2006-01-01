package Utilities 
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
  
  function evaluatePoly3_derivativeAtZero 
    "Evaluate polynomial of order 3 that passes the origin with a predefined derivative" 
    extends Modelica.Icons.Function;
    input Real x "Value for which polynomial shall be evaluated";
    input Real x1 "Abscissa value";
    input Real y1 "y1=f(x1)";
    input Real y1d "First derivative at y1";
    input Real y0d "First derivative at f(x=0)";
    output Real y;
    annotation(smoothOrder=3);
  protected 
    Real a1;
    Real a2;
    Real a3;
    Real xx;
  algorithm 
    a1 := x1*y0d;
    a2 := 3*y1 - x1*y1d - 2*a1;
    a3 := y1 - a2 - a1;
    xx := x/x1;
    y  := xx*(a1 + xx*(a2 + xx*a3));
  end evaluatePoly3_derivativeAtZero;
  
  function regRoot2 
    "Anti-symmetric square root approximation with finite derivative in the origin" 
    
    extends Modelica.Icons.Function;
    input Real x;
    input Real x_small(min=0)=0.01 
      "approximation of function for |x| <= x_small";
    input Real k1(min=0)=1 "y = if x>=0 then sqrt(k1*x) else -sqrt(k2*|x|)";
    input Real k2(min=0)=1 "y = if x>=0 then sqrt(k1*x) else -sqrt(k2*|x|)";
    input Boolean use_yd0 = false "= true, if yd0 shall be used";
    input Real yd0(min=0)=1 "Desired derivative at x=0: dy/dx = yd0";
    output Real y;
    annotation(smoothOrder=2, Documentation(info="<html>
<p>
Approximates the function
</p>
<pre>
   y = <b>if</b> x &ge; 0 <b>then</b> <b>sqrt</b>(k1*x) <b>else</b> -<b>sqrt</b>(k2*<b>abs</b>(x)), with k1, k2 > 0
</pre>
<p>
in such a way that within the region -x_small &le; x &le; x_small, 
the function is described by two polynomials of third order
(one in the region -x_small .. 0 and one within the region 0 .. x_small)
such that the 
</p>
<ul>
<li> derviative at x=0 is finite, </li>
<li> the overall function is continuous with a first
     continuous derivative everywhere, and<li>
<li> the first and second derivative at x=0 of the
     two polynomials are identical.
</ul>
<p>
Typical screenshots for two different configurations
are shown below. The first one with k1=k2=1
</p>
<p>
<img src=\"../Images/Components/regRoot2_a.png\">
</p>
<p>
and the second one with k1=1 and k2=3. In the last
figure the (smooth) derivative of the function with
k1=1, k2=3 is shown:
</p>
<p>
<img src=\"../Images/Components/regRoot2_b.png\">
</p>
<p>
<img src=\"../Images/Components/regRoot2_c.png\">
</p>
</html>"));
  protected 
    encapsulated function regRoot2_utility 
      "Interpolating with two 3-order polynomials with a prescribed derivative at x=0" 
      import 
        Modelica_Fluid.WorkInProgress.Utilities.evaluatePoly3_derivativeAtZero;
       input Real x;
       input Real x1 "approximation of function abs(x) < x1";
       input Real k1 "y = if x>=0 then sqrt(k1*x) else -sqrt(k2*|x|); k1 >= k2";
       input Real k2 "y = if x>=0 then sqrt(k1*x) else -sqrt(k2*|x|))";
       input Boolean use_yd0 "= true, if yd0 shall be used";
       input Real yd0(min=0) "Desired derivative at x=0: dy/dx = yd0";
       output Real y;
       annotation(smoothOrder=2);
    protected 
       Real x2;
       Real xsqrt1;
       Real xsqrt2;
       Real y1;
       Real y2;
       Real y1d;
       Real y2d;
       Real w;
       Real y0d;
       Real w1;
       Real w2;
    algorithm 
       //x2 :=-x1*(k2/k1);
       x2 :=-x1;
       if x <= x2 then
          y := -sqrt(k2*abs(x));
       else
          y1 :=sqrt(k1*x1);
          y2 :=-sqrt(k2*abs(x2));
          y1d :=sqrt(k1/x1)/2;
          y2d :=sqrt(k2/abs(x2))/2;
        
          if use_yd0 then
             y0d :=yd0;
          else
             /* Determine derivative, such that first and second derivative
              of left and right polynomial are identical at x=0:
           _
           Basic equations:
              y_right = a1*(x/x1) + a2*(x/x1)^2 + a3*(x/x1)^3
              y_left  = b1*(x/x2) + b2*(x/x2)^2 + b3*(x/x2)^3
              yd_right*x1 = a1 + 2*a2*(x/x1) + 3*a3*(x/x1)^2
              yd_left *x2 = b1 + 2*b2*(x/x2) + 3*b3*(x/x2)^2
              ydd_right*x1^2 = 2*a2 + 6*a3*(x/x1)
              ydd_left *x2^2 = 2*b2 + 6*b3*(x/x2)
           _
           Conditions (6 equations for 6 unknowns):
                     y1 = a1 + a2 + a3
                     y2 = b1 + b2 + b3
                 y1d*x1 = a1 + 2*a2 + 3*a3
                 y2d*x2 = b1 + 2*b2 + 3*b3
                    y0d = a1/x1 = b1/x2
                   y0dd = 2*a2/x1^2 = 2*b2/x2^2
           _
           Derived equations:
              b1 = a1*x2/x1
              b2 = a2*(x2/x1)^2
              b3 = y2 - b1 - b2
                 = y2 - a1*(x2/x1) - a2*(x2/x1)^2
              a3 = y1 - a1 - a2
           _
           Remaining equations
              y1d*x1 = a1 + 2*a2 + 3*(y1 - a1 - a2)
                     = 3*y1 - 2*a1 - a2
              y2d*x2 = a1*(x2/x1) + 2*a2*(x2/x1)^2 +
                       3*(y2 - a1*(x2/x1) - a2*(x2/x1)^2)
                     = 3*y2 - 2*a1*(x2/x1) - a2*(x2/x1)^2
              y0d    = a1/x1
           _
           Solving these equations results in y0d below
           (note, the denominator "(1-w)" is always non-zero, because w is negative) 
           */
             w :=x2/x1;
             y0d := ( (3*y2 - x2*y2d)/w - (3*y1 - x1*y1d)*w) /(2*x1*(1 - w));
          end if;
        
          /* Modify derivative y0d, such that the polynomial is 
           monotonically increasing. A sufficient condition is
             0 <= y0d <= sqrt(8.75*k_i/|x_i|)
        */
          w1 :=sqrt(8.75*k1/x1);
          w2 :=sqrt(8.75*k2/abs(x2));
          y0d :=min(y0d, 0.9*min(w1, w2));
        
          /* Perform interpolation in scaled polynomial:
           y_new = y/y1
           x_new = x/x1
        */
          y := y1*(if x >= 0 then evaluatePoly3_derivativeAtZero(x/x1,1,1,y1d*x1/y1,y0d*x1/y1) else 
                                  evaluatePoly3_derivativeAtZero(x/x1,x2/x1,y2/y1,y2d*x1/y1,y0d*x1/y1));
       end if;
    end regRoot2_utility;
  algorithm 
    y := smooth(2,if x >= x_small then sqrt(k1*x) else 
                  if x <= -x_small then -sqrt(k2*abs(x)) else 
                  if k1 >= k2 then regRoot2_utility(x,x_small,k1,k2,use_yd0,yd0) else 
                                  -regRoot2_utility(-x,x_small,k2,k1,use_yd0,yd0));
  end regRoot2;
  
  function regSquare2 
    "Anti-symmetric square approximation with non-zero derivative in the origin" 
    extends Modelica.Icons.Function;
    input Real x;
    input Real x_small(min=0)=0.01 
      "approximation of function for |x| <= x_small";
    input Real k1(min=0)=1 "y = (if x>=0 then k1 else k2)*x*|x|";
    input Real k2(min=0)=1 "y = (if x>=0 then k1 else k2)*x*|x|";
    input Boolean use_yd0 = false "= true, if yd0 shall be used";
    input Real yd0(min=0)=1 "Desired derivative at x=0: dy/dx = yd0";
    output Real y;
    annotation(smoothOrder=2, Documentation(info="<html>
<p>
Approximates the function
</p>
<pre>
   y = <b>if</b> x &ge; 0 <b>then</b> k1*x*x <b>else</b> -k2*x*x, with k1, k2 > 0
</pre>
<p>
in such a way that within the region -x_small &le; x &le; x_small, 
the function is described by two polynomials of third order
(one in the region -x_small .. 0 and one within the region 0 .. x_small)
such that the 
</p>
<ul>
<li> derviative at x=0 is not zero (in order that the
     inverse of the function does not have an infinite derivative), </li>
<li> the overall function is continuous with a first
     continuous derivative everywhere, and<li>
<li> the first and second derivative at x=0 of the
     two polynomials are identical.
</ul>
<p>
Typical screenshots for k1=1, k2=3 are shown below.
In the last
figure the (smooth) derivative of the function with
k1=1, k2=3 is shown:
</p>
<p>
<img src=\"../Images/Components/regSquare2_b.png\">
</p>
<p>
<img src=\"../Images/Components/regSquare2_c.png\">
</p>
</html>"));
  protected 
    encapsulated function regSquare2_utility 
      "Interpolating with two 3-order polynomials with a prescribed derivative at x=0" 
      import 
        Modelica_Fluid.WorkInProgress.Utilities.evaluatePoly3_derivativeAtZero;
       input Real x;
       input Real x1 "approximation of function abs(x) < x1";
       input Real k1 "y = (if x>=0 then k1 else -k2)*x*|x|; k1 >= k2";
       input Real k2 "y = (if x>=0 then k1 else -k2)*x*|x|";
       input Boolean use_yd0 = false "= true, if yd0 shall be used";
       input Real yd0(min=0)=1 "Desired derivative at x=0: dy/dx = yd0";
       output Real y;
       annotation(smoothOrder=2);
    protected 
       Real x2;
       Real y1;
       Real y2;
       Real y1d;
       Real y2d;
       Real w;
       Real w1;
       Real w2;
       Real y0d;
    algorithm 
       // x2 :=-x1*(k2/k1)^2;
       x2 := -x1;
       if x <= x2 then
          y := -k2*x^2;
       else
           y1 := k1*x1^2;
           y2 :=-k2*x2^2;
          y1d := k1*2*x1;
          y2d :=-k2*2*x2;
          if use_yd0 then
             y0d :=yd0;
          else
             /* Determine derivative, such that first and second derivative
              of left and right polynomial are identical at x=0:
              see derivation in function regRoot2
           */
             w :=x2/x1;
             y0d := ( (3*y2 - x2*y2d)/w - (3*y1 - x1*y1d)*w) /(2*x1*(1 - w));
          end if;
        
          /* Modify derivative y0d, such that the polynomial is 
           monotonically increasing. A sufficient condition is
             0 <= y0d <= sqrt(5)*k_i*|x_i|
        */
          w1 :=sqrt(5)*k1*x1;
          w2 :=sqrt(5)*k2*abs(x2);
          y0d :=min(y0d, 0.9*min(w1, w2));
        
          y := if x >= 0 then evaluatePoly3_derivativeAtZero(x,x1,y1,y1d,y0d) else 
                              evaluatePoly3_derivativeAtZero(x,x2,y2,y2d,y0d);
       end if;
    end regSquare2_utility;
  algorithm 
    y := smooth(2,if x >= x_small then k1*x^2 else 
                  if x <= -x_small then -k2*x^2 else 
                  if k1 >= k2 then regSquare2_utility(x,x_small,k1,k2,use_yd0,yd0) else 
                                  -regSquare2_utility(-x,x_small,k2,k1,use_yd0,yd0));
  end regSquare2;
  
  model test_regRoot_derivatives "Test whether regRoot2 can be differentiated" 
    extends Modelica.Icons.Example;
    parameter Real x_small = 0.01;
    Real x;
    Real y;
    Real yd;
    Real k1;
    Real k2;
    annotation (experiment(StopTime=2, NumberOfIntervals=5000),
                                        experimentSetupOutput);
  equation 
    x = time - 1;
    k1 = 1 + 0.1*time;
    k2 = 2 + 0.2*time;
    
    y = (if x >= 0 then k1 else k2)*Modelica_Fluid.Utilities.regRoot(x,x_small);
    yd = der(y);
  end test_regRoot_derivatives;
  
  model test_regRoot2_derivatives "Test whether regRoot2 can be differentiated" 
    extends Modelica.Icons.Example;
    parameter Real x_small = 0.01;
    parameter Real k1a = 1;
    parameter Real k1b = 0.2;
    parameter Real k2a = 2;
    parameter Real k2b = 0.2;
    Real x;
    Real y;
    Real yd;
    Real ydd;
    Real k1;
    Real k2;
    Real y2;
    Real y2d;
    Real y2dd;
    annotation (experiment(StopTime=2, NumberOfIntervals=5000),
                                        experimentSetupOutput);
  equation 
    x = time - 1;
    
    k1 = k1a + k1b*time;
    k2 = k2a + k2b*time;
    
    y = regRoot2(x,x_small, k1, k2);
    yd = der(y);
    ydd = der(yd);
    
    y2 = regRoot2(x,x_small, k1, k2, true, 10);
    y2d = der(y2);
    y2dd = der(y2d);
  end test_regRoot2_derivatives;
  
  model test_regSquare2_derivatives 
    "Test whether regSquare2 can be differentiated" 
    extends Modelica.Icons.Example;
    parameter Real x_small = 0.1;
    parameter Real k1a = 1;
    parameter Real k1b = 0.2;
    parameter Real k2a = 2;
    parameter Real k2b = 0.2;
    Real x;
    Real y;
    Real yd;
    Real ydd;
    Real k1;
    Real k2;
    Real y2;
    Real y2d;
    Real y2dd;
    annotation (experiment(StopTime=2, NumberOfIntervals=5000),
                                        experimentSetupOutput);
  equation 
    x = time - 1;
    
  /*
  k1 = k1a + k1b*time;
  k2 = k2a + k2b*time;
*/
    k1=1;
    k2=3;
    
    y = regSquare2(x,x_small, k1, k2);
    yd = der(y);
    ydd = der(yd);
    
    y2 = regSquare2(x,x_small, k1, k2, true, 0.02);
    y2d = der(y2);
    y2dd = der(y2d);
  end test_regSquare2_derivatives;
  
    model FiniteVolume 
    "One dimensional volume according to the finite volume method with 1 mass, 1 energy and 2 momentum balances" 
    
      import SI = Modelica.SIunits;
      import Modelica.Math;
    
      replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
          annotation (choicesAllMatching = true);
    
      Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare model Medium = Medium) 
        annotation(extent=[-120, -10; -100, 10]);
      Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare model Medium =  Medium) 
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
  import Modelica_Fluid.Types.InitTypes.*;
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
  if initOption == NoInit then
    // no initial equations
  elseif initOption == InitialValues then
    if not Medium.singleState then
      medium.p = p_start;
    end if;
    if use_T_start then
      medium.T = T_start;
    else
      medium.h = h_start;
    end if;
    medium.Xi = X_start[1:Medium.nXi];
  elseif initOption == SteadyState then
    if not Medium.singleState then
       der(medium.p) = 0;
    end if;
    der(medium.h) = 0;
    der(medium.Xi) = zeros(Medium.nXi);
  elseif initOption == SteadyStateHydraulic then
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
  
  package PipeFriction 
    partial model PartialPipeFriction 
      parameter Boolean from_dp=true 
        " = true, use m_flow = f(dp) otherwise use dp = f(m_flow), i.e., inverse equation"
        annotation (Evaluate=true, Dialog(tab="Advanced"));
      parameter Boolean lumped_dp=true "=true, use lumped pressure drop" annotation(Evaluate=true, Dialog(tab="No input", enable=false));
      parameter Integer n(min=1) = 1 "Pipe segmentation" annotation(Dialog(tab="No input", enable=false));
      parameter Integer np(min=1) = 1 "Number of friction loss components" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Length length "Length of pipe" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Diameter d_h "Hydraulic diameter"          annotation(Dialog(tab="No input", enable=false));
      SI.Pressure[np] dp "Pressure drop due to friction loss";
      SI.MassFlowRate[np] m_flow "Mass flow rate";
      input SI.Density[np] d_1 
        "if m_flow>=0 then upstream else downstream density";
      input SI.Density[np] d_2 
        "if m_flow>=0 then downstream else upstream density";
      replaceable package Medium = PackageMedium extends 
        Modelica.Media.Interfaces.PartialMedium annotation(Dialog(tab="No input", enable=false));
      annotation (Icon(Ellipse(extent=[-60,60; 60,-60], style(
              color=0,
              rgbcolor={0,0,0},
              gradient=3,
              fillColor=69,
              rgbfillColor={0,128,255})),
                                        Text(
            extent=[-38,22; 40,-18],
            string="%name",
            style(
              color=0,
              rgbcolor={0,0,0},
              gradient=3,
              fillColor=69,
              rgbfillColor={0,128,255},
              fillPattern=7))));
    end PartialPipeFriction;
    
    model GenericFrictionLoss "Work around for partial" 
      extends Modelica_Fluid.WorkInProgress.Interfaces.PartialPressureLoss;
    annotation(structurallyIncomplete);
    end GenericFrictionLoss;
    
    model PipeFriction_RoughSurface 
      extends PartialPipeFriction;
      parameter Boolean use_Re = false 
        "= true, if turbulent region is defined by Re, otherwise by p_small or m_flow_small"
        annotation(Evaluate=true, Dialog(tab="Advanced"));
      parameter SI.AbsolutePressure dp_small = 1 
        "Turbulent flow if |dp| >= dp_small" 
        annotation(Dialog(tab="Advanced", enable=not use_Re and from_dp));
      parameter SI.MassFlowRate m_flow_small = 0.0001 
        "Turbulent flow if |m_flow| >= m_flow_small" 
        annotation(Dialog(tab="Advanced", enable=not use_Re and not from_dp));
      
      final parameter SI.Length[np] dx=if lumped_dp then ones(np)*length else {if i > 1 and 
          i < np then length/n else length/n/2 for i in 1:np} 
        "Length of pipe segment, two halves on each end if lumped_dp=false";
      
      parameter SI.Length roughness(min=1e-10) = 0.05 
        "Absolute roughness of pipe (> 0 required)";
      outer Medium.ThermodynamicState[n] state;
      outer Medium.ThermodynamicState state_a;
      outer Medium.ThermodynamicState state_b;
      SI.DynamicViscosity[np] eta_1 
        "Upstream dynamic viscosity if m_flow >= 0, dummy variable if use_Re=false";
      SI.DynamicViscosity[np] eta_2 
        "Downstream dynamic viscosity if m_flow >= 0, dummy variable if use_Re=false";
      GenericFrictionLoss[np] ploss(
         m_flow = m_flow,
         dp = dp,
         d_a = d_1,
         d_b = d_2,
         eta_a = eta_1,
         eta_b = eta_2,
         each dp_small=dp_small,
         each m_flow_small=m_flow_small,
         each from_dp=from_dp,
         each use_Re = use_Re,
         final lossFactors=
         Modelica_Fluid.WorkInProgress.Utilities.PressureLossFactors.wallFriction(
         length=dx, diameter=d_h, roughness=roughness));
      
    annotation(structurallyIncomplete, Documentation(info="<html>
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
    equation 
      for i in 1:np loop
        if use_Re then
          eta_1[i] = if i > 1 then Medium.dynamicViscosity(state[i - 1]) else 
            Medium.dynamicViscosity(state_a);
          eta_2[i] = if i < np then Medium.dynamicViscosity(state[i]) else 
            Medium.dynamicViscosity(state_b);
        else
          eta_1[i] = 0;
          eta_2[i] = 0;
        end if;
      end for;
    end PipeFriction_RoughSurface;
    
    model PipeFriction_SimpleLinear 
      extends PartialPipeFriction;
      parameter SI.MassFlowRate m_flow_nominal "Nominal mass flow rate";
      parameter SI.Density d_nominal "Nominal density";
      parameter SI.Pressure dp_nominal "Nominal pressure drop";
    annotation(structurallyIncomplete);
    equation 
      if lumped_dp then
        dp[1] = dp_nominal*d_nominal/m_flow_nominal*m_flow[1]/(if from_dp then (if 
          noEvent(dp[1] >= 0) then d_1[1] else d_2[1]) else (if noEvent(m_flow[1] >= 0) then 
                d_1[1] else d_2[1]));
      else
        for i in 1:np loop
          dp[i] = dp_nominal/(if (i > 1 and i < np) then n else 2*n)*d_nominal/m_flow_nominal*m_flow[i]/(if from_dp then (
            if noEvent(dp[i] >= 0) then d_1[i] else d_2[i]) else (if noEvent(m_flow[
            i] >= 0) then d_1[i] else d_2[i]));
        end for;
      end if;
    end PipeFriction_SimpleLinear;
  end PipeFriction;
  
  package PipeHeatTransfer 
    
    partial model PartialPipeHeatTransfer 
      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(Dialog(tab="No input", enable=false));
      parameter Integer n(min=1)=1 "Number of pipe segments" annotation(Dialog(tab="No input", enable=false));
      SI.HeatFlowRate[n] Q_flow "Heat flow rates";
      parameter SI.Area A_h "Total heat transfer area" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Length d_h "Hydraulic diameter" annotation(Dialog(tab="No input", enable=false));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
        "Thermal port" 
        annotation (extent=[-20,60; 20,80]);
      SI.Temperature[n] T;
    equation 
      
      annotation (Icon(Ellipse(extent=[-60,64; 60,-56], style(
              color=42,
              rgbcolor={127,0,0},
              gradient=3,
              fillColor=1,
              rgbfillColor={232,0,0})), Text(
            extent=[-38,26; 40,-14],
            style(
              color=42,
              rgbcolor={127,0,0},
              gradient=3,
              fillColor=1,
              rgbfillColor={232,0,0},
              fillPattern=7),
            string="%name")));
    end PartialPipeHeatTransfer;
    
    partial model PipeHT_constAlpha 
      extends PartialPipeHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0=200;
    equation 
      for i in 1:n loop
        thermalPort[i].Q_flow=alpha0*A_h/n*(thermalPort[i].T-T[i]);
      end for;
      thermalPort.Q_flow=Q_flow;
    end PipeHT_constAlpha;
  end PipeHeatTransfer;
end Utilities;
