package Components "Basic components for fluids" 
  extends Modelica.Icons.Library;
  
  annotation (Documentation(info="<html>
<p>
This package contains examples for basic component models of 
the fluid library.
</p>
</html>"));
  
  model Orifice 
    "Simple orifice with constant loss factor zeta and small laminar region" 
    
    import SI = Modelica.SIunits;
    import Modelica.SIunits.Conversions.*;
    extends Interfaces.PartialTwoPortTransport;
    parameter Real zeta=1 "Loss factor from: dp = 0.5*rho*zeto*v^2";
    parameter SI.Area A=0.01 "Area of orifice";
    parameter Real C(unit="m^3/(Pa.s)") = 1 
      "Flow conductance for small pressure drops";
    parameter Boolean dp_given=true 
      "|Advanced|| True, if m_dot is computed as function of pressure drop dp (otherwise, use inverse function to avoid nonlinear loop)";
    annotation (
      Diagram, 
      Icon(
        Text(
          extent=[-126, -76; 130, -110], 
          style(color=0), 
          string="zeta=%zeta"), 
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
This component should model an orifice with constant loss
factor zeta. This means for turbulent flow the pressure loss dp = port_a.p - port_b.p
is computed as function of the upstream velocity v and the
upstream density rho as:
</p>
<pre>
     dp = 0.5*rho*zeta*v^2
</pre>
<p>
A small laminar region is included to treat reversing and
small mass flow rate in an appropriate way.
</p>
<p>
Note, this element will be changed. It is provided here
just in order to build up benchmarks. In particular parameter \"C\"
should be removed and the loss factor zeta should be
optionally internally computed from technological parameters
(using tables from the books of Miller or Idelchik).
</p>
<p>
In the 'Advanced' menu it is possible to select which
form of the equation should be used (direct or inverse equation).
In some situations this may avoid non-linear systems of equations.
At the last Modelica meeting it was discussed to introduce
a language element for alternatives, such that the tool can
figure out what alternative to use. If this would be available,
parameter dp_given could be removed and the equations would
be written as:
</p>
<pre>
  alternative
     m_dot = f1(dp)
  or
     dp = f2(m_dot)
  end alternative;
</pre>
<p>
The tool has then \"somehow\" to select the better alternative.
If you have examples where this is an issue, please provide
them, in order that it is possible to experiment with.
</p>
</html>"));
    SI.Pressure dp "Pressure loss";
  protected 
    parameter Real k=C*zeta/(2*A*A);
  equation 
    /* Function is a second order polynomial in the mass flow rate,
     i.e., it is continuous and differentiable upto any order.
         rho*C*dp = m_dot + k*m_dot^2
     At low mass flow rates the "k*m_dot^2" term can be neglected
     and we have laminar flow:
         m_dot = rho*C*dp
     At high mass flow rates the linear "m_dot" term can be neglected
     and we have turbulent flow:
         rho*C*dp = k*m_dot^2, i.e.,
               dp = k/(rho*C) * m_dot^2
     on the other hand we have for turbulent flow:
         dp = 0.5*rho*zeta*v^2
            = 0.5*rho*zeta*(m_dot/(rho*A))^2
            = 0.5*zeta/(rho*A^2) * m_dot^2
     Comparision results in
         k = C*zeta/(2*A^2)
  */
    dp = port_a.p - port_b.p;
    if dp_given then
      port_a.m_dot = noEvent((if dp >= 0 then 1 else -1)*(-1 + sqrt(1 + 4*k*C*
        abs(dp)))/(2*k));
    else
      C*dp = port_a.m_dot + k*port_a.m_dot^2*noEvent(if port_a.m_dot >= 0 then 
        +1 else -1);
    end if;
  end Orifice;
  
  model ShortPipe 
    "Short pipe where mass flow rate is a function of pressure drop (only transport, no storage of mass or energy)"
     
    
    /* This currently gives not a nice layout
      Images(Parameters(name="|frictionType = DetailedFriction|", source=
            "../../Images/PipeFriction1_small.png")),
*/
    
    import SI = Modelica.SIunits;
    import Modelica.Math;
    import Modelica.SIunits.Conversions.*;
    import Modelica_Fluid.Examples.Types.FrictionTypes;
    import Modelica_Fluid.Examples.Types.CrossSectionTypes;
    import FT = Modelica_Fluid.Examples.Types.FrictionTypes;
    extends Modelica_Fluid.Interfaces.PartialTwoPortTransport;
    SI.Pressure dp "Pressure loss due to friction";
    Real zero=port_a.p - port_b.p - dp "momentum balance (may be modified)";
    
    parameter FT.Temp frictionType=FT.ConstantTurbulent 
      "Type of friction to determine pressure loss";
    parameter Medium.AbsolutePressure dp_nominal(min=1.e-10) = from_bar(1.0) 
      "|frictionType = ConstantLaminar or ConstantTurbulent| Nominal pressure drop";
    
    parameter Medium.MassFlowRate m_dot_nominal(min=1.e-10) = 1 
      "|frictionType = ConstantLaminar or ConstantTurbulent| Nominal mass flow rate at nominal pressure drop";
    parameter Types.Length_mm length=1000 
      "|frictionType = DetailedFriction| Length of pipe";
    parameter Types.Length_mm roughness=0 
      "|frictionType = DetailedFriction| Roughness of pipe";
    parameter Modelica_Fluid.Examples.Types.CrossSectionTypes.Temp crossSectionType=
        Modelica_Fluid.Examples.Types.CrossSectionTypes.Circular 
      "|frictionType = DetailedFriction| Type of cross section of pipe";
    parameter Types.Length_mm diameter=100 
      "|crossSectionType = circular| Inner diameter of pipe";
    parameter Types.Length_mm width=50 
      "|crossSectionType = rectangular| Inner width of pipe";
    parameter Types.Length_mm height=20 
      "|crossSectionType = rectangular| Inner height of pipe";
    parameter SI.Area area=0.01 
      "|crossSectionType = general| Cross sectional area of pipe";
    parameter SI.Length perimeter=0.1 
      "|crossSectionType = general| Wetted perimeter of cross sectional area";
    parameter Boolean from_dp=true 
      "|Advanced|| = true, use m_dot = f(dp) otherwise use dp = f(m_dot), i.e., inverse equation"
      annotation (Evaluate=true);
    parameter SI.Pressure p_small(min=1.e-10) = 1 
      "|Advanced|Only for frictionType = ConstantTurbulent| A small laminar region is introduced around p_small";
    annotation (
      Diagram, 
      Icon(
        Rectangle(extent=[-100, 60; 100, -60], style(
            color=0, 
            gradient=2, 
            fillColor=8)), 
        Rectangle(extent=[-100, 34; 100, -36], style(
            color=69, 
            gradient=2, 
            fillColor=69)), 
        Text(
          extent=[-120, 130; 116, 64], 
          string="%name", 
          style(gradient=2, fillColor=69)), 
        Text(
          extent=[-132, -64; 140, -94], 
          style(color=0), 
          string="%m_dot_nominal"), 
        Text(
          extent=[-130, -96; 142, -126], 
          style(color=0), 
          string="%dp_nominal")), 
      Documentation(info="<html>
<p>
This component models the pressure loss in a short pipe
under the assumptions that no mass or energy is stored in the
component and that the flow is quasi steady state (i.e., the
mass flow rate varies only slowly). Three loss models can be selected via
parameter <b>frictionType</b>:
</p>
<pre>
   frictionType = <b>ConstantLaminar</b>  :  dp =  k*m_dot
                = <b>ConstantTurbulent</b>:  dp =  k*m_dot^2  if m_dot &gt; 0
                                         = -k*m_dot^2  if m_dot &lt; 0
                = <b>DetailedFriction</b> :  dp = lambda(Re,Delta)*(L*rho/D)*v^2/2
                                         = lambda2(Re,Delta)*L*eta^2/(2*D^3*rho^3)
</pre>
<p>
where dp = \"port_a.p - port_b.p\" is the pressure loss and
m_dot is the mass flow rate from port_a to port_b.
</p>
<h3>ConstantLaminar and ConstantTurbulent</h3>
<p>
The pressure loss factor \"k\" is computed by providing the
mass flow rate \"m_dot_nominal\" and the corresponding
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
   from_dp = <b>true</b>:   m_dot = f1(dp)
           = <b>false</b>:  dp    = f2(m_dot)
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
    // m_dot = f1(dp);
  or
    // dp = f2(m_dot);
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
    parameter Real k=if frictionType == FrictionTypes.ConstantLaminar then 
        dp_nominal/m_dot_nominal else (if frictionType == FrictionTypes.
        ConstantTurbulent then dp_nominal/m_dot_nominal^2 else L/(2*D*D*D)) 
      "Pressure loss coefficient (dp = k*f(m_dot))";
    parameter Real delta=if from_dp then p_small else sqrt(dp_nominal/k);
    parameter Real C1=if from_dp then 0.5/sqrt(delta) - 3.0*C3*delta^2 else 0.5
        *delta "Coefficient 1 of cubic polynomial in the laminar region";
    parameter Real C3=if from_dp then -0.25/(sqrt(delta)*delta^2) else 0.5/
        delta "Coefficient 3 of cubic polynomial in the laminar region";
    
    // Auxiliary variables for DetailedFriction model
    parameter SI.Length L=length/1000 "Length of pipe in SI units";
    parameter SI.Diameter D=if crossSectionType == CrossSectionTypes.Circular
         then diameter/1000 else (if crossSectionType == CrossSectionTypes.
        Rectangular then 4*(width*height/1.e6)/(2*width*height) else 4*area/
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
    parameter Real yd2=if from_dp then 0.5 + (2.51/Math.log(10))/(aux5*aux4)
         else 2 + 4*aux1/(aux2*aux3*(Re2)^0.9);
    
    // Constants: Cubic polynomial between lg(Re1) and lg(Re2)
    parameter Real diff_x=x2 - x1;
    parameter Real m=(y2 - y1)/diff_x;
    parameter Real c2=(3*m - 2*yd1 - yd2)/diff_x;
    parameter Real c3=(yd1 + yd2 - 2*m)/(diff_x*diff_x);
    parameter Real lambda2_1=64*Re1;
    constant Real pi=Modelica.Constants.pi;
    Real dx;
    Real aux7;
    Medium.Density d 
      "Dummy, or density used for detailed friction model (from upstream port)";
    Medium.DynamicViscosity eta 
      "Dummy, or viscosity used for detailed friction model (from upstream port)";
  equation 
    zero = 0;
    
    if frictionType <> FrictionTypes.DetailedFriction then
      // Assign dummy values for auxiliary variables
      Re = 0;
      dx = 0;
      lambda = 0;
      lambda2 = 0;
      aux7 = 0;
      d = 0;
      eta = 0;
    else
      lambda = noEvent(if Re < 64 then 1 else lambda2/(Re*Re));
      
      // Use d and eta from the upstream port
      d = if port_a.m_dot > 0 then medium_a.d else medium_b.d;
      eta = if port_a.m_dot > 0 then Medium.dynamicViscosity(medium_a) else 
        Medium.dynamicViscosity(medium_b);
    end if;
    
    if from_dp then
      // equations in the form m_dot = m_dot(dp)
      if frictionType == FrictionTypes.ConstantLaminar then
        m_dot = dp/k;
      elseif frictionType == FrictionTypes.ConstantTurbulent then
        m_dot = noEvent(if dp > delta then sqrt(dp) else (if dp < -delta then -
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
        m_dot = noEvent((pi*D/4)*eta*Re*(if dp >= 0 then +1 else -1));
      end if;
    else
      // equations in the form dp = dp(m_dot)
      if frictionType == FrictionTypes.ConstantLaminar then
        dp = k*m_dot;
      elseif frictionType == FrictionTypes.ConstantTurbulent then
        dp = k*noEvent(if m_dot > delta then m_dot*m_dot else (if m_dot < -
          delta then -m_dot*m_dot else (C1 + C3*m_dot*m_dot)*m_dot));
      else
        Re = noEvent((4/pi)*abs(m_dot)/(D*eta));
        dx = noEvent(if Re < Re1 or Re > Re2 then 0 else Math.log10(Re/Re1));
        lambda2 = noEvent(if Re <= Re1 then 64*Re else (if Re >= Re2 then 0.25*
          (Re/Math.log10(Delta/3.7 + 5.74/Re^0.9))^2 else 64*Re1*(Re/Re1)^(1 + 
          dx*(c2 + dx*c3))));
        aux7 = 0;
        dp = noEvent(k*lambda2*eta*eta/d*(if m_dot >= 0 then 1 else -1));
      end if;
    end if;
  end ShortPipe;
  
  model LongPipe 
    "(will be soon replaced by a much better version) Pipe discretized according to the finite volume method (currently there is just one volume for easier discussion)"
     
    
    import SI = Modelica.SIunits;
    import Modelica.SIunits.Conversions.*;
    extends Interfaces.PartialTwoPort;
    
    constant Real pi=Modelica.Constants.pi;
    parameter Integer n(
      min=1, 
      max=1) = 1 "Number of internal volumes (currently only one)";
    parameter SI.Diameter diameter "Pipe diameter";
    parameter SI.Length length "Pipe length";
    final parameter SI.Volume V=length*pi*diameter^2/4 "Pipe volume";
    
    parameter Modelica_Fluid.Examples.Types.FrictionTypes.Temp frictionType=
        Modelica_Fluid.Examples.Types.FrictionTypes.ConstantTurbulent 
      "Type of friction to determine pressure loss";
    parameter Medium.AbsolutePressure dp_nominal(min=1.e-10) = from_bar(1.0) 
      "|frictionType = ConstantLaminar or ConstantTurbulent| Nominal pressure drop";
    
    parameter Medium.MassFlowRate m_dot_nominal(min=1.e-10) = 1 
      "|frictionType = ConstantLaminar or ConstantTurbulent| Nominal mass flow rate at nominal pressure drop";
    parameter SI.Pressure p_small(min=1.e-10) = 1 
      "|Advanced|Only for frictionType = ConstantTurbulent| A small laminar region is introduced around p_small";
    
    Interfaces.JunctionVolume volume(
      V=V, 
      redeclare package Medium = Medium, 
      initType=initType, 
      init_p=init_p, 
      p_start=p_start, 
      d_start=d_start, 
      init_T=init_T, 
      T_start=T_start, 
      h_start=h_start, 
      X_start=X_start) annotation (extent=[-10, 10; 10, -10], rotation=90);
    ShortPipe shortPipe_a(
      redeclare package Medium = Medium, 
      frictionType=frictionType, 
      dp_nominal=dp_nominal/2, 
      m_dot_nominal=m_dot_nominal, 
      p_small=p_small, 
      init_p=init_p, 
      p_start=p_start, 
      d_start=d_start, 
      init_T=init_T, 
      T_start=T_start, 
      h_start=h_start, 
      X_start=X_start) annotation (extent=[-60, -10; -40, 10]);
    ShortPipe shortPipe_b(
      redeclare package Medium = Medium, 
      frictionType=frictionType, 
      dp_nominal=dp_nominal/2, 
      m_dot_nominal=m_dot_nominal, 
      p_small=p_small, 
      init_p=init_p, 
      p_start=p_start, 
      d_start=d_start, 
      init_T=init_T, 
      T_start=T_start, 
      h_start=h_start, 
      X_start=X_start) annotation (extent=[40, -10; 60, 10]);
    annotation (
      Diagram, 
      Icon(
        Rectangle(extent=[-100, 60; 100, -60], style(
            color=0, 
            gradient=2, 
            fillColor=8)), 
        Rectangle(extent=[-100, 34; 100, -36], style(
            color=69, 
            gradient=2, 
            fillColor=69)), 
        Text(
          extent=[-120, 130; 116, 64], 
          string="%name", 
          style(gradient=2, fillColor=69)), 
        Ellipse(extent=[-60, 14; -30, -14], style(fillColor=0)), 
        Ellipse(extent=[14, 14; 44, -14], style(fillColor=0)), 
        Text(
          extent=[-134, -64; 138, -94], 
          style(color=0), 
          string="%m_dot_nominal"), 
        Text(
          extent=[-132, -96; 140, -126], 
          style(color=0), 
          string="%dp_nominal")), 
      Documentation(info="<html>
<p>
This should be a model of a pipe with discretized
balance equations. In order to ease discussion,
this model is build up from two ShortPipes and one
JunctionVolume. A much improved LongPipe model is
near its completion. It will replace this one.
</p>
</html>"));
  equation 
    connect(shortPipe_a.port_a, port_a)
      annotation (points=[-61, 0; -110, 0], style(color=69));
    connect(shortPipe_a.port_b, volume.port)
      annotation (points=[-39, 0; 0, 0], style(color=69));
    connect(shortPipe_b.port_b, port_b)
      annotation (points=[61, 0; 110, 0], style(color=69));
    connect(shortPipe_b.port_a, volume.port)
      annotation (points=[39, 0; 0, 0], style(color=69));
  end LongPipe;
  
  model Tank "Tank with one bottom inlet/outlet" 
    import Modelica.SIunits.Conversions.*;
    
    replaceable package Medium = PackageMedium extends 
      Modelica_Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
    
    Interfaces.FluidPort_b port(redeclare package Medium = Medium)
      annotation (extent=[-10, -120; 10, -100], rotation=90);
    Medium.BaseProperties medium(
      preferedMediumStates=true, 
      final p_start=p_ambient, 
      final T_start=T_start, 
      final X_start=X_start);
    
    parameter Modelica.SIunits.Area area "Tank area";
    parameter Medium.AbsolutePressure p_ambient=101325 "Tank surface pressure";
    parameter Modelica.SIunits.Height level_start(min=0) 
      "|Initialization| Initial tank level";
    parameter Medium.Temperature T_start=from_degC(20) 
      "|Initialization| Initial tank temperature";
    parameter Medium.MassFraction X_start[Medium.nX](quantity=Medium.
          substanceNames) = zeros(Medium.nX) 
      "|Initialization (only for multi-substance flow)| Initial independent tank mass fractions m_i/m";
    constant Modelica.SIunits.Acceleration g=Modelica.Constants.G_EARTH;
    Modelica.SIunits.Height level(stateSelect=StateSelect.prefer, min=0) 
      "Level height of tank";
    Modelica.SIunits.Energy U "Internal energy of tank volume";
    Modelica.SIunits.Volume V(stateSelect=StateSelect.never) 
      "Actual tank volume";
    Real m(quantity=Medium.mediumName, unit="kg") "Mass of tank volume";
    Real mX[Medium.nX](quantity=Medium.substanceNames, each unit="kg") 
      "Component masses of the independent substances";
  initial equation 
    if not Medium.incompressible then
      mX = m*X_start;
    end if;
    level = level_start;
    medium.T = T_start;
    medium.X = X_start;
  equation 
    port.p = medium.p;
    
    /* Handle reverse and zero flow */
    port.H_dot = semiLinear(port.m_dot, port.h, medium.h);
    port.mX_dot = semiLinear(port.m_dot, port.X, medium.X);
    
    /*
  More precise equations (test later):
  Momentum balance
  (integrated momentum equation for frictionless fluid with density that is
   independent of the level, i.e., the unsteady Bernoulli equation for incompressible fluid)
  v_level = der(level);
  v = -port.m_dot/(rho*A_outlet);
  level*der(v_level) + (v^2 - v_level^2)/2 - g*level + (p - p_ambient)/rho = 0;
  Energy balance
  Potential energy: E_pot = integ(dm*g*s)
                          = g*integ(rho*A*s*ds)
                          = g*rho*A*z^2/2
  Kinetic energy  : E_kin = integ(dm*v^2/2)
                          = integ(rho*A*v^2/2*ds)
                          = rho*A*v^2/2*integ(ds)
                          = rho*A*v^2/2*z
                          = M*v^2/2
  E = U + M*g*z/2 + M*v_level^2/2
  der(E) = port.H_dot + port.m_dot*v^2/2 - p_ambient*area*der(level)
*/
    
    V = area*level;
    m = V*medium.d;
    mX = m*medium.X;
    U = m*medium.u;
    
    // Mass balance
    der(m) = port.m_dot;
    der(mX) = port.mX_dot;
    
    // Momentum balance
    medium.p = m*g/area + p_ambient;
    
    // Energy balance
    der(U) = port.H_dot - p_ambient*der(V);
    
    annotation (
      Icon(
        Rectangle(extent=[-100, 90; 100, 26], style(color=7, fillColor=7)), 
        Rectangle(extent=[-100, 26; 100, -100], style(
            color=69, 
            fillColor=69, 
            fillPattern=1)), 
        Line(points=[-100, 100; -100, -100; 100, -100; 100, 100], style(
            color=0, 
            fillColor=69, 
            fillPattern=1)), 
        Text(
          extent=[-112, 162; 122, 102], 
          string="%name", 
          style(fillColor=69, fillPattern=1)), 
        Text(
          extent=[-86, -38; 94, -78], 
          style(color=0), 
          string="%level_start"), 
        Text(
          extent=[-94, 78; 94, 38], 
          style(color=0), 
          string="%p_ambient"), 
        Text(
          extent=[-94, 14; 90, -2], 
          string="level_start =", 
          style(color=0))), 
      Documentation(info="<HTML>
<p>
This is a simplified model of a tank. The top part is open to the environment.
The tank is filled with a single or multiple-substance liquid.
The whole tank is assumed to have uniform temperature and mass fractions.
</p>
</HTML>"), 
      Diagram);
    
  end Tank;
end Components;
