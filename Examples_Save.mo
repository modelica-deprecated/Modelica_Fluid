package Examples_Save 
  "Example models to show how to use the Fluid package (contains the old examples that are not yet converted to new libraries; examples do not run for this reason" 
  import Modelica.Icons;
  import Modelica.Constants;
  import SI = Modelica.SIunits;
  
  extends Icons.Library;
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
        "|Advanced|| True, if m_flow is computed as function of pressure drop dp (otherwise, use inverse function to avoid nonlinear loop)";
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
     m_flow = f1(dp)
  or
     dp = f2(m_flow)
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
         rho*C*dp = m_flow + k*m_flow^2
     At low mass flow rates the "k*m_flow^2" term can be neglected
     and we have laminar flow:
         m_flow = rho*C*dp
     At high mass flow rates the linear "m_flow" term can be neglected
     and we have turbulent flow:
         rho*C*dp = k*m_flow^2, i.e.,
               dp = k/(rho*C) * m_flow^2
     on the other hand we have for turbulent flow:
         dp = 0.5*rho*zeta*v^2
            = 0.5*rho*zeta*(m_flow/(rho*A))^2
            = 0.5*zeta/(rho*A^2) * m_flow^2
     Comparision results in
         k = C*zeta/(2*A^2)
  */
      dp = port_a.p - port_b.p;
      if dp_given then
        port_a.m_flow = noEvent((if dp >= 0 then 1 else -1)*(-1 + sqrt(1 + 4*k*C*
          abs(dp)))/(2*k));
      else
        C*dp = port_a.m_flow + k*port_a.m_flow^2*noEvent(if port_a.m_flow >= 0 then 
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
      import Modelica_Fluid.Examples_Save.Types.FrictionTypes;
      import Modelica_Fluid.Examples_Save.Types.CrossSectionTypes;
      import FT = Modelica_Fluid.Examples_Save.Types.FrictionTypes;
      extends Modelica_Fluid.Interfaces.PartialTwoPortTransport;
      SI.Pressure dp "Pressure loss due to friction";
      Real zero=port_a.p - port_b.p - dp "momentum balance (may be modified)";
      
      parameter FT.Temp frictionType=FT.ConstantTurbulent 
        "Type of friction to determine pressure loss";
      parameter Medium.AbsolutePressure dp_nominal(min=1.e-10) = from_bar(1.0) 
        "|frictionType = ConstantLaminar or ConstantTurbulent| Nominal pressure drop";
      
      parameter Medium.MassFlowRate m_flow_nominal(min=1.e-10) = 1 
        "|frictionType = ConstantLaminar or ConstantTurbulent| Nominal mass flow rate at nominal pressure drop";
      parameter Types.Length_mm length=1000 
        "|frictionType = DetailedFriction| Length of pipe";
      parameter Types.Length_mm roughness=0 
        "|frictionType = DetailedFriction| Roughness of pipe";
      parameter Modelica_Fluid.Examples_Save.Types.CrossSectionTypes.Temp 
        crossSectionType=CrossSectionTypes.Circular 
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
        "|Advanced|| = true, use m_flow = f(dp) otherwise use dp = f(m_flow), i.e., inverse equation"
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
            string="%m_flow_nominal"),
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
          dp_nominal/m_flow_nominal else (if frictionType == FrictionTypes.
          ConstantTurbulent then dp_nominal/m_flow_nominal^2 else L/(2*D*D*D)) 
        "Pressure loss coefficient (dp = k*f(m_flow))";
      parameter Real delta=if from_dp then p_small else sqrt(dp_nominal/k);
      parameter Real C1=if from_dp then 0.5/sqrt(delta) - 3.0*C3*delta^2 else 0.5
          *delta "Coefficient 1 of cubic polynomial in the laminar region";
      parameter Real C3=if from_dp then -0.25/(sqrt(delta)*delta^2) else 0.5/
          delta "Coefficient 3 of cubic polynomial in the laminar region";
      
      // Auxiliary variables for DetailedFriction model
      parameter SI.Length L=length/1000 "Length of pipe in SI units";
      parameter SI.Diameter D=if crossSectionType == CrossSectionTypes.Circular then 
                diameter/1000 else (if crossSectionType == CrossSectionTypes.
          Rectangular then 4*(width*height/1.e6)/(2*width*height) else 4*area/
          perimeter) "Diameter of pipe in SI units";
      parameter SI.ReynoldsNumber Re1=(745*exp(if Delta <= 0.0065 then 1 else 
          0.0065/Delta))^(if from_dp then 0.97 else 1) 
        "Re leaving laminar curve";
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
        d = if port_a.m_flow > 0 then medium_a.d else medium_b.d;
        eta = if port_a.m_flow > 0 then Medium.dynamicViscosity(medium_a) else 
          Medium.dynamicViscosity(medium_b);
      end if;
      
      if from_dp then
        // equations in the form m_flow = m_flow(dp)
        if frictionType == FrictionTypes.ConstantLaminar then
          m_flow = dp/k;
        elseif frictionType == FrictionTypes.ConstantTurbulent then
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
        if frictionType == FrictionTypes.ConstantLaminar then
          dp = k*m_flow;
        elseif frictionType == FrictionTypes.ConstantTurbulent then
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
      
      parameter Modelica_Fluid.Examples_Save.Types.FrictionTypes.Temp 
        frictionType=Types.FrictionTypes.ConstantTurbulent 
        "Type of friction to determine pressure loss";
      parameter Medium.AbsolutePressure dp_nominal(min=1.e-10) = from_bar(1.0) 
        "|frictionType = ConstantLaminar or ConstantTurbulent| Nominal pressure drop";
      
      parameter Medium.MassFlowRate m_flow_nominal(min=1.e-10) = 1 
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
        m_flow_nominal=m_flow_nominal,
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
        m_flow_nominal=m_flow_nominal,
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
            string="%m_flow_nominal"),
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
      parameter Medium.AbsolutePressure p_ambient=101325 
        "Tank surface pressure";
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
      port.H_flow = semiLinear(port.m_flow, port.h, medium.h);
      port.mX_flow = semiLinear(port.m_flow, port.X, medium.X);
      
      /*
  More precise equations (test later):
  Momentum balance
  (integrated momentum equation for frictionless fluid with density that is
   independent of the level, i.e., the unsteady Bernoulli equation for incompressible fluid)
  v_level = der(level);
  v = -port.m_flow/(rho*A_outlet);
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
  der(E) = port.H_flow + port.m_flow*v^2/2 - p_ambient*area*der(level)
*/
      
      V = area*level;
      m = V*medium.d;
      mX = m*medium.X;
      U = m*medium.u;
      
      // Mass balance
      der(m) = port.m_flow;
      der(mX) = port.mX_flow;
      
      // Momentum balance
      medium.p = m*g/area + p_ambient;
      
      // Energy balance
      der(U) = port.H_flow - p_ambient*der(V);
      
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
  
  package DrumBoiler 
    "Drum boiler example, see Franke, Rode, Krueger: On-line Optimization of Drum Boiler Startup, 3rd International Modelica Conference, Linkoping, 2003" 
    
    model DrumBoilerSimulation "Simulate start-up of DrumBoiler" 
      extends Modelica.Icons.Example;
      Components.DrumBoiler drumBoiler(pressureSensor(
            measurementUnit=Modelica_Fluid.Types.PressureUnits.Bar)) 
                                       annotation (extent=[-20, -40; 40, 20]);
      Modelica.Blocks.Sources.TimeTable q_F_Tab(table=[0, 0; 3600, 400; 7210,
            400]) annotation (extent=[-80, 2; -60, 22]);
      Modelica.Blocks.Sources.TimeTable Y_Valve_Tab(table=[0, 1; 3600, 1; 7210,
             1]) annotation (extent=[-80, -40; -60, -20]);
      annotation (
        Diagram,
        experiment(StopTime=7200),
        Documentation(info="<HTML>
<p>
Apply a ramp to fuel input and hold outlet valve open.
Simulate for 7200 seconds.
</p>
</HTML>"));
    equation 
      connect(q_F_Tab.y,       drumBoiler.q_F) annotation (points=[-59,12; -40,12;
            -40,-16; -25.7,-16],        style(color=3));
      connect(Y_Valve_Tab.y,       drumBoiler.Y_Valve) annotation (points=[-59,
            -30; -44,-30; -44,-34; -25.7,-34],     style(
          color=3,
          fillColor=7,
          fillPattern=1));
    end DrumBoilerSimulation;
    
    package Components 
      
      model Evaporator 
        "Simple Evaporator with two states, see Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378" 
        
        import Modelica_Fluid.Interfaces.*;
        import Modelica.SIunits.Conversions.*;
        import SI = Modelica.SIunits;
        
        // property and interface declarations
        package Medium = WaterPhaseBoundaryIF97;
        Medium.BaseProperties medium_a(region=1, p=port_a.p) "Medium in port_a";
        Medium.BaseProperties medium_b(region=2, p=port_b.p) "Medium in port_b";
        FluidPort_a port_a(redeclare package Medium = Medium) 
          annotation (extent=[-120, -10; -100, 10]);
        FluidPort_b port_b(redeclare package Medium = Medium) 
          annotation (extent=[120, -10; 100, 10]);
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
          annotation (extent=[-10, -120; 10, -100]);
        Modelica.Blocks.Interfaces.RealOutput V(
                                             redeclare type SignalType = 
              SI.Volume) "liquid volume (level)" 
          annotation (extent=[30, 100; 50, 120], rotation=90);
        Modelica.Blocks.Interfaces.RealOutput sigma_D "Thermal stress in metal"
          annotation (extent=[100, 40; 120, 60]);
        annotation (
          Coordsys(grid=[1, 1], component=[20, 20]),
          Diagram,
          Icon(
            Rectangle(extent=[-100, 59; 100, -61], style(
                color=0,
                gradient=2,
                fillColor=8)),
            Rectangle(extent=[-100, 34; 100, -36], style(
                color=69,
                gradient=2,
                fillColor=69)),
            Ellipse(extent=[18, 0; 48, -29], style(pattern=0, fillColor=7)),
            Ellipse(extent=[-1, 29; 29, 0], style(pattern=0, fillColor=7)),
            Ellipse(extent=[48, 34; 78, 5], style(pattern=0, fillColor=7)),
            Ellipse(extent=[-31, 1; -1, -28], style(pattern=0, fillColor=7)),
            Ellipse(extent=[47, 14; 77, -15], style(pattern=0, fillColor=7)),
            Ellipse(extent=[-72, 25; -42, -4], style(pattern=0, fillColor=7)),
            Ellipse(extent=[71, 0; 101, -29], style(pattern=0, fillColor=7)),
            Ellipse(extent=[74, 14; 104, -15], style(pattern=0, fillColor=7)),
            Ellipse(extent=[71, 29; 101, 0], style(pattern=0, fillColor=7)),
            Text(
              extent=[-120, 117; 116, 51],
              string="%name",
              style(gradient=2, fillColor=69)),
            Line(points=[0, -60; 0, -100], style(color=42)),
            Line(points=[40, 99; 40, 60])));
        
        // public parameters
        parameter SI.Mass m_D=300e3 "mass of surrounding drum metal";
        parameter SI.SpecificHeatCapacity cp_D=500 
          "specific heat capacity of drum metal";
        parameter SI.Volume V_t=100 "total volume inside drum";
        parameter SI.Pressure p_start=from_bar(1) "initial pressure";
        parameter SI.Volume V_start=67 "initial liquid volume";
        
      protected 
        SI.Pressure p(start=p_start, stateSelect=StateSelect.prefer) 
          "pressure inside drum boiler";
        SI.Volume V_v "volume of vapour phase";
        SI.Volume V_l(start=V_start, stateSelect=StateSelect.prefer) 
          "volumes of liquid phase";
        SI.SpecificEnthalpy h_v=medium_b.h "specific enthalpy of vapour";
        SI.SpecificEnthalpy h_l=medium_a.h "specific enthalpy of liquid";
        SI.Density rho_v=medium_b.d "density in vapour phase";
        SI.Density rho_l=medium_a.d "density in liquid phase";
        SI.Mass m "total mass of drum boiler";
        SI.Energy U "internal energy";
        SI.Temperature T_D=heatPort.T "temperature of drum";
        SI.HeatFlowRate q_F=heatPort.Q_flow "heat flow rate from furnace";
        SI.SpecificEnthalpy h_W=port_a.h "feed water enthalpy";
        SI.SpecificEnthalpy h_S=medium_b.h "steam enthalpy";
        SI.MassFlowRate qm_W=port_a.m_flow "feed water mass flow rate";
        SI.MassFlowRate qm_S=port_b.m_flow "steam mass flow rate";
      equation 
        
        // balance equations  
        m = rho_v*V_v + rho_l*V_l + m_D;
        U = rho_v*V_v*h_v + rho_l*V_l*h_l - p*V_t + m_D*cp_D*T_D;
        der(m) = qm_W + qm_S;
        der(U) = q_F + qm_W*h_W + qm_S*h_S;
        T_D = medium_a.T;
        // ideal heat transfer between metal and water
        V_t = V_l + V_v;
        
        // pressure and specific total enthalpies at ports
        port_a.p = p;
        port_b.p = p;
        port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, h_v);
        port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, h_l);
        
        // thermal stress
        sigma_D           = 60*der(T_D);
        
        // liquid level 
        V           = V_l;
      end Evaporator;
      
      model Valve 
        "Simple controlled valve with linear pressure drop coefficient" 
        
        import SI = Modelica.SIunits;
        import Modelica.SIunits.Conversions.*;
        import Modelica_Fluid.*;
        extends Interfaces.PartialTwoPortTransport;
        SI.Pressure dp "Pressure loss due to friction";
        Real residue=port_a.p - port_b.p - dp 
          "momentum balance (may be modified)";
        
        parameter Real k=1e-5 "linear valve coefficient";
        
        Modelica.Blocks.Interfaces.RealInput Y "Valve position" 
          annotation (extent=[-10, -80; 10, -60], rotation=90);
        
        annotation (Icon(
            Text(
              extent=[-126, -76; 130, -110],
              style(color=0),
              string=""),
            Text(
              extent=[-120, 130; 116, 64],
              string="%name",
              style(gradient=2, fillColor=69)),
            Line(points=[-60, -50; -60, 50; 60, -50; 60, 50; -60, -50], style(
                  color=0, thickness=2)),
            Line(points=[-60, 0; -100, 0], style(color=69)),
            Line(points=[60, 0; 100, 0], style(color=69)),
            Line(points=[0, 0; 0, -72])));
      equation 
        residue = 0;
        port_a.m_flow =Y           *k*dp;
      end Valve;
      
      model MassFlowSource_h 
        "Ideal pump that produces a mass flow rate from a large reservoir defined by input signal" 
        
        import Modelica_Fluid.*;
        import Modelica_Fluid.Interfaces.*;
        import Modelica_Media.Interfaces.*;
        
        import SI = Modelica.SIunits;
        import Modelica.SIunits.Conversions.*;
        
        FluidPort_b port(redeclare model Medium = Medium) 
          annotation (extent=[100, -10; 120, 10], rotation=0);
        replaceable package Medium = PartialMedium extends PartialMedium 
          "Medium in the component"   annotation (choicesAllMatching=true);
        
        parameter SI.SpecificEnthalpy h_ambient=5e5 "Ambient enthalphy";
        parameter SI.Pressure p_start=from_bar(1.0) 
          "|Initialization|| Initial pressure";
        annotation (
          Coordsys(
            extent=[-100, -100; 100, 100],
            grid=[2, 2],
            component=[20, 20]),
          Icon(
            Rectangle(extent=[20, 60; 100, -60], style(
                color=0,
                gradient=2,
                fillColor=8)),
            Rectangle(extent=[38, 40; 100, -40], style(
                color=69,
                gradient=2,
                fillColor=69)),
            Ellipse(extent=[-100, 80; 60, -80], style(fillColor=7)),
            Polygon(points=[-60, 70; 60, 0; -60, -68; -60, 70], style(color=73,
                   fillColor=73)),
            Text(
              extent=[-54, 32; 16, -30],
              string="m'",
              style(color=41, fillColor=41)),
            Text(extent=[-142, 142; 156, 88], string="%name"),
            Text(
              extent=[-126, -86; 146, -116],
              style(color=0),
              string="%h_ambient")),
          Window(
            x=0.45,
            y=0.01,
            width=0.44,
            height=0.65),
          Diagram);
        Modelica.Blocks.Interfaces.RealInput m_flow(
                                                redeclare type SignalType = 
              SI.MassFlowRate) 
          "Mass flow rate from an infinite reservoir in to the port as signal" 
          annotation (extent=[-140, -20; -100, 20]);
      equation 
        port.m_flow = -noEvent(max(m_flow,           0));
        port.H_flow = semiLinear(port.m_flow, port.h, h_ambient);
      end MassFlowSource_h;
      
      model DrumBoiler 
        "Complete drum boiler model, including evaporator and supplementary components" 
        
        import Modelica.SIunits.Conversions.*;
        
        Evaporator evaporator annotation (extent=[-50, -20; -30, 0]);
        annotation (
          uses(Modelica_Fluid(version="0.72")),
          Diagram,
          Coordsys(
            extent=[-100, -100; 100, 100],
            grid=[1, 1],
            component=[20, 20]),
          Icon(
            Rectangle(extent=[-100, 100; 100, -100], style(fillColor=7)),
            Text(
              extent=[-151, 165; 138, 102],
              style(fillColor=7, fillPattern=1),
              string="%name"),
            Text(
              extent=[-79, 67; 67, 21],
              string="drum",
              style(
                color=0,
                fillColor=7,
                fillPattern=1)),
            Text(
              extent=[-90, -14; 88, -64],
              string="boiler",
              style(
                color=0,
                fillColor=7,
                fillPattern=1))));
        Modelica.Thermal.HeatTransfer.PrescribedHeatFlow furnace 
          annotation (extent=[-50, -50; -30, -30], rotation=90);
        Modelica.Blocks.Interfaces.RealInput q_F 
          annotation (extent=[-139, 0; -99, -40]);
        Modelica.Blocks.Interfaces.RealInput Y_Valve 
          annotation (extent=[-139, -100; -99, -60]);
        Valve valve(redeclare package Medium = Modelica_Media.Water.WaterIF97, k=
              1.5e-5) annotation (extent=[44, -20; 64, 0]);
        Modelica_Fluid.Sources.FixedAmbient sink(redeclare package Medium = 
              Modelica_Media.Water.WaterIF97) 
          annotation (extent=[80, -20; 100, 0], rotation=180);
        Modelica_Fluid.Sensors.MassFlowSensor massFlowSensor(redeclare package 
            Medium = Modelica_Media.Water.WaterIF97) 
          annotation (extent=[10, -20; 30, 0]);
        Modelica_Fluid.Sensors.Temperature temperature(redeclare package Medium
            =                Modelica_Media.Water.WaterIF97, measurementUnit="degC") 
          annotation (extent=[10, 60; 30, 80]);
        Modelica_Fluid.Sensors.PressureSensor pressureSensor(redeclare package 
            Medium = Modelica_Media.Water.WaterIF97) 
          annotation (extent=[10, 20; 30, 40]);
        Modelica.Blocks.Continuous.PI controller(T=120, k=10) 
          annotation (extent=[-60, 30; -80, 50]);
        MassFlowSource_h pump(redeclare package Medium = 
              Modelica_Media.Water.WaterIF97) 
          annotation (extent=[-80, -20; -60, 0]);
        Modelica.Blocks.Math.Feedback feedback 
          annotation (extent=[-26, 30; -46, 50]);
        Modelica.Blocks.Sources.Constant levelSetPoint(k=67) 
          annotation (extent=[-46, 60; -26, 80]);
      protected 
        Modelica.Blocks.Interfaces.RealOutput T_S 
          annotation (extent=[34, 66; 42, 74]);
        Modelica.Blocks.Interfaces.RealOutput p_S 
          annotation (extent=[34, 26; 42, 34]);
        Modelica.Blocks.Interfaces.RealOutput qm_S 
          annotation (extent=[34, -34; 42, -26], rotation=0);
        Modelica.Blocks.Interfaces.RealOutput sigma_D 
          annotation (extent=[-24, 6; -16, 14]);
        Modelica.Blocks.Interfaces.RealOutput V_l 
          annotation (extent=[-24, 24; -16, 32]);
      public 
        Modelica.Blocks.Math.Gain MW2W(k=1e6) 
          annotation (extent=[-70, -69; -50, -50]);
      equation 
        connect(furnace.port, evaporator.heatPort) 
          annotation (points=[-40, -30; -40, -21], style(color=42));
        connect(Y_Valve, valve.Y) 
          annotation (points=[-119, -80; 54, -80; 54, -17], style(color=3));
        connect(evaporator.port_b, temperature.port) annotation (points=[
              -29, -10; -2, -10; -2, 50; 20, 50; 20, 59], style(color=69));
        connect(evaporator.port_b, pressureSensor.port) annotation (points=[-29,
               -10; -2, -10; -2, 10; 20, 10; 20, 19], style(color=69));
        connect(evaporator.port_b, massFlowSensor.port_a) 
          annotation (points=[-29, -10; 9, -10], style(color=69));
        connect(massFlowSensor.port_b, valve.port_a) 
          annotation (points=[31, -10; 43, -10], style(color=69));
        connect(valve.port_b, sink.port) annotation (points=[65,-10; 72,-10; 72,
              -10; 79,-10],      style(color=69));
        connect(pump.port, evaporator.port_a) 
          annotation (points=[-59, -10; -51, -10], style(color=69));
        connect(controller.u,feedback.y) 
          annotation (points=[-58, 40; -45, 40], style(color=3));
        connect(feedback.u2,      evaporator.V) 
          annotation (points=[-36, 32; -36, 1], style(color=3));
        connect(levelSetPoint.y,feedback.u1)             annotation (points=[-25,
               70; -20, 70; -20, 40; -28, 40], style(color=3));
        connect(pressureSensor.signal, p_S) 
          annotation (points=[30, 30; 38, 30], style(color=3));
        connect(temperature.T, T_S) 
          annotation (points=[31,70; 38,70],   style(color=3));
        connect(massFlowSensor.signal, qm_S) 
          annotation (points=[20, -20; 20, -30; 38, -30], style(color=3));
        connect(evaporator.sigma_D, sigma_D) annotation (points=[-29, -5; -26,
              -5; -26, 10; -20, 10], style(color=3));
        connect(evaporator.V, V_l) 
          annotation (points=[-36, 1; -36, 28; -20, 28], style(color=3));
        connect(controller.y,       pump.m_flow) annotation (points=[-81, 40; -90,
               40; -90, -10; -82, -10], style(color=3));
        connect(q_F,MW2W.u)       annotation (points=[-119, -20; -90, -20; -90,
               -59.5; -72, -59.5], style(color=3));
        connect(MW2W.y,furnace.Q_flow)       annotation (points=[-49, -59.5; -40,
               -59.5; -40, -50], style(color=3));
      end DrumBoiler;
      
      package WaterPhaseBoundaryIF97 
        "Physical properties for water at phase boundary at boiling and dew curves" 
        
        extends Modelica_Media.Interfaces.PartialMedium(
          mediumName="WaterIF97",
          substanceNames=fill("", 0),
          incompressible=false,
          reducedX=true,
          MassFlowRate(quantity="MassFlowRate.WaterIF97"));
        
        redeclare model extends BaseProperties 
          
        annotation(structurallyIncomplete);
          parameter Integer region=0 "specify region 1 (liquid) or 2 (vapour)";
        equation 
          
          assert(region == 1 or region == 2,
            "WaterPhaseBoundaryIF97 medium model only valid for regions 1 and 2");
          T = Modelica_Media.Water.IF97.BaseIF97.Basic.tsat(p);
          if region == 1 then
            d = Modelica_Media.Water.IF97.BaseIF97.Regions.rhol_p(p);
            h = Modelica_Media.Water.IF97.BaseIF97.Regions.hl_p(p);
          else
            d = Modelica_Media.Water.IF97.BaseIF97.Regions.rhov_p(p);
            h = Modelica_Media.Water.IF97.BaseIF97.Regions.hv_p(p);
          end if;
          u = h - p/d;
          R = 287.0; // data.R // Modelica.Constants.R/data.MM;
        end BaseProperties;
      end WaterPhaseBoundaryIF97;
    end Components;
  end DrumBoiler;
  
  package Elementary 
    "Elementary examples to demonstrate various features of the fluid library" 
    
    extends Modelica.Icons.Library;
    model SimpleMixing 
      "This example shows the difference of a JunctionVolume and of a FixedComponentVolume" 
      
      Interfaces.JunctionVolume junctionVolume(
        V=1.e-4,
        T_start=from_degC(50.0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-10, 40; 10, 20], rotation=180);
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      Sources.FixedAmbient fixedAmbient1(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-80, 20; -60, 40]);
      Components.ShortPipe shortPipe1(
        m_flow_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, 20; -20, 40]);
      annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
      Components.ShortPipe shortPipe3(
        m_flow_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[20, 20; 40, 40]);
      Sources.FixedAmbient fixedAmbient3(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[80, 20; 60, 40]);
      Sources.FixedAmbient fixedAmbient2(
        p_ambient=from_bar(1.0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        T_ambient=from_degC(10)) annotation (extent=[-80, 60; -60, 80]);
      Components.ShortPipe shortPipe2(
        m_flow_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, 60; -20, 80]);
      Sources.FixedAmbient fixedAmbient4(
        p_ambient=from_bar(1.0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        T_ambient=from_degC(1)) annotation (extent=[-80, -80; -60, -60]);
      Components.ShortPipe shortPipe4(
        m_flow_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, -80; -20, -60]);
      Components.ShortPipe shortPipe6(
        m_flow_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[20, -80; 40, -60]);
      Sources.FixedAmbient fixedAmbient6(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[80, -80; 60, -60]);
      Sources.FixedAmbient fixedAmbient5(
        p_ambient=from_bar(1.0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        T_ambient=from_degC(10)) annotation (extent=[-80, -40; -60, -20]);
      Components.ShortPipe shortPipe5(
        m_flow_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, -40; -20, -20]);
    equation 
      connect(fixedAmbient1.port, shortPipe1.port_a) 
        annotation (points=[-59, 30; -41, 30], style(color=69));
      connect(shortPipe3.port_b, fixedAmbient3.port) 
        annotation (points=[41, 30; 59, 30], style(color=69));
      connect(shortPipe3.port_a, junctionVolume.port) 
        annotation (points=[19, 30; 0, 30], style(color=69));
      connect(fixedAmbient2.port, shortPipe2.port_a) 
        annotation (points=[-59, 70; -41, 70], style(color=69));
      connect(fixedAmbient4.port, shortPipe4.port_a) 
        annotation (points=[-59, -70; -41, -70], style(color=69));
      connect(shortPipe6.port_b, fixedAmbient6.port) 
        annotation (points=[41, -70; 59, -70], style(color=69));
      connect(fixedAmbient5.port, shortPipe5.port_a) 
        annotation (points=[-59, -30; -41, -30], style(color=69));
      connect(shortPipe5.port_b, shortPipe6.port_a) annotation (points=[-19, -30;
             0, -30; 0, -70; 19, -70], style(color=69));
      connect(shortPipe4.port_b, shortPipe6.port_a) 
        annotation (points=[-19, -70; 19, -70], style(color=69));
      connect(shortPipe2.port_b, junctionVolume.port) 
        annotation (points=[-19, 70; 0, 70; 0, 30], style(color=69));
      connect(shortPipe1.port_b, junctionVolume.port) 
        annotation (points=[-19, 30; 0, 30], style(color=69));
    end SimpleMixing;
    
  end Elementary;
  
  package Tanks "Examples with Tanks" 
    extends Icons.Library;
    
    model ThreeTanksOneLiquid 
      import Modelica.SIunits.Conversions.*;
      extends Modelica.Icons.Example;
      annotation (
        Diagram,
        Coordsys(grid=[1, 1], component=[20, 20]),
        experiment(StopTime=5),
        Documentation(info="<html>
<p>
This example demonstrates the mixing of a single substance flow
between three tanks with different temperature.
</p>
<p>
When comparing with Examples.Tanks.ThreeTanks, it is demonstrated
that it is easy in this case to switch between single and multiple
substance flow. The only differences between the examples ThreeTanks
and ThreeTanksOneLiquid are that the Medium is different and that
for the model ThreeTanksOneLiquid the default values are
used for the initial mass fractions.
</p>
</html>"));
      
      Components.Tank Tank1(
        area=1,
        T_start=from_degC(50),
        level_start=3,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-90, 20; -70, 40]);
      Components.Tank Tank2(
        area=1,
        T_start=from_degC(100),
        level_start=1,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-10, 20; 10, 40]);
      Components.ShortPipe shortPipe1(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-50, -30; -30, -10]);
      
      Components.Tank Tank3(
        area=1,
        T_start=from_degC(20),
        level_start=2,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[70, 20; 90, 40]);
      Components.ShortPipe shortPipe3(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[30, -30; 50, -10]);
      Components.ShortPipe shortPipe2(
        m_flow_nominal=1000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-10, -10; 10, 10], rotation=-90);
    equation 
      connect(Tank1.port, shortPipe1.port_a) 
        annotation (points=[-80, 19; -80, -20; -51, -20], style(color=69));
      connect(shortPipe3.port_b, Tank3.port) 
        annotation (points=[51, -20; 80, -20; 80, 19], style(color=69));
      connect(Tank2.port, shortPipe2.port_a) 
        annotation (points=[0, 19; 0, 11; -6.73533e-016, 11], style(color=69));
      connect(shortPipe1.port_b, shortPipe3.port_a) 
        annotation (points=[-29, -20; 29, -20], style(color=69));
      connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
            6.73533e-016, -11; 0, -11; 0, -20; 29, -20], style(color=69));
    end ThreeTanksOneLiquid;
    
    model ThreeTanksWithPortVolume 
      import Modelica.SIunits.Conversions.*;
      extends Modelica.Icons.Example;
      annotation (
        Diagram,
        Coordsys(grid=[1, 1], component=[20, 20]),
        experiment(StopTime=5),
        Documentation(info="<html>
<p>
This example demonstrates the mixing of a single substance flow
between three tanks with different temperatures. The difference to
example \"ThreeTanksOneLiquid\" is that the port where the three
pipes are connected together contains a volume now, in order that
the mixing of the pipe flows is modelled more realistically.
</p>
</html>"));
      
      Components.Tank Tank1(
        area=1,
        T_start=from_degC(50),
        level_start=3,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-90, 20; -70, 40]);
      Components.Tank Tank2(
        area=1,
        T_start=from_degC(100),
        level_start=1,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-10, 20; 10, 40]);
      Components.ShortPipe shortPipe1(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-50, -40; -30, -20]);
      
      Components.Tank Tank3(
        area=1,
        T_start=from_degC(20),
        level_start=2,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[70, 20; 90, 40]);
      Components.ShortPipe shortPipe3(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[29, -40; 49, -20]);
      Components.ShortPipe shortPipe2(
        m_flow_nominal=1000,
        dp_nominal=from_bar(0.1),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-10, -10; 10, 10], rotation=-90);
      Interfaces.JunctionVolume junctionVolume(
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        V=1.e-4,
        T_start=from_degC(50.0),
        initType=Modelica_Media.Interfaces.PartialMedium.Choices.Init.
            InitialStates) annotation (extent=[-10, -40; 10, -20]);
    equation 
      connect(Tank1.port, shortPipe1.port_a) 
        annotation (points=[-80, 19; -80, -30; -51, -30], style(color=69));
      connect(shortPipe3.port_b, Tank3.port) 
        annotation (points=[50, -30; 80, -30; 80, 19], style(color=69));
      connect(Tank2.port, shortPipe2.port_a) 
        annotation (points=[0, 19; 0, 11; -6.73533e-016, 11], style(color=69));
      connect(shortPipe2.port_b, junctionVolume.port) annotation (points=[
            6.73533e-016, -11; 0, -11; 0, -30], style(color=69));
      connect(shortPipe1.port_b, junctionVolume.port) 
        annotation (points=[-29, -30; 0, -30], style(color=69));
      connect(shortPipe3.port_a, junctionVolume.port) 
        annotation (points=[28, -30; 0, -30], style(color=69));
    end ThreeTanksWithPortVolume;
    
    model ThreeTanksIF97 
      import Modelica.SIunits.Conversions.*;
      extends Modelica.Icons.Example;
      annotation (
        Diagram,
        Coordsys(grid=[1, 1], component=[20, 20]),
        experiment(StopTime=5),
        Documentation(info="<html>
<p>
This example demonstrates the mixing of a single substance flow
between three tanks with different temperature.
</p>
<p>
When comparing with Examples.Tanks.ThreeTanks, it is demonstrated
that it is easy in this case to switch between single and multiple
substance flow. The only differences between the examples ThreeTanks
and ThreeTanksOneLiquid are that the Medium is different and that
for the model ThreeTanksOneLiquid the default values are
used for the initial mass fractions.
</p>
</html>"));
      
      Components.Tank Tank1(
        area=1,
        T_start=from_degC(50),
        level_start=3,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-90, 20; -70, 40]);
      Components.Tank Tank2(
        area=1,
        T_start=from_degC(100),
        level_start=1,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-10, 20; 10, 40]);
      Components.ShortPipe shortPipe1(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-50, -30; -30, -10]);
      
      Components.Tank Tank3(
        area=1,
        T_start=from_degC(20),
        level_start=2,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[70, 20; 90, 40]);
      Components.ShortPipe shortPipe3(
        m_flow_nominal=2000,
        dp_nominal=from_bar(0.1),
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[30, -30; 50, -10]);
      Components.ShortPipe shortPipe2(
        m_flow_nominal=1000,
        dp_nominal=from_bar(0.1),
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-10, -10; 10, 10], rotation=-90);
    equation 
      connect(Tank1.port, shortPipe1.port_a) 
        annotation (points=[-80, 19; -80, -20; -51, -20], style(color=69));
      connect(shortPipe3.port_b, Tank3.port) 
        annotation (points=[51, -20; 80, -20; 80, 19], style(color=69));
      connect(Tank2.port, shortPipe2.port_a) 
        annotation (points=[0, 19; 0, 11; -6.73533e-016, 11], style(color=69));
      connect(shortPipe1.port_b, shortPipe3.port_a) 
        annotation (points=[-29, -20; 29, -20], style(color=69));
      connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
            6.73533e-016, -11; 0, -11; 0, -20; 29, -20], style(color=69));
    end ThreeTanksIF97;
  end Tanks;
  
  package TestComponents 
    "Test components (this package will be removed for the final version)" 
    
    extends Icons.Library;
    
    model TestLongPipe "Test ShortPipe componet" 
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      Modelica_Fluid.Sources.FixedAmbient ambient(T_ambient=from_degC(15),
          redeclare package Medium = Modelica_Media.Air.SimpleAir) 
        annotation (extent=[60, 0; 40, 20]);
      annotation (Diagram, experiment(StopTime=4));
      Modelica_Fluid.Components.LongPipe LongPipe(
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.SimpleAir,
        dp_nominal=from_bar(0.01),
        diameter=0.05,
        length=1) annotation (extent=[0, 0; 20, 20]);
      
      Modelica_Fluid.Sources.PrescribedMassFlowRate_TX MassFlowSource1(T_ambient=from_degC(30),
          redeclare package Medium = Modelica_Media.Air.SimpleAir) 
        annotation (extent=[-40, 0; -20, 20]);
      Modelica.Blocks.Sources.Ramp ramp(
        duration=3,
        height=6,
        offset=-3)   annotation (extent=[-80, 0; -60, 20]);
    equation 
      connect(LongPipe.port_b, ambient.port) 
        annotation (points=[21, 10; 39, 10], style(color=69));
      connect(MassFlowSource1.port, LongPipe.port_a) 
        annotation (points=[-19, 10; -1, 10], style(color=69));
      connect(ramp.y,       MassFlowSource1.m_flow) 
        annotation (points=[-59, 10; -42, 10], style(color=3));
    end TestLongPipe;
    
    model TestCheckStateSelection 
      "Check whether for the choosen medium the expected states are selected" 
      
      import Modelica.SIunits.Conversions.*;
      
      Modelica_Fluid.Interfaces.JunctionVolume junctionVolume1(
        V=1.e-4,
        T_start=from_degC(50.0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[-10, 90; 10, 70], rotation=180);
      extends Modelica.Icons.Example;
      Modelica_Fluid.Sources.FixedAmbient fixedAmbient1a(
        p_ambient=from_bar(1.0),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        T_ambient=from_degC(10)) annotation (extent=[-80, 70; -60, 90]);
      Modelica_Fluid.Examples_Save.Components.ShortPipe shortPipe1a(
        m_flow_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, 70; -20, 90]);
      annotation (
        Diagram(
          Text(
            extent=[-100, 102; -52, 92],
            string="SimpleLiquidWater",
            style(color=58)),
          Text(
            extent=[-113, 55; -54, 47],
            style(color=58),
            string="SimpleAir"),
          Text(
            extent=[-113, 5; -54, -3],
            style(color=58),
            string="DetailedAir")),
        Documentation(info="<html>
<p>
This model is used to select a medium model and then translate.
The log should give information about the selected states.
They should be the ones choosen as StateSelect.prefer in the
corresponding medium model.
</p>
</html>
"),     Coordsys(grid=[1, 1], component=[20, 20]));
      
      Modelica_Fluid.Examples_Save.Components.ShortPipe shortPipe1b(
        m_flow_nominal=10,
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[20, 70; 40, 90]);
      Modelica_Fluid.Sources.FixedAmbient fixedAmbient1b(
        p_ambient=from_bar(0.95),
        T_ambient=from_degC(30),
        redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
        annotation (extent=[80, 70; 60, 90]);
      
      Modelica_Fluid.Interfaces.JunctionVolume junctionVolume2(
        redeclare package Medium = Modelica_Media.Air.SimpleAir,
        p_start=from_bar(0.96),
        T_start=from_degC(50.0),
        V=0.1) annotation (extent=[-11, 40; 9, 20], rotation=180);
      
      Modelica_Fluid.Sources.FixedAmbient fixedAmbient2a(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Air.SimpleAir) 
        annotation (extent=[-80, 20; -60, 40]);
      Modelica_Fluid.Examples_Save.Components.ShortPipe shortPipe2a(
        redeclare package Medium = Modelica_Media.Air.SimpleAir,
        dp_nominal=from_bar(0.01),
        m_flow_nominal=1,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, 20; -20, 40]);
      
      Modelica_Fluid.Examples_Save.Components.ShortPipe shortPipe2b(
        redeclare package Medium = Modelica_Media.Air.SimpleAir,
        dp_nominal=from_bar(0.01),
        m_flow_nominal=1,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[20, 20; 40, 40]);
      
      Modelica_Fluid.Sources.FixedAmbient fixedAmbient2b(
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Air.SimpleAir,
        T_ambient=from_degC(10)) annotation (extent=[80, 20; 60, 40]);
      
      Modelica_Fluid.Interfaces.JunctionVolume junctionVolume3(
        p_start=from_bar(0.96),
        T_start=from_degC(50.0),
        V=0.1,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-11, -10; 9, -30], rotation=180);
      Modelica_Fluid.Sources.FixedAmbient fixedAmbient3a(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-80, -30; -60, -10]);
      Modelica_Fluid.Examples_Save.Components.ShortPipe shortPipe3a(
        dp_nominal=from_bar(0.01),
        m_flow_nominal=1,
        redeclare package Medium = Modelica_Media.Air.DetailedAir,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[-40, -30; -20, -10]);
      Modelica_Fluid.Examples_Save.Components.ShortPipe shortPipe3b(
        dp_nominal=from_bar(0.01),
        m_flow_nominal=1,
        redeclare package Medium = Modelica_Media.Air.DetailedAir,
        frictionType=Types.FrictionTypes.ConstantLaminar) 
        annotation (extent=[20, -30; 40, -10]);
      Modelica_Fluid.Sources.FixedAmbient fixedAmbient3b(
        p_ambient=from_bar(0.95),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[80, -30; 60, -10]);
    equation 
      connect(fixedAmbient1a.port, shortPipe1a.port_a) 
        annotation (points=[-59, 80; -41, 80], style(color=69));
      connect(shortPipe1b.port_b, fixedAmbient1b.port) 
        annotation (points=[41, 80; 59, 80], style(color=69));
      connect(shortPipe1b.port_a, junctionVolume1.port) 
        annotation (points=[19, 80; 0, 80], style(color=69));
      connect(shortPipe1a.port_b, junctionVolume1.port) 
        annotation (points=[-19, 80; 0, 80], style(color=69));
      connect(fixedAmbient2a.port, shortPipe2a.port_a) 
        annotation (points=[-59, 30; -41, 30], style(color=69));
      connect(shortPipe2b.port_b, fixedAmbient2b.port) 
        annotation (points=[41, 30; 59, 30], style(color=69));
      connect(shortPipe2b.port_a, junctionVolume2.port) 
        annotation (points=[19, 30; -1, 30], style(color=69));
      connect(shortPipe2a.port_b, junctionVolume2.port) 
        annotation (points=[-19, 30; -1, 30], style(color=69));
      connect(fixedAmbient3a.port, shortPipe3a.port_a) 
        annotation (points=[-59, -20; -41, -20], style(color=69));
      connect(shortPipe3b.port_b, fixedAmbient3b.port) 
        annotation (points=[41, -20; 59, -20], style(color=69));
      connect(shortPipe3b.port_a, junctionVolume3.port) 
        annotation (points=[19, -20; -1, -20], style(color=69));
      connect(shortPipe3a.port_b, junctionVolume3.port) 
        annotation (points=[-19, -20; -1, -20], style(color=69));
    end TestCheckStateSelection;
    
    model TwoVolumesAir 
      "This example shows the difference of a JunctionVolume and of a FixedComponentVolume" 
      
      import SI = Modelica.SIunits;
      parameter SI.Volume V=0.05 "Size of volume";
      
      Interfaces.JunctionVolume junctionVolume1(
        T_start=from_degC(50.0),
        V=V/2,
        redeclare package Medium = Modelica_Media.Air.DetailedAir,
        initType=Modelica_Media.Interfaces.PartialMedium.Choices.Init.
            InitialStates) annotation (extent=[-30, 40; -10, 20], rotation=180);
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      Sources.FixedAmbient fixedAmbient1(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-100, 20; -80, 40]);
      Components.ShortPipe shortPipe1(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-60, 20; -40, 40]);
      annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
      Components.ShortPipe shortPipe3(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[40, 20; 60, 40]);
      Sources.FixedAmbient fixedAmbient3(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[100, 20; 80, 40]);
      Sources.FixedAmbient fixedAmbient2(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(20),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-100, 60; -80, 80]);
      Components.ShortPipe shortPipe2(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-60, 60; -40, 80]);
      Interfaces.JunctionVolume junctionVolume2(
        T_start=from_degC(50.0),
        V=V/2,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[0, 40; 20, 20], rotation=180);
      Interfaces.JunctionVolume junctionVolume3(
        T_start=from_degC(50.0),
        V=V,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-30, -50; -10, -70], rotation=180);
      Sources.FixedAmbient fixedAmbient4(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-100, -70; -80, -50]);
      Components.ShortPipe shortPipe4(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-60, -70; -40, -50]);
      Components.ShortPipe shortPipe5(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[0, -70; 20, -50]);
      Sources.FixedAmbient fixedAmbient5(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[60, -70; 40, -50]);
      Sources.FixedAmbient fixedAmbient6(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(20),
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-100, -30; -80, -10]);
      Components.ShortPipe shortPipe6(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Air.DetailedAir) 
        annotation (extent=[-60, -30; -40, -10]);
    equation 
      connect(fixedAmbient1.port, shortPipe1.port_a) 
        annotation (points=[-79, 30; -61, 30], style(color=69));
      connect(shortPipe3.port_b, fixedAmbient3.port) 
        annotation (points=[61, 30; 79, 30], style(color=69));
      connect(fixedAmbient2.port, shortPipe2.port_a) 
        annotation (points=[-79, 70; -61, 70], style(color=69));
      connect(shortPipe2.port_b, junctionVolume1.port) 
        annotation (points=[-39, 70; -20, 70; -20, 30], style(color=69));
      connect(shortPipe1.port_b, junctionVolume1.port) 
        annotation (points=[-39, 30; -20, 30], style(color=69));
      connect(junctionVolume2.port, shortPipe3.port_a) 
        annotation (points=[10, 30; 39, 30], style(color=69));
      connect(junctionVolume1.port, junctionVolume2.port) 
        annotation (points=[-20, 30; 10, 30], style(color=69));
      connect(fixedAmbient4.port, shortPipe4.port_a) 
        annotation (points=[-79, -60; -61, -60], style(color=69));
      connect(shortPipe5.port_b, fixedAmbient5.port) 
        annotation (points=[21, -60; 39, -60], style(color=69));
      connect(shortPipe5.port_a, junctionVolume3.port) 
        annotation (points=[-1, -60; -20, -60], style(color=69));
      connect(fixedAmbient6.port, shortPipe6.port_a) 
        annotation (points=[-79, -20; -61, -20], style(color=69));
      connect(shortPipe6.port_b, junctionVolume3.port) 
        annotation (points=[-39, -20; -20, -20; -20, -60], style(color=69));
      connect(shortPipe4.port_b, junctionVolume3.port) 
        annotation (points=[-39, -60; -20, -60], style(color=69));
    end TwoVolumesAir;
    
    model TwoVolumesDetailedWater 
      "This example shows the difference of a JunctionVolume and of a FixedComponentVolume" 
      
      import SI = Modelica.SIunits;
      parameter SI.Volume V=1.e-4 "Size of volume";
      
      Interfaces.JunctionVolume junctionVolume1(
        T_start=from_degC(50.0),
        V=V/2,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-30, 40; -10, 20], rotation=180);
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      Sources.FixedAmbient fixedAmbient1(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-100, 20; -80, 40]);
      Components.ShortPipe shortPipe1(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-60, 20; -40, 40]);
      annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
      Components.ShortPipe shortPipe3(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[40, 20; 60, 40]);
      Sources.FixedAmbient fixedAmbient3(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[100, 20; 80, 40]);
      Sources.FixedAmbient fixedAmbient2(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(20),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-100, 60; -80, 80]);
      Components.ShortPipe shortPipe2(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-60, 60; -40, 80]);
      Interfaces.JunctionVolume junctionVolume2(
        T_start=from_degC(50.0),
        V=V/2,
        redeclare package Medium = Modelica_Media.Water.IF97_ph,
        initType=Modelica_Media.Interfaces.PartialMedium.Choices.Init.
            InitialStates) annotation (extent=[0, 40; 20, 20], rotation=180);
      
      Interfaces.JunctionVolume junctionVolume3(
        T_start=from_degC(50.0),
        V=V,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-30, -50; -10, -70], rotation=180);
      Sources.FixedAmbient fixedAmbient4(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-100, -70; -80, -50]);
      Components.ShortPipe shortPipe4(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-60, -70; -40, -50]);
      Components.ShortPipe shortPipe5(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[0, -70; 20, -50]);
      Sources.FixedAmbient fixedAmbient5(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[60, -70; 40, -50]);
      Sources.FixedAmbient fixedAmbient6(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(20),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-100, -30; -80, -10]);
      Components.ShortPipe shortPipe6(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-60, -30; -40, -10]);
    equation 
      connect(fixedAmbient1.port, shortPipe1.port_a) 
        annotation (points=[-79, 30; -61, 30], style(color=69));
      connect(shortPipe3.port_b, fixedAmbient3.port) 
        annotation (points=[61, 30; 79, 30], style(color=69));
      connect(fixedAmbient2.port, shortPipe2.port_a) 
        annotation (points=[-79, 70; -61, 70], style(color=69));
      connect(shortPipe2.port_b, junctionVolume1.port) 
        annotation (points=[-39, 70; -20, 70; -20, 30], style(color=69));
      connect(shortPipe1.port_b, junctionVolume1.port) 
        annotation (points=[-39, 30; -20, 30], style(color=69));
      connect(junctionVolume2.port, shortPipe3.port_a) 
        annotation (points=[10, 30; 39, 30], style(color=69));
      connect(junctionVolume1.port, junctionVolume2.port) 
        annotation (points=[-20, 30; 10, 30], style(color=69));
      connect(fixedAmbient4.port, shortPipe4.port_a) 
        annotation (points=[-79, -60; -61, -60], style(color=69));
      connect(shortPipe5.port_b, fixedAmbient5.port) 
        annotation (points=[21, -60; 39, -60], style(color=69));
      connect(shortPipe5.port_a, junctionVolume3.port) 
        annotation (points=[-1, -60; -20, -60], style(color=69));
      connect(fixedAmbient6.port, shortPipe6.port_a) 
        annotation (points=[-79, -20; -61, -20], style(color=69));
      connect(shortPipe6.port_b, junctionVolume3.port) 
        annotation (points=[-39, -20; -20, -20; -20, -60], style(color=69));
      connect(shortPipe4.port_b, junctionVolume3.port) 
        annotation (points=[-39, -60; -20, -60], style(color=69));
    end TwoVolumesDetailedWater;
    
    model OneVolumeDetailedWater 
      "This example shows the difference of a JunctionVolume and of a FixedComponentVolume" 
      
      import SI = Modelica.SIunits;
      parameter SI.Volume V=1.e-4 "Size of volume";
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
      Interfaces.JunctionVolume junctionVolume3(
        T_start=from_degC(50.0),
        V=V,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-10, 0; 10, -20], rotation=180);
      Sources.FixedAmbient fixedAmbient4(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(10),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-80, -20; -60, 0]);
      Components.ShortPipe shortPipe4(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-40, -20; -20, 0]);
      Components.ShortPipe shortPipe5(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[20, -20; 40, 0]);
      Sources.FixedAmbient fixedAmbient5(
        T_ambient=from_degC(20),
        p_ambient=from_bar(0.95),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[80, -20; 60, 0]);
      Sources.FixedAmbient fixedAmbient6(
        p_ambient=from_bar(1.0),
        T_ambient=from_degC(20),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-80, 20; -60, 40]);
      Components.ShortPipe shortPipe6(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
        annotation (extent=[-40, 20; -20, 40]);
    equation 
      connect(fixedAmbient4.port, shortPipe4.port_a) 
        annotation (points=[-59, -10; -41, -10], style(color=69));
      connect(shortPipe5.port_b, fixedAmbient5.port) 
        annotation (points=[41, -10; 59, -10], style(color=69));
      connect(shortPipe5.port_a, junctionVolume3.port) 
        annotation (points=[19, -10; 0, -10], style(color=69));
      connect(fixedAmbient6.port, shortPipe6.port_a) 
        annotation (points=[-59, 30; -41, 30], style(color=69));
      connect(shortPipe6.port_b, junctionVolume3.port) 
        annotation (points=[-19, 30; 0, 30; 0, -10], style(color=69));
      connect(shortPipe4.port_b, junctionVolume3.port) 
        annotation (points=[-19, -10; 0, -10], style(color=69));
    end OneVolumeDetailedWater;
    
    model TwoVolumesApproximationWater 
      "This example shows the difference of a JunctionVolume and of a FixedComponentVolume" 
      
      import SI = Modelica.SIunits;
      parameter SI.Volume V=1.e-4 "Size of volume";
      
      Interfaces.JunctionVolume junctionVolume1(
        T_start=from_degC(50.0),
        V=V/2,
        p_start=from_bar(100.0),
        initType=Modelica_Media.Interfaces.PartialMedium.Choices.Init.
            InitialStates,
      redeclare package Medium = 
            Modelica_Media.Water.IF97LocalApproximation_ph) 
        annotation (extent=[-12,40; 8,20],     rotation=180);
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      Sources.FixedAmbient fixedAmbient1(
        p_ambient=from_bar(100.0),
        T_ambient=from_degC(110),
      redeclare package Medium = 
            Modelica_Media.Water.IF97LocalApproximation_ph) 
        annotation (extent=[-100, 20; -80, 40]);
      
      Components.ShortPipe shortPipe1(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        T_start=400.0,
      redeclare package Medium = 
            Modelica_Media.Water.IF97LocalApproximation_ph) 
        annotation (extent=[-60, 20; -40, 40]);
      
      annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
      Components.ShortPipe shortPipe3(
        m_flow_nominal=10,
        frictionType=Types.FrictionTypes.ConstantLaminar,
        T_start=400,
      redeclare package Medium = 
            Modelica_Media.Water.IF97LocalApproximation_ph) 
        annotation (extent=[40, 20; 60, 40]);
      
      Sources.FixedAmbient fixedAmbient3(
        p_ambient=from_bar(98),
        T_ambient=from_degC(120),
      redeclare package Medium = 
            Modelica_Media.Water.IF97LocalApproximation_ph) 
        annotation (extent=[100, 20; 80, 40]);
      
    equation 
      connect(fixedAmbient1.port, shortPipe1.port_a) 
        annotation (points=[-79, 30; -61, 30], style(color=69));
      connect(shortPipe3.port_b, fixedAmbient3.port) 
        annotation (points=[61, 30; 79, 30], style(color=69));
      connect(shortPipe1.port_b, junctionVolume1.port) 
        annotation (points=[-39,30; -2,30],    style(color=69));
      connect(junctionVolume1.port, shortPipe3.port_a) annotation (points=[-2,
            30; 39,30], style(color=69, rgbcolor={0,127,255}));
    end TwoVolumesApproximationWater;
    
  end TestComponents;
  
  package Types 
    import SI = Modelica.SIunits;
    annotation (preferedView="text");
    
    type Length_mm = Real (
        quantity="Length",
        unit="mm",
        min=0);
    
      // constant SI.MassFlowRate ResidualFlow=1E-10 "Used to take care of zero mass flow";
    
    package FrictionTypes 
      "Type, constants and menu choices to define the pressure loss equations due to friction, as temporary solution until enumerations are available" 
      
      annotation (preferedView="text");
      
      extends Modelica.Icons.Library;
      constant Integer ConstantLaminar=1;
      constant Integer ConstantTurbulent=2;
      constant Integer DetailedFriction=3;
      type Temp 
        "Temporary type of FrictionTypes with choices for menus (until enumerations are available)" 
        
        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar 
              "ConstantLaminar \"dp = k*m_flow\"",
            choice=Modelica_Fluid.Types.FrictionTypes.ConstantTurbulent 
              "ConstantTurbulent \"dp = k*m_flow^2\"",
            choice=Modelica_Fluid.Types.FrictionTypes.DetailedFriction 
              "DetailedFriction \"dp = f(Re,delta,rho,L,D,nu)\""));
      end Temp;
    end FrictionTypes;
    
    package CrossSectionTypes 
      "Type, constants and menu choices to define the geometric cross of pipes, as temporary solution until enumerations are available" 
      
      annotation (preferedView="text");
      
      extends Modelica.Icons.Library;
      constant Integer Circular=1;
      constant Integer Rectangular=2;
      constant Integer General=3;
      type Temp 
        "Temporary type of CrossSectionTypes with choices for menus (until enumerations are available)" 
        
        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica_Fluid.Types.CrossSectionTypes.Circular "Circular",
            choice=Modelica_Fluid.Types.CrossSectionTypes.Rectangular 
              "Rectangular",
            choice=Modelica_Fluid.Types.CrossSectionTypes.General "General"));
      end Temp;
    end CrossSectionTypes;
    
    package InitTypes 
      "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available" 
      
      annotation (preferedView="text");
      
      extends Modelica.Icons.Library;
      constant Integer SteadyState=1;
      constant Integer SteadyPressure=2;
      constant Integer InitialStates=3;
      constant Integer NoDefaultInit=4;
      type Temp 
        "Temporary type of InitializationTypes with choices for menus (until enumerations are available)" 
        
        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica_Fluid.Types.InitializationTypes.SteadyState 
              "SteadyState (initialize in steady state)",
            choice=Modelica_Fluid.Types.InitializationTypes.SteadyPressure 
              "SteadyPressure (initialize pressure in steady state)",
            choice=Modelica_Fluid.Types.InitializationTypes.InitialStates 
              "InitialStates (initialize medium states)",
            choice=Modelica_Fluid.Types.InitializationTypes.NoDefaultInit 
              "NoDefaultInit (no default initialization)"));
      end Temp;
    end InitTypes;
    
    package FixedDensityTypes 
      "Type, constants and menu choices to define fixed density definition, as temporary solution until enumerations are available" 
      
      annotation (preferedView="text");
      
      extends Modelica.Icons.Library;
      constant Integer NoFixedDensity=1;
      constant Integer FixedDensity=2;
      constant Integer FixedDensity_via_pT=3;
      type Temp 
        "Temporary type of FixedDensityTypes with choices for menus (until enumerations are available)" 
        
        extends Integer;
        annotation (Evaluate=true, choices(
            choice=Modelica_Fluid.Types.FixedDensityTypes.NoFixedDensity 
              "NoFixedDensity (no fixed density for compressible media)",
            choice=Modelica_Fluid.Types.FixedDensityTypes.FixedDensity 
              "FixedDensity (fixed density for compressible media)",
            choice=Modelica_Fluid.Types.FixedDensityTypes.FixedDensity_via_pT 
              "FixedDensity_via_pT (fixed density defined by p, T)"));
      end Temp;
    end FixedDensityTypes;
    
  end Types;
  
  package PipeWithShockWaves 
    model ShockTube 
      parameter Modelica.SIunits.Pressure p0=5E5;
      parameter Modelica.SIunits.Temperature T0=1200;
      Modelica_Fluid.Components.IsolatedPipe IsolatedPipe1(
        dynamicMomentumBalance=true,
        L=1,
        A_a=0.01,
        k=0,
        includeThermalConductance=false,
        shortPipe(p_start=cat(1, fill(p0, div(IsolatedPipe1.n, 2)), fill(1E5, div(
              IsolatedPipe1.n, 2))), T_start=cat(1, fill(T0, div(IsolatedPipe1.n, 2)), fill(
              300, div(IsolatedPipe1.n, 2)))),
        m_flow_nominal=1E-3,
        includeKineticTerm=true,
        redeclare package Medium = Modelica_Media.Air.DryAirNasa,
        n=100,
        dp_nominal=0.001,
        viscosityFactor2=0.9,
        viscosityFactor1=0.5) 
               annotation (extent=[-20,60; 0,80]);
      UserInteraction.Outputs.SpatialPlot SpatialPlot1(
        x=linspace(0, 1, IsolatedPipe1.n),
          y=IsolatedPipe1.shortPipe.medium.T,
          minY=300,
          maxY=1200)           annotation(extent=[-100,-100; 100,-20]);
        annotation(experiment(StopTime=0.0005, Tolerance=1e-006),
            experimentSetupOutput,
          Diagram);
      UserInteraction.Outputs.SpatialPlot SpatialPlot2(
        y=IsolatedPipe1.shortPipe.medium.p,
        x=linspace(0, 1, IsolatedPipe1.n),
          minY=1E5,
          maxY=5E5)            annotation(extent=[-100,-20; 100,60]);
    end ShockTube;
    
    model TestIsolatedPipeTemperature "Test ShortPipe componet" 
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      FiniteVolume.Sources.FixedAmbient ambient(
        T_ambient=300,
        p_ambient=1E5,
        port(h(start=10000)),
        redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
      annotation (extent=[40,60; 20,80]);
      
    annotation (
      Diagram,
      experiment(
            StopTime=10,
            NumberOfIntervals=5000,
            Tolerance=1e-006,
            fixedstepsize=1e-006,
            Algorithm=""),
      experimentSetupOutput,
      Commands(file="Simulate and plot pressure.mos", file=
         "Simulate and Plot Temperature.mos"));
      
      Modelica_Fluid.Components.IsolatedPipe IsolatedPipe1(
        dp_nominal=500,
        m_flow_nominal=1,
        A_a=0.01,
        L=0.025,
        n=20,
        k=0.001,
        includeThermalConductance=true,
        dynamicMomentumBalance=false,
        redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
               annotation (extent=[-20,60; 0,80]);
      
    //  Real T[IsolatedPipe1.n](start=cat(1,fill(310, div(IsolatedPipe1.n,2)),fill(300, IsolatedPipe1.n-div(IsolatedPipe1.n,2))),
    //    fixed=fill(true,IsolatedPipe1.n));
      
      FiniteVolume.Sources.FixedAmbient ambient1(
        p_ambient=1E5,
        port(h(start=10000)),
        T_ambient=310,
        redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
      annotation (extent=[-40,80; -60,60], rotation=180);
      UserInteraction.Outputs.SpatialPlot SpatialPlot1(
        minY=295,
        maxY=315,
        x=linspace(0, 1, IsolatedPipe1.n),
        y=IsolatedPipe1.shortPipe.medium.T) 
                annotation(extent=[-100,-100; 100,0]);
    equation 
      
      connect(IsolatedPipe1.port_b, ambient.port) 
      annotation (points=[1,70; 19,70],  style(color=69));
      
      connect(ambient1.port, IsolatedPipe1.port_a) 
                                             annotation(points=[-39,70; -21,70],
          style(color=69, rgbcolor={0,127,255}));
      
    end TestIsolatedPipeTemperature;
    
    model TestThreeIsolatedPipesPressure "Test ShortPipe componet" 
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      parameter Real pressurePulsHeight=1E4;
      parameter Real pressurePulsWidth=1E-3;
      parameter Real pressurePulsStart=0.1E-3;
      parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient=300;
      parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient1=310;
      
    annotation (
      Diagram,
      experiment(
            StopTime=0.001,
            NumberOfIntervals=5000,
            Tolerance=1e-006,
            fixedstepsize=1e-006,
            Algorithm="Dassl"),
      experimentSetupOutput,
      Commands(file="Simulate and plot pressure.mos", file=
         "Simulate and Plot Temperature.mos"));
      
      Modelica_Fluid.Components.IsolatedPipe IsolatedPipe1(
        m_flow_nominal=1,
        A_a=0.01,
        dp_nominal=50,
        dynamicMomentumBalance=true,
        nVolumes=100,
        L=0.5,
        redeclare package Medium = Modelica_Media.Air.DryAirNasa,
        includeKineticTerm=true,
        includeViscosity=true) 
               annotation (extent=[12,70; 32,90]);
      Modelica_Fluid.Sources.PrescribedAmbient_pT prescribedAmbient(
          redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
                      annotation (extent=[-22,70; -2,90]);
      Modelica.Blocks.Sources.Ramp Ramp1(
          duration=scalar({pressurePulsWidth/2}),
          height=scalar({pressurePulsHeight}),
          offset=1E5,
          startTime=scalar({pressurePulsStart})) 
                     annotation (extent=[-100,52; -80,72]);
      Modelica.Blocks.Sources.Ramp Ramp2(
          duration=scalar({pressurePulsWidth/2}),
          height=scalar({-pressurePulsHeight}),
          offset=0,
          startTime=scalar({pressurePulsStart + pressurePulsWidth/2})) 
                      annotation (extent=[-100,106; -80,86]);
      Modelica.Blocks.Math.Add Add1 
                                  annotation (extent=[-60,76; -40,96]);
      Modelica_Fluid.Sources.FixedAmbient ambient1(
        p_ambient=1E5,
        port(h(start=10000)),
        T_ambient=T_ambient,
          redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
      annotation (extent=[100,70; 80,90]);
      Modelica_Fluid.Components.IsolatedPipe IsolatedPipe2(
        m_flow_nominal=1,
        dp_nominal=50,
        A_a=0.01/2,
        dynamicMomentumBalance=true,
        nVolumes=100,
        L=0.5,
        redeclare package Medium = Modelica_Media.Air.DryAirNasa,
        includeKineticTerm=true,
        includeViscosity=true) 
               annotation (extent=[48,70; 68,90]);
      UserInteraction.Outputs.SpatialPlot SpatialPlot1(
        maxY=1.1e5,
        minY=0.9e5,
          x=linspace(0, 1, IsolatedPipe1.nVolumes),
        y=IsolatedPipe1.pipeSegment.medium.p) 
                            annotation(extent=[-100,-40; 0,20]);
      UserInteraction.Outputs.SpatialPlot SpatialPlot2(
        maxY=1.1e5,
        minY=0.9e5,
          x=linspace(0, 1, IsolatedPipe2.nVolumes),
        y=IsolatedPipe2.pipeSegment.medium.p) 
                            annotation(extent=[0,-40; 100,20]);
      Modelica_Fluid.Sources.FixedAmbient ambient2(
        p_ambient=1E5,
        port(h(start=10000)),
        T_ambient=T_ambient,
          redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
      annotation (extent=[100,30; 80,50]);
      Modelica_Fluid.Components.IsolatedPipe IsolatedPipe3(
        m_flow_nominal=1,
        dp_nominal=25,
        A_a=0.01/2,
        dynamicMomentumBalance=true,
        nVolumes=50,
        L=0.25,
        redeclare package Medium = Modelica_Media.Air.DryAirNasa,
        includeKineticTerm=true,
        includeViscosity=true) 
               annotation (extent=[48,30; 68,50]);
      UserInteraction.Outputs.SpatialPlot SpatialPlot3(
        maxY=1.1e5,
        minY=0.9e5,
          x=linspace(0, 1, IsolatedPipe3.nVolumes),
        y=IsolatedPipe3.pipeSegment.medium.p) 
                            annotation(extent=[6,-100; 54,-40]);
      Modelica.Blocks.Sources.Constant Constant1(k=T_ambient1) 
        annotation (extent=[-60,46; -40,66]);
    equation 
      
      connect(prescribedAmbient.port, IsolatedPipe1.port_a) 
      annotation (points=[-1,80; 11,80], style(color=69));
      connect(Ramp1.y,Add1.u2)           annotation (points=[-79,62; -72,62; -72,80;
            -62,80],          style(color=3));
      connect(Ramp2.y,Add1.u1) 
      annotation (points=[-79,96; -70,96; -70,92; -62,92],   style(color=3));
      connect(Add1.y, prescribedAmbient.p_in) 
      annotation (points=[-39,86; -26,86], style(color=3));
      connect(IsolatedPipe2.port_a, IsolatedPipe1.port_b) 
                                                annotation(points=[47,80; 33,80],
          style(color=69, rgbcolor={0,127,255}));
      connect(IsolatedPipe2.port_b, ambient1.port) 
                                             annotation(points=[69,80; 79,80],
          style(color=69, rgbcolor={0,127,255}));
      connect(IsolatedPipe3.port_a, IsolatedPipe1.port_b) 
                                                annotation(points=[47,40; 40,40;
            40,80; 33,80], style(color=69, rgbcolor={0,127,255}));
      connect(IsolatedPipe3.port_b,ambient2. port) 
                                             annotation(points=[69,40; 79,40],
          style(color=69, rgbcolor={0,127,255}));
      
      connect(Constant1.y, prescribedAmbient.T_ambient)       annotation (points=[
            -39,56; -32,56; -32,74; -24,74], style(color=3, rgbcolor={0,0,255}));
    end TestThreeIsolatedPipesPressure;
    
    model TestThreeIsolatedPipesTemperature "Test ShortPipe componet" 
      
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      
      parameter Integer n=20;
      parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient=300;
      parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient1=310;
      
      FiniteVolume.Sources.FixedAmbient ambient(
        T_ambient=T_ambient1,
        p_ambient=1E5,
        port(h(start=10000)),
        redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
      annotation (extent=[-74,70; -54,90]);
      
    annotation (
      Diagram,
      experiment(
            StopTime=100,
            NumberOfIntervals=5000,
            Tolerance=1e-006,
            fixedstepsize=1e-006,
            Algorithm=""),
      experimentSetupOutput,
      Commands(file="Simulate and plot pressure.mos", file=
         "Simulate and Plot Temperature.mos"));
      
      Modelica_Fluid.Components.IsolatedPipe IsolatedPipe1(
        m_flow_nominal=1,
        A_a=0.01,
        dp_nominal=50,
        L=0.05,
        n=20,
        includeThermalConductance=true,
        k=0.1,
        dynamicMomentumBalance=false,
        redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
               annotation (extent=[-40,70; -20,90]);
      FiniteVolume.Sources.FixedAmbient ambient1(
        p_ambient=1E5,
        port(h(start=10000)),
        T_ambient=T_ambient1,
        redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
      annotation (extent=[48,70; 28,90]);
      Modelica_Fluid.Components.IsolatedPipe IsolatedPipe2(
        m_flow_nominal=1,
        dp_nominal=50,
        A_a=0.01/2,
        L=0.05,
        n=20,
        includeThermalConductance=true,
        k=0.1,
        dynamicMomentumBalance=false,
        redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
               annotation (extent=[-4,70; 16,90]);
      FiniteVolume.Sources.FixedAmbient ambient2(
        p_ambient=1E5,
        port(h(start=10000)),
        T_ambient=T_ambient,
        redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
      annotation (extent=[48,30; 28,50]);
      Modelica_Fluid.Components.IsolatedPipe IsolatedPipe3(
        m_flow_nominal=1,
        dp_nominal=25,
        A_a=0.01/2,
        L=0.025,
        n=20,
        includeThermalConductance=true,
        k=0.1,
        dynamicMomentumBalance=false,
        redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
               annotation (extent=[-4,30; 16,50]);
      UserInteraction.Outputs.SpatialPlot SpatialPlot1(
        minY=295,
        maxY=315,
          maxX=IsolatedPipe1.L,
          y=IsolatedPipe1.shortPipe.medium.T,
          x=linspace(0, IsolatedPipe1.L, IsolatedPipe1.n)) 
                                       annotation(extent=[-100,-40; 0,20]);
      UserInteraction.Outputs.SpatialPlot SpatialPlot2(
        minY=295,
        maxY=315,
          maxX=IsolatedPipe2.L,
          x=linspace(0, IsolatedPipe2.L, IsolatedPipe2.n),
          y=IsolatedPipe2.shortPipe.medium.T) 
                                       annotation(extent=[0,-40; 100,20]);
      UserInteraction.Outputs.SpatialPlot SpatialPlot3(
        minY=295,
        maxY=315,
          maxX=IsolatedPipe3.L,
          x=linspace(0, IsolatedPipe3.L, IsolatedPipe3.n),
          y=IsolatedPipe3.shortPipe.medium.T) 
                                       annotation(extent=[6,-100; 54,-40]);
    equation 
      
      connect(IsolatedPipe2.port_a,IsolatedPipe1. port_b) 
                                                annotation(points=[-5,80; -19,80],
                 style(color=69, rgbcolor={0,127,255}));
      
      connect(IsolatedPipe2.port_b,ambient1. port) 
                                             annotation(points=[17,80; 27,80],
          style(color=69, rgbcolor={0,127,255}));
      
      connect(IsolatedPipe3.port_a,IsolatedPipe1. port_b) 
                                                annotation(points=[-5,40; -12,40;
              -12,80; -19,80],
                            style(color=69, rgbcolor={0,127,255}));
      
      connect(IsolatedPipe3.port_b,ambient2. port) 
                                             annotation(points=[17,40; 27,40],
          style(color=69, rgbcolor={0,127,255}));
      
      connect(ambient.port, IsolatedPipe1.port_a) 
                                            annotation(points=[-53,80; -41,80],
          style(
          color=69,
          rgbcolor={0,127,255},
          gradient=2,
          fillColor=76,
          rgbfillColor={170,170,255}));
      
    end TestThreeIsolatedPipesTemperature;
  end PipeWithShockWaves;
  
  package Examples_Francesco 
    
    model SourceP "Pressure source for water/steam flows" 
      extends Icons.SourceP;
      replaceable package Medium = Modelica_Media.Interfaces.PartialMedium 
        "Medium model" annotation(choicesAllMatching = true);
      parameter Modelica.SIunits.Pressure p0=1.01325e5 "Nominal pressure";
      parameter Real R=0 "Hydraulic resistance";
      parameter Modelica.SIunits.SpecificEnthalpy h=1e5 
        "Nominal specific enthalpy";
      Modelica.SIunits.Pressure p "Actual pressure";
      Modelica_Fluid.Interfaces.FluidPort_b flange(redeclare package Medium = 
            Medium)  annotation (extent=[80, -20; 120, 20]);
      Modelica.Blocks.Interfaces.RealInput in_p0(
                                              n=1) 
        annotation (extent=[-60, 72; -20, 112], rotation=-90);
      Modelica.Blocks.Interfaces.RealInput in_h(
                                             n=1) 
        annotation (extent=[20, 70; 60, 110], rotation=-90);
    equation 
      if R == 0 then
        flange.p = p;
      else
        flange.p = p + flange.m_flow*R;
      end if;
      if cardinality(in_p0)==0 then
        p = p0;
        in_p0.signal[1] = 0;
      else
        p = in_p0.signal[1];
      end if;
      if cardinality(in_h)==0 then
        flange.H_flow = semiLinear(flange.m_flow,flange.h,h);
        in_h.signal[1] = 0;
      else
        flange.H_flow = semiLinear(flange.m_flow,flange.h,in_h.signal[1]);
      end if;
      annotation (
        Diagram,
        Icon(Text(extent=[-106, 90; -52, 50], string="p0"), Text(extent=[66, 90;
                 98, 52], string="h")),
        Documentation(info="<HTML>
<p><b>Modelling options</b></p>
<p>If <tt>R</tt> is set to zero, the pressure source is ideal; otherwise, the outlet pressure decreases proportionally to the outgoing flowrate.</p>
<p>If the <tt>in_p0</tt> connector is wired, then the source pressure is given by the corresponding signal, otherwise it is fixed to <tt>p0</tt>.</p>
<p>If the <tt>in_h</tt> connector is wired, then the source pressure is given by the corresponding signal, otherwise it is fixed to <tt>h</tt>.</p>
</HTML>", revisions="<html>
<ul>
<li><i>18 Jun 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Removed <tt>p0_fix</tt> and <tt>hfix</tt>; the connection of external signals is now detected automatically.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"));
    end SourceP;
    
    model SinkP "Pressure sink for water/steam flows" 
      extends Icons.SourceP;
      replaceable package Medium = Modelica_Media.Interfaces.PartialMedium 
        "Medium model" annotation(choicesAllMatching = true);
      parameter Modelica.SIunits.Pressure p0=1.01325e5 "Nominal pressure";
      parameter Real R=0 "Hydraulic resistance";
      parameter Modelica.SIunits.SpecificEnthalpy h=1e5 
        "Nominal specific enthalpy";
      Modelica.SIunits.Pressure p;
      Modelica_Fluid.Interfaces.FluidPort_a flange(redeclare package Medium = 
            Medium)  annotation (extent=[-120, -20; -80, 20]);
      Modelica.Blocks.Interfaces.RealInput in_p0(
                                              n=1) 
        annotation (extent=[-60, 68; -20, 108], rotation=-90);
      Modelica.Blocks.Interfaces.RealInput in_h(
                                             n=1) 
        annotation (extent=[20, 68; 60, 108], rotation=-90);
    equation 
      if R == 0 then
        flange.p = p;
      else
        flange.p = p + flange.m_flow*R;
      end if;
      if cardinality(in_p0)==0 then
        p = p0;
        in_p0.signal[1] = 0;
      else
        p = in_p0.signal[1];
      end if;
      if cardinality(in_h)==0 then
        flange.H_flow = semiLinear(flange.m_flow,flange.h,h);
        in_h.signal[1] = 0;
      else
        flange.H_flow = semiLinear(flange.m_flow,flange.h,in_h.signal[1]);
      end if;
      annotation (
        Icon(Text(extent=[-106, 92; -56, 50], string="p0"), Text(extent=[54, 94;
                 112, 52], string="h")),
        Diagram,
        Documentation(info="<HTML>
<p><b>Modelling options</b></p>
<p>If <tt>R</tt> is set to zero, the pressure sink is ideal; otherwise, the inlet pressure increases proportionally to the incoming flowrate.</p>
<p>If the <tt>in_p0</tt> connector is wired, then the source pressure is given by the corresponding signal, otherwise it is fixed to <tt>p0</tt>.</p>
<p>If the <tt>in_h</tt> connector is wired, then the source pressure is given by the corresponding signal, otherwise it is fixed to <tt>h</tt>.</p>
</HTML>", revisions="<html>
<ul>
<li><i>18 Jun 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Removed <tt>p0_fix</tt> and <tt>hfix</tt>; the connection of external signals is now detected automatically.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"));
    end SinkP;
    
    partial model ValveBase "Base model for valves" 
      extends Icons.Valve;
      replaceable package Medium = Modelica_Media.Interfaces.PartialMedium 
        "Medium model" 
                     annotation(choicesAllMatching= true);
      Medium.BaseProperties fluid;
      parameter Integer CvData "(0: Av | 1: Kv | 2: Cv | 3: OpPoint)";
      parameter Modelica.SIunits.Area Avnom=0 "Av (metric) flow coefficient";
      parameter Real Kvnom(unit="m^3/h")=0 "Kv (metric) flow coefficient";
      parameter Real Cvnom(unit="USG/min")=0 "Cv (US) flow coefficient";
      parameter Modelica.SIunits.Pressure dpnom "Nominal pressure drop";
      parameter Modelica.SIunits.MassFlowRate wnom=0 "Nominal mass flowrate";
      parameter ThermoPowerMedia.LiquidDensity rhonom=0 "Nominal inlet density";
      parameter Boolean CheckValve=false "Reverse flow stopped";
      parameter Real b=0.01 "Regularisation factor";
      replaceable function FlowChar = linear "Flow characteristic";
      Modelica.SIunits.MassFlowRate w "Mass flowrate";
      Real Av;
      Medium.Density rho "Inlet density";
      Medium.Temperature Tin;
      Modelica_Fluid.Interfaces.FluidPort_a inlet(redeclare package Medium = 
            Medium) annotation (extent=[-120, -20; -80, 20]);
      Modelica_Fluid.Interfaces.FluidPort_b outlet(redeclare package Medium = 
            Medium)  annotation (extent=[80, -20; 120, 20]);
      Modelica.Blocks.Interfaces.RealInput theta(
                                              n=1) 
        annotation (extent=[-20, 60; 20, 100], rotation=-90);
      annotation (
        Icon(Text(extent=[-100, -40; 100, -80], string="%name")),
        Diagram,
        Documentation(info="<HTML>
<p>This is the base model for the <tt>ValveLiq</tt>, <tt>ValveLiqChoked</tt>, and <tt>ValveVap</tt> valve models. The model is based on the IEC 534 / ISA S.75 standards for valve sizing.
<p>The model optionally supports reverse flow conditions (assuming symmetrical behaviour) or check valve operation, and has been suitably modified to avoid numerical singularities at zero pressure drop. 
<p>The flow characteristic can be customised.
<p><b>Modelling options</b></p>
<p>The following options are available to specify the valve flow coefficient in fully open conditions:
<ul><li><tt>CvData = 0</tt>: the flow coefficient is given by the metric Av coefficient <tt>Avnom</tt> (m^2).
<li><tt>CvData = 1</tt>: the flow coefficient is given by the metric Kv coefficient <tt>Kvnom</tt> (m^3/h).
<li><tt>CvData = 2</tt>: the flow coefficient is given by the US Cv coefficient <tt>Cvnom</tt> (USG/min).
<li><tt>CvData = 3</tt>: the flow coefficient is given implicitly by the operating point (<tt>wnom</tt>,<tt>dpnom</tt>, <tt>rhonom<tt/>).
</ul>
<p>The nominal pressure drop <tt>dpnom</tt> must always be specified; to avoid numerical singularities, the flow characteristic is modified for pressure drops less than <tt>b*dpnom</tt> (the default value is 1% of the nominal pressure drop). Increase this parameter if numerical instabilities occur in valves with very low pressure drops.
<p>If <tt>CheckValve</tt> is true, then the flow is stopped when the outlet pressure is higher than the inlet pressure; otherwise, reverse flow takes place.
<p>The default flow characteristic <tt>FlowChar</tt> is linear; this can be replaced by any user-defined function (e.g. equal percentage, quick opening, etc.).
</HTML>", revisions="<html>
<ul>
<li><i>1 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Valve models restructured using inheritance. <br>
       Adapted to Modelica_Media.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"));
    // Default functions for valve characteristics
    function linear 
      input Real x;
      output Real y;
      annotation (derivative=der_linear);
    algorithm 
          y := x;
    end linear;
      
    function one 
      input Real x;
      output Real y;
      annotation (derivative=der_one);
    algorithm 
        y := 1;
    end one;
      
    function der_linear 
      input Real x;
      input Real der_x;
      output Real der_y;
    algorithm 
      der_y := der_x;
    end der_linear;
      
    function der_one 
      input Real x;
      input Real der_x;
      output Real der_y;
    algorithm 
      der_y := 0;
    end der_one;
      
    equation 
      if CvData == 0 then
        Av = Avnom;
      elseif CvData == 1 then
        Av = 2.7778e-5*Kvnom;
      elseif CvData == 2 then
        Av = 2.4027e-5*Cvnom;
      elseif CvData == 3 then
        Av = wnom/sqrt(rhonom*dpnom);
      end if;
      inlet.m_flow + outlet.m_flow = 0;
      w = inlet.m_flow;
      inlet.H_flow=semiLinear(inlet.m_flow,inlet.h,fluid.h);
      outlet.H_flow=semiLinear(outlet.m_flow,outlet.h,fluid.h);
      inlet.H_flow+outlet.H_flow=0;
      fluid.p=inlet.p;
      Tin=fluid.T;
      rho=fluid.d;
      
    end ValveBase;
    
    model ValveLiq "Valve for liquid water flow" 
      extends ValveBase;
      Real z "Normalized pressure drop";
      Real sqrtz;
      annotation (
        Icon(Text(extent=[-100, -40; 100, -80], string="%name")),
        Diagram,
        Documentation(info="<HTML>
<p>Liquid water valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluid. <p>
Extends the <tt>ValveBase</tt> model (see the corresponding documentation for common valve features).
</html>", revisions="<html>
<ul>
<li><i>1 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Valve model restructured using inheritance. <br>
       Adapted to Modelica_Media.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</HTML>"));
    equation 
      z = (inlet.p - outlet.p)/dpnom;
      if CheckValve then
        sqrtz = (if z >= 0 then z/sqrt(z + b) else 0);
      else
        sqrtz = noEvent(z/sqrt(abs(z) + b));
      end if;
      w = FlowChar(theta.signal[1])*Av*sqrt(rho*dpnom)*sqrtz;
    end ValveLiq;
    
    partial model PumpBase "Base model for centrifugal pumps" 
      extends Icons.Pump;
      import Modelica.SIunits.Conversions.NonSIunits.*;
      replaceable package Medium = Modelica_Media.Interfaces.PartialMedium 
        "Medium model" 
                     annotation(choicesAllMatching= true);
      Medium.BaseProperties fluid "Fluid propertie";
      replaceable package SatMedium = 
          Modelica_Media.Interfaces.PartialTwoPhaseMedium 
        "Saturated medium model (required only for NPSH computation)" 
                     annotation(choicesAllMatching= true);
      parameter Integer Np0(min=1) = 1 "Nominal number of pumps in parallel";
      parameter Modelica.SIunits.Pressure pin_start 
        "Inlet Pressure Start Value";
      parameter Modelica.SIunits.Pressure pout_start 
        "Outlet Pressure Start Value";
      parameter Modelica.SIunits.SpecificEnthalpy hstart=1e5 
        "Fluid Specific Enthalpy Start Value";
      parameter Medium.Density rho0=1000 "Nominal Liquid Density";
      parameter AngularVelocity_rpm n0=1500 "Nominal rotational speed";
      parameter Modelica.SIunits.Volume V "Pump Internal Volume";
      parameter Real etaMech(
        min=0,
        max=1) = 0.98 "Mechanical Efficiency";
      parameter Integer CharData(
        min=0,
        max=1) = 1 "0: (A,B,C,D,E,F) - 1: OpPoints";
      parameter Real A0=0;
      parameter Real B0=0;
      parameter Real C0=0;
      parameter Real D0=0;
      parameter Real E0=0;
      parameter Real F0=0;
      parameter Modelica.SIunits.Height head_nom[3] 
        "Pump head for three operating points";
      parameter Modelica.SIunits.VolumeFlowRate q_nom[3] 
        "Volume flow rate for three operating points (single pump)";
      parameter Modelica.SIunits.Power P_cons[3] 
        "Power consumption for three operating points (single pump)";
      parameter Boolean CheckValve=false "Reverse flow stopped";
      parameter Boolean ThermalCapacity=false 
        "Fluid heat capacity accounted for";
      parameter Boolean ComputeNPSHa=false 
        "Compute NPSH Available at the inlet";
      constant Modelica.SIunits.Acceleration g=Modelica.Constants.g_n;
      Medium.Density rho "Liquid density";
      Medium.Temperature Tin "Liquid inlet temperature";
      Modelica.SIunits.MassFlowRate w "Mass flowrate (single pump)";
      Modelica.SIunits.SpecificEnthalpy h(start=hstart) 
        "Fluid specific enthalpy";
      AngularVelocity_rpm n "Shaft r.p.m.";
      Integer Np(min=1) "Number of pumps in parallel";
      
      Modelica.SIunits.Power P "Power Consumption (single pump)";
      Modelica.SIunits.Power Ptot "Power Consumption (total)";
      constant Modelica.SIunits.Power P_eps=1e-8 
        "Small coefficient to avoid numerical singularities";
      Modelica.SIunits.Power Phyd "Hydraulic power (single pump)";
      Real eta "Global Efficiency";
      Modelica.SIunits.Length NPSHa "Net Positive Suction Head available";
      Modelica.SIunits.Pressure pv "Saturated liquid pressure";
      Boolean FlowOn(start=true);
      Real s "Auxiliary Variable";
      Real A;
      Real B;
      Real C;
      Real D;
      Real E;
      Real F;
      Modelica_Fluid.Interfaces.FluidPort_a infl(redeclare package Medium = 
            Medium, p(start=pin_start)) 
        annotation (extent=[-100, 2; -60, 42]);
      Modelica_Fluid.Interfaces.FluidPort_b outfl(redeclare package Medium = 
            Medium, p(start=pout_start)) 
        annotation (extent=[40,52; 80,92]);
      Modelica.Blocks.Interfaces.RealInput in_Np(
                                              n=1) "Number of  parallel pumps" 
        annotation (extent=[18, 70; 38, 90], rotation=-90);
    equation 
      if cardinality(in_Np)==0 then
        Np = Np0;
        in_Np.signal[1] = 0;
      else
        Np = in_Np.signal[1];
      end if;
      
      if CharData == 0 then
        A = A0;
        B = B0;
        C = C0;
        D = D0;
        E = E0;
        F = F0;
      else
        if CharData == 1 then
          head_nom[1]*g = -A*q_nom[1]^2 + B*q_nom[1] + C;
          head_nom[2]*g = -A*q_nom[2]^2 + B*q_nom[2] + C;
          head_nom[3]*g = -A*q_nom[3]^2 + B*q_nom[3] + C;
          P_cons[1] = D*(n0^2)*q_nom[1] - E*n0*(q_nom[1]^2) + F*(n0^2);
          P_cons[2] = D*(n0^2)*q_nom[2] - E*n0*(q_nom[2]^2) + F*(n0^2);
          P_cons[3] = D*(n0^2)*q_nom[3] - E*n0*(q_nom[3]^2) + F*(n0^2);
        end if;
      end if;
      
      infl.m_flow + outfl.m_flow = 0;
      w = infl.m_flow/Np;
      FlowOn = s > 0;
      
      if (FlowOn or (not CheckValve)) then
        w = s;
        (outfl.p - infl.p)/rho = -A*(w/rho)^2 + B*(n/n0)*w/rho + C*(n/n0)^2;
      else
        (outfl.p - infl.p)/rho = C*(n/n0)^2 - s*1e3;
        w = 0;
      end if;
      
      P = D*(n^2)*(w/rho) - E*n*((w/rho)^2) + F*(n^2);
      Ptot = P*Np;
      
      infl.H_flow=semiLinear(infl.m_flow,infl.h,fluid.h);
      outfl.H_flow=semiLinear(outfl.m_flow,outfl.h,fluid.h);
      fluid.p=infl.p;
      
      rho = fluid.d;
      Tin = fluid.T;
      
      h = fluid.h;
      
      if ThermalCapacity then
        (rho*V*der(h)) = (outfl.m_flow/Np)*outfl.h + (infl.m_flow/Np)*infl.h + Phyd;
      else
        0 = (outfl.m_flow/Np)*outfl.h + (infl.m_flow/Np)*infl.h + Phyd;
      end if;
      Phyd = P*etaMech;
      eta = ((outfl.p - infl.p)*w/rho)/(P + P_eps);
      
      if ComputeNPSHa then
        pv=SatMedium.saturationPressure(fluid.T);
        NPSHa=(infl.p-pv)/(rho*Modelica.Constants.g_n);
      else
        pv=0;
        NPSHa=0;
      end if;
      
      annotation (
        Icon,
        Diagram,
        Documentation(info="<HTML>
<p>This is the base model for the <tt>Pump</tt> and <tt>
PumpMech</tt> pump models.
<p>The model describes a centrifugal pump, or a group of <tt>Np</tt> identical pumps in parallel. The hydraulic characteristic (head vs. flowrate) is represented, as well as the pump power consumption.
<p>In order to avoid singularities in the computation of the outlet enthalpy at zero flowrate, the thermal capacity of the fluid inside the pump body can be taken into account.
<p>The model can either support reverse flow conditions or include a built-in check valve to avoid flow reversal.
<p><b>Modelling options</b></p>
<p>The following options are available to specify the pump characteristics:
<ul><li><tt>CharData = 0</tt>: the coefficients of the characteristics (<tt>A,B,C,D,E,F</tt>) are provided directly
<li><tt>CharData = 1</tt>: the characteristics are specified by providing a vector of three operating points (in terms of heads <tt>head[3]</tt>, volume flow rate <tt>q[3]</tt>, power consumption <tt>P_cons[3]</tt>, nominal fluid density <tt>rho0</tt>, and nominal rotational speed <tt>n0</tt>) for a single pump.
</ul>
<p>If the <tt>in_Np</tt> input connector is wired, it provides the number of pumps in parallel; otherwise,  <tt>Np0</tt> parallel pumps are assumed.</p>
<p>If <tt>ThermalCapacity</tt> is set to true, the heat capacity of the fluid inside the pump is taken into account: this is necessary to avoid singularities in the computation of the outlet enthalpy in case of zero flowrate. If zero flowrate conditions are always avoided, this effect can be neglected by setting <tt>ThermalCapacity</tt> to false, thus avoiding a fast state variable in the model.
<p>The <tt>CheckValve</tt> parameter determines whether the pump has a built-in check valve or not.
</HTML>", revisions="<html>
<ul>
<li><i>2 Aug 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       NPSHa computation added. Changed parameter names</li>
<li><i>5 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model restructured by using inheritance. Adapted to Modelica_Media.</li>
<li><i>15 Jan 2004</i>
    by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br>
       <tt>ThermalCapacity</tt> and <tt>CheckValve</tt> added.</li>
<li><i>15 Dec 2003</i>
    by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br>
       First release.</li>
</ul>
</html>"));
    end PumpBase;
    
    model Pump "Centrifugal pump with ideally controlled speed" 
      extends PumpBase;
      import Modelica.SIunits.Conversions.NonSIunits.*;
      parameter AngularVelocity_rpm n_const=n0 "Constant rotational speed";
      Modelica.Blocks.Interfaces.RealInput in_n(
                                             n=1) "RPM" 
        annotation (extent=[-36, 70; -16, 90], rotation=-90);
    equation 
      if cardinality(in_n)==0 then
        n = n_const;
        in_n.signal[1] = 0;
      else
        n = in_n.signal[1];
      end if;
      
      annotation (
        Icon(
          Text(extent=[-58,94; -30,74], string="n"),
          Text(extent=[-10,102; 18,82], string="Np")),
        Diagram,
        Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>Np</tt> pumps in parallel) with controlled speed, either fixed or provided by an external signal.
<p>The model extends <tt>PumpBase</tt>
<p>If the <tt>in_n</tt> input connector is wired, it provides rotational speed of the pumps (rpm); otherwise, a constant rotational speed equal to <tt>n_const</tt> (which can be different from <tt>n0</tt>) is assumed.</p>
</HTML>", revisions="<html>
<ul>
<li><i>5 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model restructured by using inheritance. Adapted to Modelica_Media.</li>
<li><i>15 Jan 2004</i>
    by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br>
       <tt>ThermalCapacity</tt> and <tt>CheckValve</tt> added.</li>
<li><i>15 Dec 2003</i>
    by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br>
       First release.</li>
</ul>
</html>"));
    end Pump;
    
    model PumpMech "Centrifugal pump with mechanical connector for the shaft" 
      extends PumpBase;
      Modelica.SIunits.Angle phi "Shaft angle";
      Modelica.SIunits.AngularVelocity omega "Shaft angular velocity";
      Modelica.Mechanics.Rotational.Interfaces.Flange_a MechPort 
        annotation (extent=[80,4; 110,32]);
    equation 
      phi = MechPort.phi;
      omega = der(phi);
      n = Modelica.SIunits.Conversions.to_rpm(omega);
      P = omega*MechPort.tau;
      annotation (
        Icon(
          Text(extent=[-10,104; 18,84], string="Np")),
        Diagram,
        Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>Np</tt> pumps in parallel) with a mechanical rotational connector for the shaft, to be used when the pump drive has to be modelled explicitly. In the case of <tt>Np</tt> pumps in parallel, the mechanical connector is relative to a single pump.
<p>The model extends <tt>PumpBase</tt>
 </HTML>", revisions="<html>
<ul>
<li><i>5 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model restructured by using inheritance. Adapted to Modelica_Media.</li>
<li><i>15 Jan 2004</i>
    by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br>
       <tt>ThermalCapacity</tt> and <tt>CheckValve</tt> added.</li>
<li><i>15 Dec 2003</i>
    by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br>
       First release.</li>
</ul>
</html>"));
    end PumpMech;
    
    model ValveLin "Valve for water/steam flows with linear pressure drop" 
      extends Icons.Valve;
      replaceable package Medium = Modelica_Media.Interfaces.PartialMedium 
        "Medium model" 
                   annotation(choicesAllMatching= true);
      parameter ThermoPowerMedia.HydraulicConductance Kv 
        "Nominal hydraulic conductance";
      Modelica.SIunits.MassFlowRate w "Mass flowrate";
      Modelica.SIunits.SpecificEnthalpy h "Fluid specific enthalpy";
      Modelica_Fluid.Interfaces.FluidPort_a inlet(redeclare package Medium=Medium) 
                    annotation (extent=[-120, -20; -80, 20]);
      Modelica_Fluid.Interfaces.FluidPort_b outlet(redeclare package Medium=Medium) 
                     annotation (extent=[80, -20; 120, 20]);
      Modelica.Blocks.Interfaces.RealInput cmd(
                                            n=1) 
        annotation (extent=[-20, 60; 20, 100], rotation=-90);
    equation 
      inlet.m_flow + outlet.m_flow = 0;
      w = inlet.m_flow;
      inlet.H_flow=semiLinear(inlet.m_flow,inlet.h,h);
      outlet.H_flow=semiLinear(outlet.m_flow,outlet.h,h);
      inlet.H_flow+outlet.H_flow=0;
      w = Kv*cmd.signal[1]*(inlet.p - outlet.p);
      annotation (
        Icon(Text(extent=[-100, -40; 100, -74], string="%name")),
        Diagram,
        Documentation(info="<HTML>
<p>This very simple model provides a pressure drop which is proportional to the flowrate and to the <tt>cmd</tt> signal, without computing any fluid property.</p>
</HTML>", revisions="<html>
<ul>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"));
    end ValveLin;
    annotation (uses(Modelica(version="2.1 Beta1")));
    
    model TestWaterPump "Test case for WaterPump" 
      annotation (
        Diagram,
        experiment(StopTime=10, Tolerance=1e-006),
        Documentation(info="<HTML>
<p>This model tests the <tt>Pump</tt> model with the check valve option active.
<p>The sink pressure is varied sinusoidally with a period of 10 s, so as to operate the pump in all the possible working conditions, including stopped flow.
<p>
Simulation Interval = [0...10] sec <br> 
Integration Algorithm = DASSL <br>
Algorithm Tolerance = 1e-6 
<p><b>Revision history:</b></p>
<ul>
<li><i>5 Feb 2004</i>
    by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco
Schiavo</a>:<br>
       First release.</li>
</ul>
</HTML>"));
      SourceP Source(p0=1e5, redeclare package Medium = 
            Modelica_Media.Water.WaterIF97) 
        annotation (extent=[-86, 24; -66, 44]);
      ValveLin ValveLin1(Kv=1e-5, redeclare package Medium = 
            Modelica_Media.Water.WaterIF97) 
        annotation (extent=[8, 22; 28, 42]);
      SinkP SinkP1(p0=3e5, redeclare package Medium = 
            Modelica_Media.Water.WaterIF97) 
        annotation (extent=[42, 22; 62, 42]);
      Pump Pump1(
        rho0=1000,
        CharData=1,
        pin_start=1e5,
        pout_start=4e5,
        hstart=1e5,
        ThermalCapacity=true,
        V=0.01,
        P_cons={800,1800,2000},
        CheckValve=true,
        head_nom={60,30,0},
        q_nom={0,0.001,0.0015},
      redeclare package Medium = Modelica_Media.Water.WaterIF97,
      redeclare package SatMedium = Modelica_Media.Water.WaterIF97,
        ComputeNPSHa=true)  annotation (extent=[-52, 26; -32, 46]);
      Modelica.Blocks.Sources.Constant Constant1 
        annotation (extent=[-72, 64; -52, 84]);
      Modelica.Blocks.Sources.Sine Sine1(
        amplitude=5e5,
        freqHz=0.1,
        offset=3e5,
        phase=0,
        startTime=0)   annotation (extent=[78, 66; 58, 86]);
    equation 
      connect(ValveLin1.outlet, SinkP1.flange) 
        annotation (points=[28, 32; 42, 32]);
      connect(Source.flange, Pump1.infl) 
        annotation (points=[-66, 34; -58, 34; -58, 38.2; -50, 38.2]);
      connect(Constant1.y,       ValveLin1.cmd) annotation (points=[-51, 74; -16,
             74; -16, 72; 18, 72; 18, 40], style(color=3));
      connect(Sine1.y,       SinkP1.in_p0) annotation (points=[57, 76; 52, 76;
            52, 40.8; 48, 40.8], style(color=3));
      connect(Pump1.outfl, ValveLin1.inlet) 
        annotation (points=[-36,43.2; -14,43.2; -14,32; 8,32]);
    end TestWaterPump;
    
    model TestValve "Test case for valves" 
      package Medium = Modelica_Media.Water.WaterIF97;
      SourceP SourceP1(p0=10e5,
      redeclare package Medium = Medium) 
        annotation (extent=[-100,30; -80,50]);
      SourceP SourceP2(p0=8e5,
      redeclare package Medium = Medium) 
        annotation (extent=[-100, -50; -80, -30]);
      SinkP SinkP1(p0=1e5,
      redeclare package Medium = Medium) 
        annotation (extent=[64, -4; 84, 16]);
      ValveLiq V1(
        rhonom=1000,
        CvData=3,
        dpnom=9e5,
        wnom=1.5,
      redeclare package Medium = Medium) 
                  annotation (extent=[-50, 58; -30, 78]);
      ValveLiq V2(
        dpnom=5e5,
        rhonom=1000,
        CvData=3,
        wnom=1.2,
      redeclare package Medium = Medium) 
                  annotation (extent=[-38, 26; -18, 46]);
      ValveLiq V3(
        dpnom=3e5,
        rhonom=1000,
        CvData=3,
        wnom=1.1,
      redeclare package Medium = Medium) 
                  annotation (extent=[-38, -38; -18, -18]);
      ValveLiq V4(
        dpnom=8e5,
        rhonom=1000,
        CvData=3,
        wnom=1.3,
      redeclare package Medium = Medium) 
                  annotation (extent=[-38, -78; -18, -58]);
      ValveLiq V5(
        dpnom=4e5,
        wnom=2,
        rhonom=1000,
        CvData=3,
      redeclare package Medium = Medium) 
                  annotation (extent=[28, -4; 48, 16]);
      annotation (
        Diagram,
        experiment(StopTime=4, Tolerance=1e-006),
        Documentation(info="<HTML>
<p>This model tests the <tt>ValveLiq</tt> model zero or reverse flow conditions.
<p>Simulate the model for 4 s. At t = 1 s the V5 valve closes in 1 s, the V2 and V3 valves close in 2 s and the V1 and V4 valves open in 2 s. The flow in valve V3 reverses between t = 1.83 and t = 1.93.
<p><b>Revision history:</b></p>
<ul>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco
Casella</a>:<br>
       First release.</li>
</ul>
</HTML>"));
      SinkP SinkP2(p0=1e5,
      redeclare package Medium = Medium) 
        annotation (extent=[-16, 58; 4, 78]);
      SinkP SinkP3(p0=1e5, redeclare package Medium = Medium) 
        annotation (extent=[-2, -78; 18, -58]);
      Modelica.Blocks.Sources.Ramp CloseLoad(
        duration=1,
        height=-0.99,
        offset=1,
        startTime=1)    annotation (extent=[8, 28; 28, 48]);
      Modelica.Blocks.Sources.Ramp OpenRelief(
        duration=2,
        height=1,
        offset=0,
        startTime=1) 
                    annotation (extent=[-92, 74; -72, 94]);
      Modelica.Blocks.Sources.Ramp CloseValves(
        duration=2,
        height=-1,
        offset=1,
        startTime=1) 
                    annotation (extent=[-96, -12; -76, 8]);
    equation 
      connect(V1.outlet, SinkP2.flange) annotation (points=[-30, 68; -16, 68]);
      connect(V5.outlet, SinkP1.flange) annotation (points=[48, 6; 64, 6]);
      connect(V4.outlet, SinkP3.flange) annotation (points=[-18, -68; -2, -68]);
      connect(CloseLoad.y,       V5.theta) 
        annotation (points=[29, 38; 38, 38; 38, 14], style(color=3));
      connect(OpenRelief.y,       V1.theta) 
        annotation (points=[-71, 84; -40, 84; -40, 76], style(color=3));
      connect(OpenRelief.y,       V4.theta) annotation (points=[-71, 84; -70, -6;
             -42, -6; -42, -50; -28, -50; -28, -60], style(color=3));
      connect(CloseValves.y,       V3.theta) 
        annotation (points=[-75, -2; -28, -2; -28, -20], style(color=3));
      connect(CloseValves.y,       V2.theta) annotation (points=[-75, -2; -42,
            -2; -42, 54; -28, 54; -28, 44], style(color=3));
      connect(SourceP1.flange, V1.inlet) annotation (points=[-80,40; -66,40; -66,68;
            -50,68], style(color=69, rgbcolor={0,127,255}));
      connect(SourceP1.flange, V2.inlet) annotation (points=[-80,40; -60,40; -60,36;
            -38,36], style(color=69, rgbcolor={0,127,255}));
      connect(V2.outlet, V5.inlet) annotation (points=[-18,36; 5,36; 5,6; 28,6],
          style(color=69, rgbcolor={0,127,255}));
      connect(V3.outlet, V5.inlet) annotation (points=[-18,-28; 6,-28; 6,6; 28,6],
          style(color=69, rgbcolor={0,127,255}));
      connect(SourceP2.flange, V4.inlet) annotation (points=[-80,-40; -60,-40; -60,
            -68; -38,-68], style(color=69, rgbcolor={0,127,255}));
      connect(SourceP2.flange, V3.inlet) annotation (points=[-80,-40; -58,-40; -58,
            -28; -38,-28], style(color=69, rgbcolor={0,127,255}));
    end TestValve;
    
    package Icons 
      partial model SourceP 
        annotation (Icon(
            Ellipse(extent=[-80, 80; 80, -80], style(color=0, fillPattern=1)),
            Text(
              extent=[-20, 34; 28, -26],
              string="P",
              style(color=7, fillPattern=1)),
            Text(extent=[-100, -78; 100, -106], string="%name")));
      end SourceP;
      
      partial model Valve 
        annotation (Icon(
            Line(points=[0, 40; 0, 0], style(
                color=0,
                thickness=2,
                fillPattern=1)),
            Polygon(points=[-80, 40; -80, -40; 0, 0; -80, 40], style(
                color=0,
                thickness=2,
                fillPattern=1)),
            Polygon(points=[80, 40; 0, 0; 80, -40; 80, 40], style(
                color=0,
                thickness=2,
                fillPattern=1)),
            Rectangle(extent=[-20, 60; 20, 40], style(
                color=0,
                fillColor=0,
                fillPattern=1))), Diagram);
      end Valve;
      
      partial model Pump 
        annotation (Icon(
            Polygon(points=[-40, -24; -60, -60; 60, -60; 40, -24; -40, -24],
                style(pattern=0, fillColor=74)),
            Ellipse(extent=[-60, 80; 60, -40], style(gradient=3)),
            Polygon(points=[-30, 52; -30, -8; 48, 20; -30, 52], style(
                pattern=0,
                gradient=2,
                fillColor=7)),
            Text(extent=[-100, -64; 100, -90], string="%name")));
      end Pump;
    end Icons;
  end Examples_Francesco;
end Examples_Save;
