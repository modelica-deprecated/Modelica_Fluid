package Utilities 
  "Utility models to construct fluid components (should not be used directly) " 
  model PortVolume 
    "Fixed volume associated with a port by the finite volume method (used to build up physical components; fulfills mass and energy balance)" 
    import Modelica_Fluid.Types.InitTypes.*;
    extends Interfaces.PartialInitializationParameters;
    
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Medium model" 
       annotation (choicesAllMatching=true);
    parameter SI.Volume V "Volume";
    
    Interfaces.FluidPort_a port(
      redeclare package Medium = Medium) "Fluid port" 
      annotation (extent=[-10, -10; 10, 10], rotation=0);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
      "Thermal port" 
      annotation (extent=[-20,88; 20,108]);
    
    Medium.BaseProperties medium(preferredMediumStates=true,
                p(start=p_start), T(start=T_start),
                h(start=h_start), Xi(start= X_start[1:Medium.nXi]));
    SI.Energy U "Internal energy of fluid";
    SI.Mass m "Mass of fluid";
    SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
    annotation (
     Icon(
        Ellipse(extent=[-100, 100; 100, -100], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=3,
            fillColor=68,
            rgbfillColor={170,213,255})),
        Text(extent=[-144, 178; 146, 116], string="%name"),
        Text(
          extent=[-130, -108; 144, -150],
          style(color=0),
          string="V=%V")), Documentation(info="<html>
<p>
This component models the <b>volume</b> of <b>fixed size</b> that is
associated with the <b>fluid port</b> to which it is connected.
This means that all medium properties inside the volume, are identical
to the port medium properties. In particular, the specific enthalpy
inside the volume (= medium.h) is always identical to the specific enthalpy
in the port (port.h = medium.h). Usually, this model is used when
discretizing a component according to the finite volume method into
volumes in internal ports that only store energy and mass and into
transport elements that just transport energy, mass and momentum
between the internal ports without storing these quantities during the
transport. This splitting is only possible under certain assumptions.
</p>
</html>"),
      Diagram);
  equation 
    // medium properties set by port values
      port.p = medium.p;
      port.h = medium.h;
      port.Xi = medium.Xi;
      thermalPort.T = medium.T;
    
    // Total quantities
       m    = V*medium.d "Total Mass";
       mXi = m*medium.Xi "Independent component masses";
       U    = m*medium.u "Internal energy";
    
    // Mass and energy balance
       der(m)    = port.m_flow "Total mass balance";
       der(mXi)  = port.mXi_flow "Independent component mass balance";
       der(U)    = port.H_flow + thermalPort.Q_flow "Energy balance";
    
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
  end PortVolume;
  extends Modelica.Icons.Library;
  
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
    SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
    
    parameter Modelica_Fluid.Types.FrictionTypes.Temp frictionType=Modelica_Fluid.Types.
        FrictionTypes.ConstantTurbulent 
      "Type of friction to determine pressure loss";
    parameter SI.AbsolutePressure dp_nominal(min=1.e-10)=
      Modelica.SIunits.Conversions.from_bar(1.0) " Nominal pressure drop" 
      annotation (Dialog(enable=frictionType==FT.ConstantLaminar or frictionType==FT.ConstantTurbulent, group=
            "frictionType = ConstantLaminar or ConstantTurbulent"));
    
    parameter SI.MassFlowRate m_flow_nominal(min=1.e-10) = 1 
      " Nominal mass flow rate at nominal pressure drop" annotation (Dialog(
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
    parameter SI.Length perimeter=0.1 
      " Wetted perimeter of cross sectional area" 
      annotation (Dialog(enable=frictionType==FT.DetailedFriction and crossSectionType==CT.General, group="frictionType = DetailedFriction"));
    parameter Boolean from_dp=true 
      " = true, use m_flow = f(dp) otherwise use dp = f(m_flow), i.e., inverse equation"
      annotation (Evaluate=true, Dialog(tab="Advanced"));
    parameter SI.Pressure p_small(min=1.e-10) = 1 
      " A small laminar region is introduced around p_small" annotation (Dialog(
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
</html>",   revisions="<html>
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
  
  function checkAmbient "Check whether ambient definition is correct" 
    
    extends Modelica.Icons.Function;
    input String mediumName;
    input String substanceNames[:] "Names of substances";
    input Boolean singleState;
    input Boolean define_p;
    input Real X_ambient[:];
    input String modelName = "??? ambient ???";
  protected 
    Integer nX = size(X_ambient,1);
    String X_str;
  algorithm 
    assert(not singleState or singleState and define_p, "
Wrong value of parameter define_p (= false) in model \""   + modelName + "\":
The selected medium \"" + mediumName + "\" has Medium.singleState=true.
Therefore, an ambient density cannot be defined and
define_p = true is required.
");
    
    for i in 1:nX loop
      assert(X_ambient[i] >= 0.0, "
Wrong ambient mass fractions in medium \""
  + mediumName + "\" in model \"" + modelName + "\":
The ambient value X_ambient(" + String(i) + ") = " + String(
        X_ambient[i]) + "
is negative. It must be positive.
");
    end for;
    
    if nX > 0 and abs(sum(X_ambient) - 1.0) > 1.e-10 then
       X_str :="";
       for i in 1:nX loop
          X_str :=X_str + "   X_ambient[" + String(i) + "] = " + String(X_ambient[
          i]) + " \"" + substanceNames[i] + "\"\n";
       end for;
       Modelica.Utilities.Streams.error(
          "The ambient mass fractions in medium \"" + mediumName + "\" in model \"" + modelName + "\"\n" +
          "do not sum up to 1. Instead, sum(X_ambient) = " + String(sum(X_ambient)) + ":\n"
          + X_str);
    end if;
  end checkAmbient;
  
  function regRoot 
    "Anti-symmetric square root approximation with finite derivative in the origin" 
    extends Modelica.Icons.Function;
    input Real x;
    input Real delta=0.01 
      "Range of significant deviation from sqrt(abs(x))*sgn(x)";
    output Real y;
    annotation(derivative(zeroDerivative=delta)=Utilities.regRoot_der,
      Documentation(info="<html>
This function approximates sqrt(abs(x))*sgn(x), such that the derivative is finite and smooth in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = regRoot(x)</td><td>y ~= sqrt(abs(x))*sgn(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = regRoot(x)</td><td>y ~= x/sqrt(delta)</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
<p>
With the default value of delta=0.01, the difference between sqrt(x) and regRoot(x) is 16% around x=0.01, 0.25% around x=0.1 and 0.0025% around x=1.
</p> 
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  algorithm 
    y := x/(x*x+delta*delta)^0.25;
    
  annotation (Documentation(info="<html>
This function approximates sqrt(x)*sign(x), such that the derivative is finite and smooth in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= sqrt(abs(x))*sign(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= x/delta</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
<p>
With the default value of delta=0.01, the difference between sqrt(x) and sqrtReg(x) is 0.5% around x=0.1 and 0.005% around x=1.
</p> 
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  end regRoot;
  
  function regRoot_der "Derivative of regRoot" 
    extends Modelica.Icons.Function;
    input Real x;
    input Real delta=0.01 "Range of significant deviation from sqrt(x)";
    input Real dx "Derivative of x";
    output Real dy;
  algorithm 
    dy := dx*0.5*(x*x+2*delta*delta)/((x*x+delta*delta)^1.25);
    
  annotation (Documentation(info="<html>
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  end regRoot_der;
  
  function regSquare 
    "Anti-symmetric square approximation with non-zero derivative in the origin" 
    extends Modelica.Icons.Function;
    input Real x;
    input Real delta=0.01 "Range of significant deviation from x^2*sgn(x)";
    output Real y;
    annotation(Documentation(info="<html>
This function approximates x^2*sgn(x), such that the derivative is non-zero in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = regSquare(x)</td><td>y ~= x^2*sgn(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = regSquare(x)</td><td>y ~= x*delta</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
<p>
With the default value of delta=0.01, the difference between x^2 and regSquare(x) is 41% around x=0.01, 0.4% around x=0.1 and 0.005% around x=1.
</p> 
</p> 
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  algorithm 
    y := x*sqrt(x*x+delta*delta);
    
  annotation (Documentation(info="<html>
This function approximates sqrt(x)*sign(x), such that the derivative is finite and smooth in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= sqrt(abs(x))*sign(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= x/delta</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
<p>
With the default value of delta=0.01, the difference between sqrt(x) and sqrtReg(x) is 0.5% around x=0.1 and 0.005% around x=1.
</p> 
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  end regSquare;
  
  function regPow 
    "Anti-symmetric power approximation with non-zero derivative in the origin" 
    extends Modelica.Icons.Function;
    input Real x;
    input Real a;
    input Real delta=0.01 "Range of significant deviation from x^a*sgn(x)";
    output Real y;
    annotation(Documentation(info="<html>
This function approximates abs(x)^a*sign(x), such that the derivative is positive, finite and smooth in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = regPow(x)</td><td>y ~= abs(x)^a*sgn(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = regPow(x)</td><td>y ~= x*delta^(a-1)</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  algorithm 
    y := x*(x*x+delta*delta)^((a-1)/2);
    
  annotation (Documentation(info="<html>
This function approximates sqrt(x)*sign(x), such that the derivative is finite and smooth in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= sqrt(abs(x))*sign(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= x/delta</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
<p>
With the default value of delta=0.01, the difference between sqrt(x) and sqrtReg(x) is 0.5% around x=0.1 and 0.005% around x=1.
</p> 
</html>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
  end regPow;
end Utilities;
