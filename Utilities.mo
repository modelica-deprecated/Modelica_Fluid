package Utilities 
  "Utility models to construct fluid component (should not be used directly) " 
  extends Modelica.Icons.Library;
    model FiniteVolume 
    "One dimensional volume according to the finite volume method with 1 mass, 1 energy and 2 momentum balances" 
    
      import SI = Modelica.SIunits;
      import Modelica.Math;
    
      replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component"    annotation (
          choicesAllMatching =                                                                            true);
    
      Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare model Medium = Medium) 
        annotation(extent=[-120, -10; -100, 10]);
      Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare model Medium = Medium) 
        annotation(extent=[120, -10; 100, 10]);
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
        annotation(extent=[-10, 60; 10, 80], rotation=-90);
      Medium.BaseProperties medium(preferredMediumStates=true) 
      "Medium properties in the middle of the finite volume";
      SI.Mass M "Total mass in volume";
      SI.Mass[Medium.nX_i] MXi "Independent component masses";
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
    
      input Real dp_a 
      "Pressure loss due to pipe friction between port_a and middle of pipe";
      input Real dp_b 
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
        M   = medium.d*A_m*dx;
        MXi = M*medium.X_i;
        U   = M*medium.u;
    
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
      port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.X_i, medium.X_i);
      port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.X_i, medium.X_i);
    
      // Heat port has the medium temperature
      heatPort.T = medium.T;
    
    end FiniteVolume;
  
model PipeSegment 
    "One segment of a pipe with 1 mass, 1 energy, 2 momementum balances and pipe friction" 
    
  import SI = Modelica.SIunits;
  extends Modelica_Fluid.Utilities.FiniteVolume(medium(
             p(start=p_start), d(start=d_start),
             T(start=T_start), h(start=h_start), X_i(start=X_start[1:Medium.nX_i])));
  extends Modelica_Fluid.Interfaces.PartialMenuInitialization;
    
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
  if initType == Modelica_Fluid.Types.InitTypes.InitialStates then
    if not Medium.singleState then
      if init_p then
         medium.p = p_start;
      else
         medium.d = d_start;
      end if;
    end if;
      
    if init_T then
       medium.T = T_start;
    else
       medium.h = h_start;
    end if;
      
    medium.X_i = X_start[1:Medium.nX_i];
      
  elseif initType == Modelica_Fluid.Types.InitTypes.SteadyState then
    if not Medium.singleState then
      der(medium.d) = 0;
    end if;
    der(medium.u) = 0;
    der(medium.X_i) = zeros(Medium.nX_i);
      
  elseif initType == Modelica_Fluid.Types.InitTypes.SteadyMass then
    if not Medium.singleState then
      der(medium.d) = 0;
    end if;
      
    if init_T then
       medium.T = T_start;
    else
       medium.h = h_start;
    end if;
      
    der(medium.X_i) = zeros(Medium.nX_i);
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
  dp_a/dp_nominal = if linearPressureDrop then m_flow_a/m_flow_nominal else abs(
    m_flow_a)*m_flow_a/m_flow_nominal^2;
  dp_b/dp_nominal = if linearPressureDrop then -m_flow_b/m_flow_nominal else abs(
    -m_flow_b)*(-m_flow_b)/m_flow_nominal^2;
    
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
    SI.Pressure dp "Pressure loss due to pipe friction";
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
    input Boolean singleState;
    input Boolean define_p;
    input Real X_ambient[:];
  protected 
    Integer nX = size(X_ambient,1);
  algorithm 
    assert(not singleState or singleState and define_p, "
Wrong value of parameter define_p (= false) in ambient source component:
The selected medium \"" + mediumName + "\" has Medium.singleState=true.
Therefore, an ambient density cannot be defined and
define_p = true is required.
");
    
    for i in 1:nX loop
      assert(X_ambient[i] >= 0.0, "
Wrong ambient mass fractions in medium \"" + mediumName + "\":
The ambient value X_ambient(" + String(i) + ") = " + String(
        X_ambient[i]) + "
is negative. It must be positive.
");
    end for;
    
    assert(nX==0 or (nX>0 and abs(sum(X_ambient) - 1.0) < 1.e-10), "
The ambient mass fractions in medium \""   + mediumName + "\"
do not sum up to 1. Instead, sum(X_ambient) = "
                                         + String(sum(X_ambient)) + ".
");
    
  end checkAmbient;
end Utilities;
