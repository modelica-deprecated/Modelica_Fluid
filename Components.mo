package Components "Basic components for fluid models" 
  model ShortPipe 
    "Simple pipe model with pressure loss (no storage of mass and energy in pipe)" 
    extends Modelica_Fluid.Interfaces.PartialTwoPortTransport;
    extends Modelica_Fluid.Utilities.PipeFriction;
    annotation (Icon(
        Rectangle(extent=[-100,60; 100,-60],   style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[-100,34; 100,-36],   style(
            color=69,
            gradient=2,
            fillColor=69)),
        Text(
          extent=[-150,140; 150,80],
          string="%name",
          style(gradient=2, fillColor=69))), Documentation(info="<html>
<p>
Model <b>ShortPipe</b> defines a simple pipe model 
with pressure loss due to friction. It is assumed that
no mass or energy is stored in the pipe. 
The details of the pipe friction model are described
<a href=\"Modelica://Modelica_Fluid.Utilities.PipeFriction\">here</a>.
</p>
</html>"));
  equation 
    dp = port_a.p - port_b.p;
    if frictionType == Modelica_Fluid.Types.FrictionTypes.DetailedFriction then
       d = if port_a.m_flow > 0 then medium_a.d else medium_b.d;
       eta = if port_a.m_flow > 0 then Medium.dynamicViscosity(medium_a) else 
                                      Medium.dynamicViscosity(medium_b);
    else
      // Assign dummy values for auxiliary variables
       d = 0;
       eta = 0;
    end if;
  end ShortPipe;
  extends Modelica.Icons.Library;
  
  annotation (preferedView="info",
              Documentation(info="<html>
<p>
This package will contain basic component models of the fluid library.
It is currently empty as all components are being evaluated together 
with examples first (see sub-package Examples). 
</p>
</html>"));
  
model IsolatedPipe 
    "Model of an isolated pipe consisting of n pipe segments/FiniteVolumes" 
  import SI = Modelica.SIunits;
    
  replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
                                                                       annotation (
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

  model Tank "Tank with one bottom inlet/outlet" 
    import Modelica.SIunits.Conversions.*;
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
    
    Interfaces.FluidPort_b port(redeclare package Medium = Medium) 
      annotation (extent=[-10, -120; 10, -100], rotation=90);
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      Xi(start=X_start[1:Medium.nXi]));
    
    parameter Modelica.SIunits.Area area "Tank area";
    parameter Medium.AbsolutePressure p_ambient=101325 "Tank surface pressure";
    parameter Modelica.SIunits.Height level_start(min=0) 
      "|Initialization| Initial tank level";
    parameter Medium.Temperature T_start=from_degC(20) 
      "|Initialization| Initial tank temperature";
    parameter Medium.MassFraction X_start[Medium.nX](quantity=Medium.
          substanceNames) = zeros(Medium.nX) 
      "|Initialization (only for multi-substance flow)| Initial independent tank mass fractions m_i/m";
    constant Modelica.SIunits.Acceleration g=Modelica.Constants.g_n;
    Modelica.SIunits.Height level(stateSelect=StateSelect.prefer, min=0) 
      "Level height of tank";
    Modelica.SIunits.Energy U "Internal energy of tank volume";
    Modelica.SIunits.Volume V(stateSelect=StateSelect.never) 
      "Actual tank volume";
    Real m(quantity=Medium.mediumName, unit="kg") "Mass of tank volume";
    Real mX[Medium.nX](quantity=Medium.substanceNames, each unit="kg") 
      "Component masses of the independent substances";
  initial equation 
    if not Medium.singleState then
      mX = m*X_start[1:Medium.nXi];
    end if;
    level = level_start;
    medium.T = T_start;
    medium.Xi = X_start[1:Medium.nXi];
  equation 
    port.p = medium.p;
    
    /* Handle reverse and zero flow */
    port.H_flow = semiLinear(port.m_flow, port.h, medium.h);
    port.mXi_flow = semiLinear(port.m_flow, port.Xi, medium.X);
    
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
    mX = m*medium.Xi;
    U = m*medium.u;
    
    // Mass balance
    der(m) = port.m_flow;
    der(mX) = port.mXi_flow;
    
    // Momentum balance
    medium.p = m*g/area + p_ambient;
    
    // Energy balance
    der(U) = port.H_flow - 0.5*(p_ambient+medium.p)*der(V);
    
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
