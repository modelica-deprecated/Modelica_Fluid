package Components "Basic components for fluid models" 
  model PortVolume 
    "Fixed volume associated with a port by the finite volume method (used to build up physical components; fulfills mass and energy balance)" 
    import Modelica_Fluid.Types.InitTypes.*;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Medium model" 
       annotation (choicesAllMatching=true);
    parameter SI.Volume V=1e-3 "Volume";
    parameter Types.InitTypes.Temp initOption = NoInit "Initialization option" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
      "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_start = Medium.reference_p 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.Temperature T_start=
      if use_T_start then 293.15 else Medium.T_phX(p_start,h_start,X_start) 
      "Start value of temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.h_pTX(p_start, T_start, X_start) else 1e4 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.reference_X 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    
    Modelica_Fluid.Interfaces.FluidPort_a port(
      redeclare package Medium = Medium) "Fluid port" 
      annotation (extent=[-10, -10; 10, 10], rotation=0);
    Interfaces.HeatPort_a thermalPort "Thermal port" 
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
      medium.p = p_start;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initOption == SteadyState then
      der(medium.p) = 0;
      der(medium.h) = 0;
      der(medium.Xi) = zeros(Medium.nXi);
    elseif initOption == SteadyStateHydraulic then
      der(medium.p) = 0;
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
  
  model MixingVolume 
    "Mixing volume with inlet and outlet ports (flow reversal is allowed)" 
    import Modelica_Fluid.Types.InitTypes.*;
    import Modelica.Constants.*;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Medium model" 
       annotation (choicesAllMatching=true);
    parameter SI.Volume V=1e-3 "Volume";
    parameter Boolean allowFlowReversal = true 
      "Flow reversal at the ports is allowed by the equations";
    parameter Types.InitTypes.Temp initOption = NoInit "Initialization option" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
      "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_start = Medium.reference_p 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.Temperature T_start=
      if use_T_start then 293.15 else Medium.T_phX(p_start,h_start,X_start) 
      "Start value of temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.h_pTX(p_start, T_start, X_start) else 1e4 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.reference_X 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    Interfaces.FluidPort_a port_a(redeclare package Medium = Medium,
                                  m_flow(min=if allowFlowReversal then -inf else 0)) 
      "Fluid inlet port" annotation (extent=[-112,-10; -92,10]);
    Interfaces.FluidPort_b port_b(redeclare package Medium = Medium,
                                  m_flow(max=if allowFlowReversal then +inf else 0)) 
      "Fluid outlet port" annotation (extent=[90,-10; 110,10]);
    Interfaces.HeatPort_a thermalPort "Thermal port" 
      annotation (extent=[-20,88; 20,108]);
    Medium.BaseProperties medium(preferredMediumStates=true,
                 p(start=p_start), h(start=h_start),
                 T(start=T_start), Xi(start=X_start[1:Medium.nXi]));
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
</html>"),
      Diagram);
    
  equation 
    // boundary conditions
    port_a.p = medium.p;
    port_b.p = medium.p;
    port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium.h);
    port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium.h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium.Xi);
    port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium.Xi);
    thermalPort.T = medium.T;
    
    // Total quantities
    m    = V*medium.d;
    mXi = m*medium.Xi;
    U    = m*medium.u;
    
    // Mass and energy balance
    der(m)   = port_a.m_flow + port_b.m_flow;
    der(mXi) = port_a.mXi_flow + port_b.mXi_flow;
    der(U)   = port_a.H_flow + port_b.H_flow + thermalPort.Q_flow;
    
  initial equation 
    // Initial conditions
    if initOption == NoInit then
      // no initial equations
    elseif initOption == InitialValues then
      medium.p = p_start;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initOption == SteadyState then
      der(medium.p) = 0;
      der(medium.h) = 0;
      der(medium.Xi) = zeros(Medium.nXi);
    elseif initOption == SteadyStateHydraulic then
      der(medium.p) = 0;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    else
      assert(false, "Unsupported initialization option");
    end if;
  end MixingVolume;
  
  model PressureDropPipe 
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
    if frictionType == Modelica_Fluid.Types.FrictionTypes.DetailedFriction then
       d = if port_a.m_flow > 0 then medium_a.d else medium_b.d;
       eta = if port_a.m_flow > 0 then Medium.dynamicViscosity(medium_a) else 
                                      Medium.dynamicViscosity(medium_b);
    else
      // Assign dummy values for auxiliary variables
       d = 0;
       eta = 0;
    end if;
  end PressureDropPipe;
  extends Modelica.Icons.Library;
  
  annotation (preferedView="info",
              Documentation(info="<html>
<p>
This package will contain basic component models of the fluid library.
It is currently empty as all components are being evaluated together 
with examples first (see sub-package Examples). 
</p>
</html>"));
  
  model Tank "Tank with one bottom inlet/outlet port" 
    import Modelica.SIunits.Conversions.*;
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
    parameter SI.Area area "Tank area";
    parameter SI.Volume V0 "Volume of the liquid when the level is zero";
    parameter Medium.AbsolutePressure p_ambient=101325 "Tank surface pressure";
    parameter SI.Height level_start(min=0) "Initial tank level" 
      annotation(Dialog(group=Initialization));
    parameter Medium.Temperature T_start=from_degC(20) 
      "Initial tank temperature" 
      annotation(Dialog(group=Initialization));
    parameter Medium.MassFraction X_start[Medium.nX](quantity=Medium.
          substanceNames) = zeros(Medium.nX) 
      "Initial independent tank mass fractions m_i/m" 
      annotation(Dialog(group=Initialization),enable=(Medium.nXi>0));
    
    Interfaces.FluidPort_b port(redeclare package Medium = Medium) 
      annotation (extent=[-10, -120; 10, -100], rotation=90);
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      Xi(start=X_start[1:Medium.nXi]));
    
    constant Modelica.SIunits.Acceleration g=Modelica.Constants.g_n;
    SI.Height level(start=level_start,stateSelect=StateSelect.prefer) 
      "Level height of tank";
    SI.Energy U "Internal energy of tank volume";
    SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
    Real m(quantity=Medium.mediumName, unit="kg") "Mass of tank volume";
    Real mX[Medium.nX](quantity=Medium.substanceNames, each unit="kg") 
      "Component masses of the independent substances";
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
    
    V = area*level+V0;
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
  /*
initial equation 
  if not Medium.singleState then
    mX = m*X_start[1:Medium.nXi];
  end if;
  level = level_start;
  medium.T = T_start;
  medium.Xi = X_start[1:Medium.nXi];
*/
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
  
  partial model ValveBase "Base model for valves" 
    extends Interfaces.PartialTwoPortTransport;
    import Modelica_Fluid.Types.CvTypes;
    parameter CvTypes.Temp CvData = CvTypes.Av "Selection of flow coefficient" 
       annotation(Dialog(group = "Flow Coefficient"));
    parameter SI.Area Av(fixed = if CvData==CvTypes.Av then true else false,
                         start = m_flow_nom/(sqrt(rhonom*dpnom))*
                                             flowCharacteristic(stemPosition_nom)) = 0 
      "Av (metric) flow coefficient" 
       annotation(Dialog(group = "Flow Coefficient",
                         enable = (CvData==CvTypes.Av)));
    parameter Real Kv(unit="m3/h")=0 "Kv (metric) flow coefficient" 
      annotation(Dialog(group = "Flow Coefficient",
                        enable = (CvData==CvTypes.Kv)));
    parameter Real Cv(unit="USG/min")=0 "Cv (US) flow coefficient" 
      annotation(Dialog(group = "Flow Coefficient",
                        enable = (CvData==CvTypes.Cv)));
    parameter Medium.AbsolutePressure pnom = 1e5 "Nominal inlet pressure" 
      annotation(Dialog(group="Nominal operating point"));
    parameter SI.Pressure dpnom = 1e5 "Nominal pressure drop" 
      annotation(Dialog(group="Nominal operating point"));
    parameter Medium.MassFlowRate m_flow_nom = 1 "Nominal mass flowrate" 
      annotation(Dialog(group="Nominal operating point"));
    parameter Medium.Density rhonom = 1000 "Nominal inlet density" 
      annotation(Dialog(group="Nominal operating point"));
    parameter Real stemPosition_nom = 1 "Nominal stem position" 
      annotation(Dialog(group="Nominal operating point"));
    parameter Boolean CheckValve=false "Reverse flow stopped";
    parameter Real delta=0.01 "Regularisation factor";
    replaceable function flowCharacteristic = Utilities.linear 
      "Inherent flow characteristic";
    Modelica.Blocks.Interfaces.RealInput stemPosition 
      "Stem position in the range 0-1" annotation (extent=[-10,70; 10,90],    rotation=-90);
    
    Medium.Density rho "Density at port a";
    Medium.Temperature T "Temperature at port a";
  protected 
    function sqrtR = Utilities.regRoot(delta = delta*dpnom);
    annotation (
      Icon(Text(extent=[-100, -40; 100, -80], string="%name"),
        Line(points=[0,60; 0,0],   style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Polygon(points=[-100,52; -100,-54; 0,0; -100,52],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Polygon(points=[100,52; 0,0; 100,-54; 100,52],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Rectangle(extent=[-20,70; 20,50],   style(
            color=0,
            fillColor=0,
            fillPattern=1))),
      Diagram,
      Documentation(info="<HTML>
<p>This is the base model for the <tt>ValveLiq</tt>, <tt>ValveLiqChoked</tt>, and <tt>ValveVap</tt> valve models. The model is based on the IEC 534 / ISA S.75 standards for valve sizing.
<p>The model optionally supports reverse flow conditions (assuming symmetrical behaviour) or check valve operation, and has been suitably modified to avoid numerical singularities at zero pressure drop. 
<p>The flow characteristic can be customised.
<p><b>Modelling options</b></p>
<p>The following options are available to specify the valve flow coefficient in fully open conditions:
<ul><li><tt>CvData = ThermoPower.Water.ValveBase.CvTypes.Av</tt>: the flow coefficient is given by the metric <tt>Av</tt> coefficient (m^2).
<li><tt>CvData = ThermoPower.Water.ValveBase.CvTypes.Kv</tt>: the flow coefficient is given by the metric <tt>Kv</tt> coefficient (m^3/h).
<li><tt>CvData = ThermoPower.Water.ValveBase.CvTypes.Cv</tt>: the flow coefficient is given by the US <tt>Cv</tt> coefficient (USG/min).
<li><tt>CvData = ThermoPower.Water.ValveBase.CvTypes.OpPoint</tt>: the flow coefficient must be specified by an additional initial equation (e.g. <tt>w=0.5</tt>); the start value given by <tt>Av</tt> is used to initialise the numerical solution of the equation.
</ul>
<p>The nominal pressure drop <tt>dpnom</tt> must always be specified; to avoid numerical singularities, the flow characteristic is modified for pressure drops less than <tt>b*dpnom</tt> (the default value is 1% of the nominal pressure drop). Increase this parameter if numerical instabilities occur in valves with very low pressure drops.
<p>If <tt>CheckValve</tt> is true, then the flow is stopped when the outlet pressure is higher than the inlet pressure; otherwise, reverse flow takes place.
<p>The default flow characteristic <tt>flowCharacteristic</tt> is linear; this can be replaced by any user-defined function (e.g. equal percentage, quick opening, etc.).
</HTML>",
        revisions="<html>
<ul>
<li><i>6 Apr 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Enumeration-type choice of CvData.</li>
<li><i>16 Dec 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Standard medium definition added.</li>
<li><i>18 Nov 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       <tt>Avnom</tt> removed; <tt>Av</tt> can now be set directly. <tt>Kvnom</tt> and <tt>Cvnom</tt> renamed to <tt>Kv</tt> and <tt>Cv</tt>.<br>
<tt>CvData=3</tt> no longer uses <tt>dpnom</tt>, <tt>wnom</tt> and <tt>rhonom</tt>, and requires an additional initial equation to set the flow coefficient based on the initial working conditions.
<li><i>1 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Valve models restructured using inheritance. <br>
       Adapted to Modelica.Media.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"));
  initial equation 
    if CvData == CvTypes.Kv then
      Av = 2.7778e-5*Kv "Unit conversion";
    elseif CvData == CvTypes.Cv then
      Av = 2.4027e-5*Cv "Unit conversion";
    end if;
    assert(CvData>=0 and CvData<=3, "Invalid CvData");
  equation 
    T = medium_a.T;
    rho = medium_a.d;
  end ValveBase;
  
  model ValveIncompressible "Valve for (almost) incompressible fluids" 
    extends ValveBase;
    import Modelica_Fluid.Types.CvTypes;
  annotation (
    Icon(Text(extent=[-100, -40; 100, -80], string="%name")),
    Diagram,
    Documentation(info="<HTML>
<p>Liquid water valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluid. <p>
Extends the <tt>ValveBase</tt> model (see the corresponding documentation for common valve features).
</html>",
      revisions="<html>
<ul>
<li><i>1 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Valve model restructured using inheritance. <br>
       Adapted to Modelica.Media.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</HTML>"));
  initial equation 
    if CvData == CvTypes.OpPoint then
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*sqrt(rhonom)*sqrtR(dpnom) 
        "Determination of Av by the operating point";
    end if;
    
  equation 
    if CheckValve then
      m_flow = flowCharacteristic(stemPosition)*Av*sqrt(rho)*
               smooth(0,if dp>=0 then sqrtR(dp) else 0);
    else
      m_flow = flowCharacteristic(stemPosition)*Av*sqrt(rho)*sqrtR(dp);
    end if;
  end ValveIncompressible;
  
  model ValveVaporizing 
    "Valve for possibly vaporizing (almost) incompressible fluids, accounts for choked flow conditions" 
    import Modelica_Fluid.Types.CvTypes;
    extends ValveBase(
      redeclare replaceable package Medium = 
      Modelica.Media.Interfaces.PartialTwoPhaseMedium);
    parameter Real Flnom=0.9 "Liquid pressure recovery factor";
    replaceable function Flfun = Utilities.one 
      "Pressure recovery characteristic";
    Real Ff "Ff coefficient (see IEC/ISA standard)";
    Real Fl "Pressure recovery coefficient Fl (see IEC/ISA standard)";
    Medium.AbsolutePressure pv "Saturation pressure";
    SI.Pressure dpEff "Effective pressure drop";
    Medium.AbsolutePressure pin "Inlet pressure";
    Medium.AbsolutePressure pout "Outlet pressure";
    annotation (
      Icon(Text(extent=[-100, -40; 100, -80], string="%name")),
      Diagram,
      Documentation(info="<HTML>
<p>Liquid water valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluid, with possible choked flow conditions. <p>
Extends the <tt>ValveBase</tt> model (see the corresponding documentation for common valve features).<p>
The model operating range includes choked flow operation, which takes place for low outlet pressures due to flashing in the vena contracta; otherwise, non-choking conditions are assumed.
<p>The default liquid pressure recovery coefficient <tt>Fl</tt> is constant and given by the parameter <tt>Flnom</tt>. The relative change (per unit) of the recovery coefficient can be specified as a given function of the valve opening by customising the <tt>Flfun</tt> function.
</HTML>",
        revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Rewritten with sqrtReg.</li>
<<li><i>16 Dec 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Standard medium definition added.</li>
li><i>1 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Valve model restructured using inheritance. <br>
       Adapted to Modelica.Media.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</HTML>"));
  initial equation 
    if CvData == CvTypes.OpPoint then
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*sqrt(rhonom)*sqrtR(dpnom) 
        "Determination of Av by the operating point";
    end if;
  equation 
    pin = port_a.p;
    pout = port_b.p;
    pv = Medium.saturationPressure(T);
    Ff = 0.96 - 0.28*sqrt(pv/Medium.fluidConstants[1].criticalPressure);
    Fl = Flnom*FlCharacteristic(stemPosition);
    dpEff = if pout < (1 - Fl^2)*pin + Ff*Fl^2*pv then 
              Fl^2*(pin - Ff*pv) else dp 
      "Effective pressure drop, accounting for possible choked conditions";
    if CheckValve then
       m_flow = flowCharacteristic(stemPosition)*Av*sqrt(rho)*
           (if dpEff>=0 then sqrtR(dpEff) else 0);
     else
       m_flow = flowCharacteristic(stemPosition)*Av*sqrt(rho)*sqrtR(dpEff);
    end if;
  end ValveVaporizing;
  
  model ValveCompressible 
    "Valve for compressible fluids, accounts for choked flow conditions" 
    extends ValveBase;
    import Modelica_Fluid.Types.CvTypes;
    parameter Real Fxt_full=0.5 "Fk*xt critical ratio at full opening";
    replaceable function xtCharacteristic = Utilities.one 
      "Critical ratio characteristic";
    Real Fxt;
    Real x "Pressure drop ratio";
    Real xs "Saturated pressure drop ratio";
    Real Y "Compressibility factor";
    Medium.AbsolutePressure p "Inlet pressure";
  protected 
    parameter Real Fxt_nom(fixed=false) "Nominal Fxt";
    parameter Real x_nom(fixed=false) "Nominal pressure drop ratio";
    parameter Real xs_nom(fixed=false) "Nominal saturated pressure drop ratio";
    parameter Real Y_nom(fixed=false) "Nominal compressibility factor";
    
    annotation (
    Icon(Text(extent=[-100, -40; 100, -80], string="%name")),
    Diagram,
    Documentation(info="<HTML>
<p>Liquid water valve model according to the IEC 534/ISA S.75 standards for valve sizing, compressible fluid. <p>
Extends the <tt>ValveBase</tt> model (see the corresponding documentation for common valve features).
<p>The product Fk*xt is given by the parameter <tt>Fxtnom</tt>, and is assumed constant by default. The relative change (per unit) of the xt coefficient with the valve opening can be specified by customising the <tt>xtCharacteristic</tt> function.
</HTML>",
      revisions="<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Rewritten with sqrtReg.</li>
<li><i>16 Dec 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Standard medium definition added.</li>
<li><i>1 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Valve model restructured using inheritance. <br>
       Adapted to Modelica.Media.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</HTML>"));
  initial equation 
    if CvData == CvTypes.OpPoint then
      // Determination of Av by the nominal operating point conditions
      Fxt_nom = Fxt_full*xtCharacteristic(stemPosition_nom);
      x_nom = dpnom/pnom;
      xs_nom = smooth(0, if x_nom > Fxt_nom then Fxt_nom else x_nom);
      Y_nom = 1 - abs(xs_nom)/(3*Fxt_nom);
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*Y_nom*sqrt(rhonom)*sqrtR(pnom*xs_nom);
    else
      // Dummy values
      Fxt_nom = 0;
      x_nom = 0;
      xs_nom = 0;
      Y_nom = 0;
    end if;
    
  equation 
    p = noEvent(if dp>=0 then port_a.p else port_b.p);
    Fxt = Fxt_full*xtCharacteristic(stemPosition);
    x = dp/p;
    xs = smooth(0, if x < -Fxt then -Fxt else if x > Fxt then Fxt else x);
    Y = 1 - abs(xs)/(3*Fxt);
    if CheckValve then
      m_flow = flowCharacteristic(stemPosition)*Av*Y*sqrt(rho)*
        smooth(0,if xs>=0 then sqrtR(p*xs) else 0);
    else
      m_flow = flowCharacteristic(stemPosition)*Av*Y*sqrt(rho)*sqrtR(p*xs);
    end if;
  end ValveCompressible;
  
  partial model PumpBase "Base model for centrifugal pumps" 
    import Modelica.SIunits.Conversions.NonSIunits.*;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Medium model";
    Medium.BaseProperties fluid(p(start=pin_start),h(start=hstart)) 
      "Fluid properties at the inlet";
    replaceable package SatMedium = 
        Modelica.Media.Interfaces.PartialTwoPhaseMedium 
      "Saturated medium model (required only for NPSH computation)";
    replaceable function flowCharacteristic = 
        PumpCharacteristics.BaseFlowCharacteristic 
      "Head vs. q_flow characteristic at nominal speed";
    replaceable function powerCharacteristic = 
        PumpCharacteristics.BasePowerCharacteristic 
      "Power consumption characteristic at nominal speed";
    parameter Integer Np0(min=1) = 1 "Nominal number of pumps in parallel";
    parameter Real etaMech(
      min=0,
      max=1) = 0.98 "Mechanical Efficiency";
    parameter SI.Mass M = 0 "Fluid mass inside the pump";
    parameter Boolean CheckValve=false "Reverse flow stopped";
    parameter Boolean ComputeNPSHa=false "Compute NPSH Available at the inlet";
    parameter Medium.AbsolutePressure pin_start "Inlet Pressure Start Value" 
      annotation(Dialog(tab="Initialization"));
    parameter Medium.AbsolutePressure pout_start "Outlet Pressure Start Value" 
      annotation(Dialog(tab="Initialization"));
    parameter Medium.SpecificEnthalpy hstart=1e5 
      "Fluid Specific Enthalpy Start Value" 
      annotation(Dialog(tab="Initialization"));
    parameter SI.MassFlowRate m_flow_start = 0 
      "Start value of mass flow rate (single pump)" 
      annotation(Dialog(tab="Initialization"));
    constant SI.Acceleration g=Modelica.Constants.g_n;
  //  parameter Choices.Init.Options.Temp initOpt=Choices.Init.Options.noInit 
  //    "Initialisation option";
    Modelica_Fluid.Interfaces.FluidPort_a inlet(redeclare package Medium = Medium,
        p(start=pin_start)) 
    annotation (extent=[-100,-40; -60,0]);
    Modelica_Fluid.Interfaces.FluidPort_b outlet(redeclare package Medium = 
          Medium, p(start=pout_start)) 
    annotation (extent=[40,12; 80,52]);
    SI.Pressure dp = outlet.p - inlet.p "Pressure increase";
    SI.Height head = dp/(rho*g) "Pump head";
    Medium.Density rho "Liquid density";
    Medium.SpecificEnthalpy h_out 
      "Enthalpy of the liquid flowing out of the pump";
    Medium.Temperature Tin "Liquid inlet temperature";
    SI.MassFlowRate m_flow = inlet.m_flow "Mass flow rate (total)";
    SI.MassFlowRate m_flow_single = m_flow/Np "Mass flow rate (single pump)";
    SI.VolumeFlowRate q_flow = m_flow/rho "Volume flow rate (total)";
    SI.VolumeFlowRate q_flow_single = q_flow/Np 
      "Volume flow rate (single pump)";
    AngularVelocity_rpm n "Shaft r.p.m.";
    Integer Np(min=1) "Number of pumps in parallel";
    SI.Power P_single "Power Consumption (single pump)";
    SI.Power P_tot = P_single*Np "Power Consumption (total)";
    constant SI.Power P_eps=1e-8 
      "Small coefficient to avoid numerical singularities in efficiency computations";
    SI.Power Phyd_single "Hydraulic power (single pump)";
    Real eta "Global Efficiency";
    SI.Length NPSHa "Net Positive Suction Head available";
    Medium.AbsolutePressure pv "Saturation pressure of inlet liquid";
    Real s(start = m_flow_start) 
      "Curvilinear abscissa for the flow curve in parametric form";
    Modelica.Blocks.Interfaces.IntegerInput in_Np 
      annotation (extent=[16,34; 36,54], rotation=-90);
  equation 
    // Number of pumps in parallel
    Np = in_Np;
    if cardinality(in_Np)==0 then
      in_Np = Np0 "Number of pumps selected by parameter";
    end if;
    
    // Flow equations
    if noEvent(s > 0 or (not CheckValve)) then
      // Flow characteristics when check valve is open
      m_flow_single = s;
      head = flowCharacteristic(q_flow_single, n);
    else
      // Flow characteristics when check valve is closed
      head = flowCharacteristic(0,n) - s;
      q_flow_single = 0;
    end if;
    
    // Power consumption  
    P_single = powerCharacteristic(q_flow_single, head, dp, n) 
      "Power consumption (single pump)";
    // Hydraulic power
    Phyd_single = P_single*etaMech 
      "Hydraulic power transferred to the fluid (single pump)";
    eta = (dp*q_flow_single)/(P_single + P_eps) "Hydraulic efficiency";
    
    // Fluid properties
    fluid.p = inlet.p;
    fluid.h = inlet.h;
    fluid.Xi = inlet.Xi;
    rho = fluid.d;
    Tin = fluid.T;
    
    // Mass and energy balances
    inlet.m_flow + outlet.m_flow = 0 "Mass balance";
    inlet.mXi_flow + outlet.mXi_flow = zeros(Medium.nXi) 
      "Substance mass balance";
    inlet.H_flow=semiLinear(inlet.m_flow,inlet.h,h_out) 
      "Enthalpy flow at the inlet";
    outlet.H_flow=semiLinear(outlet.m_flow,outlet.h,h_out) 
      "Enthalpy flow at the outlet";
    if M > 0 then
      M * der(h_out) = m_flow_single*(inlet.h - outlet.h) + Phyd_single 
        "Dynamic energy balance (density variations neglected)";
    else
      inlet.H_flow + outlet.H_flow + Phyd_single*Np = 0 "Static energy balance";
    end if;
    
    // NPSH computations
    if ComputeNPSHa then
      pv=SatMedium.saturationPressure(fluid.T);
      NPSHa=(inlet.p-pv)/(rho*Modelica.Constants.g_n);
    else
      pv=0;
      NPSHa=0;
    end if;
  /*
initial equation 
  if initOpt == Choices.Init.Options.noInit then
    // do nothing
  elseif initOpt == Choices.Init.Options.steadyState then
    if ThermalCapacity then
      der(h)=0;
    end if;
  else
    assert(false, "Unsupported initialisation option");
  end if;
*/
    annotation (
      Icon(
        Polygon(points=[-40,-64; -60,-100; 60,-100; 40,-64; -40,-64],
            style(pattern=0, fillColor=74)),
        Ellipse(extent=[-60,40; 60,-80],   style(gradient=3)),
        Polygon(points=[-30,12; -30,-48; 48,-20; -30,12],   style(
            pattern=0,
            gradient=2,
            fillColor=7)),
        Text(extent=[-100,-110; 100,-136], string="%name"),
        Text(extent=[-10,60; 18,40],  string="Np")),
      Diagram,
      Documentation(info="<HTML>
<p>This is the base model for the <tt>Pump</tt> and <tt>
PumpMech</tt> pump models.
<p>The model describes a centrifugal pump, or a group of <tt>Np</tt> identical pumps in parallel. The hydraulic characteristic (head vs. flowrate) is represented, as well as the pump power consumption.
<p>In order to avoid singularities in the computation of the outlet enthalpy at zero flowrate, the thermal capacity of the fluid inside the pump body can be taken into account.
<p>The model can either support reverse flow conditions or include a built-in check valve to avoid flow reversal.
<p><b>Modelling options</b></p>
<p>The following options are available to specify the pump characteristics:
<ul><li><tt>OpPoints = false</tt>: the coefficients of the characteristics must be specified by 6 additional initial equations
<li><tt>OpPoints = true</tt>: the characteristics are specified by providing a vector of three operating points (in terms of heads <tt>head[3]</tt>, volume flow rate <tt>q[3]</tt>, power consumption <tt>P_cons[3]</tt>, nominal fluid density <tt>rho0</tt>, and nominal rotational speed <tt>n0</tt>) for a single pump.
</ul>
<p>If the <tt>in_Np</tt> input connector is wired, it provides the number of pumps in parallel; otherwise,  <tt>Np0</tt> parallel pumps are assumed.</p>
<p>If <tt>ThermalCapacity</tt> is set to true, the heat capacity of the fluid inside the pump is taken into account: this is necessary to avoid singularities in the computation of the outlet enthalpy in case of zero flowrate. If zero flowrate conditions are always avoided, this effect can be neglected by setting <tt>ThermalCapacity</tt> to false, thus avoiding a fast state variable in the model.
<p>The <tt>CheckValve</tt> parameter determines whether the pump has a built-in check valve or not.
</HTML>",
        revisions="<html>
<ul>
<li><i>6 Apr 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       <tt>CharData</tt> substituted by <tt>OpPoints</tt></li>
<li><i>16 Dec 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Standard medium definition added.</li>
<li><i>2 Aug 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Optional NPSHa computation added. Changed parameter names</li>
<li><i>5 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model restructured by using inheritance. Adapted to Modelica.Media.</li>
<li><i>15 Jan 2004</i>
    by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br>
       <tt>ThermalCapacity</tt> and <tt>CheckValve</tt> added.</li>
<li><i>15 Dec 2003</i>
    by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br>
       First release.</li>
</ul>
</html>"));
    
  end PumpBase;
  
  package PumpCharacteristics 
    import NonSI = Modelica.SIunits.Conversions.NonSIunits;
    
    partial function BaseFlowCharacteristic 
      "Base class for pump flow characteristics" 
      extends Modelica.Icons.Function;
      input SI.VolumeFlowRate q_flow "Volumetric flow rate";
      input NonSI.AngularVelocity_rpm n "Rotational speed";
      output SI.Height head "Pump head";
    end BaseFlowCharacteristic;
    
    partial function BasePowerCharacteristic 
      "Base class for pump power consumption characteristics" 
      extends Modelica.Icons.Function;
      input SI.VolumeFlowRate q_flow "Volumetric flow rate";
      input SI.Height head "Pump head";
      input SI.Pressure dp "Outlet pressure minus inlet pressure";
      input NonSI.AngularVelocity_rpm n "Rotational speed";
      output SI.Power consumption "Power consumption";
    end BasePowerCharacteristic;
    
    function LinearFlowCharacteristic 
      extends BaseFlowCharacteristic;
      input SI.Height head_nom[2] "Pump head for two operating points"     annotation(Dialog);
      input SI.VolumeFlowRate q_nom[2] 
        "Volume flow rate for two operating points (single pump)" annotation(Dialog);
      input NonSI.AngularVelocity_rpm n0 "Nominal rotational speed"       annotation(Dialog);
    protected 
      constant Real g = Modelica.Constants.g_n;
      /* Linear system to determine the coefficients:
  h_nom[1]*g = q_nom[1]*(n0/n0)*c[1] + (nO/n0)^2*c[2];
  h_nom[2]*g = q_nom[2]*(n0/n0)*c[1] + (nO/n0)^2*c[2];
  */
      Real c[2] = Modelica.Math.Matrices.solve([q_nom,ones(2)],head_nom*g) 
        "Coefficients of linear head curve";
    algorithm 
      // head * g = q*(n/n0)*c[1] + (n/n0)^2*c[2];
      head := 1/g * (q_flow*(n/n0)*c[1] + (n/n0)^2*c[2]);
    end LinearFlowCharacteristic;
    
    function QuadraticFlowCharacteristic 
      extends BaseFlowCharacteristic;
      input SI.Height head_nom[3] "Pump head for three operating points" annotation(Dialog);
      input SI.VolumeFlowRate q_nom[3] 
        "Volume flow rate for three operating points (single pump)" annotation(Dialog);
      input NonSI.AngularVelocity_rpm n0 "Nominal rotational speed"       annotation(Dialog);
    protected 
      constant Real g = Modelica.Constants.g_n;
      Real q_nom2[3] = {q_nom[1]^2,q_nom[2]^2, q_nom[3]^2} 
        "Squared nominal flow rates";
      /* Linear system to determine the coefficients:
  h_nom[1]*g = q_nom[1]^2*c[1] + q_nom[1]*(n0/n0)*c[2] + (nO/n0)^2*c[3];
  h_nom[2]*g = q_nom[2]^2*c[1] + q_nom[2]*(n0/n0)*c[2] + (nO/n0)^2*c[3];
  h_nom[3]*g = q_nom[3]^2*c[1] + q_nom[3]*(n0/n0)*c[2] + (nO/n0)^2*c[3];
  */
      Real c[3] = Modelica.Math.Matrices.solve([q_nom2, q_nom,ones(3)],head_nom*g) 
        "Coefficients of quadratic head curve";
    algorithm 
      head := 1/g * (q_flow^2*c[1] + q_flow*(n/n0)*c[2] + (n/n0)^2*c[3]);
    end QuadraticFlowCharacteristic;
    
    function ConstantEfficiencyPowerCharacteristic 
       extends BasePowerCharacteristic;
       input Real efficiency "Hydraulic efficiency";
    algorithm 
      consumption := q_flow*dp/efficiency;
    end ConstantEfficiencyPowerCharacteristic;
    
    function QuadraticPowerCharacteristic 
      extends BasePowerCharacteristic;
      input SI.Height P_nom[3] "Power consumption for three operating points" 
                                                                         annotation(Dialog);
      input SI.VolumeFlowRate q_nom[3] 
        "Volume flow rate for three operating points (single pump)" annotation(Dialog);
      input NonSI.AngularVelocity_rpm n0 "Nominal rotational speed"       annotation(Dialog);
    protected 
      Real q_nom2[3] = {q_nom[1]^2,q_nom[2]^2, q_nom[3]^2} 
        "Squared nominal flow rates";
      /* Linear system to determine the coefficients:
  P_nom[1] = (n0/n0)^2*q_nom[1]*c[1] + (n0/n0)*q_nom[1]^2*c[2] + (nO/n0)^2*c[3];
  P_nom[2] = (n0/n0)^2*q_nom[2]*c[1] + (n0/n0)*q_nom[2]^2*c[2] + (nO/n0)^2*c[3];
  P_nom[3] = (n0/n0)^2*q_nom[3]*c[1] + (n0/n0)*q_nom[3]^2*c[2] + (nO/n0)^2*c[3];
  */
      Real c[3] = Modelica.Math.Matrices.solve([q_nom,q_nom2,ones(3)],P_nom) 
        "Coefficients of quadratic power consumption curve";
    algorithm 
      consumption :=c[1]*(n/n0)^2*q_flow + c[2]*(n/n0)*q_flow^2 + c[3]*(n/n0)^2;
    end QuadraticPowerCharacteristic;
    
  end PumpCharacteristics;
  
  model Pump "Centrifugal pump with ideally controlled speed" 
    extends PumpBase;
    import Modelica.SIunits.Conversions.NonSIunits.*;
    parameter AngularVelocity_rpm n_const "Constant rotational speed";
    Modelica.Blocks.Interfaces.RealInput in_n "RPM" 
      annotation (extent=[-36,34; -16,54],   rotation=-90);
  equation 
      n = in_n "Rotational speed";
    if cardinality(in_n)==0 then
      in_n = n_const "Rotational speed provided by parameter";
    end if;
    annotation (
      Icon(
        Text(extent=[-58,58; -30,38], string="n")),
      Diagram,
      Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>Np</tt> pumps in parallel) with controlled speed, either fixed or provided by an external signal.
<p>The model extends <tt>PumpBase</tt>
<p>If the <tt>in_n</tt> input connector is wired, it provides rotational speed of the pumps (rpm); otherwise, a constant rotational speed equal to <tt>n_const</tt> (which can be different from <tt>n0</tt>) is assumed.</p>
</HTML>",
        revisions="<html>
<ul>
<li><i>5 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model restructured by using inheritance. Adapted to Modelica.Media.</li>
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
    SI.Angle phi "Shaft angle";
    SI.AngularVelocity omega "Shaft angular velocity";
    Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft 
    annotation (extent=[80,4; 110,32]);
  equation 
    phi = shaft.phi;
    omega = der(phi);
    n = Modelica.SIunits.Conversions.to_rpm(omega);
    P_single = omega*shaft.tau;
  annotation (
    Icon(Rectangle(extent=[60,26; 84,12], style(
            color=10,
            rgbcolor={95,95,95},
            gradient=2,
            fillColor=10,
            rgbfillColor={95,95,95}))),
    Diagram,
    Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>Np</tt> pumps in parallel) with a mechanical rotational connector for the shaft, to be used when the pump drive has to be modelled explicitly. In the case of <tt>Np</tt> pumps in parallel, the mechanical connector is relative to a single pump.
<p>The model extends <tt>PumpBase</tt>
 </HTML>",
       revisions="<html>
<ul>
<li><i>5 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model restructured by using inheritance. Adapted to Modelica.Media.</li>
<li><i>15 Jan 2004</i>
    by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br>
       <tt>ThermalCapacity</tt> and <tt>CheckValve</tt> added.</li>
<li><i>15 Dec 2003</i>
    by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco Schiavo</a>:<br>
       First release.</li>
</ul>
</html>"));
  end PumpMech;
  
  model ValveLinear "Valve for water/steam flows with linear pressure drop" 
    extends Interfaces.PartialTwoPortTransport;
    parameter Types.HydraulicConductance Kv = 1/1e5 
      "Hydraulic conductance at full opening";
    Modelica.Blocks.Interfaces.RealInput opening 
    annotation (extent=[-20, 60; 20, 100], rotation=-90);
  equation 
    m_flow = Kv*opening*dp;
    
  annotation (
    Icon(Text(extent=[-100,-66; 100,-100],  string="%name"),
        Line(points=[0,54; 0,-10], style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Polygon(points=[-100,42; -100,-64; 0,-10; -100,42],style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Polygon(points=[100,42; 0,-10; 100,-64; 100,42],style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Rectangle(extent=[-20,74; 20,54],   style(
            color=0,
            fillColor=0,
            fillPattern=1))),
    Diagram,
    Documentation(info="<HTML>
<p>This very simple model provides a pressure drop which is proportional to the flowrate and to the <tt>cmd</tt> signal, without computing any fluid property.</p>
</HTML>",
      revisions="<html>
<ul>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"));
  end ValveLinear;
  
  model Evaporator 
    "Simple Evaporator with two states, see Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378" 
    import Modelica_Fluid.Interfaces.*;
    import Modelica.SIunits.Conversions.*;
    // property and interface declarations
    replaceable package Medium = 
        Modelica.Media.Interfaces.PartialTwoPhaseMedium 
      extends Modelica.Media.Interfaces.PartialTwoPhaseMedium "Medium model" 
        annotation (choicesAllMatching=true);
    Medium.SaturationProperties sat 
      "State vector to compute saturation properties";
    FluidPort_a feedwater(redeclare package Medium = Medium) 
      annotation (extent=[-120, -10; -100, 10]);
    FluidPort_b steam(redeclare package Medium = Medium) 
      annotation (extent=[120, -10; 100, 10]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
      annotation (extent=[-10, -120; 10, -100]);
    Modelica.Blocks.Interfaces.RealOutput V(
      redeclare type SignalType = SI.Volume) "liquid volume" 
      annotation (extent=[30, 100; 50, 120], rotation=90);
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
    parameter Medium.SpecificHeatCapacity cp_D=500 
      "specific heat capacity of drum metal";
    parameter SI.Volume V_t=100 "total volume inside drum";
    parameter Medium.AbsolutePressure p_start=from_bar(1) "initial pressure";
    parameter SI.Volume V_start=67 "initial liquid volume";
    
    Medium.AbsolutePressure p(start=p_start, fixed = true, stateSelect=StateSelect.prefer) 
      "pressure inside drum boiler";
    Medium.Temperature T "temperature inside drum boiler";
    SI.Volume V_v "volume of vapour phase";
    SI.Volume V_l(start=V_start, fixed = true, stateSelect=StateSelect.prefer) 
      "volumes of liquid phase";
    Medium.SpecificEnthalpy h_v=Medium.dewEnthalpy(sat) 
      "specific enthalpy of vapour";
    Medium.SpecificEnthalpy h_l=Medium.bubbleEnthalpy(sat) 
      "specific enthalpy of liquid";
    Medium.Density rho_v=Medium.dewDensity(sat) "density in vapour phase";
    Medium.Density rho_l=Medium.bubbleDensity(sat) "density in liquid phase";
    SI.Mass m "total mass of drum boiler";
    SI.Energy U "internal energy";
    Medium.Temperature T_D=heatPort.T "temperature of drum";
    SI.HeatFlowRate q_F=heatPort.Q_flow "heat flow rate from furnace";
    Medium.SpecificEnthalpy h_W=feedwater.h "feed water enthalpy";
    Medium.SpecificEnthalpy h_S=steam.h "steam enthalpy";
    SI.MassFlowRate qm_W=feedwater.m_flow "feed water mass flow rate";
    SI.MassFlowRate qm_S=steam.m_flow "steam mass flow rate";
  equation 
    // balance equations  
    m = rho_v*V_v + rho_l*V_l + m_D "Total mass";
    U = rho_v*V_v*h_v + rho_l*V_l*h_l - p*V_t + m_D*cp_D*T_D "Total energy";
    der(m) = qm_W + qm_S "Mass balance";
    der(U) = q_F + qm_W*h_W + qm_S*h_S "Energy balance";
    V_t = V_l + V_v;
    
    // Properties of saturated liquid and steam
    sat.psat = p;
    sat.Tsat = T;
    sat.Tsat = Medium.saturationTemperature(p);
    
    // ideal heat transfer between metal and water
    T_D = T;
    
    // boundary conditions at the ports
    feedwater.p = p;
    feedwater.H_flow = semiLinear(feedwater.m_flow, feedwater.h, h_l);
    steam.p = p;
    steam.H_flow = semiLinear(steam.m_flow, steam.h, h_v);
    
    // liquid volume 
    V = V_l;
  end Evaporator;
  
  model Flow1D   
  end Flow1D;
  
  model Pipe   
  end Pipe;
end Components;
