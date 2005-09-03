package Components "Basic components for fluid models" 
  model PortVolume 
    "Fixed volume associated with a port by the finite volume method (used to build up physical components; fulfills mass and energy balance)" 
    import SI = Modelica.SIunits;
    
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Medium model" 
       annotation (choicesAllMatching=true);
    
    parameter SI.Volume V=1e-6 "Volume";
    
    parameter Boolean use_p_start=true "select p_start or d_start" 
      annotation (Evaluate=true, Dialog(group="Initial pressure or initial density"));
    parameter Medium.AbsolutePressure p_start = Medium.reference_p 
      "Initial pressure" 
      annotation (Dialog(group="Initial pressure or initial density", enable=use_p_start));
    parameter Medium.Density d_start=1 "Initial density" 
      annotation (Dialog(group="Initial pressure or initial density", enable=not use_p_start));
    parameter Boolean use_T_start=true "select T_start or h_start" 
      annotation (Evaluate=true, Dialog(group="Initial temperature or initial specific enthalpy"));
    parameter Medium.Temperature T_start = Modelica.SIunits.Conversions.from_degC(20) 
      "Initial temperature" 
      annotation (Dialog(group="Initial temperature or initial specific enthalpy", enable=use_T_start));
    parameter Medium.SpecificEnthalpy h_start = 1.e4 
      "Initial specific enthalpy" 
      annotation (Dialog(group="Initial temperature or initial specific enthalpy", enable=not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.reference_X 
      "Initial mass fractions m_i/m" 
      annotation (Dialog(group="Only for multi-substance flow", enable=Medium.nX > 0));
    
    Modelica_Fluid.Interfaces.FluidPort_a port(redeclare package Medium = 
                 Medium)                                annotation (extent=[-10, -10; 10, 10], rotation=0);
    Medium.BaseProperties medium(preferredMediumStates=true);
    SI.Energy U "Internal energy of port volume";
    SI.Mass m "Mass of junction volume";
    SI.Mass mXi[Medium.nXi] "Independent substance masses of port volume";
    
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
    
  /*
initial equation 
  if not Medium.singleState then
    if use_p_start then
       medium.p = p_start;
    else
       medium.d = d_start;
    end if;
  end if;
  
  if use_T_start then
     medium.T = T_start;
  else
     medium.h = h_start;
  end if;  
  medium.Xi = X_start[1:Medium.nXi];
*/
    
    Interfaces.HeatPort_a thermalPort annotation (extent=[-20,88; 20,108]);
  equation 
    // port = medium properties
      port.p = medium.p;
      port.h = medium.h;
      port.Xi = medium.Xi;
      thermalPort.T = medium.T;
    
    // Total quantities
       m    = V*medium.d;
       mXi = m*medium.Xi;
       U    = m*medium.u;
    
    // Mass and energy balance
       der(m)    = port.m_flow;
       der(mXi)  = port.mXi_flow;
       der(U)    = port.H_flow + thermalPort.Q_flow;
  end PortVolume;

  model MixingVolume 
    "Mixing volume with inlet and outlet ports (flow reversal is allowed)" 
    import SI = Modelica.SIunits;
    
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Medium model" 
       annotation (choicesAllMatching=true);
    
    parameter SI.Volume V=1e-6 "Fixed size of junction volume";
    
    parameter Boolean use_p_start=true "select p_start or d_start" 
      annotation (Evaluate=true, Dialog(group="Initial pressure or initial density"));
    parameter Medium.AbsolutePressure p_start = Medium.reference_p 
      "Initial pressure" 
      annotation (Dialog(group="Initial pressure or initial density", enable=use_p_start));
    parameter Medium.Density d_start=1 "Initial density" 
      annotation (Dialog(group="Initial pressure or initial density", enable=not use_p_start));
    parameter Boolean use_T_start=true "select T_start or h_start" 
      annotation (Evaluate=true, Dialog(group="Initial temperature or initial specific enthalpy"));
    parameter Medium.Temperature T_start = Modelica.SIunits.Conversions.from_degC(20) 
      "Initial temperature" 
      annotation (Dialog(group="Initial temperature or initial specific enthalpy", enable=use_T_start));
    parameter Medium.SpecificEnthalpy h_start = 1.e4 
      "Initial specific enthalpy" 
      annotation (Dialog(group="Initial temperature or initial specific enthalpy", enable=not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.reference_X 
      "Initial mass fractions m_i/m" 
      annotation (Dialog(group="Only for multi-substance flow", enable=Medium.nX > 0));
    
    Modelica_Fluid.Interfaces.FluidPort_a port_a(
      redeclare package Medium = Medium) annotation (extent=[-112,-4; -92,16]);
    Modelica_Fluid.Interfaces.FluidPort_b port_b(
      redeclare package Medium = Medium) annotation (extent=[92,-8; 112,12]);
    Medium.BaseProperties medium(preferredMediumStates=true);
    SI.Energy U "Internal energy of port volume";
    SI.Mass m "Mass of junction volume";
    SI.Mass mXi[Medium.nXi] "Independent substance masses of port volume";
    
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
    
  /*
initial equation 
  if not Medium.singleState then
    if use_p_start then
       medium.p = p_start;
    else
       medium.d = d_start;
    end if;
  end if;
  
  if use_T_start then
     medium.T = T_start;
  else
     medium.h = h_start;
  end if;
  
  medium.Xi = X_start[1:Medium.nXi];
*/
    Interfaces.HeatPort_a thermalPort annotation (extent=[-20,88; 20,108]);
  equation 
    // boundary condition
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
  end MixingVolume;

  model PressureDrop 
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
  end PressureDrop;
  import Modelica.SIunits.*;
  extends Modelica.Icons.Library;
  
  annotation (preferedView="info",
              Documentation(info="<html>
<p>
This package will contain basic component models of the fluid library.
It is currently empty as all components are being evaluated together 
with examples first (see sub-package Examples). 
</p>
</html>"));
  
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
    parameter Modelica.SIunits.Height level_start(min=0) "Initial tank level" 
      annotation(Dialog(group=Initialization));
    parameter Medium.Temperature T_start=from_degC(20) 
      "Initial tank temperature" 
      annotation(Dialog(group=Initialization));
    parameter Medium.MassFraction X_start[Medium.nX](quantity=Medium.
          substanceNames) = zeros(Medium.nX) 
      "Initial independent tank mass fractions m_i/m" 
      annotation(Dialog(group=Initialization),enable=(Medium.nXi>0));
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
  
  partial model ValveBase "Base model for valves" 
    extends Interfaces.PartialTwoPortTransport;
    import Modelica_Fluid.Types.CvTypes;
    parameter CvTypes.Temp CvData "Selection of flow coefficient";
    parameter Area Av(fixed = if CvData==CvTypes.Av then true else false,
                      start = m_flow_nom/(sqrt(rhonom*dpnom))) = 0 
      "Av (metric) flow coefficient";
    parameter Real Kv(unit="m^3/h")=0 "Kv (metric) flow coefficient";
    parameter Real Cv(unit="USG/min")=0 "Cv (US) flow coefficient";
    parameter Real pnom "Nominal inlet pressure";
    parameter Pressure dpnom "Nominal pressure drop";
    parameter Pressure rhonom = 1000 "Nominal inlet density";
    parameter Real stemPosition_nom = 1 "Nominal stem position";
    parameter MassFlowRate m_flow_nom "Nominal mass flowrate";
    parameter Boolean CheckValve=false "Reverse flow stopped";
    parameter Real b=0.01 "Regularisation factor";
    replaceable function FlowChar = Utilities.linear 
      "Inherent flow characteristic";
    Medium.Density rho "Density at port a";
    Medium.Temperature T "Temperature at port a";
  protected 
    function sqrtR = Utilities.sqrtReg(delta = b*dpnom);
  public 
    Modelica.Blocks.Interfaces.RealInput stemPosition 
      "Stem position in the range 0-1" annotation (extent=[-10,70; 10,90],    rotation=-90);
    
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
<p>The default flow characteristic <tt>FlowChar</tt> is linear; this can be replaced by any user-defined function (e.g. equal percentage, quick opening, etc.).
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
      m_flow_nom = FlowChar(stemPosition_nom)*Av*sqrt(rhonom)*sqrtR(dp) 
        "Determination of Av by the operating point";
    end if;
    
  equation 
    if CheckValve then
      m_flow = FlowChar(stemPosition)*Av*sqrt(rho)*
               smooth(0,if dp>=0 then sqrtR(dp) else 0);
    else
      m_flow = FlowChar(stemPosition)*Av*sqrt(rho)*sqrtR(dp);
    end if;
  end ValveIncompressible;
  
  model ValveIncompressibleChoked 
    "Valve for (almost) incompressible fluids, accounts for choked flow conditions" 
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
    Pressure dpEff "Effective pressure drop";
    Pressure pin "Inlet pressure";
    Pressure pout "Outlet pressure";
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
      m_flow_nom = FlowChar(stemPosition_nom)*Av*sqrt(rhonom)*sqrtR(dp) 
        "Determination of Av by the operating point";
    end if;
  equation 
    pin = port_a.p;
    pout = port_b.p;
    pv = Medium.saturationPressure(T);
    Ff = 0.96 - 0.28*sqrt(pv/Medium.fluidConstants[1].criticalPressure);
    Fl = Flnom*Flfun(stemPosition);
    dpEff = if pout < (1 - Fl^2)*pin + Ff*Fl^2*pv then 
              Fl^2*(pin - Ff*pv) else dp 
      "Effective pressure drop, accounting for possible choked conditions";
    if CheckValve then
       m_flow = FlowChar(stemPosition)*Av*sqrt(rho)*
           (if dpEff>=0 then sqrtR(dpEff) else 0);
     else
       m_flow = FlowChar(stemPosition)*Av*sqrt(rho)*sqrtR(dpEff);
    end if;
  end ValveIncompressibleChoked;
  
  model ValveCompressible 
    "Valve for compressible fluids, accounts for choked flow conditions" 
    extends ValveBase;
    import Modelica_Fluid.Types.CvTypes;
    parameter Real Fxt_full=0.5 "Fk*xt critical ratio at full opening";
    replaceable function xtfun = Utilities.one "Critical ratio characteristic";
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
<p>The product Fk*xt is given by the parameter <tt>Fxtnom</tt>, and is assumed constant by default. The relative change (per unit) of the xt coefficient with the valve opening can be specified by customising the <tt>xtfun</tt> function.
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
      Fxt_nom = Fxt_full*xtfun(stemPosition_nom);
      x_nom = dpnom/pnom;
      xs_nom = smooth(0, if x_nom > Fxt_nom then Fxt_nom else x_nom);
      Y_nom = 1 - abs(xs_nom)/(3*Fxt_nom);
      m_flow_nom = FlowChar(stemPosition_nom)*Av*Y_nom*sqrt(rhonom)*sqrtR(pnom*xs_nom);
    else
      // Dummy values
      Fxt_nom = 0;
      x_nom = 0;
      xs_nom = 0;
      Y_nom = 0;
    end if;
    
  equation 
    p = noEvent(if dp>=0 then port_a.p else port_b.p);
    Fxt = Fxt_full*xtfun(stemPosition);
    x = dp/p;
    xs = smooth(0, if x < -Fxt then -Fxt else if x > Fxt then Fxt else x);
    Y = 1 - abs(xs)/(3*Fxt);
    if CheckValve then
      m_flow = FlowChar(stemPosition)*Av*Y*sqrt(rho)*
        smooth(0,if xs>=0 then sqrtR(p*xs) else 0);
    else
      m_flow = FlowChar(stemPosition)*Av*Y*sqrt(rho)*sqrtR(p*xs);
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
    parameter Integer Np0(min=1) = 1 "Nominal number of pumps in parallel";
    parameter Pressure pin_start "Inlet Pressure Start Value";
    parameter Pressure pout_start "Outlet Pressure Start Value";
    parameter SpecificEnthalpy hstart=1e5 "Fluid Specific Enthalpy Start Value";
    parameter Density rho0=1000 "Nominal Liquid Density";
    parameter AngularVelocity_rpm n0=1500 "Nominal rotational speed";
    parameter Mass M=1 "Fluid mass inside the pump";
    parameter Real etaMech(
      min=0,
      max=1) = 0.98 "Mechanical Efficiency";
    parameter Boolean OpPoints = true 
      "Characteristic curves defined by operation points";
    parameter Height head_nom[3] "Pump head for three operating points";
    parameter VolumeFlowRate q_nom[3] 
      "Volume flow rate for three operating points (single pump)";
    parameter Power P_cons[3] 
      "Power consumption for three operating points (single pump)";
    parameter Boolean CheckValve=false "Reverse flow stopped";
    parameter Boolean ComputeNPSHa=false "Compute NPSH Available at the inlet";
  //  parameter Choices.Init.Options.Temp initOpt=Choices.Init.Options.noInit 
  //    "Initialisation option";
    constant Acceleration g=Modelica.Constants.g_n;
    Modelica_Fluid.Interfaces.FluidPort_a inlet(redeclare package Medium = Medium,
        p(start=pin_start)) 
    annotation (extent=[-100, 2; -60, 42]);
    Modelica_Fluid.Interfaces.FluidPort_b outlet(redeclare package Medium = 
          Medium, p(start=pout_start)) 
    annotation (extent=[40,54; 80,94]);
    Modelica.Blocks.Interfaces.IntegerInput in_Np "Number of  parallel pumps" 
    annotation (extent=[16,76; 36,96], rotation=-90);
    
    Density rho "Liquid density";
    SpecificEnthalpy h "Enthalpy of liquid inside the pump";
    Medium.Temperature Tin "Liquid inlet temperature";
    MassFlowRate m_flow_1 "Mass flowrate (single pump)";
    AngularVelocity_rpm n "Shaft r.p.m.";
    Integer Np(min=1) "Number of pumps in parallel";
    
    Power P "Power Consumption (single pump)";
    Power Ptot "Power Consumption (total)";
    constant Power P_eps=1e-8 
      "Small coefficient to avoid numerical singularities";
    Power Phyd "Hydraulic power (single pump)";
    Real eta "Global Efficiency";
    Length NPSHa "Net Positive Suction Head available";
    Medium.AbsolutePressure pv "Saturated liquid pressure";
    Boolean FlowOn(start=true);
    Real s "Auxiliary Variable";
    parameter Real A(fixed=false) "Flow characteristics coefficient";
    parameter Real B(fixed=false) "Flow characteristics coefficient";
    parameter Real C(fixed=false) "Flow characteristics coefficient";
    parameter Real D(fixed=false) "Consumption characteristics coefficient";
    parameter Real E(fixed=false) "Consumption characteristics coefficient";
    parameter Real F(fixed=false) "Consumption characteristics coefficient";
  initial equation 
    // Computation of flow characteristics coefficients from the nominal
    // operating points
    if OpPoints then
      head_nom[1]*g = -A*q_nom[1]^2 + B*q_nom[1] + C;
      head_nom[2]*g = -A*q_nom[2]^2 + B*q_nom[2] + C;
      head_nom[3]*g = -A*q_nom[3]^2 + B*q_nom[3] + C;
      P_cons[1] = D*(n0^2)*q_nom[1] - E*n0*(q_nom[1]^2) + F*(n0^2);
      P_cons[2] = D*(n0^2)*q_nom[2] - E*n0*(q_nom[2]^2) + F*(n0^2);
      P_cons[3] = D*(n0^2)*q_nom[3] - E*n0*(q_nom[3]^2) + F*(n0^2);
    end if;
  equation 
    // Number of pumps in parallel
    Np = in_Np;
    if cardinality(in_Np)==0 then
      in_Np = Np0 "Number of pumps selected by parameter";
    end if;
    
    // Flow equations
    m_flow_1 = inlet.m_flow/Np "Single pump flow rate";
    FlowOn = s > 0;
    if (FlowOn or (not CheckValve)) then
      // Flow characteristics when check valve is open
      m_flow_1 = s;
      (outlet.p - inlet.p)/rho = -A*(m_flow_1/rho)^2 + B*(n/n0)*m_flow_1/rho + C*(n/n0)^2;
    else
      // Flow characteristics when check valve is closed
      (outlet.p - inlet.p)/rho = C*(n/n0)^2 - s*1e3;
      m_flow_1 = 0;
    end if;
    
    // Power consumption  
    P = D*(n^2)*(m_flow_1/rho) - E*n*((m_flow_1/rho)^2) + F*(n^2) 
      "Power consumption (single pump)";
    Ptot = P*Np "Power consumption (total)";
    // Hydraulic power
    Phyd = P*etaMech "Hydraulic power transferred to the fluid (single pump)";
    eta = ((outlet.p - inlet.p)*m_flow_1/rho)/(P + P_eps) 
      "Hydraulic efficiency";
    
    // Fluid properties
    fluid.p = inlet.p;
    rho = fluid.d;
    Tin = fluid.T;
    
    // Mass and energy balances
    inlet.m_flow + outlet.m_flow = 0 "Mass balance";
    inlet.mXi_flow + outlet.mXi_flow = zeros(Medium.nXi) 
      "Substance mass balance";
    if M > 0 then
      M * der(h) = m_flow_1*(inlet.h - outlet.h) + Phyd 
        "Energy balance (density variations neglected)";
    else
      0 = m_flow_1*(inlet.h - outlet.h) + Phyd "Energy balance";
    end if;
    h = fluid.h;
    inlet.H_flow=semiLinear(inlet.m_flow,inlet.h,fluid.h);
    outlet.H_flow=semiLinear(outlet.m_flow,outlet.h,fluid.h);
    
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
        Polygon(points=[-40,-22; -60,-58; 60,-58; 40,-22; -40,-22],
            style(pattern=0, fillColor=74)),
        Ellipse(extent=[-60,82; 60,-38],   style(gradient=3)),
        Polygon(points=[-30,54; -30,-6; 48,22; -30,54],     style(
            pattern=0,
            gradient=2,
            fillColor=7)),
        Text(extent=[-100,-62; 100,-88],   string="%name")),
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
  
  model Pump "Centrifugal pump with ideally controlled speed" 
    extends PumpBase;
    import Modelica.SIunits.Conversions.NonSIunits.*;
    parameter AngularVelocity_rpm n_const=n0 "Constant rotational speed";
    Modelica.Blocks.Interfaces.RealInput in_n "RPM" 
      annotation (extent=[-36, 70; -16, 90], rotation=-90);
  equation 
      n = in_n "Rotational speed";
    if cardinality(in_n)==0 then
      in_n = n_const "Rotational speed provided by parameter";
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
    extends Icons.Water.PumpMech;
    Angle phi "Shaft angle";
    AngularVelocity omega "Shaft angular velocity";
    Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft 
    annotation (extent=[80,4; 110,32]);
  equation 
    phi = shaft.phi;
    omega = der(phi);
    n = Modelica.SIunits.Conversions.to_rpm(omega);
    P = omega*shaft.tau;
    
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
    parameter Real Kv(unit="kg/(s.Pa)") "Hydraulic conductance at full opening";
    Modelica.Blocks.Interfaces.RealInput opening 
    annotation (extent=[-20, 60; 20, 100], rotation=-90);
  equation 
    m_flow = Kv*opening*dp;
    
  annotation (
    Icon(Text(extent=[-100, -40; 100, -74], string="%name"),
        Line(points=[0,54; 0,-6],  style(
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
    /* Do we need them?
  Medium.BaseProperties medium_a(h=feedwater.h, p=feedwater.p) 
    "Medium in feedwater";
  Medium.BaseProperties medium_b(h=steam.h, p=steam.p) "Medium in steam";
  */
    FluidPort_a feedwater(redeclare package Medium = Medium) 
      annotation (extent=[-120, -10; -100, 10]);
    FluidPort_b steam(redeclare package Medium = Medium) 
      annotation (extent=[120, -10; 100, 10]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
      annotation (extent=[-10, -120; 10, -100]);
    Modelica.Blocks.Interfaces.RealOutput V(
      redeclare type SignalType = Volume) "liquid volume (level)" 
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
    parameter Mass m_D=300e3 "mass of surrounding drum metal";
    parameter SpecificHeatCapacity cp_D=500 
      "specific heat capacity of drum metal";
    parameter Volume V_t=100 "total volume inside drum";
    parameter Pressure p_start=from_bar(1) "initial pressure";
    parameter Volume V_start=67 "initial liquid volume";
    
  /* why protected ?
protected 
*/
    
    Medium.AbsolutePressure p(start=p_start, fixed = true, stateSelect=StateSelect.prefer) 
      "pressure inside drum boiler";
    Medium.Temperature T "temperature inside drum boiler";
    Volume V_v "volume of vapour phase";
    Volume V_l(start=V_start, fixed = true, stateSelect=StateSelect.prefer) 
      "volumes of liquid phase";
    Medium.SpecificEnthalpy h_v=Medium.dewEnthalpy(sat) 
      "specific enthalpy of vapour";
    Medium.SpecificEnthalpy h_l=Medium.bubbleEnthalpy(sat) 
      "specific enthalpy of liquid";
    Medium.Density rho_v=Medium.dewDensity(sat) "density in vapour phase";
    Medium.Density rho_l=Medium.bubbleDensity(sat) "density in liquid phase";
    Mass m "total mass of drum boiler";
    Energy U "internal energy";
    Medium.Temperature T_D=heatPort.T "temperature of drum";
    HeatFlowRate q_F=heatPort.Q_flow "heat flow rate from furnace";
    SpecificEnthalpy h_W=feedwater.h "feed water enthalpy";
    SpecificEnthalpy h_S=steam.h "steam enthalpy";
    MassFlowRate qm_W=feedwater.m_flow "feed water mass flow rate";
    MassFlowRate qm_S=steam.m_flow "steam mass flow rate";
  equation 
    // balance equations  
    m = rho_v*V_v + rho_l*V_l; // why m_D? The balance equation is only applied to the fluid!
    U = rho_v*V_v*h_v + rho_l*V_l*h_l - p*V_t + m_D*cp_D*T_D;
    der(m) = qm_W + qm_S;
    der(U) = q_F + qm_W*h_W + qm_S*h_S;
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
    
    // thermal stress
    // Why 60? Parameter needed!!!
    sigma_D = -60*der(T_D);
    
    // liquid level 
    V = V_l;
  end Evaporator;
  
  
end Components;
