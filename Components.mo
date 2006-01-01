package Components "Basic components for fluid models" 
  
  model MixingVolume 
    "Mixing volume with inlet and outlet ports (flow reversal is allowed)" 
    import Modelica_Fluid.Types.InitTypes.*;
    import Modelica.Constants.*;
    extends Interfaces.PartialInitializationParameters;
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
    parameter SI.Volume V "Volume";
    parameter Boolean allowFlowReversal = true 
      "Flow reversal at the ports is allowed by the equations" annotation(Dialog(tab="Advanced"));
    Interfaces.FluidPort_a port_a(redeclare package Medium = Medium,
                                  m_flow(min=if allowFlowReversal then -inf else 0)) 
      "Fluid inlet port" annotation (extent=[-112,-10; -92,10]);
    Interfaces.FluidPort_b port_b(redeclare package Medium = Medium,
                                  m_flow(max=if allowFlowReversal then +inf else 0)) 
      "Fluid outlet port" annotation (extent=[90,-10; 110,10]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
      "Thermal port" 
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
This model describes pressure losses due to friction in a pipe. It is assumed that no mass or energy is stored in the pipe. 
The details of the pipe friction model are described
<a href=\"Modelica://Modelica_Fluid.Utilities.PipeFriction\">here</a>.
</p>
</html>"));
  equation 
    if frictionType == Modelica_Fluid.Types.FrictionTypes.DetailedFriction then
       if from_dp then
          d = if dp >= 0 then medium_a.d else medium_b.d;
          eta = if dp >= 0 then Medium.dynamicViscosity(medium_a) else 
                               Medium.dynamicViscosity(medium_b);
       else
          d = if m_flow >= 0 then medium_a.d else medium_b.d;
          eta = if m_flow >= 0 then Medium.dynamicViscosity(medium_a) else 
                               Medium.dynamicViscosity(medium_b);
       end if;
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
  
  model StaticHead 
    "Models the static head between two ports at different heights" 
    extends Modelica_Fluid.Interfaces.PartialTwoPortTransport;
    parameter SI.Height H_b_a "Height of port b over port a";
    parameter SI.Acceleration g = Modelica.Constants.g_n "Gravity acceleration"
                                                                                 annotation(Dialog(tab="Advanced"));
    Medium.Density d "Fluid density";
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
This model describes the static head due to the relative height between the two connectors. No mass, energy and momentum storage, and no pressure drop due to friction are considered.
</p>
</html>", revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Added to Modelica_Fluid</li>
</ul>
</html>"));
  equation 
   d = if dp > 0 then medium_a.d else medium_b.d;
   port_a.p = port_b.p + H_b_a*g*d;
  end StaticHead;
  
  model Tank "Tank with one bottom inlet/outlet port" 
    import Modelica_Fluid.Types.InitTypes.*;
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
    parameter SI.Area area "Tank area";
    parameter SI.Area pipeArea "Area of outlet pipe";
    parameter SI.Volume V0 = 0 "Volume of the liquid when the level is zero";
    parameter SI.Height H0 = 0 
      "Height of zero level reference over the bottom port";
    parameter SI.Acceleration g = Modelica.Constants.g_n 
      "Acceleration of gravity";
    parameter Medium.AbsolutePressure p_ambient=101325 "Tank surface pressure";
    parameter Types.InitTypes.Temp initOption = NoInit "Initialization option" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
      "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate = true);
    parameter Medium.Temperature T_start=
      if use_T_start then 293.15 else Medium.T_phX(p_ambient,h_start,X_start) 
      "Start value of temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.h_pTX(p_ambient, T_start, X_start[1:Medium.nXi]) else 1e4 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.reference_X 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter SI.Height level_start(min=0) "Initial tank level" 
      annotation(Dialog(tab="Initialization"));
    
    Interfaces.FluidPort_b port(redeclare package Medium = Medium,
                                m_flow(start=0), mXi_flow(each start=0)) 
      annotation (extent=[-10, -120; 10, -100], rotation=90);
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      Xi(start=X_start[1:Medium.nXi]));
    
    SI.Height level(start=level_start,stateSelect=StateSelect.prefer) 
      "Height of tank level over the zero reference";
    Medium.AbsolutePressure p_bottom "Pressure at bottom of tank";
    SI.Energy U "Internal energy of tank volume";
    SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
    Real m(quantity=Medium.mediumName, unit="kg", stateSelect=StateSelect.never) 
      "Mass of tank volume";
    Real mX[Medium.nX](quantity=Medium.substanceNames, each unit="kg") 
      "Component masses of the independent substances";
  equation 
    medium.p = p_ambient;
    V = area*level+V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mX = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";
    
    // Mass balance
    der(m) = port.m_flow;
    der(mX) = port.mXi_flow;
    
    // Momentum balance
    p_bottom = (medium.d*g*(level+H0)) + p_ambient;
    port.p = p_bottom - smooth(2,noEvent(if port.m_flow < 0 then port.m_flow^2/(2*medium.d*pipeArea^2) else 0));
    
    // Energy balance
    if Medium.singleState then
      der(U) = port.H_flow "Mechanical work is neglected";
    else
      der(U) = port.H_flow - p_ambient*der(V);
    end if;
    
    /* Handle reverse and zero flow */
    port.H_flow = semiLinear(port.m_flow, port.h, medium.h);
    port.mXi_flow = semiLinear(port.m_flow, port.Xi, medium.X);
    
  initial equation 
    if initOption == NoInit then
      // no initial equations
    elseif initOption == InitialValues then
      level = level_start;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initOption == SteadyStateHydraulic then
      der(level) = 0;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    else
      assert(false, "Unsupported initialization option");
    end if;
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
          style(color=0),
          string="level_start")),
      Documentation(info="<HTML>
<p>
This is a simplified model of a tank. The top part is open to the environment at the fixed pressure <tt>p_ambient</tt>. Heat transfer to the environment and to the tank walls is neglected.
The tank is filled with a single or multiple-substance liquid, assumed to have uniform temperature and mass fractions.
<p>The geometry of the tank is specified by the following parameters: <tt>V0</tt> is the volume of the liquid when the level is at the zero reference; <tt>area</tt> is the cross-sectional area of the tank; <tt>H0</tt> is the height of the zero-reference level plane over the port connector. It is thus possible to model rounded-bottom tanks, as long as they have a cylindrical shape in the range of operating levels.
<p>The tank can be initialized with the following options:
<ul>
<li>NoInit: no explicit initial conditions
<li>InitialValues: initial values of temperature (or specific enthalpy), composition and level are specified
<li>SteadyStateHydraulic: initial values of temperature (or specific enthalpy) and composition are specified; the initial level is determined so that levels and pressure are at steady state.
</ul>
Full steady state initialization is not supported, because the corresponding intial equations for temperature/enthalpy are undetermined (the flow rate through the port at steady state is zero). 
</p>
</HTML>",   revisions="<html>
<ul>
<li><i>1 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from a previous version of Modelica_Fluid</li>
</ul>
</html>"),
      Diagram);
  end Tank;
  
  model ValveIncompressible "Valve for (almost) incompressible fluids" 
    extends Interfaces.PartialValve;
    import Modelica_Fluid.Types.CvTypes;
  annotation (
    Icon(Text(extent=[-100, -40; 100, -80], string="%name")),
    Diagram,
    Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluids. <p>
Extends the <tt>Interfaces.PartialValve</tt> model (see the corresponding documentation for common valve features).
<p>This model can be used with any low compressibility fluids, such as liquids or gases at very low pressure drops.
</html>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
  initial equation 
    if CvData == CvTypes.OpPoint then
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*sqrt(d_nom)*sqrtR(dp_nom) 
        "Determination of Av by the operating point";
    end if;
    
  equation 
    if CheckValve then
      m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*
               smooth(0,if dp>=0 then sqrtR(dp) else 0);
    else
      m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*sqrtR(dp);
    end if;
  end ValveIncompressible;
  
  model ValveVaporizing 
    "Valve for possibly vaporizing (almost) incompressible fluids, accounts for choked flow conditions" 
    import Modelica_Fluid.Types.CvTypes;
    extends Interfaces.PartialValve(
      redeclare replaceable package Medium = 
      Modelica.Media.Interfaces.PartialTwoPhaseMedium);
    parameter Real Fl_nom=0.9 "Liquid pressure recovery factor";
    replaceable function FlCharacteristic = 
        Modelica_Fluid.Types.ValveCharacteristics.one 
      extends Modelica_Fluid.Types.ValveCharacteristics.baseFun 
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
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluid, with possible choked flow conditions. <p>
Extends the <tt>Interfaces.PartialValve</tt> model (see the corresponding documentation for common valve features).<p>
The model operating range includes choked flow operation, which takes place for low outlet pressures due to flashing in the vena contracta; otherwise, non-choking conditions are assumed.
<p>This model can be used with two-phase medium models, to describe the liquid and (possible) two-phase conditions.
<p>The default liquid pressure recovery coefficient <tt>Fl</tt> is constant and given by the parameter <tt>Fl_nom</tt>. The relative change (per unit) of the recovery coefficient can be specified as a given function of the valve opening by replacing the <tt>FlCharacteristic</tt> function.
</HTML>",
        revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
  initial equation 
    if CvData == CvTypes.OpPoint then
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*sqrt(d_nom)*sqrtR(dp_nom) 
        "Determination of Av by the operating point";
    end if;
  equation 
    pin = port_a.p;
    pout = port_b.p;
    pv = Medium.saturationPressure(T);
    Ff = 0.96 - 0.28*sqrt(pv/Medium.fluidConstants[1].criticalPressure);
    Fl = Fl_nom*FlCharacteristic(stemPosition);
    dpEff = if pout < (1 - Fl^2)*pin + Ff*Fl^2*pv then 
              Fl^2*(pin - Ff*pv) else dp 
      "Effective pressure drop, accounting for possible choked conditions";
    if CheckValve then
       m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*
           (if dpEff>=0 then sqrtR(dpEff) else 0);
     else
       m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*sqrtR(dpEff);
    end if;
  end ValveVaporizing;
  
  model ValveCompressible 
    "Valve for compressible fluids, accounts for choked flow conditions" 
    extends Interfaces.PartialValve;
    import Modelica_Fluid.Types.CvTypes;
    parameter Real Fxt_full=0.5 "Fk*xt critical ratio at full opening";
    replaceable function xtCharacteristic = 
        Modelica_Fluid.Types.ValveCharacteristics.one 
      extends Modelica_Fluid.Types.ValveCharacteristics.baseFun 
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
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, compressible fluid. <p>
Extends the <tt>Interfaces.PartialValve</tt> model (see the corresponding documentation for common valve features).
<p>The product Fk*xt is given by the parameter <tt>Fxt_full</tt>, and is assumed constant by default. The relative change (per unit) of the xt coefficient with the valve opening can be specified by replacing the <tt>xtCharacteristic</tt> function.
</HTML>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
  initial equation 
    if CvData == CvTypes.OpPoint then
      // Determination of Av by the nominal operating point conditions
      Fxt_nom = Fxt_full*xtCharacteristic(stemPosition_nom);
      x_nom = dp_nom/p_nom;
      xs_nom = smooth(0, if x_nom > Fxt_nom then Fxt_nom else x_nom);
      Y_nom = 1 - abs(xs_nom)/(3*Fxt_nom);
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*Y_nom*sqrt(d_nom)*sqrtR(p_nom*xs_nom);
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
      m_flow = flowCharacteristic(stemPosition)*Av*Y*sqrt(d)*
        smooth(0,if xs>=0 then sqrtR(p*xs) else 0);
    else
      m_flow = flowCharacteristic(stemPosition)*Av*Y*sqrt(d)*sqrtR(p*xs);
    end if;
  end ValveCompressible;
  
  model ValveLinear "Valve for water/steam flows with linear pressure drop" 
    extends Interfaces.PartialTwoPortTransport;
    parameter Types.HydraulicConductance Kv 
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
<p>This very simple model provides a pressure drop which is proportional to the flowrate and to the <tt>opening</tt> signal, without computing any fluid property.
<p>A medium model must be nevertheless be specified, so that the fluid ports can be connected to other components using the same medium model.
</HTML>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
  end ValveLinear;
  
  model Pump "Centrifugal pump with ideally controlled speed" 
    extends Interfaces.PartialPump;
    import Modelica.SIunits.Conversions.NonSIunits.*;
    parameter AngularVelocity_rpm N_const = N_nom "Constant rotational speed";
    Modelica.Blocks.Interfaces.RealInput N_in "Prescribed rotational speed" 
      annotation (extent=[-36,34; -16,54],   rotation=-90);
  equation 
      N = N_in "Rotational speed";
    if cardinality(N_in)==0 then
      N_in = N_const "Rotational speed provided by parameter";
    end if;
    annotation (
      Icon(
        Text(extent=[-58,58; -30,38], string="n")),
      Diagram,
      Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>Np</tt> pumps in parallel) with controlled speed, either fixed or provided by an external signal.
<p>The model extends <tt>PartialPump</tt>
<p>If the <tt>N_in</tt> input connector is wired, it provides rotational speed of the pumps (rpm); otherwise, a constant rotational speed equal to <tt>n_const</tt> (which can be different from <tt>N_nom</tt>) is assumed.</p>
</HTML>",
        revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end Pump;
  
  model PumpShaft "Centrifugal pump with mechanical connector for the shaft" 
    extends Interfaces.PartialPump;
    SI.Angle phi "Shaft angle";
    SI.AngularVelocity omega "Shaft angular velocity";
    Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft 
    annotation (extent=[80,4; 110,32]);
  equation 
    phi = shaft.phi;
    omega = der(phi);
    N = Modelica.SIunits.Conversions.to_rpm(omega);
    W_single = omega*shaft.tau;
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
<p>The model extends <tt>PartialPump</tt>
 </HTML>",
       revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end PumpShaft;
  
  model Evaporator 
    "Simple Evaporator with two states, see Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378" 
    import Modelica.SIunits.Conversions.*;
    import Modelica_Fluid.Types.InitTypes.*;
    import Modelica.Constants.*;
    replaceable package Medium = 
        Modelica.Media.Interfaces.PartialTwoPhaseMedium 
      extends Modelica.Media.Interfaces.PartialTwoPhaseMedium "Medium model" 
        annotation (choicesAllMatching=true);
    parameter SI.Mass m_D "mass of surrounding drum metal";
    parameter Medium.SpecificHeatCapacity cp_D 
      "specific heat capacity of drum metal";
    parameter SI.Volume V_t "total volume inside drum";
    parameter Types.InitTypes.Temp initOption = NoInit "Initialization option" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_start = Medium.reference_p 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization"));
    parameter SI.Volume V_l_start = V_t/2 
      "Start value of liquid volumeStart value of volume" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean allowFlowReversal = false 
      "Flow reversal at the ports is allowed by the equations" annotation(Dialog(tab="Advanced"));
    
    Interfaces.FluidPort_a feedwater(redeclare package Medium = Medium,
                             m_flow(min = if allowFlowReversal then -inf else 0)) 
      annotation (extent=[-120, -10; -100, 10]);
    Interfaces.FluidPort_b steam(redeclare package Medium = Medium,
                             m_flow(max = if allowFlowReversal then +inf else 0)) 
      annotation (extent=[120, -10; 100, 10]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
      annotation (extent=[-10, -120; 10, -100]);
    Modelica.Blocks.Interfaces.RealOutput V(
      redeclare type SignalType = SI.Volume) "liquid volume" 
      annotation (extent=[30, 100; 50, 120], rotation=90);
    
    Medium.SaturationProperties sat 
      "State vector to compute saturation properties";
    Medium.AbsolutePressure p(start=p_start, stateSelect=StateSelect.prefer) 
      "pressure inside drum boiler";
    Medium.Temperature T "temperature inside drum boiler";
    SI.Volume V_v "volume of vapour phase";
    SI.Volume V_l(start=V_l_start, stateSelect=StateSelect.prefer) 
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
    
    // Check that two-phase equilibrium is actually possible
    assert(p<Medium.fluidConstants[1].criticalPressure-10000,
           "Evaporator model requires subcritical pressure");
  initial equation 
    // Initial conditions
    if initOption == NoInit then
      // no initial equations
    elseif initOption == InitialValues then
     p = p_start;
     V_l = V_l_start;
    elseif initOption == SteadyState then
      der(p) = 0;
      der(V_l) = 0;
    else
      assert(false, "Unsupported initialization option");
    end if;
    
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
        Line(points=[40, 99; 40, 60])),
      Documentation(revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     Initialization options fixed</li>
<li><i>6 Sep 2005</i><br>
    Model by Ruediger Franke modified after the 45th Design Meeting</li>
</ul>
</html>", info="<html>
Model of a simple evaporator with two states. The model assumes two-phase equilibrium inside the component; saturated steam goes out of the steam outlet.
<p>
References: Astroem, Bell: Drum-boiler dynamics, Automatica 36, 2000, pp.363-378
</html>"));
  equation 
    
  end Evaporator;
  
end Components;
