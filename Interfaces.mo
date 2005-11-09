package Interfaces 
  "Interfaces for steady state and unsteady, mixed-phase, multi-substance, incompressible and compressible flow" 
  
  annotation (Documentation(info="<html>
</html>", revisions="<html>
<ul>
<li><i>Nov. 2, 2005</i>
       by Francesco Casella: restructured after 45th Design Meeting.</li>
<li><i>Nov. 20-21, 2002</i>
       by Hilding Elmqvist, Mike Tiller, Allan Watson, John Batteh, Chuck Newman,
       Jonas Eborn: Improved at the 32nd Modelica Design Meeting.
<li><i>Nov. 11, 2002</i>
       by Hilding Elmqvist, Martin Otter: improved version.</li>
<li><i>Nov. 6, 2002</i>
       by Hilding Elmqvist: first version.</li>
<li><i>Aug. 11, 2002</i>
       by Martin Otter: Improved according to discussion with Hilding
       Elmqvist and Hubertus Tummescheit.<br>
       The PortVicinity model is manually
       expanded in the base models.<br>
       The Volume used for components is renamed
       PartialComponentVolume.<br>
       A new volume model \"Fluid.Components.PortVolume\"
       introduced that has the medium properties of the port to which it is
       connected.<br>
       Fluid.Interfaces.PartialTwoPortTransport is a component
       for elementary two port transport elements, whereas PartialTwoPort
       is a component for a container component.</li>
</li>
</ul>
</html>"));
  
  extends Modelica.Icons.Library;
  import SI = Modelica.SIunits;
  
  connector FluidPort 
    "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)" 
    
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Medium model" annotation (choicesAllMatching=true);
    
    Medium.AbsolutePressure p "Pressure in the connection point";
    flow Medium.MassFlowRate m_flow 
      "Mass flow rate from the connection point into the component";
    
    Medium.SpecificEnthalpy h 
      "Specific mixture enthalpy in the connection point";
    flow Medium.EnthalpyFlowRate H_flow 
      "Enthalpy flow rate into the component (if m_flow > 0, H_flow = m_flow*h)";
    
    Medium.MassFraction Xi[Medium.nXi] 
      "Independent mixture mass fractions m_i/m in the connection point";
    flow Medium.MassFlowRate mXi_flow[Medium.nXi] 
      "Mass flow rates of the independent substances from the connection point into the component (if m_flow > 0, mXi_flow = m_flow*Xi)";
    
    Medium.ExtraProperty C[Medium.nC] 
      "properties c_i/m in the connection point";
    flow Medium.ExtraPropertyFlowRate mC_flow[Medium.nC] 
      "Flow rates of auxiliary properties from the connection point into the component (if m_flow > 0, mC_flow = m_flow*C)";
    
  end FluidPort;
  
  connector FluidPort_a "Fluid connector with filled icon" 
    extends FluidPort;
    annotation (Diagram(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
               fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
               fillColor=69)), Text(extent=[-88, 206; 112, 112], string="%name")),
         Icon(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
              fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
              fillColor=69)), Text(
          extent=[-126, 160; 130, 104],
          string="%name",
          style(
            color=0,
            fillColor=69,
            fillPattern=1))));
  end FluidPort_a;
  
  connector FluidPort_b "Fluid connector with outlined icon" 
    extends FluidPort;
    annotation (Diagram(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
               fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
               fillColor=69)), Ellipse(extent=[-80, 80; 80, -80], style(color=69,
               fillColor=7)), Text(extent=[-88, 192; 112, 98], string="%name")),
         Icon(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
              fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
              fillColor=69)), Ellipse(extent=[-80, 80; 80, -80], style(color=69,
               fillColor=7)), Text(
          extent=[-126, 160; 130, 104],
          string="%name",
          style(
            color=0,
            fillColor=69,
            fillPattern=1))));
  end FluidPort_b;
  
  connector HeatPort = Modelica.Thermal.HeatTransfer.Interfaces.HeatPort 
    "Thermal port for 1-dim. heat transfer";
  
  connector HeatPort_a = Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a 
    "Thermal port for 1-dim. heat transfer (filled rectangular icon)" 
    annotation (Icon);
  
  connector HeatPort_b = Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b 
    "Thermal port for 1-dim. heat transfer (unfilled rectangular icon)" 
    annotation (Icon(
                    Text(
          extent=[-98,196; 102,102],
          string="%name",
          style(color=42))));
  
  partial model PartialInitializationParameters 
    "Define parameter menu to initialize medium in component that has one medium model" 
    import Modelica_Fluid.Types.InitTypes.*;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Medium model" 
       annotation (choicesAllMatching=true);
    parameter Types.InitTypes.Temp initOption = NoInit "Initialization option" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_start = Medium.reference_p 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
      "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate=true);
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
  end PartialInitializationParameters;
  
  partial model PartialSource 
    "Partial component source with one fluid connector" 
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium model within the source" 
       annotation (choicesAllMatching=true);
    FluidPort_b port(redeclare package Medium = Medium) 
      annotation (extent=[100, -10; 120, 10], rotation=0);
    Medium.BaseProperties medium "Medium in the source";
  equation 
    port.p = medium.p;
    port.H_flow = semiLinear(port.m_flow, port.h, medium.h);
    port.mXi_flow = semiLinear(port.m_flow, port.Xi, medium.Xi);
    annotation (Documentation(info="<html>
<p>
Partial component to model the <b>volume interface</b> of a <b>source</b>
component, such as a mass flow source. The essential
features are:
</p>
<ul>
<li> The pressure in the connection port (= port.p) is identical to the
     pressure in the volume (= medium.p).</li>
<li> The enthalpy flow rate (= port.H_flow) and the mass flow rates of the
     substances (= port.mX_flow) depend on the direction of the mass flow rate.</li>
</ul>
</html>"));
  end PartialSource;
  
  partial model PartialTwoPortTransport 
    "Partial element transporting fluid between two ports without storing mass or energy" 
    import Modelica.SIunits.*;
    import Modelica.Constants.*;
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component"  annotation (
        choicesAllMatching =                                                                            true);
    parameter Boolean allowFlowReversal = true 
      "Flow reversal at the ports is allowed by the equations"  annotation(Dialog(tab="Advanced"));
    
    FluidPort_a port_a(redeclare package Medium = Medium,
                       m_flow(min=if allowFlowReversal then -inf else 0)) 
      annotation (extent=[-120, -10; -100, 10]);
    FluidPort_b port_b(redeclare package Medium = Medium,
                       m_flow(max=if allowFlowReversal then +inf else 0)) 
      annotation (extent=[120, -10; 100, 10]);
    Medium.BaseProperties medium_a "Medium properties in port_a";
    Medium.BaseProperties medium_b "Medium properties in port_b";
    Medium.MassFlowRate m_flow 
      "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
    Pressure dp(start=0) "Pressure difference between port_a and port_b";
    
    annotation (
      Coordsys(grid=[1, 1], component=[20, 20]),
      Diagram,
      Documentation(info="<html>
<p>
This component transports fluid between its two ports, without
storing mass or energy. Reversal and zero mass flow rate is taken
care of, for details see definition of built-in operator semiLinear().
<p>
When using this partial component, an equation for the momentum
balance has to be added by specifying a relationship
between the pressure drop <tt>dp</tt> and the mass flow rate <tt>m_flow</tt>.
</p>
</html>"));
  equation 
    // Properties in the ports
    port_a.p   = medium_a.p;
    port_a.h   = medium_a.h;
    port_a.Xi = medium_a.Xi;
    port_b.p   = medium_b.p;
    port_b.h   = medium_b.h;
    port_b.Xi = medium_b.Xi;
    
    /* Handle reverse and zero flow */
    port_a.H_flow   = semiLinear(port_a.m_flow, port_a.h,  port_b.h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, port_b.Xi);
    
    /* Energy, mass and substance mass balance */
    port_a.H_flow + port_b.H_flow = 0;
    port_a.m_flow + port_b.m_flow = 0;
    port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
    
    // Design direction of mass flow rate
    m_flow = port_a.m_flow;
    
    // Pressure difference between ports
    dp = port_a.p - port_b.p;
    
  end PartialTwoPortTransport;
  
  partial model PartialAbsoluteSensor 
    "Partial component to model a sensor that measures a potential variable" 
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor" annotation (
        choicesAllMatching =                                                                        true);
    FluidPort_a port(redeclare package Medium = Medium) 
      annotation (extent=[-10, -120; 10, -100], rotation=90);
    
    annotation (Documentation(info="<html>
<p>
Partial component to model an <b>absolute sensor</b>. Can be used for pressure sensor models.
Use for other properties such as temperature or density is discouraged, because the enthalpy at the connector can have different meanings, depending on the connection topology. Use <tt>PartialFlowSensor</tt> instead.
as signal.
</p>
</html>"));
  equation 
    port.m_flow = 0;
    port.H_flow = 0;
    port.mXi_flow = zeros(Medium.nXi);
  end PartialAbsoluteSensor;
  
  partial model PartialFlowSensor 
    "Partial component to model sensors that measure flow properties" 
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
        choicesAllMatching = true);
    Medium.SpecificEnthalpy h "enthalpy in flow";
    Medium.MassFraction[Medium.nXi] Xi "flow composition";
    
    FluidPort_a port_a(redeclare package Medium = Medium) 
      annotation (extent=[-120, -10; -100, 10]);
    FluidPort_b port_b(redeclare package Medium = Medium) 
      annotation (extent=[120, -10; 100, 10]);
    
    annotation (Documentation(info="<html>
<p>
Partial component to model a <b>sensor</b> that measures any intensive properties
of a flow, e.g., to get temperature or density in the flow
between fluid connectors.<br>
The model includes zero-volume balance equations. Sensor models inheriting from
this partial class should add a medium instance to calculate the measured property.
</p>
</html>"));
  equation 
    port_a.p   = port_b.p;
    // Local *zero-volume* enthalpy and composition
    port_a.H_flow = semiLinear(port_a.m_flow,port_a.h,h);
    port_b.H_flow = semiLinear(port_b.m_flow,port_b.h,h);
    port_a.mXi_flow = semiLinear(port_a.m_flow,port_a.Xi,Xi);
    port_b.mXi_flow = semiLinear(port_b.m_flow,port_b.Xi,Xi);
    // Static balances
    0 = port_a.m_flow + port_b.m_flow;
    0 = port_a.H_flow + port_b.H_flow;
    zeros(Medium.nXi) = port_a.mXi_flow + port_b.mXi_flow;
  end PartialFlowSensor;
  
protected 
  partial model PartialRelativeSensor 
    "Partial component to model a sensor that measures the difference of effort variables at two ports" 
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
        choicesAllMatching =                                                                         true);
    
    FluidPort_a port_a(redeclare package Medium = Medium) 
      annotation (extent=[-120, -10; -100, 10]);
    FluidPort_b port_b(redeclare package Medium = Medium) 
      annotation (extent=[120, -10; 100, 10]);
    
    annotation (Documentation(info="<html>
<p>
Partial component to model a <b>sensor</b> that measures
the <b>difference between two effort variables</b>, e.g. to obtain the temperature difference 
between fluid connectors.
</p>
</html>"));
  equation 
    port_a.m_flow = 0;
    port_a.H_flow = 0;
    port_a.mXi_flow = zeros(Medium.nXi);
    
    port_b.m_flow = 0;
    port_b.H_flow = 0;
    port_b.mXi_flow = zeros(Medium.nXi);
  end PartialRelativeSensor;
  
public 
  partial model PartialValve "Base model for valves" 
    import Modelica_Fluid.Types.CvTypes;
    
    parameter Medium.AbsolutePressure pin_start = p_nom 
      "Start value of inlet pressure" 
      annotation(Dialog(tab = "Initialization"));
    
    extends Interfaces.PartialTwoPortTransport(
      medium_a(p(start=pin_start), T(start=T_start),
               h(start=h_start),   Xi(start=X_start[1:Medium.nXi])),
      medium_b(p(start=pout_start), T(start=T_start),
               h(start=h_start),   Xi(start=X_start[1:Medium.nXi])));
    
    parameter CvTypes.Temp CvData = CvTypes.Av "Selection of flow coefficient" 
       annotation(Dialog(group = "Flow Coefficient"));
    parameter SI.Area Av(fixed = if CvData==CvTypes.Av then true else false,
                         start = m_flow_nom/(sqrt(d_nom*dp_nom))*
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
    parameter Medium.AbsolutePressure p_nom "Nominal inlet pressure" 
      annotation(Dialog(group="Nominal operating point"));
    parameter SI.Pressure dp_nom "Nominal pressure drop" 
      annotation(Dialog(group="Nominal operating point"));
    parameter Medium.MassFlowRate m_flow_nom "Nominal mass flowrate" 
      annotation(Dialog(group="Nominal operating point"));
    parameter Medium.Density d_nom = 1000 "Nominal inlet density" 
      annotation(Dialog(group="Nominal operating point"));
    parameter Real stemPosition_nom = 1 "Nominal stem position" 
      annotation(Dialog(group="Nominal operating point"));
    parameter Boolean CheckValve=false "Reverse flow stopped";
    
    replaceable function flowCharacteristic = 
        Modelica_Fluid.Types.ValveCharacteristics.linear 
      extends Modelica_Fluid.Types.ValveCharacteristics.baseFun 
      "Inherent flow characteristic" 
      annotation(choicesAllMatching=true);
    parameter Medium.AbsolutePressure pout_start = p_nom-dp_nom 
      "Start value of outlet pressure" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
      "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate = true);
    parameter Medium.Temperature T_start=
      if use_T_start then 293.15 else Medium.T_phX(pin_start,h_start,X_start) 
      "Start value of inlet temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.h_pTX(pin_start, T_start, X_start[1:Medium.nXi]) else 1e4 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.reference_X 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    
    parameter Real delta=0.01 "Regularisation factor" annotation(Dialog(tab="Advanced"));
    
    Modelica.Blocks.Interfaces.RealInput stemPosition 
      "Stem position in the range 0-1" annotation (extent=[-10,70; 10,90],    rotation=-90);
    
    Medium.Density d "Density at port a";
    Medium.Temperature T "Temperature at port a";
  protected 
    function sqrtR = Utilities.regRoot(delta = delta*dp_nom);
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
<p>This is the base model for the <tt>ValveIncompressible</tt>, <tt>ValveVaporizing</tt>, and <tt>ValveCompressible</tt> valve models. The model is based on the IEC 534 / ISA S.75 standards for valve sizing.
<p>The model optionally supports reverse flow conditions (assuming symmetrical behaviour) or check valve operation, and has been suitably modified to avoid numerical singularities at zero pressure drop. 
<p><b>Modelling options</b></p>
<p>The following options are available to specify the valve flow coefficient in fully open conditions:
<ul><li><tt>CvData = Modelica_Fluid.Types.CvTypes.Av</tt>: the flow coefficient is given by the metric <tt>Av</tt> coefficient (m^2).
<li><tt>CvData = Modelica_Fluid.Types.CvTypes.Kv</tt>: the flow coefficient is given by the metric <tt>Kv</tt> coefficient (m^3/h).
<li><tt>CvData = Modelica_Fluid.Types.CvTypes.Cv</tt>: the flow coefficient is given by the US <tt>Cv</tt> coefficient (USG/min).
<li><tt>CvData = Modelica_Fluid.Types.CvTypes.OpPoint</tt>: the flow is computed from the nominal operating point specified by <tt>p_nom</tt>, <tt>dp_nom</tt>, <tt>m_flow_nom</tt>, <tt>d_nom</tt>, <tt>stemPosition_nom</tt>.
</ul>
<p>The nominal pressure drop <tt>dp_nom</tt> must always be specified; to avoid numerical singularities, the flow characteristic is modified for pressure drops less than <tt>b*dp_nom</tt> (the default value is 1% of the nominal pressure drop). Increase this parameter if numerical problems occur in valves with very low pressure drops.
<p>If <tt>CheckValve</tt> is true, then the flow is stopped when the outlet pressure is higher than the inlet pressure; otherwise, reverse flow takes place.
<p>The inherent flow characteristic <tt>flowCharacteristic</tt>, linear by default, can be replaced by any user-defined function (e.g. equal percentage, quick opening, etc.).
</HTML>",
        revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
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
    d = medium_a.d;
  end PartialValve;
  
  partial model PartialPump "Base model for centrifugal pumps" 
    import Modelica.SIunits.Conversions.NonSIunits.*;
    import Modelica.Constants.*;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Medium model" annotation(choicesAllMatching=true);
    Medium.BaseProperties fluid(p(start=pin_start),h(start=h_start)) 
      "Fluid properties at the inlet";
    replaceable package SatMedium = 
        Modelica.Media.Interfaces.PartialTwoPhaseMedium 
      "Saturated medium model (required only for NPSH computation)";
    replaceable function flowCharacteristic = 
        Modelica_Fluid.Types.PumpCharacteristics.baseFlow 
      "Head vs. q_flow characteristic at nominal speed and density" 
      annotation(Dialog(group="Characteristics"), choicesAllMatching=true);
    parameter Boolean usePowerCharacteristic = false 
      "Use powerCharacteristic (vs. efficiencyCharacteristic)" 
       annotation(Dialog(group="Characteristics"));
    replaceable function powerCharacteristic = 
      Modelica_Fluid.Types.PumpCharacteristics.basePower 
      "Power consumption vs. q_flow at nominal speed and density" 
      annotation(Dialog(group="Characteristics", enable = usePowerCharacteristic),
                 choicesAllMatching=true);
    replaceable function efficiencyCharacteristic = 
      Modelica_Fluid.Types.PumpCharacteristics.constantEfficiency(eta_nom = 0.8) 
      extends Modelica_Fluid.Types.PumpCharacteristics.baseEfficiency 
      "Efficiency vs. q_flow at nominal speed and density" 
      annotation(Dialog(group="Characteristics",enable = not usePowerCharacteristic),
                 choicesAllMatching=true);
    parameter AngularVelocity_rpm N_nom = 1500 "Nominal rotational speed" 
      annotation(Dialog(group="Characteristics"));
    parameter Medium.Density d_nom = 1000 "Nominal fluid density" 
      annotation(Dialog(group="Characteristics"));
    parameter Integer Np_nom(min=1) = 1 "Nominal number of pumps in parallel";
    parameter SI.Mass M = 0 "Fluid mass inside the pump";
    parameter Boolean checkValve=true "Reverse flow stopped";
    parameter Boolean allowFlowReversal = true 
      "Flow reversal at the ports is allowed by the equations";
    parameter Boolean computeNPSHa=false "Compute NPSH Available at the inlet";
    parameter Medium.AbsolutePressure pin_start "Inlet Pressure Start Value" 
      annotation(Dialog(tab="Initialization"));
    parameter Medium.AbsolutePressure pout_start "Outlet Pressure Start Value" 
      annotation(Dialog(tab="Initialization"));
    parameter Boolean use_T_start = true 
      "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate = true);
    parameter Medium.Temperature T_start=
      if use_T_start then 293.15 else Medium.T_phX(pin_start,h_start,Medium.reference_X[1:Medium.nXi]) 
      "Start value of temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.h_pTX(pin_start, T_start, Medium.reference_X[1:Medium.nXi]) else 1e4 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter SI.MassFlowRate m_flow_start = 0 
      "Start value of mass flow rate (total)" 
      annotation(Dialog(tab="Initialization"));
    constant SI.Acceleration g=Modelica.Constants.g_n;
  //  parameter Choices.Init.Options.Temp initOpt=Choices.Init.Options.noInit 
  //    "Initialisation option";
    Modelica_Fluid.Interfaces.FluidPort_a inlet(redeclare package Medium = Medium,
        p(start=pin_start),
        m_flow(start = m_flow_start,
               min = if allowFlowReversal and not checkValve then -inf else 0)) 
    annotation (extent=[-100,-40; -60,0]);
    Modelica_Fluid.Interfaces.FluidPort_b outlet(redeclare package Medium = Medium,
        p(start=pout_start),
        m_flow(start = -m_flow_start,
               max = if allowFlowReversal and not checkValve then +inf else 0)) 
    annotation (extent=[40,12; 80,52]);
    SI.Pressure dp = outlet.p - inlet.p "Pressure increase";
    SI.Height head = dp/(d*g) "Pump head";
    Medium.Density d "Liquid density at the inlet";
    Medium.SpecificEnthalpy h_out(start=h_start) 
      "Enthalpy of the liquid flowing out of the pump";
    Medium.Temperature Tin "Liquid inlet temperature";
    SI.MassFlowRate m_flow = inlet.m_flow "Mass flow rate (total)";
    SI.MassFlowRate m_flow_single = m_flow/Np "Mass flow rate (single pump)";
    SI.VolumeFlowRate q_flow = m_flow/d "Volume flow rate (total)";
    SI.VolumeFlowRate q_flow_single = q_flow/Np 
      "Volume flow rate (single pump)";
    AngularVelocity_rpm N "Shaft rotational speed";
    Integer Np(min=1) "Number of pumps in parallel";
    SI.Power W_single "Power Consumption (single pump)";
    SI.Power W_tot = W_single*Np "Power Consumption (total)";
    constant SI.Power W_eps=1e-8 
      "Small coefficient to avoid numerical singularities in efficiency computations";
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
      in_Np = Np_nom "Number of pumps selected by parameter";
    end if;
    
    // Flow equations
    if noEvent(s > 0 or (not checkValve)) then
      // Flow characteristics when check valve is open
      q_flow_single = s;
      head = (N/N_nom)^2*flowCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-6)));
    else
      // Flow characteristics when check valve is closed
      head = (N/N_nom)^2*flowCharacteristic(0) - s;
      q_flow_single = 0;
    end if;
    
    // Power consumption  
    if usePowerCharacteristic then
      W_single = (N/N_nom)^3*(d/d_nom)*powerCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-6))) 
        "Power consumption (single pump)";
      eta = (dp*q_flow_single)/(W_single + W_eps) "Hydraulic efficiency";
    else
      eta = efficiencyCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-10)));
      W_single = dp*q_flow/eta;
    end if;
    // Fluid properties
    fluid.p = inlet.p;
    fluid.h = inlet.h;
    fluid.Xi = inlet.Xi;
    d = fluid.d;
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
      M * der(h_out) = m_flow_single*(inlet.h - outlet.h) + W_single 
        "Dynamic energy balance (density variations neglected)";
    else
      inlet.H_flow + outlet.H_flow + W_single*Np = 0 "Static energy balance";
    end if;
    
    // NPSH computations
    if computeNPSHa then
      pv=SatMedium.saturationPressure(fluid.T);
      NPSHa=(inlet.p-pv)/(d*Modelica.Constants.g_n);
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
<p>The model describes a centrifugal pump, or a group of <tt>Np</tt> identical pumps in parallel. The pump model is based on the theory of kinematic similarity: the pump characteristics are given for nominal operating conditions (rotational speed and fluid density), and then adapted to actual operating condition, according to the similarity equations. 
<p><b>Modelling options</b></p>
<p> The nominal hydraulic characteristic (head vs. volume flow rate) is given by the the replaceable function <tt>flowCharacteristic</tt>. 
<p> The pump energy balance can be specified in two alternative ways:
<ul>
<li><tt>usePowerCharacteristic = false</tt> (default option): the replaceable function <tt>efficiencyCharacteristic</tt> (efficiency vs. volume flow rate in nominal conditions) is used to determine the efficiency, and then the power consumption. The default is a constant efficiency of 0.8.
<li><tt>usePowerCharacteristic = true</tt>: the replaceable function <tt>powerCharacteristic</tt> (power consumption vs. volume flow rate in nominal conditions) is used to determine the power consumption, and then the efficiency.
</ul>
<p>
Several functions are provided in the package <tt>PumpCharacteristics</tt> to specify the characteristics as a function of some operating points at nominal conditions.
<p>Depending on the value of the <tt>checkValve</tt> parameter, the model either supports reverse flow conditions, or includes a built-in check valve to avoid flow reversal.
<p>If the <tt>in_Np</tt> input connector is wired, it provides the number of pumps in parallel; otherwise,  <tt>Np_n</tt> parallel pumps are assumed.</p>
<p>It is possible to take into account the heat capacity of the fluid inside the pump by specifying its mass <tt>M</tt> at nominal conditions; this is necessary to avoid singularities in the computation of the outlet enthalpy in case of zero flow rate. If zero flow rate conditions are always avoided, this dynamic effect can be neglected by leaving the default value <tt>M = 0</tt>, thus avoiding a fast state variable in the model.
<p>If <tt>computeNPSHa = true</tt>, the available net positive suction head is also computed; this requires a two-phase medium model to provide the fluid saturation pressure.
</HTML>",
        revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
    
  end PartialPump;
end Interfaces;
