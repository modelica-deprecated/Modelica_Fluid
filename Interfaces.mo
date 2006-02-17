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
    annotation (defaultComponentName="port_a",
                Diagram(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
               fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
               fillColor=69)), Text(extent=[-88, 206; 112, 112], string="%name")),
         Icon(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
              fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
              fillColor=69))));
  end FluidPort_a;
  
  connector FluidPort_b "Fluid connector with outlined icon" 
    extends FluidPort;
    annotation (defaultComponentName="port_b",
                Diagram(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
               fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
               fillColor=69)), Ellipse(extent=[-80, 80; 80, -80], style(color=69,
               fillColor=7)), Text(extent=[-88, 192; 112, 98], string="%name")),
         Icon(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
              fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
              fillColor=69)), Ellipse(extent=[-80, 80; 80, -80], style(color=69,
               fillColor=7))));
  end FluidPort_b;
  
  connector FluidPort_ArrayIcon 
    "Fluid connector with icon suited for an array of FluidPorts" 
    extends FluidPort;
    annotation (defaultComponentName="ports",
                Diagram(Rectangle(extent=[-100, 100; 100, -100], style(color=69,
               fillColor=69)), Rectangle(extent=[-100, 100; 100, -100], style(color=16,
               fillColor=69)), Text(extent=[-88, 206; 112, 112], string="%name")),
         Icon(Rectangle(extent=[-100, 100; 100, -100], style(color=69,
              fillColor=69)), Rectangle(extent=[-100, 100; 100, -100], style(color=16,
              fillColor=69))));
  end FluidPort_ArrayIcon;
  
  partial model PartialInitializationParameters 
    "Define parameter menu to initialize medium in component that has one medium model" 
    import Modelica_Fluid.Types;
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
    parameter Modelica_Fluid.Types.InitWithGlobalDefault.Temp initOption=
              Modelica_Fluid.Types.InitWithGlobalDefault.UseGlobalFluidOption 
      "Initialization option" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_start = Medium.p_default 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
      "= true, use T_start, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate=true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
      "Start value of temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
  protected 
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
    parameter Types.Init.Temp initOption2=
        if initOption == Modelica_Fluid.Types.InitWithGlobalDefault.UseGlobalFluidOption then 
             fluidOptions.default_initOption else initOption 
        annotation(Evaluate=true, Hide=true);
  end PartialInitializationParameters;
  
  partial model PartialGuessValueParameters 
    "Define parameter menu to initialize guess values of medium in component that has one medium model" 
    import Modelica_Fluid.Types;
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
    parameter Medium.AbsolutePressure p_start = Medium.p_default 
      "Guess value of pressure" 
      annotation(Dialog(tab = "Guess Value Initialization"));
    parameter Boolean use_T_start = true 
      "= true, use T_start, otherwise h_start" 
      annotation(Dialog(tab = "Guess Value Initialization"), Evaluate=true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
      "Guess value of temperature" 
      annotation(Dialog(tab = "Guess Value Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default 
      "Guess value of specific enthalpy" 
      annotation(Dialog(tab = "Guess Value Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Guess value of mass fractions m_i/m" 
      annotation (Dialog(tab="Guess Value Initialization", enable=Medium.nXi > 0));
  end PartialGuessValueParameters;
  
  partial model PartialSource 
    "Partial component source with one fluid connector" 
    import Modelica.Constants;
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium model within the source" 
       annotation (choicesAllMatching=true);
    FluidPort_b port(redeclare package Medium = Medium,
                     m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
      annotation (extent=[90,-10; 110,10],    rotation=0);
    Medium.BaseProperties medium "Medium in the source";
    parameter Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Temp 
      flowDirection=
              Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption 
      "Unidirectional (out of port_b) or bidirectional flow component" 
                                                                annotation(Dialog(tab="Advanced"));
  protected 
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
    parameter Boolean allowFlowReversal=
       flowDirection == Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Bidirectional
       or flowDirection == Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption
       and fluidOptions.default_flowDirection ==Modelica_Fluid.Types.FlowDirection.Bidirectional 
      "= false, if flow only out of port_b, otherwise reversing flow allowed" 
       annotation(Evaluate=true, Hide=true);
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
</html>"),
      Diagram,
      Coordsys(grid=[1,1], scale=0));
  end PartialSource;
  
  partial model PartialTwoPortTransport 
    "Partial element transporting fluid between two ports without storing mass or energy" 
    import SI = Modelica.SIunits;
    import Modelica.Constants;
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component"  annotation (
        choicesAllMatching =                                                                            true);
    
    parameter Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Temp 
      flowDirection= Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
    
    FluidPort_a port_a(redeclare package Medium = Medium,
                       m_flow(start=0,min=if allowFlowReversal then -Constants.inf else 0)) 
      "Fluid connector a (positive design flow direction is from port_a to port_b)"
      annotation (extent=[-110,-10; -90,10]);
    FluidPort_b port_b(redeclare package Medium = Medium,
                       m_flow(start=0,max=if allowFlowReversal then +Constants.inf else 0)) 
      "Fluid connector b (positive design flow direction is from port_a to port_b)"
      annotation (extent=[110,-10; 90,10]);
    Medium.BaseProperties medium_a "Medium properties in port_a";
    Medium.BaseProperties medium_b "Medium properties in port_b";
    Medium.MassFlowRate m_flow(start=0) 
      "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
    SI.VolumeFlowRate V_flow_a = port_a.m_flow/medium_a.d 
      "Volume flow rate near port_a";
    SI.Pressure dp(start=0) "Pressure difference between port_a and port_b";
    
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
</html>"),
      Icon);
  protected 
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
    parameter Boolean allowFlowReversal=
       flowDirection == Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Bidirectional
       or flowDirection == Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption
       and fluidOptions.default_flowDirection ==Modelica_Fluid.Types.FlowDirection.Bidirectional 
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
       annotation(Evaluate=true, Hide=true);
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
      annotation (extent=[-10,-110; 10,-90],    rotation=90);
    
    annotation (Documentation(info="<html>
<p>
Partial component to model an <b>absolute sensor</b>. Can be used for pressure sensor models.
Use for other properties such as temperature or density is discouraged, because the enthalpy at the connector can have different meanings, depending on the connection topology. Use <tt>PartialFlowSensor</tt> instead.
as signal.
</p>
</html>"),
      Diagram,
      Coordsys(grid=[1,1], scale=0));
  equation 
    port.m_flow = 0;
    port.H_flow = 0;
    port.mXi_flow = zeros(Medium.nXi);
  end PartialAbsoluteSensor;
  
  partial model PartialFlowSensor 
    "Partial component to model sensors that measure flow properties" 
    
    import Modelica.Constants;
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
        choicesAllMatching = true);
    Medium.SpecificEnthalpy h "enthalpy in flow";
    Medium.MassFraction[Medium.nXi] Xi "flow composition";
    
    FluidPort_a port_a(redeclare package Medium = Medium,
                       m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
      annotation (extent=[-110,-10; -90,10]);
    FluidPort_b port_b(redeclare package Medium = Medium,
                       m_flow(max=if allowFlowReversal then +Constants.inf else 0)) 
      annotation (extent=[110,-10; 90,10]);
    
    parameter Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Temp 
      flowDirection=
              Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption 
      "Unidirectional (port_a -> port_b) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
    
    annotation (Documentation(info="<html>
<p>
Partial component to model a <b>sensor</b> that measures any intensive properties
of a flow, e.g., to get temperature or density in the flow
between fluid connectors.<br>
The model includes zero-volume balance equations. Sensor models inheriting from
this partial class should add a medium instance to calculate the measured property.
</p>
</html>"),
      Diagram,
      Coordsys(grid=[1,1], scale=0));
  protected 
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
    parameter Boolean allowFlowReversal=
       flowDirection == Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Bidirectional
       or flowDirection == Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption
       and fluidOptions.default_flowDirection ==Modelica_Fluid.Types.FlowDirection.Bidirectional 
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
       annotation(Evaluate=true, Hide=true);
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
      if use_T_start then Medium.T_default else Medium.temperature_phX(pin_start,h_start,X_start) 
      "Start value of inlet temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(pin_start, T_start, X_start[1:Medium.nXi]) else Medium.h_default 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    
    parameter Real delta=0.01 "Regularisation factor" annotation(Dialog(tab="Advanced"));
    
    Modelica.Blocks.Interfaces.RealInput stemPosition 
      "Stem position in the range 0-1" annotation (extent=[-20,70; 20,110],   rotation=-90);
    
    Medium.Density d "Density at port a";
    Medium.Temperature T "Temperature at port a";
  protected 
    function sqrtR = Utilities.regRoot(delta = delta*dp_nom);
    annotation (
      Icon(Text(extent=[-143,-66; 148,-106],  string="%name"),
        Line(points=[0,60; 0,0],   style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Polygon(points=[-100,50; -100,-50; 0,0; -100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Polygon(points=[100,50; 0,0; 100,-50; 100,50],  style(
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
</html>"),
      Coordsys(grid=[2,2], scale=0));
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
    import Modelica.Constants;
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
    parameter Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Temp 
      flowDirection=
              Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption 
      "Unidirectional (inlet -> outlet) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
    parameter Boolean computeNPSHa=false "Compute NPSH Available at the inlet";
    parameter Medium.AbsolutePressure pin_start 
      "Guess value for inlet pressure" 
      annotation(Dialog(tab="Initialization"));
    parameter Medium.AbsolutePressure pout_start 
      "Guess value for outlet pressure" 
      annotation(Dialog(tab="Initialization"));
    parameter Boolean use_T_start = true 
      "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate = true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(pin_start,h_start,X_start) 
      "Guess value for temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(pin_start, T_start, X_start) else Medium.h_default 
      "Guess value for specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Guess value for mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter SI.MassFlowRate m_flow_start = 0 
      "Guess value for mass flow rate (total)" 
      annotation(Dialog(tab="Initialization"));
    constant SI.Acceleration g=Modelica.Constants.g_n;
  //  parameter Choices.Init.Options.Temp initOpt=Choices.Init.Options.noInit 
  //    "Initialisation option";
    Modelica_Fluid.Interfaces.FluidPort_a inlet(redeclare package Medium = Medium,
        p(start=pin_start),
        m_flow(start = m_flow_start,
               min = if allowFlowReversal and not checkValve then -Constants.inf else 0)) 
    annotation (extent=[-100,-40; -60,0]);
    Modelica_Fluid.Interfaces.FluidPort_b outlet(redeclare package Medium = Medium,
        p(start=pout_start),
        m_flow(start = -m_flow_start,
               max = if allowFlowReversal and not checkValve then +Constants.inf else 0)) 
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
  protected 
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
   parameter Boolean allowFlowReversal=
       flowDirection == Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Bidirectional
       or flowDirection == Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.UseGlobalFluidOption
       and fluidOptions.default_flowDirection ==Modelica_Fluid.Types.FlowDirection.Bidirectional 
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
       annotation(Evaluate=true, Hide=true);
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
      // head = (N/N_nom)^2*flowCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-6)));
      head = noEvent((((if abs(N) > 1e-6 then N else 1e-6))/N_nom)^2*flowCharacteristic(q_flow_single*N_nom/((if abs(N) > 1e-6 then N else 1e-6))));
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
  
  partial model PartialTwoPort 
    "Partial component with Medium definition and two fluid ports" 
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching = true);
    
    FluidPort_a port_a(redeclare package Medium = Medium) 
      "Fluid connector a (positive design flow direction is from port_a to port_b)"
      annotation (extent=[-110,-10; -90,10]);
    FluidPort_b port_b(redeclare package Medium = Medium) 
      "Fluid connector b (positive design flow direction is from port_a to port_b)"
      annotation (extent=[110,-10; 90,10]);
    annotation (Documentation(info="<html>
 
</html>"));
  end PartialTwoPort;
  
partial model Flow1D 
    import Modelica_Fluid.Types;
    import Modelica.Constants.*;
  replaceable package Medium = PackageMedium 
    extends Modelica.Media.Interfaces.PartialMedium "Fluid medium model" 
   annotation (choicesAllMatching=true);
    
//Discretization
  parameter Integer n(min=1) "Number of pipe segments";
  final parameter Integer np=if lumped_dp then 1 else n + 1 
      "Number of momentum balances"                                                     annotation(Dialog(tab="Advanced"),Evaluate=true);
    
//Advanced model options
  parameter Boolean allowFlowReversal=true 
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
                                                                     annotation(Dialog(tab="Advanced", group="Mass and energy balances"));
  parameter Boolean static=true "= true, no mass or energy is stored" 
                                annotation(Dialog(tab="Advanced", group="Mass and energy balances", enable=not dynamicTerm),Evaluate=true);
  parameter Boolean lumped_dp=true 
      " = true, lumped pressure drop, reduces number of pressure states to one"
                                                                              annotation(Dialog(tab="Advanced", group="Momentum balance", enable=not dynamicTerm),Evaluate=true);
  parameter Boolean kineticTerm=false " = true, include kinetic term" 
                                              annotation(Dialog(tab="Advanced", group="Momentum balance", enable=not lumped_dp),Evaluate=true);
  parameter Boolean dynamicTerm=false 
      " = true, include dynamic term, only if not lumped_dp and not static"                               annotation(Dialog(tab="Advanced", group="Momentum balance", enable=(not static and not lumped_dp)),Evaluate=true);
    
//Initialization
    parameter Types.InitWithGlobalDefault.Temp initOption=Types.
        InitWithGlobalDefault.UseGlobalFluidOption "Initialization option" 
    annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start=true "Use T_start if true, otherwise h_start"
    annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_a_start=Medium.p_default 
      "Start value of pressure at port a" 
    annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_b_start=Medium.p_default 
      "Start value of pressure at port b" 
    annotation(Dialog(tab = "Initialization"));
    final parameter Medium.AbsolutePressure[n] p_start=if n > 1 then linspace(
        p_a_start - (p_a_start - p_b_start)/(2*n), 
        p_b_start + (p_a_start - p_b_start)/(2*n), 
        n) else {(p_a_start + p_b_start)/2} "Start value of pressure";
    parameter Medium.Temperature T_start=if use_T_start then Medium.T_default
         else Medium.temperature_phX(
        (p_a_start + p_b_start)/2, 
        h_start, 
        X_start) "Start value of temperature" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=if use_T_start then 
        Medium.specificEnthalpy_pTX(
        (p_a_start + p_b_start)/2, 
        T_start, 
        X_start) else Medium.h_default "Start value of specific enthalpy" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter Medium.MassFlowRate mflow_start "Start value for mass flow rate" 
                                                                             annotation(Evaluate=true, Dialog(tab = "Initialization"));
    final parameter SI.Pressure dp_start=p_a_start - p_b_start;
    
//Geometry parameters
  parameter Modelica_Fluid.Types.CrossSectionTypes.Temp crossSectionType=
      Modelica_Fluid.Types.CrossSectionTypes.Circular 
      "Type of cross section of pipe" 
    annotation (Dialog(tab="General", group="Pipe geometry"));
  parameter SI.Length length "Length"   annotation(Dialog(tab="General", group="Pipe geometry"));
  parameter SI.Diameter d_inner "Inner diameter of circular pipe" annotation(Dialog(group="Pipe geometry", enable=crossSectionType==1));
  parameter SI.Length h_inner "Inner height of rectangular pipe"    annotation(Dialog(group="Pipe geometry", enable=crossSectionType==2));
  parameter SI.Length w_inner "Inner width of rectangular pipe"    annotation(Dialog(group="Pipe geometry", enable=crossSectionType==2));
  parameter SI.Length P_inner=if crossSectionType == 1 then Modelica.Constants.pi*d_inner else if crossSectionType == 2 then 2*h_inner + 2*
      w_inner else 1 "Inner perimeter" 
                                      annotation(Dialog(tab="General", group="Pipe geometry", enable=crossSectionType==3));
  inner parameter SI.Area A_inner=if crossSectionType == 1 then Modelica.Constants.pi*d_inner*d_inner/4 else if crossSectionType
       == 2 then h_inner*w_inner else 1 "Inner cross section area" 
                                          annotation(Dialog(tab="General", group="Pipe geometry", enable=crossSectionType==3));
  final parameter SI.Volume V=A_inner*length "Volume" 
                                                     annotation(Dialog(tab="General", group="Pipe geometry"));
    
//Pressure Drop
  replaceable package WallFriction = 
      Modelica_Fluid.PressureLosses.Utilities.WallFriction.QuadraticTurbulent 
    extends 
      Modelica_Fluid.PressureLosses.Utilities.WallFriction.PartialWallFriction 
      "Characteristic of wall friction" 
                                       annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
  parameter SI.Diameter d_h=4*A_inner/P_inner "Hydraulic diameter" annotation(Dialog(tab="General", group="Pressure loss"));
  parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                     annotation(Evaluate=true);
  parameter SI.Length roughness(min=0) = 2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" 
      annotation(Dialog(tab="General", group="Pressure loss",enable=WallFriction.use_roughness));
  parameter Boolean use_eta_nominal=false 
      "= true, if eta_nominal is used, otherwise computed from medium"                            annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
  parameter Boolean use_d_nominal=false 
      "= true, if d_nominal is used, otherwise computed from medium"                              annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
  parameter SI.DynamicViscosity eta_nominal=0.01 
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                            annotation(Dialog(tab="General", group="Pressure loss",enable=use_nominal));
  parameter SI.Density d_nominal=0.01 
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                             annotation(Dialog(tab="General", group="Pressure loss",enable=use_nominal));
  parameter Boolean show_Re=false 
      "= true, if Reynolds number is included for plotting" 
     annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
  parameter Boolean from_dp=true 
      " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
    annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
  parameter SI.AbsolutePressure dp_small=1 
      "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
    annotation(Dialog(tab="Advanced",group="Pressure loss", enable=from_dp and WallFriction.use_dp_small));
  parameter SI.MassFlowRate m_flow_small=0.01 
      "Turbulent flow if |m_flow| >= m_flow_small (only used if WallFriction=QuadraticTurbulent)"
    annotation(Dialog(tab="Advanced",group="Pressure loss", enable=not from_dp and WallFriction.use_m_flow_small));
  SI.ReynoldsNumber[n] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
      m_flow,
      (eta_a + eta_b)/2,
      diameter) if                                                                                               show_Re 
      "Reynolds number of pipe flow";
    
//Total quantities
  SI.Energy[n] U "Internal energy of fluid";
  SI.Mass[n] m "Mass of fluid";
  SI.Mass[n,Medium.nXi] mXi "Masses of independent components in the fluid";
    
//Flow quantities
  inner Medium.MassFlowRate[n + 1] m_flow(each min=if allowFlowReversal then -inf else 
              0, each start=mflow_start, each fixed=false) 
      "Mass flow rates of fluid across segment boundaries";
  SI.Velocity[n+1] v "velocity at volume boundaries (for display purposes)";
  Medium.MassFlowRate[n + 1,Medium.nXi] mXi_flow 
      "Independent mass flow rates across segment boundaries";
  Medium.EnthalpyFlowRate[n + 1] H_flow 
      "Enthalpy flow rates of fluid across segment boundaries";
  Real[n+1] I_flow;
  SI.Pressure[np] dp(start=dp0) 
      "Pressure drop due to friction loss and gravity";
    
//Source terms
  Medium.MassFlowRate[n] ms_flow "Mass flow rate, source or sink";
  Medium.MassFlowRate[n,Medium.nXi] msXi_flow 
      "Independent mass flow rates, source or sink";
  SI.HeatFlowRate[n] Qs_flow "Heat flow rate, source or sink";
    
//Fluid ports
  Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = 
        Medium, m_flow(min=if allowFlowReversal and not static then -inf else 0)) 
      "Fluid inlet port" 
                       annotation (extent=[-112,-10; -92,10]);
  Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = 
        Medium, m_flow(max=if allowFlowReversal and not static then +inf else 0)) 
      "Fluid outlet port" 
                        annotation (extent=[90,-10; 110,10]);
  Medium.BaseProperties[n] medium(
    each preferredMediumStates=if static then false else true,
    p(start=p_start),
    each h(start=h_start),
    each T(start=T_start),
    each Xi(start=X_start[1:Medium.nXi]));
 /* Medium.BaseProperties medium_a(
    preferredMediumStates=false,
    p(start=p_a_start),
    h(start=h_start),
    T(start=T_start),
    Xi(start=X_start[1:Medium.nXi]));
  Medium.BaseProperties medium_b(
    preferredMediumStates=false,
    p(start=p_b_start),
    h(start=h_start),
    T(start=T_start),
    Xi(start=X_start[1:Medium.nXi]));*/
  annotation (Diagram, Icon(Rectangle(extent=[-100,40; 100,-40], style(
          color=69,
          gradient=2,
          fillColor=69))),
      Documentation(info="<html>
<p>
From Katrins email, Nov. 28, 2005:
</p>
 
<p>
Distributed volume model, properties and flow variables are arrays, no components as in the isolated pipe. Momentum and energy balances on a staggered grid, half a momentum balance on each end of the pipe. The medium properties in the ports are those of the upstream volume. I am strongly in favour with not using the energy balance 2 (the one where the momentum balance has been substracted) here, because you are loosing all the benefits of a staggered grid. You need twice as many momentum balances to calculate the algebraic pressures at the volume boundary which appear now in the energy balance. (And I am not sure if this can be  properly handled by the tool). The pressure drop is then also part of the energy balance and needs to be in accordance with the chosen grid. I agree with you that neglecting potential and kinetic energy in this case might not comply with a highly accurate formulation for teaching purposes, but for most applications it is more than sufficient. However, velocity and gravity can play a significant role in the momentum balance, which should have the option to include those terms (-> dynamic pressure). Not intertwining the two balances has also the advantage to be able to neglect specific terms in one of the balances and not in both.
</p>
 
<p>
The model contains source terms in mass and energy balances, which are not determined here. Therefore it is a partial model and could also be used for reactions or partial condensing gases with neglectable liquid volume (-> i.e. moist air).
</p>
 
<pre>
Modelling options (via boolean flags) are:
- static or dynamic (mass and energy) balances
- lumped pressure drop
- lumped composition (not sure yet if that is feasible)
- including velocity and gravity term in momentum balance, perhaps also the dynamic term.
</pre>
 
<p>
One issue not solved yet: For pressure drop and velocity term the densities at the ports are required. A medium function computing density from p and h would be most convenient. 
</p>
</html>"));
    
  protected 
   SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else eta[1]);//approximation
  SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else eta[n]);//approximation
  SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else 
            (if use_eta_nominal then fill(eta_nominal, n) else 
      Medium.dynamicViscosity(medium));
  SI.Density[n] d=if use_d_nominal then ones(n)*d_nominal else medium.d;
  SI.Density d_a=if use_d_nominal then d_nominal else medium[1].d;//approximation
  SI.Density d_b=if use_d_nominal then d_nominal else medium[n].d;//approximation
  /*SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else noEvent(if (dp[1] >= 0) then 
      Medium.dynamicViscosity(medium_a) else Medium.dynamicViscosity(medium[1])));
  SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else noEvent(if (dp[np] >= 0) then 
      Medium.dynamicViscosity(medium[n]) else Medium.dynamicViscosity(medium_b)));
  SI.Density d_a=if use_d_nominal then d_nominal else noEvent(if (dp[1] >= 0) then 
            medium_a.d else medium[1].d);
  SI.Density d_b=if use_d_nominal then d_nominal else noEvent(if (dp[np] >= 0) then 
           medium[n].d else medium_b.d);
 The outlet port medium model (medium_a or medium_b) describes the medium properties just after the outlet connection
  point which would give wrong densities at the pipe outlet for pressure drop calculations
  if fluid flows with substantially different densities are mixed*/
  outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
  parameter Types.Init.Temp initOption2=if initOption == Types.
      InitWithGlobalDefault.UseGlobalFluidOption then fluidOptions.
      default_initOption else initOption 
      annotation(Evaluate=true, Hide=true);
  parameter SI.Pressure[np] dp0={if lumped_dp then dp_start else dp_start/(if i
       > 1 and i < np then n else 2*n) for i in 1:np};
    
//Momentum balance terms
  //SI.Force[np] DI_flow "Delta momentum flow across flow grid boundaries";
  SI.Force[np] F_f "Friction force";
  SI.Force[np] F_p "Pressure forces";
    
initial equation 
  // Initial conditions
  if not static then
    if initOption2 == Types.Init.NoInit then
    // no initial equations
    elseif initOption2 == Types.Init.SteadyState then
    //steady state initialization
      if use_T_start then
        der(medium.T) = zeros(n);
      else
        der(medium.h) = zeros(n);
      end if;
      if not (lumped_dp or Medium.singleState) then
        der(medium.p) = zeros(n);
      elseif lumped_dp then
      //  der(medium[1].p) = 0;
      end if;
      if dynamicTerm then
        der(m_flow[1:n]) = zeros(n);
      end if;
      for i in 1:n loop
        der(medium[i].Xi) = zeros(Medium.nXi);
      end for;
    elseif initOption2 == Types.Init.InitialValues then
    //Initialization with initial values
      if use_T_start then
        medium.T = ones(n)*T_start;
      else
        medium.h = ones(n)*h_start;
      end if;
    elseif initOption2 == Types.Init.SteadyStateHydraulic then
    //Steady state initialization for hydraulic states (p, m_flow)
      if use_T_start then
        medium.T = ones(n)*T_start;
      else
        medium.h = ones(n)*h_start;
      end if;
      if not (lumped_dp or Medium.singleState) then
        der(medium.p) = zeros(n);
      elseif lumped_dp then
        der(medium[1].p) = 0;
      end if;
      if dynamicTerm and not Medium.singleState then
        der(m_flow[1:n])=zeros(n);
      end if;
    else
      assert(false, "Unsupported initialization option");
    end if;
  end if;
    
equation 
    
  //Port medium models
/*  port_a.p = medium_a.p;
  port_b.p = medium_b.p;
  port_a.h = medium_a.h;
  port_b.h = medium_b.h;
  port_a.Xi = medium_a.Xi;
  port_b.Xi = medium_b.Xi;*/
    
  // Boundary conditions
  port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium[1].h);
  port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium[n].h);
  port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium[1].Xi);
  port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium[n].Xi);
  port_a.m_flow = m_flow[1];
  port_b.m_flow = -m_flow[n + 1];
    
  // Distributed flow quantities
  for i in 2:n loop
    H_flow[i] = semiLinear(m_flow[i], medium[i - 1].h, medium[i].h);
    mXi_flow[i, :] = semiLinear(m_flow[i], medium[i - 1].Xi, medium[i].Xi);
    v[i] = m_flow[i]/(medium[i - 1].d + medium[i].d)*2/A_inner;
    if kineticTerm then
      I_flow[i] = Modelica_Fluid.Utilities.regSquare(m_flow[i], delta=max(0.0001,
      0.01*mflow_start))/A_inner/noEvent(if m_flow[i]>=0 then d[i-1] else d[i]);
    else
      I_flow[i] = 0;
    end if;
  end for;
  H_flow[1] = port_a.H_flow;
  H_flow[n + 1] = -port_b.H_flow;
  mXi_flow[1, :] = port_a.mXi_flow;
  mXi_flow[n + 1, :] = -port_b.mXi_flow;
  v[1] = m_flow[1]/d_a/A_inner;
  v[n + 1] = m_flow[n + 1]/d_b/A_inner;
  if kineticTerm then
    I_flow[1] = Modelica_Fluid.Utilities.regSquare(m_flow[1], delta=max(0.0001,
      0.01*mflow_start))/A_inner/d_a;
    I_flow[n + 1] = Modelica_Fluid.Utilities.regSquare(m_flow[n + 1], max(
      0.0001, 0.01*mflow_start))/d_b/A_inner;
  else
    I_flow[1] = 0;
    I_flow[n + 1] = 0;
  end if;
    
  // Total quantities
  for i in 1:n loop
    m[i] =V/n*medium[i].d;
    mXi[i, :] = m[i]*medium[i].Xi;
    U[i] = m[i]*medium[i].u;
  end for;
    
  //Mass and energy balance
  for i in 1:n loop
    if static then
      0 = m_flow[i] - m_flow[i + 1] + ms_flow[i];
      zeros(Medium.nXi) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
      0 = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
    else
      der(m[i]) = m_flow[i] - m_flow[i + 1] + ms_flow[i];
      der(mXi[i, :]) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
      der(U[i]) = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
    end if;
  end for;
  for i in 1:n loop
  assert(allowFlowReversal or (m_flow[i]>=0),"Flow reversal not allowed");
  end for;
    
//Pressure drop and gravity
 if from_dp and not WallFriction.dp_is_zero then
    if lumped_dp then
      m_flow[1] = WallFriction.massFlowRate_dp(dp[1] - height_ab*fluidOptions.g*(
        d_a + d_b)/2, d_a, d_b, eta_a, eta_b, length, d_h, roughness,
        dp_small);
    else
      m_flow[1] = WallFriction.massFlowRate_dp(dp[1] - height_ab/2/n*
        fluidOptions.g*(d_a + d[1])/2, d_a, d[1], eta_a, eta[1],
        length/n/2, d_h, roughness, dp_small);
      for i in 2:n loop
        m_flow[i] = WallFriction.massFlowRate_dp(dp[i] - height_ab/n*
          fluidOptions.g*(d[i-1] + d[i])/2, d[i-1],
          d[i], eta[i - 1], eta[i], length/n, d_h, roughness,
          dp_small);
      end for;
      m_flow[n + 1] = WallFriction.massFlowRate_dp(dp[np] - height_ab/n/2*
        fluidOptions.g*(d[n] + d_b)/2, d[n], d_b, eta[n], eta_b,
        length/n/2, d_h, roughness, dp_small);
    end if;
  else
    if lumped_dp then
      dp[1] = WallFriction.pressureLoss_m_flow(m_flow[1], d_a, d_b, eta_a,
        eta_b, length, d_h, roughness, m_flow_small) + height_ab*
        fluidOptions.g*(d_a + d_b)/2;
    else
      dp[1] = WallFriction.pressureLoss_m_flow(m_flow[1], d_a, d[1], eta_a,
        eta[1], length/n/2, d_h, roughness, m_flow_small) + height_ab/2/n*
        fluidOptions.g*(d_a + d[1])/2;
      for i in 2:n loop
        dp[i] = WallFriction.pressureLoss_m_flow(m_flow[i], d[i-1], d[i], eta[i-1],
          eta[i], length/n, d_h, roughness, m_flow_small) + height_ab/n*
          fluidOptions.g*(d[i-1] + d[i])/2;
      end for;
      dp[np] = WallFriction.pressureLoss_m_flow(m_flow[np], d[n], d_b, eta[n],
        eta_b, length/n/2, d_h, roughness, m_flow_small) + height_ab/n/2*
        fluidOptions.g*(d[n] + d_b)/2;
    end if;
  end if;
    
//Momentum Balance
if lumped_dp then
    F_p[1] = (port_a.p - port_b.p)*A_inner;
    F_f[1] = -dp[1]*A_inner;
    zeros(np) = F_p + F_f;
    medium.p = ones(n)*(port_a.p + port_b.p)/2;
  else
    F_p[1] = (port_a.p-medium[1].p)*A_inner;
    F_f[1] = -dp[1]*A_inner;
    (if dynamicTerm then der(m_flow[1])*length/n/2 else 0) = F_p[1] + F_f[1] + (I_flow[1]-I_flow[2])/2;
    for i in 2:n loop
      F_p[i] = (medium[i-1].p-medium[i].p)*A_inner;
      F_f[i] = -dp[i]*A_inner;
      (if dynamicTerm then der(m_flow[i])*length/n else 0) = F_p[i] + F_f[i] + (I_flow[i-1]-I_flow[i+1])/2;
    end for;
    F_p[np] = (medium[n].p-port_b.p)*A_inner;
    F_f[np] = -dp[np]*A_inner;
    (if dynamicTerm then der(m_flow[n + 1])*length/n/2 else 0) = F_p[np] + F_f[np] + (I_flow[n]-I_flow[n+1])/2;
  end if;
    
end Flow1D;

partial model PartialPipeWall "Wall interface" 
  parameter Integer n(min=1) "Pipe segmentation" annotation(Dialog(tab="No Input", enable=false));
  parameter SI.Diameter a_inner "Inner cross section area" annotation(Dialog(tab="No Input", enable=false));
  parameter SI.Length a_outer "Outer cross section area" annotation(Dialog(tab="No Input", enable=false));
  parameter SI.Length length "Pipe length" annotation(Dialog(tab="No Input", enable=false));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort_a 
      "Thermal port" 
    annotation (extent=[-20,40; 20,60]);
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort_b 
      "Thermal port" 
    annotation (extent=[-20,-40; 20,-60]);
    
 annotation (Diagram, Icon(Rectangle(extent=[-100,40; 100,-40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)), Text(
        extent=[-78,20; 80,-18],
        string="%name",
        style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=7))));
end PartialPipeWall;
end Interfaces;
