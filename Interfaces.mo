package Interfaces 
  "Interfaces for steady state and unsteady, mixed-phase, multi-substance, incompressible and compressible flow without momentum." 
  
  annotation (Documentation(info="<HTML>
<p>
This library provides connector definitions for the following
type of flow systems:
</p>
<ul>
<li> One-dimensional flow of a <b>single substance</b>
     or of a <b>mixture of substances</b> with optional <b>multiple phases</b>.</li>
<li> <b>Incompressible</b> and <b>compressible</b> medium.</li>
<li> <b>Steady state</b> and <b>unsteady</b> flow.</li>
<li> <b>Momentum</b> of the flow is <b>neglected</b>.</li>
<li> The <b>kinetic energy</b> in the flow is <b>neglected</b><br>
     (this is usually a good approximation, if |v| &lt; 0.1 Ma, i.e.,
     if the flow speed is less than about 10 % of the speed of sound in
     the medium).</li>
</ul>
<p>
Additionally, the medium model interfaces are defined with package
<b>PartialMedium</b>. The definitions are made in such a way that
it is not possible to connect connectors of different media together.
</p>
<dl>
<dt><b>Main Author:</b>
<dd>Hilding Elmqvist, Dynasim</dl>
<p><b>Release Notes:</b></p>
<ul>
<li><i>Nov. 6, 2002</i>
       by Hilding Elmqvist: first version.</li>
<li><i>Nov. 11, 2002</i>
       by Hilding Elmqvist, Martin Otter: improved version.</li>
<li><i>Nov. 20-21, 2002</i>
       by Hilding Elmqvist, Mike Tiller, Allan Watson, John Batteh, Chuck Newman,
       Jonas Eborn: Improved at the 32nd Modelica Design Meeting.
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
</HTML>"));
  
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
    annotation (Icon(Text(
          extent=[-98,196; 102,102],
          string="%name",
          style(color=42))));
  
  connector HeatPort_b = Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b 
    "Thermal port for 1-dim. heat transfer (unfilled rectangular icon)" 
    annotation (Icon(
                    Text(
          extent=[-98,196; 102,102],
          string="%name",
          style(color=42))));
  
  partial model PartialMenuInitialization 
    "Define Medium model and parameter menu to initialize medium in component that has states" 
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component"  annotation (
        choicesAllMatching =                                                                            true);
    
    parameter Modelica_Fluid.Types.InitTypes.Temp initType=Modelica_Fluid.Types.InitTypes.NoInit 
      "Type of initialization"                  annotation (Evaluate=true, Dialog(
          tab="Initialization"));
    parameter Boolean init_p=true 
      "= true, if p_start is used, otherwise d_start (true is required for incompressible medium)"
      annotation (Dialog(tab="Initialization", group=
            "Initialization of mass balance"));
    parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure p_start=
        Medium.reference_p "Start value of pressure p, if init_p = true" 
    annotation(Dialog(enable=init_p,
        tab="Initialization",
        group="Initialization of mass balance"));
    parameter Modelica.Media.Interfaces.PartialMedium.Density d_start=1 
      "Start value of density d, if init_p = false" annotation (Dialog(enable=not init_p, tab=
            "Initialization", group="Initialization of mass balance"));
    parameter Boolean init_T=true 
      "= true, if T_start is used, otherwise h_start" annotation (Dialog(tab=
            "Initialization", group="Initialization of energy balance"));
    parameter Modelica.Media.Interfaces.PartialMedium.Temperature T_start=
        Modelica.SIunits.Conversions.from_degC(20) 
      "Start value of temperature T, if init_T = true" 
      annotation (Dialog(enable=init_T, tab="Initialization", group=
            "Initialization of energy balance"));
    parameter Modelica.Media.Interfaces.PartialMedium.SpecificEnthalpy h_start=
        1.e4 "Start value of specific enthalpy h, if init_T = false" annotation (
       Dialog(enable=not init_T, tab="Initialization", group="Initialization of energy balance"));
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X_start[Medium.nX]= Medium.reference_X 
      "Start values of mass fractions X" 
      annotation (Dialog(tab="Initialization", group=
            "Initialization of mass fractions (only for multi-substance fluids)"));
  end PartialMenuInitialization;
  
  partial model PartialMenuStartGuesses 
    "Define Medium model and parameter menu to provide start values as guesses (for components that do not have states)" 
    
    import Modelica.SIunits.Conversions;
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
    
    parameter Modelica.Media.Interfaces.PartialMedium.AbsolutePressure p_start=
        Medium.reference_p "Guess value of pressure p" 
                                           annotation (Dialog(tab=
            "Initialization", group="Initialization of mass balance"));
    parameter Modelica.Media.Interfaces.PartialMedium.Density d_start=1.e3 
      "Guess value of density d" annotation (Dialog(tab="Initialization", group=
            "Initialization of mass balance"));
    parameter Modelica.Media.Interfaces.PartialMedium.Temperature T_start=
        Conversions.from_degC(20) "Guess value of temperature T" annotation (
        Dialog(tab="Initialization", group="Initialization of energy balance"));
    parameter Modelica.Media.Interfaces.PartialMedium.SpecificEnthalpy h_start=
        1.e4 "Guess value of specific enthalpy h" annotation (Dialog(tab=
            "Initialization", group="Initialization of energy balance"));
    parameter Modelica.Media.Interfaces.PartialMedium.MassFraction X_start[Medium.nX]=
        Medium.reference_X "Guess values of mass fractions X" 
      annotation (Dialog(tab="Initialization", group=
            "Initialization of mass fractions (only for multi-substance fluids)"));
  end PartialMenuStartGuesses;
  
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
  
  partial model PartialAbsoluteSensor 
    "Partial component to model a sensor that measures a potential variable" 
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor" annotation (
        choicesAllMatching =                                                                        true);
    FluidPort_a port(redeclare package Medium = Medium) 
      annotation (extent=[-10, -120; 10, -100], rotation=90);
    
    annotation (Documentation(info="<html>
<p>
Partial component to model an <b>absolute sensor</b>,
e.g., to get the pressure or temperature from a fluid connector
as signal.
</p>
</html>"));
  equation 
    port.m_flow = 0;
    port.H_flow = 0;
    port.mXi_flow = zeros(Medium.nXi);
  end PartialAbsoluteSensor;
  
  partial model PartialFlowRateSensor 
    "Partial component to model a sensor that measures a flow rate" 
    
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
a <b>flow rate</b>, e.g., to get the mass flow rate 
between fluid connectors.
</p>
</html>"));
  equation 
    port_a.p   = port_b.p;
    port_a.h   = port_b.h;
    port_a.Xi = port_b.Xi;
    0 = port_a.m_flow + port_b.m_flow;
    0 = port_a.H_flow + port_b.H_flow;
    zeros(Medium.nXi) = port_a.mXi_flow + port_b.mXi_flow;
  end PartialFlowRateSensor;
  
  partial model PartialTwoPortTransport 
    "Partial element transporting fluid between two ports without storing mass or energy" 
    
    replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component"  annotation (
        choicesAllMatching =                                                                            true);
    
    FluidPort_a port_a(redeclare package Medium = Medium) 
      annotation (extent=[-120, -10; -100, 10]);
    FluidPort_b port_b(redeclare package Medium = Medium) 
      annotation (extent=[120, -10; 100, 10]);
    Medium.BaseProperties medium_a(T(start = 300.0),p(start = 1.0e5)) 
      "Medium properties in port_a";
    Medium.BaseProperties medium_b(T(start = 300.0),p(start = 1.0e5)) 
      "Medium properties in port_b";
    Medium.MassFlowRate m_flow 
      "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
    annotation (
      Coordsys(grid=[1, 1], component=[20, 20]),
      Diagram,
      Documentation(info="<html>
<p>
This component transports fluid between its two ports, without
storing mass or energy. Reversal and zero mass flow rate is taken
care off, for details see definition of built-in operator semiLinear().
When using this partial component, an equation for the momentum
balance has to be added by specifying a relationship
between the pressure drop \"port_a.p - port_b.p\" and the
mass flow rate \"m_flow = port_a.m_flow\".
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
    port_a.H_flow   = semiLinear(port_a.m_flow, port_a.h,   port_b.h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, port_b.Xi);
    
    /* Energy, mass and substance mass balance */
    port_a.H_flow + port_b.H_flow = 0;
    port_a.m_flow + port_b.m_flow = 0;
    port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
    
    // Design direction of mass flow rate
    m_flow = port_a.m_flow;
  end PartialTwoPortTransport;
  
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
  
  model PortVolume 
    "Fixed volume associated with a port by the finite volume method (used to build up physical components; fulfills mass and energy balance)" 
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
  equation 
    // port = medium properties
      port.p = medium.p;
      port.h = medium.h;
      port.Xi = medium.Xi;
    
    // Total quantities
       m    = V*medium.d;
       mXi = m*medium.Xi;
       U    = m*medium.u;
    
    // Mass and energy balance
       der(m)    = port.m_flow;
       der(mXi) = port.mXi_flow;
       der(U)    = port.H_flow;
  end PortVolume;
end Interfaces;
