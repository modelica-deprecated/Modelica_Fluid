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
     
    
    replaceable package Medium = Modelica_Media.Interfaces.PartialMedium 
      "Medium model" annotation (choicesAllMatching=true);
    
    Medium.AbsolutePressure p "Pressure in the connection point";
    flow Medium.MassFlowRate m_dot 
      "Mass flow rate from the connection point into the component";
    
    Medium.SpecificEnthalpy h 
      "Specific mixture enthalpy in the connection point";
    flow Medium.EnthalpyFlowRate H_dot 
      "Enthalpy flow rate into the component (if m_dot > 0, H_dot = m_dot*h)";
    
    Medium.MassFraction X[Medium.nX](quantity=Medium.substanceNames) 
      "Independent mixture mass fractions m_i/m in the connection point";
    flow Medium.MassFlowRate mX_dot[Medium.nX](quantity=Medium.substanceNames) 
      "Mass flow rates of the independent substances from the connection point into the component (if m_dot > 0, mX_dot = m_dot*X)";
    
  end FluidPort;
  
  connector FluidPort_a "Fluid connector with filled icon" 
    extends FluidPort;
    annotation (Diagram(Rectangle(extent=[-100, 100; 100, -100], style(color=69, 
               fillColor=69)), Text(extent=[-88, 206; 112, 112], string="%name")), 
         Icon(Rectangle(extent=[-100, 100; 100, -100], style(color=69, 
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
    annotation (Diagram(Rectangle(extent=[-100, 100; 100, -100], style(color=69, 
               fillColor=7)), Text(extent=[-88, 192; 112, 98], string="%name")), 
         Icon(Rectangle(extent=[-100, 100; 100, -100], style(color=69, 
              fillColor=7)), Text(
          extent=[-126, 160; 130, 104], 
          string="%name", 
          style(
            color=0, 
            fillColor=69, 
            fillPattern=1))));
  end FluidPort_b;
  
  partial model PartialInit 
    "Define Medium model and parameter menu to initialize medium in component that has states"
     
    
    import Modelica.SIunits.Conversions.*;
    replaceable package Medium = PackageMedium extends 
      Modelica_Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
    
    parameter Medium.Choices.Init.Temp initType=Medium.Choices.Init.NoInit 
      "|Initialization||Type of initialization" annotation (Evaluate=true);
    parameter Boolean init_p=true 
      "|Initialization|Initialization of mass balance| = true, if p_start is used, otherwise d_start (true is required for incompressible medium)";
    parameter Medium.AbsolutePressure p_start=101325 
      "|Initialization|Initialization of mass balance| Start value of pressure p, if init_p = true";
    parameter Medium.Density d_start=1 
      "|Initialization|Initialization of mass balance| Start value of density d, if init_p = false";
    parameter Boolean init_T=true 
      "|Initialization|Initialization of energy balance| = true, if T_start is used, otherwise h_start";
    parameter Medium.Temperature T_start=from_degC(20) 
      "|Initialization|Initialization of energy balance| Start value of temperature T, if init_T = true";
    parameter Medium.SpecificEnthalpy h_start=1.e4 
      "|Initialization|Initialization of energy balance| Start value of specific enthalpy h, if init_T = false";
    parameter Medium.MassFraction X_start[:]=zeros(Medium.nX) 
      "|Initialization|Initialization of mass fractions (only for multi-substance fluids)| Start values of independent mass fractions X";
  end PartialInit;
  
  partial model PartialInitAlgebraic 
    "Define Medium model and parameter menu to initialize medium for algebraic equations (e.g. ShortPipe)"
     
    
    import Modelica.SIunits.Conversions.*;
    
    replaceable package Medium = PackageMedium extends 
      Modelica_Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
    
    parameter Boolean init_p=true 
      "|Initialization|Initialization of mass balance| = true, if p_start is used, otherwise d_start (true is required for incompressible medium)";
    parameter Medium.AbsolutePressure p_start=101325 
      "|Initialization|Initialization of mass balance| Start value of pressure p, if init_p = true";
    parameter Medium.Density d_start=1.e3 
      "|Initialization|Initialization of mass balance| Start value of density d, if init_p = false";
    parameter Boolean init_T=true 
      "|Initialization|Initialization of energy balance| = true, if T_start is used, otherwise h_start";
    parameter Medium.Temperature T_start=from_degC(20) 
      "|Initialization|Initialization of energy balance| Start value of temperature T, if init_T = true";
    parameter Medium.SpecificEnthalpy h_start=1.e4 
      "|Initialization|Initialization of energy balance| Start value of specific enthalpy h, if init_T = false";
    parameter Medium.MassFraction X_start[:]=zeros(Medium.nX) 
      "|Initialization|Initialization of mass fractions (only for multi-substance fluids)| Start values of independent mass fractions X";
  end PartialInitAlgebraic;
  
  partial model PartialSource 
    "Partial component source with one fluid connector" 
    
    replaceable package Medium = PackageMedium extends 
      Modelica_Media.Interfaces.PartialMedium 
      "Medium model within the source (= different from the port medium, if port.m_dot > 0)"
       annotation (choicesAllMatching=true);
    
    FluidPort_b port(redeclare package Medium = Medium)
      annotation (extent=[100, -10; 120, 10], rotation=0);
    
    Medium.BaseProperties medium "Medium in the source";
  equation 
    port.p = medium.p;
    
    
      /* Handle reverse and zero flow (for details, see Fluid.Interfaces.semiLinear) */
    port.H_dot = semiLinear(port.m_dot, port.h, medium.h);
    port.mX_dot = semiLinear(port.m_dot, port.X, medium.X);
    annotation (Documentation(info="<html>
<p>
Partial component to model the <b>volume interface</b> of a <b>source</b>
component, such as a mass flow source. The essential
features are:
</p>
<ul>
<li> The pressure in the connection port (= port.p) is identical to the
     pressure in the volume (= medium.p).</li>
<li> The enthalpy flow rate (= port.H_dot) and the mass flow rates of the
     substances (= port.mX_dot) depend on the direction of the mass flow rate,
     according to the semiLinear(..) equations.</li>
</ul>
</html>"));
  end PartialSource;
  
  partial model PartialTwoPortTransport 
    "Partial element transporting fluid between two ports without storing mass or energy"
     
    
    import SI = Modelica.SIunits;
    import Cv = Modelica.SIunits.Conversions;
    
    extends PartialInitAlgebraic;
    
    FluidPort_a port_a(redeclare package Medium = Medium)
      annotation (extent=[-120, -10; -100, 10]);
    FluidPort_b port_b(redeclare package Medium = Medium)
      annotation (extent=[120, -10; 100, 10]);
    Medium.BaseProperties medium_a(
      final init_p=init_p, 
      final p_start=p_start, 
      final d_start=d_start, 
      final init_T=init_T, 
      final T_start=T_start, 
      final h_start=h_start, 
      final X_start=X_start) "Medium properties in port_a";
    Medium.BaseProperties medium_b(
      final init_p=init_p, 
      final p_start=p_start, 
      final d_start=d_start, 
      final init_T=init_T, 
      final T_start=T_start, 
      final h_start=h_start, 
      final X_start=X_start) "Medium properties in port_b";
    Medium.MassFlowRate m_dot 
      "Mass flow rate from port_a to port_b (m_dot > 0 is design flow direction)";
    
    annotation (
      Coordsys(grid=[1, 1], component=[20, 20]), 
      Diagram, 
      Documentation(info="<html>
<p>
This component transports fluid between its two ports, without
storing mass or energy. Reversal and zero mass flow rate is taken
care off, for details see definition of built-in operator semiLinear().
When using this partial component,
the momentum equation has to be added by specifying a relationship
between the pressure drop \"dp = port_a.p - port_b.p\" and the
mass flow rate \"port_a.m_dot\".
</p>
</html>"));
  equation 
    /* Compute medium variables from port variables */
    medium_a.p = port_a.p;
    medium_a.h = port_a.h;
    medium_a.X = port_a.X;
    medium_b.p = port_b.p;
    medium_b.h = port_b.h;
    medium_b.X = port_b.X;
    
    /* Handle reverse and zero flow */
    port_a.H_dot = semiLinear(port_a.m_dot, port_a.h, port_b.h);
    port_a.mX_dot = semiLinear(port_a.m_dot, port_a.X, port_b.X);
    
    /* Energy, mass and substance mass balance */
    port_a.H_dot + port_b.H_dot = 0;
    port_a.m_dot + port_b.m_dot = 0;
    port_a.mX_dot + port_b.mX_dot = zeros(Medium.nX);
    
    // Design direction of mass flow rate
    m_dot = port_a.m_dot;
  end PartialTwoPortTransport;
  
  partial model PartialOnePort 
    "Partial fluid component with one fluid connector (to be used as container of other components)"
     
    
    extends PartialInit;
    
    FluidPort_a port(redeclare package Medium = Medium)
      annotation (extent=[-10, -120; 10, -100], rotation=90);
    annotation (
      Icon, 
      Coordsys(grid=[1, 1], component=[20, 20]), 
      Diagram, 
      Documentation(info="<html>
<p>
This partial component should be used for new <b>container</b> models
with <b>two fluid ports</b> that contain
connections of elementary fluid component models (such as Fluid.Components.Pipe or
Fluid.Components.PortVolume).
</p>
<p>
Note, whenever a new elementary component
is implemented that is defined directly by equations, reversal and zero mass
flow rate has to be taken care off. This is most easily performed by using
the built-in operator semiLinear(). When connecting elementary
components together, the flow reversal is already handeled in these components.
</p>
</html>"));
  end PartialOnePort;
  
  partial model PartialTwoPort 
    "Partial fluid component with two fluid connectors (to be used as container of other components)"
     
    
    extends PartialInit;
    
    FluidPort_a port_a(redeclare package Medium = Medium)
      annotation (extent=[-120, -10; -100, 10]);
    FluidPort_b port_b(redeclare package Medium = Medium)
      annotation (extent=[120, -10; 100, 10]);
    annotation (
      Coordsys(grid=[1, 1], component=[20, 20]), 
      Diagram, 
      Documentation(info="<html>
<p>
This partial component should be used for new <b>container</b> models with <b>two fluid ports</b>
that contain connections of elementary fluid component models (such as Fluid.Components.Pipe or
Fluid.Components.PortVolume).
</p>
<p>
Note, whenever a new elementary component
is implemented that is defined directly by equations, reversal and zero mass
flow rate has to be taken care off. This is most easily performed by using
the submodel Modelica_Fluid.Interfaces.semiLinear. When connecting elementary
components together, the flow reversal is already handeled in these components.
</p>
</html>"));
  end PartialTwoPort;
  
  model ThermalAdaptor 
    "Adaptor between a HeatPort (heat transfer port) and a FluidPort" 
    
    replaceable package Medium = PackageMedium "Medium model" annotation (
        choicesAllMatching=true);
    Medium.BaseProperties medium "Medium used in the this model";
    
    FluidPort_a fluidPort(redeclare package Medium = Medium)
      annotation (extent=[100, -10; 120, 10], rotation=0);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort
      annotation (extent=[-120, -10; -100, 10]);
    annotation (
      Coordsys(
        extent=[-100, -100; 100, 100], 
        grid=[2, 2], 
        component=[20, 20]), 
      Icon(
        Text(extent=[-134, 134; 134, 72], string="%name"), 
        Rectangle(extent=[0, 60; 60, -60], style(
            color=0, 
            gradient=1, 
            fillColor=9)), 
        Rectangle(extent=[-20, 60; 0, -60], style(color=0, fillColor=42)), 
        Rectangle(extent=[-100, 6; -20, -6], style(color=42, fillColor=42)), 
        Rectangle(extent=[60, 60; 100, -60], style(
            color=69, 
            fillColor=69, 
            fillPattern=1))), 
      Diagram, 
      Documentation(info="<html>
<p>
This component is an adaptor between a connector from the Modelica.Thermal.HeatTransfer
library and a connector from the Fluid library. This adaptor is used
to model the heat transfer of a part into a fluid by connecting the
heat transfer components via this adaptor to a connection point of
the fluid.
</p>
</html>"));
  equation 
    
    
      // Intensive quantities of the fluidPort are used to compute medium properties
    medium.p = fluidPort.p;
    medium.h = fluidPort.h;
    medium.X = fluidPort.X;
    
    // No mass flow from the heatPort to the fluidPort
    fluidPort.m_dot = 0;
    fluidPort.mX_dot = zeros(Medium.nX);
    
    // Energy balance between the two ports
    heatPort.Q_dot + fluidPort.H_dot = 0;
    
    // Boundary condition
    heatPort.T = medium.T;
  end ThermalAdaptor;
  
  model JunctionVolume 
    "Fixed volume associated with a port by the finite volume method (the medium properties of the volume are the ones of the port)"
     
    
    import SI = Modelica.SIunits;
    import Cv = Modelica.SIunits.Conversions;
    
    extends PartialInit;
    
    FluidPort_b port(redeclare package Medium = Medium)
      annotation (extent=[-10, -10; 10, 10], rotation=0);
    
    parameter SI.Volume V=1e-6 "Fixed size of junction volume";
    
    Medium.BaseProperties medium(
      preferedMediumStates=true, 
      final initType=initType, 
      final init_p=init_p, 
      final p_start=p_start, 
      final d_start=d_start, 
      final init_T=init_T, 
      final T_start=T_start, 
      final h_start=h_start, 
      final X_start=X_start);
    
    SI.Energy U "Internal energy of port volume";
    Real m(quantity=Medium.mediumName, unit="kg") "Mass of junction volume";
    Real mX[Medium.nX](quantity=Medium.substanceNames, each unit="kg") 
      "Independent substance masses of junction volume";
  equation 
    medium.p = port.p;
    medium.h = port.h;
    medium.X = port.X;
    m = V*medium.d;
    U = m*medium.u;
    mX = m*medium.X;
    der(m) = port.m_dot;
    der(U) = port.H_dot;
    der(mX) = port.mX_dot;
    annotation (Icon(
        Ellipse(extent=[-100, 100; 100, -100], style(
            color=0, 
            gradient=3, 
            fillColor=71)), 
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
inside the volume (=h) is always identical to the specific enthalpy
in the port (port.h = h). Usually, this model is used when
discretizing a component according to the finite volume method into
volumes in internal ports that only store energy and mass and into
transport elements that just transport energy, mass and momentum
between the internal ports without storing these quantities during the
transport.
</p>
</html>"));
    
  end JunctionVolume;
  
  function checkAmbient "Check whether ambient definition is correct" 
    
    extends Modelica.Icons.Function;
    input String mediumName;
    input Boolean incompressible;
    input Boolean define_p;
    input Boolean reducedX;
    input Integer nX;
    input Real X_ambient[:];
  algorithm 
    assert(not incompressible or incompressible and define_p, "
Wrong value of parameter define_p (= false) in ambient source component:
The selected medium \"" + mediumName + "\" is incompressible.
Therefore, an ambient density cannot be defined and
define_p = true is required.
");
    
    for i in 1:nX loop
      assert(X_ambient[i] >= 0.0, "
Wrong ambient mass fractions in medium \"" + mediumName + "\":
The ambient value X_ambient(" + integerString(i) + ") = " + realString(
        X_ambient[i]) + "
is negative. It must be positive.
");
    end for;
    
    assert(reducedX or not reducedX and nX > 0 and abs(sum(X_ambient) - 1.0) < 
      1.e-10, "
Wrong ambient mass fractions in medium \"" + mediumName + "\":
This medium requires that the ambient mass fractions X_ambient
sum up to 1. However, sum(X_ambient) = " + realString(sum(X_ambient)) + ".
");
    
    assert(not reducedX or reducedX and sum(X_ambient) < 1 + 1.e-10, "
Wrong ambient mass fractions in medium \"" + mediumName + "\":
This medium requires that the sum of the ambient mass fractions X_ambient
is at most 1. However, sum(X_ambient) = " + realString(sum(X_ambient)) + ".
");
  end checkAmbient;
end Interfaces;
