within Modelica_Fluid;
package Sensors 
  "Ideal sensor components (provide the variables in a fluid connector as signals)" 
  extends Modelica_Fluid.Icons.VariantLibrary;
  
  annotation (preferedView="info", Documentation(info="<html>
<p align = justify>
Package <b>Sensors</b> consists of idealized sensor components that
provide variables of a medium model and/or fluid ports as
output signals. These signals can be, e.g., further processed
with components of the Modelica.Blocks library.
Also more realistic sensor models can be built, by further
processing (e.g., by attaching block Modelica.Blocks.FirstOrder to
model the time constant of the sensor).
 
</p>
 
<p align = justify>For the thermodynamic state variables temperature, specific entalpy, specific entropy and density the fluid library provides two different types of sensors: <b>one port</b> and <b>two port</b> sensors. </p>
 
<ul>
<li>
The <b>one port</b> sensors have the advantage of easily introducing them to and removing them from a complete model, due to no connections have to be broken. But if the sensor is placed in a three or more way connection without any explicit junction model, the resulting temperature signal can be something unexpected. <a href= \"Modelica_Fluid.Test.TestComponents.Sensors.TestTemperatureSensor\">Modelica_Fluid.Test.TestComponents.Sensors.TestTemperatureSensor </a> provides a test case, which shows the arising trouble. That means, that one port sensors should only be used in combination with explicit junction models (<a href= \"Modelica_Fluid.Junctions\">Modelica_Fluid.Junctions</a>) or in the case that the user is really sure to connect the sensor only between two other components! </li> 
 
<li> The <b>two port</b> sensors are the save way for getting the expecting results. If they are connected in series the output signal corresponds to the passing fluid. To provide a temperature signal of a volume, e.g. a tank, the two port sensors can be used only with one port connected to the volume.</li>
</ul>
 
</html>",
      revisions="<html>
<ul>
<li><i>31 Oct 2007</i>
    by Carsten Heinrich<br>
       updated sensor models, included one and two port sensors for thermodynamic state variables</li>
</ul>
</html>"));
  
  model Pressure "Ideal pressure sensor" 
    extends Sensors.BaseClasses.PartialAbsoluteSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput p(unit = "Pa") "Pressure at port" 
      annotation (extent=[100,-10; 120,10]);
    
    annotation (
    Diagram(
        Line(points=[70,0; 100,0], style(rgbcolor={0,0,127})),
        Line(points=[0,-70; 0,-100], style(color=69))),
    Icon(
        Line(points=[70,0; 100,0], style(rgbcolor={0,0,127})),
        Line(points=[0,-70; 0,-100], style(color=69)),
        Text(extent=[-126,160; 138,98], string="%name"),
        Text(
          extent=[212,-51; 52,-103],
          style(color=0),
          string="p")),
      Documentation(info="<HTML>
<p>
This component monitors the absolute pressure at its fluid port. The sensor is 
ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"),   Coordsys(grid=[1,1], component=[20,20]));
  equation 
    p = port.p;
  end Pressure;
  
  model DensityOnePort "Ideal one port density sensor" 
    extends Sensors.BaseClasses.PartialAbsoluteSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput d(unit = "kg/m3") 
      "Density in port medium" 
      annotation (extent=[90,-10; 110,10],    rotation=0);
    
  annotation (defaultComponentName="density",
    Diagram(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127}))),
    Icon(
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127})),
        Text(extent=[-126,160; 138,98], string="%name"),
        Text(
          extent=[170,-53; 10,-105],
          style(color=0),
          string="d")),
    Documentation(info="<HTML>
<p>
This component monitors the density of the medium in the flow
between fluid ports. The sensor is 
ideal, i.e., it does not influence the fluid.
</p>
<p>If using the one port sensor please read the <a href = Modelica_Fluid.Sensors>Information</a>  first.</p>
 
</HTML>
"));
  equation 
    d = Medium.density(Medium.setState_phX(port.p, inflow(port.h_outflow), inflow(port.Xi_outflow)));
  end DensityOnePort;
  
  model DensityTwoPort "Ideal two port density sensor" 
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput d(unit="kg/m3") 
      "Density of the passing fluid" 
      annotation (extent=[-10,-120; 10,-100], rotation=-90);
    
  annotation (defaultComponentName="density",
    Diagram(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127}))),
    Icon(
      Text(extent=[-126,160; 138,98],   string="%name"),
        Text(
          extent=[168,-71; 8,-123],
          style(color=0),
          string="d"),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127})),
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255}))),
    Documentation(info="<HTML>
<p>
This component monitors the volume flow rate flowing from port_a to port_b. 
The sensor is ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    d = if port_a.m_flow > 0 then 
           Medium.density(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow)) else 
           Medium.density(Medium.setState_phX(port_a.p, port_a.h_outflow, port_a.Xi_outflow));
  end DensityTwoPort;
  
  model TemperatureOnePort "Ideal one port temperature sensor" 
      extends Sensors.BaseClasses.PartialAbsoluteSensor;
    
    Modelica.Blocks.Interfaces.RealOutput T(unit = "K") 
      "Temperature in port medium" 
      annotation (extent=[60,-10; 80,10],     rotation=0);
    
  annotation (defaultComponentName="temperature",
    Diagram(
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127})),
        Ellipse(extent=[-20, -98; 20, -60], style(
            color=0,
            thickness=2,
            fillColor=42)),
        Rectangle(extent=[-12, 40; 12, -68], style(color=42, fillColor=42)),
        Polygon(points=[-12, 40; -12, 80; -10, 86; -6, 88; 0, 90; 6, 88; 10, 86;
               12, 80; 12, 40; -12, 40], style(color=0, thickness=2)),
        Line(points=[-12, 40; -12, -64], style(color=0, thickness=2)),
        Line(points=[12, 40; 12, -64], style(color=0, thickness=2)),
        Line(points=[-40, -20; -12, -20], style(color=0)),
        Line(points=[-40, 20; -12, 20], style(color=0)),
        Line(points=[-40, 60; -12, 60], style(color=0))),
      Icon(
        Ellipse(extent=[-20,-88; 20,-50],   style(
            color=0,
            thickness=2,
            fillColor=42)),
        Rectangle(extent=[-12,50; 12,-58],   style(color=42, fillColor=42)),
        Polygon(points=[-12,50; -12,90; -10,96; -6,98; 0,100; 6,98; 10,96; 12,
              90; 12,50; -12,50],        style(color=0, thickness=2)),
        Line(points=[-12,50; -12,-54],   style(color=0, thickness=2)),
        Line(points=[12,50; 12,-54],   style(color=0, thickness=2)),
        Line(points=[-40,-10; -12,-10],   style(color=0)),
        Line(points=[-40,30; -12,30],   style(color=0)),
        Line(points=[-40,70; -12,70],   style(color=0)),
        Text(
          extent=[120,-40; 0,-90],
          style(color=0),
          string="T"),
        Text(extent=[-126,160; 138,98],   string="%name"),
        Line(points=[12,0; 60,0],  style(rgbcolor={0,0,127}))),
      Documentation(info="<HTML>
<p>
This component monitors the temperature of the medium in the flow
between fluid ports. The sensor is 
ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    T = Medium.temperature(Medium.setState_phX(port.p, inflow(port.h_outflow), inflow(port.Xi_outflow)));
  end TemperatureOnePort;
  
  model TemperatureTwoPort "Ideal two port temperature sensor" 
    extends Sensors.BaseClasses.PartialFlowSensor;
    
    Modelica.Blocks.Interfaces.RealOutput T(unit="K") 
      "Temperature of the passing fluid" 
      annotation (extent=[-10,-120; 10,-100], rotation=-90);
    
  annotation (defaultComponentName="temperature",
    Diagram(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127}))),
    Icon(
      Text(extent=[-126,160; 138,98],   string="%name"),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127})),
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Ellipse(extent=[-20,-88; 20,-50],   style(
            color=0,
            thickness=2,
            fillColor=42)),
        Rectangle(extent=[-12,50; 12,-58],   style(color=42, fillColor=42)),
        Polygon(points=[-12,50; -12,90; -10,96; -6,98; 0,100; 6,98; 10,96; 12,90;
              12,50; -12,50],            style(color=0, thickness=2)),
        Line(points=[-12,50; -12,-54],   style(color=0, thickness=2)),
        Line(points=[12,50; 12,-54],   style(color=0, thickness=2)),
        Line(points=[-40,-10; -12,-10],   style(color=0)),
        Line(points=[-40,30; -12,30],   style(color=0)),
        Line(points=[-40,70; -12,70],   style(color=0)),
        Text(
          extent=[120,-40; 0,-90],
          style(color=0),
          string="T"),
        Line(points=[12,0; 60,0],  style(rgbcolor={0,0,127}))),
    Documentation(info="<HTML>
<p>
This component monitors the temperature of the passing fluid. 
The sensor is ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    T = if port_a.m_flow > 0 then 
           Medium.temperature(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow)) else 
           Medium.temperature(Medium.setState_phX(port_a.p, port_a.h_outflow, port_a.Xi_outflow));
  end TemperatureTwoPort;
  
  model SpecificEnthalpyOnePort "Ideal one port specific enthalphy sensor" 
    extends Sensors.BaseClasses.PartialAbsoluteSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput h_out(unit="J/kg") 
      "Specific enthalpy in port medium" 
      annotation (extent=[90,-10; 110,10],    rotation=0);
    
  annotation (defaultComponentName="specificEnthalpy",
    Diagram(
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127}))),
    Icon(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127})),
        Text(extent=[-126,160; 138,98], string="%name"),
        Text(
          extent=[212,-51; 52,-103],
          style(color=0),
          string="h")),
    Documentation(info="<HTML>
<p>
This component monitors the specific enthalphy of the medium in the flow
between fluid ports. The sensor is ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    h_out = inflow(port.h_outflow);
  end SpecificEnthalpyOnePort;
  
  model SpecificEnthalpyTwoPort 
    "Ideal two port sensor for the specific enthalpy" 
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput h_out(unit="J/kg") 
      "Specific enthalpy of the passing fluid" 
      annotation (extent=[-10,-120; 10,-100], rotation=-90);
    
  annotation (defaultComponentName="specificEnthalpy",
    Diagram(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127}))),
    Icon(
      Text(extent=[-126,160; 138,98],   string="%name"),
        Text(
          extent=[168,-71; 8,-123],
          style(color=0),
          string="h"),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127})),
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255}))),
    Documentation(info="<HTML>
<p>
This component monitors the specific enthalpy of passing fluid. 
The sensor is ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    h_out = if port_a.m_flow > 0 then port_b.h_outflow else port_a.h_outflow;
  end SpecificEnthalpyTwoPort;
  
  model SpecificEntropyOnePort "Ideal one port specific entropy sensor" 
    extends Sensors.BaseClasses.PartialAbsoluteSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput s(unit="J/(kg.K)") 
      "Specific entropy in port medium" 
      annotation (extent=[90,-10; 110,10],    rotation=0);
    
  annotation (defaultComponentName="specificEntropy",
    Diagram(
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127}))),
    Icon(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127})),
        Text(extent=[-126,160; 138,98], string="%name"),
        Text(
          extent=[170,-55; 10,-107],
          style(color=0),
          string="s")),
    Documentation(info="<HTML>
<p>
This component monitors the specific enthalphy of the medium in the flow
between fluid ports. The sensor is ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    s = Medium.specificEntropy(Medium.setState_phX(port.p, inflow(port.h_outflow), inflow(port.Xi_outflow)));
  end SpecificEntropyOnePort;
  
  model SpecificEntropyTwoPort "Ideal two port sensor for the specific entropy" 
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput s(unit="J/(kg.K)") 
      "Specific entropy of the passing fluid" 
      annotation (extent=[-10,-120; 10,-100], rotation=-90);
    
  annotation (defaultComponentName="specificEntropy",
    Diagram(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127}))),
    Icon(
      Text(extent=[-126,160; 138,98],   string="%name"),
        Text(
          extent=[168,-71; 8,-123],
          style(color=0),
          string="s"),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127})),
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255}))),
    Documentation(info="<HTML>
<p>
This component monitors the specific entropy of the passing fluid. 
The sensor is ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    s = if port_a.m_flow > 0 then 
           Medium.specificEntropy(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow)) else 
           Medium.specificEntropy(Medium.setState_phX(port_a.p, port_a.h_outflow, port_a.Xi_outflow));
    
  end SpecificEntropyTwoPort;
  
  model MassFlowRate "Ideal sensor for mass flow rate" 
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput m_flow(unit = "kg/s") 
      "Mass flow rate from port_a to port_b" annotation (extent=[-10,-120; 10,
          -100], rotation=-90);
    
  annotation (
    Diagram(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127}))),
    Icon(
      Text(extent=[-140,147; 140,82],   string="%name"),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Text(
          extent=[178,-81; 18,-133],
          style(color=0),
          string="m_flow"),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127})),
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255}))),
    Documentation(info="<HTML>
<p>
This component monitors the mass flow rate flowing from port_a to port_b. 
The sensor is ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    m_flow = port_a.m_flow;
  end MassFlowRate;
  
  model VolumeFlowRate "Ideal sensor for volume flow rate" 
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput V_flow(unit = "m3/s") 
      "Volume flow rate from port_a to port_b" 
      annotation (extent=[-10,-120; 10,-100], rotation=-90);
    
  annotation (
    Diagram(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127}))),
    Icon(
      Text(extent=[-126,160; 138,98],   string="%name"),
        Text(
          extent=[188,-71; 28,-123],
          style(color=0),
          string="V_flow"),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127})),
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255}))),
    Documentation(info="<HTML>
<p>
This component monitors the volume flow rate flowing from port_a to port_b. 
The sensor is ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    V_flow = port_a.m_flow/(if port_a.m_flow > 0 then 
                   Medium.density(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow)) else 
                   Medium.density(Medium.setState_phX(port_a.p, port_a.h_outflow, port_a.Xi_outflow)));
  end VolumeFlowRate;
  
  model RelativePressure "Ideal relative pressure sensor" 
    extends Modelica.Icons.TranslationalSensor;
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
        choicesAllMatching = true);
    
    Modelica_Fluid.Interfaces.FluidPort_a port_a(
                                  redeclare package Medium = Medium) 
      annotation (extent=[-120, -10; -100, 10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_b(
                                  redeclare package Medium = Medium) 
      annotation (extent=[120, -10; 100, 10]);
    
    Modelica.Blocks.Interfaces.RealOutput p_rel(unit="Pa") 
      "Relative pressure signal"                                                      annotation (extent=[-10, -80; 10, -100], rotation=90);
    annotation (
      Icon(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(extent=[-140, 94; 144, 34], string="%name"),
        Text(
          extent=[92, -62; 34, -122],
          string="p_rel",
          style(color=0))),
      Diagram(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(
          extent=[64, -74; 32, -102],
          string="p_rel",
          style(color=0))),
      Documentation(info="<HTML>
<p>
The relative pressure \"port_a.p - port_b.p\" is determined between
the two ports of this component and is provided as output signal. The
sensor should be connected in parallel with other equipment, no flow
through the sensor is allowed.
</p>
</HTML>
"));
  equation 
    // Zero flow equations for connectors
    port_a.m_flow = 0;
    port_b.m_flow = 0;
    
    // No contribution of specific quantities
    port_a.h_outflow = 0;
    port_b.h_outflow = 0;
    port_a.Xi_outflow = zeros(Medium.nXi);
    port_b.Xi_outflow = zeros(Medium.nXi);
    port_a.C_outflow  = zeros(Medium.nC);
    port_b.C_outflow  = zeros(Medium.nC);
    
    // Relative pressure  
    p_rel = port_a.p - port_b.p;
  end RelativePressure;
  
  model RelativeTemperature "Ideal relative temperature sensor" 
    extends Modelica.Icons.TranslationalSensor;
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
        choicesAllMatching = true);
    Modelica_Fluid.Interfaces.FluidPort_a port_a(
                                  redeclare package Medium = Medium) 
      annotation (extent=[-120, -10; -100, 10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_b(
                                  redeclare package Medium = Medium) 
      annotation (extent=[120, -10; 100, 10]);
    
    Modelica.Blocks.Interfaces.RealOutput T_rel(unit="K") 
      "Relative temperature signal"                                                                               annotation (extent=[-10, -80; 10, -100], rotation=90);
    annotation (
      Icon(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(extent=[-140, 94; 144, 34], string="%name"),
        Text(
          extent=[92, -62; 34, -122],
          string="p_rel",
          style(color=0))),
      Diagram(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(
          extent=[64, -74; 32, -102],
          style(color=0),
          string="T_rel")),
      Documentation(info="<HTML>
<p>
The relative temperature \"port_a.T - port_b.T\" is determined between
the two ports of this component and is provided as output signal. The
sensor should be connected in parallel with other equipment, no flow
through the sensor is allowed.
</p>
</HTML>
"));
  equation 
    // Zero flow equations for connectors
    port_a.m_flow = 0;
    port_b.m_flow = 0;
    
    // No contribution of specific quantities
    port_a.h_outflow = 0;
    port_b.h_outflow = 0;
    port_a.Xi_outflow = zeros(Medium.nXi);
    port_b.Xi_outflow = zeros(Medium.nXi);
    port_a.C_outflow  = zeros(Medium.nC);
    port_b.C_outflow  = zeros(Medium.nC);
    
    // Relative temperature
    T_rel = Medium.temperature(Medium.setState_phX(port_a.p, inflow(port_a.h_outflow), inflow(port_a.Xi_outflow))) -
            Medium.temperature(Medium.setState_phX(port_b.p, inflow(port_b.h_outflow), inflow(port_b.Xi_outflow)));
  end RelativeTemperature;
  
/*
  model SpecificInternalEnergy "Ideal specific internal energy sensor" 
    extends BaseClasses.Sensors.PartialSensor;
    extends Modelica.Icons.RotationalSensor;
    Medium.BaseProperties medium;
    Modelica.Blocks.Interfaces.RealOutput u(unit = "J/kg") 
      "Specific internal energy in port medium" annotation (extent=[100,-10; 120,10]);
    
  annotation (
    Diagram(
        Line(points=[70,0; 100,0], style(rgbcolor={0,0,127})),
        Line(points=[0,-70; 0,-100], style(color=69))),
    Icon(
        Line(points=[70,0; 100,0], style(rgbcolor={0,0,127})),
        Line(points=[0,-70; 0,-100], style(color=69)),
        Text(extent=[-126,160; 138,98], string="%name"),
        Text(
          extent=[212,-51; 52,-103],
          style(color=0),
          string="u")),
    Documentation(info="<HTML>
<p>
This component monitors the specific internal energy of the medium in the fluid port. The sensor is 
ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    port.p   = medium.p;
    port.h   = medium.h;
    port.Xi = medium.Xi;
    u = medium.u;
  end SpecificInternalEnergy;
  
  model RelDensity "Ideal relative density sensor" 
    extends Interfaces.PartialRelativeSensor;
    extends Modelica.Icons.TranslationalSensor;
    Medium.BaseProperties medium_a;
    Medium.BaseProperties medium_b;
    Modelica.Blocks.Interfaces.RealOutput d_rel(redeclare type SignalType = 
          SI.Density) "Relative density signal" annotation (extent=[-10, -80; 10, -100], rotation=90);
    annotation (
      Icon(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(extent=[-140, 94; 144, 34], string="%name"),
        Text(
          extent=[92, -62; 34, -122],
          string="d_rel",
          style(color=0))),
      Diagram(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(
          extent=[64, -74; 32, -102],
          string="d_rel",
          style(color=0))),
      Documentation(info="<HTML>
<p>
The relative density \"d(port_a) - d(port_b)\" is determined between
the two ports of this component and is provided as output signal.
</p>
</HTML>
"));
  equation 
    port_a.p   = medium_a.p;
    port_a.h   = medium_a.h;
    port_a.Xi = medium_a.Xi;
    port_b.p   = medium_b.p;
    port_b.h   = medium_b.h;
    port_b.Xi = medium_b.Xi;
    
    d_rel = medium_a.d - medium_b.d;
  end RelDensity;
  
  model RelTemperature "Ideal relative temperature sensor" 
    extends Interfaces.PartialRelativeSensor;
    extends Modelica.Icons.TranslationalSensor;
    Medium.BaseProperties medium_a;
    Medium.BaseProperties medium_b;
    Modelica.Blocks.Interfaces.RealOutput T_rel(redeclare type SignalType = 
          SI.Temperature) "Relative temperature signal" annotation (extent=[-10, -80; 10, -100], rotation=90);
    annotation (
      Icon(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(extent=[-140, 94; 144, 34], string="%name"),
        Text(
          extent=[92, -62; 34, -122],
          string="T_rel",
          style(color=0))),
      Diagram(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(
          extent=[64, -74; 32, -102],
          string="T_rel",
          style(color=0))),
      Documentation(info="<HTML>
<p>
The relative temperature \"T(port_a) - T(port_b)\" is determined between
the two ports of this component and is provided as output signal.
</p>
</HTML>
"));
  equation 
    port_a.p   = medium_a.p;
    port_a.h   = medium_a.h;
    port_a.Xi = medium_a.Xi;
    port_b.p   = medium_b.p;
    port_b.h   = medium_b.h;
    port_b.Xi = medium_b.Xi;
    
    T_rel = medium_a.T - medium_b.T;
  end RelTemperature;
  
  model RelSpecificEnthalpy "Ideal relative specific enthalpy sensor" 
    extends Interfaces.PartialRelativeSensor;
    extends Modelica.Icons.TranslationalSensor;
    Modelica.Blocks.Interfaces.RealOutput h_rel(redeclare type SignalType = 
          SI.SpecificEnthalpy) "Relative specific enthalpy signal" annotation (extent=[-10, -80; 10, -100], rotation=90);
    annotation (
      Icon(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(extent=[-140, 94; 144, 34], string="%name"),
        Text(
          extent=[92, -62; 34, -122],
          string="h_rel",
          style(color=0))),
      Diagram(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(
          extent=[64, -74; 32, -102],
          string="h_rel",
          style(color=0))),
      Documentation(info="<HTML>
<p>
The relative specific enthalpy \"port_a.h - port_b.h\" is determined between
the two ports of this component and is provided as output signal.
</p>
</HTML>
"));
  equation 
    h_rel = port_a.h - port_b.h;
  end RelSpecificEnthalpy;
  
  model RelSpecificInternalEnergy 
    "Ideal relative specific internal energy sensor" 
    extends Interfaces.PartialRelativeSensor;
    extends Modelica.Icons.TranslationalSensor;
    Medium.BaseProperties medium_a;
    Medium.BaseProperties medium_b;
    Modelica.Blocks.Interfaces.RealOutput u_rel(redeclare type SignalType = 
          SI.SpecificEnergy) "Relative specific internal energy signal" annotation (extent=[-10, -80; 10, -100], rotation=90);
    annotation (
      Icon(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(extent=[-140, 94; 144, 34], string="%name"),
        Text(
          extent=[92, -62; 34, -122],
          string="u_rel",
          style(color=0))),
      Diagram(
        Line(points=[-100, 0; -70, 0], style(color=69)),
        Line(points=[70, 0; 100, 0], style(color=69)),
        Line(points=[0, -30; 0, -80], style(rgbcolor={0,0,127})),
        Text(
          extent=[64, -74; 32, -102],
          string="u_rel",
          style(color=0))),
      Documentation(info="<HTML>
<p>
The relative specific internal energy \"u(port_a) - u(port_b)\" is determined between
the two ports of this component and is provided as output signal.
</p>
</HTML>
"));
  equation 
    port_a.p   = medium_a.p;
    port_a.h   = medium_a.h;
    port_a.Xi = medium_a.Xi;
    port_b.p   = medium_b.p;
    port_b.h   = medium_b.h;
    port_b.Xi = medium_b.Xi;
    
    u_rel = medium_a.u - medium_b.u;
  end RelSpecificInternalEnergy;
*/
  
  package BaseClasses 
    extends Modelica_Fluid.Icons.BaseClassLibrary;
    
    partial model PartialAbsoluteSensor 
      "Partial component to model a sensor that measures a potential variable" 
      
      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium 
        "Medium in the sensor" 
        annotation(choicesAllMatching=true);
      
      Modelica_Fluid.Interfaces.FluidPort_a port(
        redeclare package Medium=Medium) 
        annotation (extent=[-10,-110; 10,-90],rotation=90);
      
      annotation (Documentation(info="<html>
<p>
Partial component to model an <b>absolute sensor</b>. Can be used for pressure sensor models.
Use for other properties such as temperature or density is discouraged, because the enthalpy at the connector can have different meanings, depending on the connection topology. Use <tt>PartialFlowSensor</tt> instead.
as signal.
</p>
</html>"),
        Diagram,
        Coordsys(grid=[1,1], scale=0),
        Icon);
    equation 
      port.m_flow = 0;
      port.h_outflow = 0;
      port.Xi_outflow = zeros(Medium.nXi);
      port.C_outflow = zeros(Medium.nC);
    end PartialAbsoluteSensor;
    
    partial model PartialFlowSensor 
      "Partial component to model sensors that measure flow properties" 
      import Modelica.Constants;
      import Modelica_Fluid.Types.FlowDirection;
      
      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium 
        "Medium in the sensor" 
        annotation(choicesAllMatching=true);
      
      Modelica_Fluid.Interfaces.FluidPort_a port_a(
        redeclare package Medium=Medium,
        m_flow(min=if allowFlowReversal then -Constants.inf else 0.0)) 
        annotation (extent=[-110,-10; -90,10]);
      Modelica_Fluid.Interfaces.FluidPort_b port_b(
        redeclare package Medium=Medium,
        m_flow(max=if allowFlowReversal then +Constants.inf else 0.0)) 
        annotation (extent=[110,-10; 90,10]);
      
      parameter FlowDirection.Temp flowDirection=FlowDirection.Bidirectional 
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
      parameter Boolean allowFlowReversal=
         flowDirection == Modelica_Fluid.Types.FlowDirection.Bidirectional 
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
         annotation(Evaluate=true, Hide=true);
    equation 
      // mass balance
      0 = port_a.m_flow + port_b.m_flow;
      
      // momentum equation (no pressure loss)
      port_a.p = port_b.p;
      
      // isenthalpic state transformation (no storage and no loss of energy)
      port_a.h_outflow = inflow(port_b.h_outflow);
      port_b.h_outflow = inflow(port_a.h_outflow);
      
      port_a.Xi_outflow = inflow(port_b.Xi_outflow);
      port_b.Xi_outflow = inflow(port_a.Xi_outflow);
      
      port_a.C_outflow = inflow(port_b.C_outflow);
      port_b.C_outflow = inflow(port_a.C_outflow);
    end PartialFlowSensor;
    
  protected 
    partial model PartialRelativeSensor 
      "Partial component to model a sensor that measures the difference of effort variables at two ports" 
      
      replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
          choicesAllMatching =                                                                         true);
      
      Modelica_Fluid.Interfaces.FluidPort_a port_a(
                                    redeclare package Medium = Medium) 
        annotation (extent=[-120, -10; -100, 10]);
      Modelica_Fluid.Interfaces.FluidPort_b port_b(
                                    redeclare package Medium = Medium) 
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
    
  end BaseClasses;
end Sensors;
