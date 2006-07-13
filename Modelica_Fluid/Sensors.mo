package Sensors 
  "Ideal sensor components (provide the variables in a fluid connector as signals)" 
  extends Modelica_Fluid.Icons.VariantLibrary;
  import SI = Modelica.SIunits;
  
  annotation (preferedView="info", Documentation(info="<html>
<p>
Package <b>Sensors</b> consists of idealized sensor components that
provide variables of a medium model and/or fluid ports as
output signals. These signals can be, e.g., further processed
with components of the Modelica.Blocks library.
Also more realistic sensor models can be built, by further
processing (e.g., by attaching block Modelica.Blocks.FirstOrder to
model the time constant of the sensor).
</p>
</html>"));
  
  model Pressure "Ideal pressure sensor" 
    import SI = Modelica.SIunits;
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
  
  model Density "Ideal density sensor" 
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Medium.BaseProperties medium;
    Modelica.Blocks.Interfaces.RealOutput d(unit = "kg/m3") 
      "Density in port medium" 
      annotation (extent=[-10,-120; 10,-100], rotation=-90);
    
  annotation (
    Diagram(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127}))),
    Icon(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127})),
        Text(extent=[-126,160; 138,98], string="%name"),
        Text(
          extent=[212,-51; 52,-103],
          style(color=0),
          string="d")),
    Documentation(info="<HTML>
<p>
This component monitors the density of the medium in the flow
between fluid ports. The sensor is 
ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    port_a.p = medium.p;
    h  = medium.h;
    Xi = medium.Xi;
    d  = medium.d;
  end Density;
  
  model Temperature "Ideal temperature sensor" 
    import SI = Modelica.SIunits;
    extends Sensors.BaseClasses.PartialFlowSensor;
    Medium.BaseProperties medium;
    Modelica.Blocks.Interfaces.RealOutput T(unit = "K") 
      "Temperature in port medium" 
      annotation (extent=[-10,-120; 10,-100], rotation=-90);
    
  annotation (
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
        Line(points=[0,-70; 0,-100], style(rgbcolor={0,0,127})),
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
        Line(points=[-100,0; -14,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[14,0; 100,0],   style(color=69, rgbcolor={0,128,255}))),
      Documentation(info="<HTML>
<p>
This component monitors the temperature of the medium in the flow
between fluid ports. The sensor is 
ideal, i.e., it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    port_a.p = medium.p;
    h  = medium.h;
    Xi = medium.Xi;
    T  = medium.T;
  end Temperature;
  
  model MassFlowRate "Ideal sensor for mass flow rate" 
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput m_flow(unit = "kg/s") 
      "mass flow rate from port_a to port_b" annotation (extent=[-10,-120; 10,
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
    Medium.BaseProperties medium;
    Modelica.Blocks.Interfaces.RealOutput V_flow(unit = "m3/s") 
      "volume flow rate from port_a to port_b" 
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
    port_a.p   = medium.p;
    h  = medium.h;
    Xi = medium.Xi;
    V_flow = port_a.m_flow/medium.d;
  end VolumeFlowRate;
  
  model SpecificEnthalpy "Ideal specific enthalphy sensor" 
    extends Sensors.BaseClasses.PartialFlowSensor;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput h_out(unit="J/kg") 
      "Specific enthalpy in port medium" 
      annotation (extent=[-10,-120; 10,-100], rotation=-90);
    
  annotation (
    Diagram(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
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
    h_out = h;
  end SpecificEnthalpy;
  
  model RelativePressure "Ideal relative pressure sensor" 
    extends Modelica.Icons.TranslationalSensor;
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
        choicesAllMatching = true);
    
    Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) 
      annotation (extent=[-120, -10; -100, 10]);
    Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) 
      annotation (extent=[120, -10; 100, 10]);
    
    Modelica.Blocks.Interfaces.RealOutput p_rel(redeclare type SignalType = 
          SI.Pressure) "Relative pressure signal" annotation (extent=[-10, -80; 10, -100], rotation=90);
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
    port_a.H_flow = 0;
    port_a.mXi_flow = zeros(Medium.nXi);
    port_b.m_flow = 0;
    port_b.H_flow = 0;
    port_b.mXi_flow = zeros(Medium.nXi);
    
    p_rel = port_a.p - port_b.p;
  end RelativePressure;
  
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
    partial model PartialAbsoluteSensor 
      "Partial component to model a sensor that measures a potential variable" 
      
      replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the sensor" annotation (
          choicesAllMatching =                                                                        true);
      Interfaces.FluidPort_a port(redeclare package Medium = Medium) 
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
      
      replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
          choicesAllMatching = true);
      Medium.SpecificEnthalpy h "enthalpy in flow";
      Medium.MassFraction[Medium.nXi] Xi "flow composition";
      
      Interfaces.FluidPort_a port_a(redeclare package Medium = Medium,
                         m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
        annotation (extent=[-110,-10; -90,10]);
      Interfaces.FluidPort_b port_b(redeclare package Medium = Medium,
                         m_flow(max=if allowFlowReversal then +Constants.inf else 0)) 
        annotation (extent=[110,-10; 90,10]);
      
      parameter Modelica_Fluid.Types.FlowDirection.Temp flowDirection=
                Modelica_Fluid.Types.FlowDirection.Unidirectional 
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
      
      replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
          choicesAllMatching =                                                                         true);
      
      Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) 
        annotation (extent=[-120, -10; -100, 10]);
      Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) 
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
