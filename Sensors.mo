package Sensors "Sensor components" 
  extends Modelica.Icons.Library;
  import SI = Modelica.SIunits;
  
  model TemperatureSensor "Absolute temperature sensor" 
    extends Interfaces.PartialOnePort;
    annotation (
      Diagram(
        Ellipse(extent=[-20, -98; 20, -60], style(
            color=0, 
            thickness=2, 
            fillColor=42)), 
        Rectangle(extent=[-12, 40; 12, -68], style(color=42, fillColor=42)), 
        Line(points=[12, 0; 90, 0]), 
        Polygon(points=[-12, 40; -12, 80; -10, 86; -6, 88; 0, 90; 6, 88; 10, 86; 
               12, 80; 12, 40; -12, 40], style(color=0, thickness=2)), 
        Line(points=[-12, 40; -12, -64], style(color=0, thickness=2)), 
        Line(points=[12, 40; 12, -64], style(color=0, thickness=2)), 
        Line(points=[-40, -20; -12, -20], style(color=0)), 
        Line(points=[-40, 20; -12, 20], style(color=0)), 
        Line(points=[-40, 60; -12, 60], style(color=0)), 
        Text(
          extent=[100, -40; 36, -104], 
          string="°C", 
          style(color=0))), 
      Icon(
        Ellipse(extent=[-20, -98; 20, -60], style(
            color=0, 
            thickness=2, 
            fillColor=42)), 
        Rectangle(extent=[-12, 40; 12, -68], style(color=42, fillColor=42)), 
        Line(points=[12, 0; 90, 0]), 
        Polygon(points=[-12, 40; -12, 80; -10, 86; -6, 88; 0, 90; 6, 88; 10, 86; 
               12, 80; 12, 40; -12, 40], style(color=0, thickness=2)), 
        Line(points=[-12, 40; -12, -64], style(color=0, thickness=2)), 
        Line(points=[12, 40; 12, -64], style(color=0, thickness=2)), 
        Line(points=[-40, -20; -12, -20], style(color=0)), 
        Line(points=[-40, 20; -12, 20], style(color=0)), 
        Line(points=[-40, 60; -12, 60], style(color=0)), 
        Text(
          extent=[88, -40; 24, -104], 
          string="°C", 
          style(color=0)), 
        Text(extent=[-132, 144; 108, 84], string="%name")), 
      Documentation(info="<HTML>
<p>
This model monitors the temperature at its port. The sensor is 
ideal, i.e. it does not influence the fluid. The measurement unit is 
degree centigrade.
</p>
</HTML>
"));
    Modelica.Blocks.Interfaces.OutPort signal(redeclare type SignalType = 
          SI.CelsiusTemperature) annotation (extent=[90, -10; 110, 10]);
    Medium.BaseProperties medium(p=port.p, h=port.h);
  equation 
    signal.signal[1] = SI.Conversions.to_degC(medium.T);
    port.m_dot = 0;
    port.H_dot = 0;
    port.mX_dot = zeros(Medium.nX);
  end TemperatureSensor;
  
  model PressureSensor "Pressure sensor" 
    extends Interfaces.PartialOnePort;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.OutPort signal(redeclare type SignalType = 
          SI.Conversions.NonSIunits.Pressure_bar) "Pressure at port"
      annotation (extent=[90, -10; 110, 10]);
    annotation (
      Diagram(Line(points=[0, -70; 0, -100]), Line(points=[69, 0; 90, 0], style(
              color=42))), 
      Icon(
        Text(
          extent=[-50, -6; 58, -68], 
          string="bar", 
          style(color=0)), 
        Line(points=[69, 0; 90, 0], style(color=42)), 
        Line(points=[0, -70; 0, -100]), 
        Text(extent=[-132, 144; 108, 84], string="%name")), 
      Documentation(info="<HTML>
<p>
This model monitors the pressure at its port. The sensor is 
ideal, i.e. it does not influence the fluid. The measurement unit is bar.
</p>
</HTML>
"));
  equation 
    signal.signal[1] = SI.Conversions.to_bar(port.p);
    port.m_dot = 0;
    port.H_dot = 0;
    port.mX_dot = zeros(Medium.nX);
  end PressureSensor;
  
  model MassFlowSensor "Mass flow rate sensor" 
    extends Interfaces.PartialTwoPortTransport;
    extends Modelica.Icons.RotationalSensor;
    SI.Pressure dp "Pressure loss due to friction";
    Real residue=port_a.p - port_b.p - dp "momentum balance (may be modified)";
    Modelica.Blocks.Interfaces.OutPort signal(redeclare type SignalType = 
          SI.MassFlowRate) "Mass flow from port_a -> port_b"
      annotation (extent=[-10, -110; 10, -90], rotation=270);
    annotation (
      Diagram(
        Line(points=[-70, 0; -100, 0]), 
        Line(points=[0, -70; 0, -90]), 
        Line(points=[71, 0; 100, 0])), 
      Icon(
        Text(
          extent=[21, -58; 76, -116], 
          string="m_dot", 
          style(color=0)), 
        Line(points=[-70, 0; -100, 0]), 
        Line(points=[70, 0; 100, 0]), 
        Line(points=[0, -70; 0, -90]), 
        Text(extent=[-132, 144; 108, 84], string="%name"), 
        Text(
          extent=[-42, 0; 44, -70], 
          style(color=0), 
          string="kg/s")), 
      Documentation(info="<HTML>
<p>
This model monitors the mass flow rate flowing through it. The sensor is 
ideal, i.e. it does not influence the flow. The output signal is positive 
if the fluid flows from port_a to port_b. The measurement unit is kg/s.
</p>
</HTML>
"));
  equation 
    residue = 0;
    dp = 0;
    signal.signal[1] = m_dot;
  end MassFlowSensor;
  
end Sensors;
