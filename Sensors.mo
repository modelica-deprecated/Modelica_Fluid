package Sensors "Sensor components" 
  extends Modelica.Icons.Library;
  import SI = Modelica.SIunits;
  
  
  model PressureSensor "Pressure sensor" 
    extends Interfaces.PartialOnePort;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput signal(  redeclare type SignalType = 
          SI.Conversions.NonSIunits.Pressure_bar) "Pressure at port" 
      annotation (extent=[90, -10; 110, 10]);
    parameter Types.PressureUnits.Temp measurementUnit=
      Types.PressureUnits.Pascal "Measurement unit of sensor signal" 
      annotation (Evaluate=true);
    annotation (
      Diagram(Line(points=[0, -70; 0, -100]), Line(points=[69, 0; 90, 0], style(
              color=42))),
      Icon(
        Text(
          extent=[-50, -6; 58, -68],
          string="p",
          style(color=0)),
        Line(points=[69, 0; 90, 0], style(color=42)),
        Line(points=[0, -70; 0, -100]),
        Text(extent=[-132, 144; 108, 84], string="%name")),
      Documentation(info="<HTML>
<p>
This model monitors the pressure at its port. The sensor is 
ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    signal =
      if measurementUnit == Types.PressureUnits.Pascal then port.p else 
           if measurementUnit == Types.PressureUnits.Bar then SI.Conversions.to_bar(port.p) else 
           if measurementUnit == Types.PressureUnits.KiloPascal then port.p*1e-3 else 
           if measurementUnit == Types.PressureUnits.MegaPascal then port.p*1e-6 else 
           port.p; // else should be Undefined
    port.m_flow = 0;
    port.H_flow = 0;
    port.mX_flow = zeros(Medium.nX);
  end PressureSensor;
  
  model MassFlowSensor "Mass flow rate sensor" 
    extends Interfaces.PartialTwoPortTransport;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.RealOutput signal(
                                              redeclare type SignalType = 
          SI.MassFlowRate) "Mass flow from port_a -> port_b" 
      annotation (extent=[-10, -110; 10, -90], rotation=270);
    parameter Types.MassFlowRateUnits.Temp measurementUnit=
      Types.MassFlowRateUnits.KilogramPerSecond 
      "Measurement unit of sensor signal" 
      annotation (Evaluate=true);
    SI.Pressure dp "Pressure loss due to friction";
    Real residue=port_a.p - port_b.p - dp "momentum balance (may be modified)";
    annotation (
      Diagram(
        Line(points=[-70, 0; -100, 0]),
        Line(points=[0, -70; 0, -90]),
        Line(points=[71, 0; 100, 0])),
      Icon(
        Text(
          extent=[21, -58; 76, -116],
          string="m_flow",
          style(color=0)),
        Line(points=[-70, 0; -100, 0]),
        Line(points=[70, 0; 100, 0]),
        Line(points=[0, -70; 0, -90]),
        Text(extent=[-132, 144; 108, 84], string="%name"),
        Text(
          extent=[-42, 0; 44, -70],
          style(color=0),
          string="m_flow")),
      Documentation(info="<HTML>
<p>
This model monitors the mass flow rate flowing through it. The sensor is 
ideal, i.e. it does not influence the flow. The output signal is positive 
if the fluid flows from port_a to port_b.
</p>
</HTML>
"));
  equation 
    residue = 0;
    dp = 0;
    signal           =
      if measurementUnit == Types.MassFlowRateUnits.KilogramPerSecond then m_flow else 
           if measurementUnit == Types.MassFlowRateUnits.TonPerHour then m_flow*3.6 else 
           m_flow; // else should be Undefined
  end MassFlowSensor;
  
  model Temperature "Absolute temperature sensor" 
    extends Interfaces.PartialOnePort;
    Modelica.Blocks.Interfaces.RealOutput T(unit = Unit) 
                                 annotation (extent=[100,-10; 120,10]);
    Medium.BaseProperties medium(p=port.p, h=port.h);
    parameter Modelica_Fluid.Types.TemperatureUnits Unit="K" 
      "Measurement unit of sensor signal" 
      annotation (Evaluate=true);
    annotation (
      Diagram(
        Ellipse(extent=[-20, -98; 20, -60], style(
            color=0,
            thickness=2,
            fillColor=42)),
        Rectangle(extent=[-12, 40; 12, -68], style(color=42, fillColor=42)),
        Line(points=[12,0; 100,0]),
        Polygon(points=[-12, 40; -12, 80; -10, 86; -6, 88; 0, 90; 6, 88; 10, 86;
               12, 80; 12, 40; -12, 40], style(color=0, thickness=2)),
        Line(points=[-12, 40; -12, -64], style(color=0, thickness=2)),
        Line(points=[12, 40; 12, -64], style(color=0, thickness=2)),
        Line(points=[-40, -20; -12, -20], style(color=0)),
        Line(points=[-40, 20; -12, 20], style(color=0)),
        Line(points=[-40, 60; -12, 60], style(color=0))),
      Icon(
        Ellipse(extent=[-20, -98; 20, -60], style(
            color=0,
            thickness=2,
            fillColor=42)),
        Rectangle(extent=[-12, 40; 12, -68], style(color=42, fillColor=42)),
        Line(points=[12,0; 100,0]),
        Polygon(points=[-12, 40; -12, 80; -10, 86; -6, 88; 0, 90; 6, 88; 10, 86;
               12, 80; 12, 40; -12, 40], style(color=0, thickness=2)),
        Line(points=[-12, 40; -12, -64], style(color=0, thickness=2)),
        Line(points=[12, 40; 12, -64], style(color=0, thickness=2)),
        Line(points=[-40, -20; -12, -20], style(color=0)),
        Line(points=[-40, 20; -12, 20], style(color=0)),
        Line(points=[-40, 60; -12, 60], style(color=0)),
        Text(
          extent=[180,-28; 20,-80],
          style(color=0),
          string="%Unit"),
        Text(extent=[-126,160; 138,98],   string="%name")),
      Documentation(info="<HTML>
<p>
This model monitors the temperature at its port. The sensor is 
ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    if Unit == "K" then
       T = medium.T;
    elseif Unit == "degC" then
       T = SI.Conversions.to_degC(medium.T);
    elseif Unit == "degF" then
       T = SI.Conversions.to_degF(medium.T);
    elseif Unit == "degR" then
       T = SI.Conversions.to_degRk(medium.T);
    else
       assert(false, "parameter Unit = \"" + Unit + "\" but must be " +
                     "K, degC, degF, or degR");
    end if;
    port.m_flow = 0;
    port.H_flow = 0;
    port.mX_flow = zeros(Medium.nX);
  end Temperature;

  model RelativeTemperature "Absolute temperature sensor" 
    extends Interfaces.PartialOnePort;
    Modelica.Blocks.Interfaces.RealOutput T(unit = Unit) 
                                 annotation (extent=[100,-10; 120,10]);
    Medium.BaseProperties medium(p=port.p, h=port.h);
    parameter Modelica_Fluid.Types.TemperatureUnits Unit="K" 
      "Measurement unit of sensor signal" 
      annotation (Evaluate=true);
    annotation (
      Diagram(
        Ellipse(extent=[-20, -98; 20, -60], style(
            color=0,
            thickness=2,
            fillColor=42)),
        Rectangle(extent=[-12, 40; 12, -68], style(color=42, fillColor=42)),
        Line(points=[12,0; 100,0]),
        Polygon(points=[-12, 40; -12, 80; -10, 86; -6, 88; 0, 90; 6, 88; 10, 86;
               12, 80; 12, 40; -12, 40], style(color=0, thickness=2)),
        Line(points=[-12, 40; -12, -64], style(color=0, thickness=2)),
        Line(points=[12, 40; 12, -64], style(color=0, thickness=2)),
        Line(points=[-40, -20; -12, -20], style(color=0)),
        Line(points=[-40, 20; -12, 20], style(color=0)),
        Line(points=[-40, 60; -12, 60], style(color=0))),
      Icon(
        Ellipse(extent=[-20, -98; 20, -60], style(
            color=0,
            thickness=2,
            fillColor=42)),
        Rectangle(extent=[-12, 40; 12, -68], style(color=42, fillColor=42)),
        Line(points=[12,0; 100,0]),
        Polygon(points=[-12, 40; -12, 80; -10, 86; -6, 88; 0, 90; 6, 88; 10, 86;
               12, 80; 12, 40; -12, 40], style(color=0, thickness=2)),
        Line(points=[-12, 40; -12, -64], style(color=0, thickness=2)),
        Line(points=[12, 40; 12, -64], style(color=0, thickness=2)),
        Line(points=[-40, -20; -12, -20], style(color=0)),
        Line(points=[-40, 20; -12, 20], style(color=0)),
        Line(points=[-40, 60; -12, 60], style(color=0)),
        Text(
          extent=[180,-28; 20,-80],
          style(color=0),
          string="%Unit"),
        Text(extent=[-126,160; 138,98],   string="%name")),
      Documentation(info="<HTML>
<p>
This model monitors the temperature at its port. The sensor is 
ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    if Unit == "K" then
       T = medium.T;
    elseif Unit == "degC" then
       T = SI.Conversions.to_degC(medium.T);
    elseif Unit == "degF" then
       T = SI.Conversions.to_degF(medium.T);
    elseif Unit == "degR" then
       T = SI.Conversions.to_degRk(medium.T);
    else
       assert(false, "parameter Unit = \"" + Unit + "\" but must be " +
                     "K, degC, degF, or degR");
    end if;
    port.m_flow = 0;
    port.H_flow = 0;
    port.mX_flow = zeros(Medium.nX);
  end RelativeTemperature;
end Sensors;
