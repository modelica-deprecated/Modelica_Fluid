package Sensors "Sensor components" 
  extends Modelica.Icons.Library;
  import SI = Modelica.SIunits;
  
  model TemperatureSensor "Absolute temperature sensor" 
    extends Interfaces.PartialOnePort;
    Modelica.Blocks.Interfaces.OutPort signal(redeclare type SignalType = 
          SI.CelsiusTemperature) annotation (extent=[90, -10; 110, 10]);
    Medium.BaseProperties medium(p=port.p, h=port.h);
    parameter TemperatureUnits.Temp measurementUnit =
      TemperatureUnits.Kelvin "Measurement unit of sensor signal" 
      annotation( Evaluate=true);
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
          string="T", 
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
          string="T", 
          style(color=0)), 
        Text(extent=[-132, 144; 108, 84], string="%name")), 
      Documentation(info="<HTML>
<p>
This model monitors the temperature at its port. The sensor is 
ideal, i.e. it does not influence the fluid.
</p>
</HTML>
"));
  equation 
    signal.signal[1] =
      if measurementUnit == TemperatureUnits.Kelvin then medium.T
      else if measurementUnit == TemperatureUnits.Celsius then SI.Conversions.to_degC(medium.T)
      else if measurementUnit == TemperatureUnits.Fahrenheit then SI.Conversions.to_degF(medium.T)
      else if measurementUnit == TemperatureUnits.Rankine then SI.Conversions.to_degRk(medium.T)
      else medium.T; // else should be Undefined
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
    parameter PressureUnits.Temp measurementUnit =
      PressureUnits.Pascal "Measurement unit of sensor signal" 
      annotation( Evaluate=true);
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
    signal.signal[1] =
      if measurementUnit == PressureUnits.Pascal then port.p
      else if measurementUnit == PressureUnits.Bar then SI.Conversions.to_bar(port.p)
      else if measurementUnit == PressureUnits.KiloPascal then port.p*1e-3
      else if measurementUnit == PressureUnits.MegaPascal then port.p*1e-6
      else port.p; // else should be Undefined
    port.m_dot = 0;
    port.H_dot = 0;
    port.mX_dot = zeros(Medium.nX);
  end PressureSensor;
  
  model MassFlowSensor "Mass flow rate sensor" 
    extends Interfaces.PartialTwoPortTransport;
    extends Modelica.Icons.RotationalSensor;
    Modelica.Blocks.Interfaces.OutPort signal(redeclare type SignalType = 
          SI.MassFlowRate) "Mass flow from port_a -> port_b"
      annotation (extent=[-10, -110; 10, -90], rotation=270);
    parameter MassFlowRateUnits.Temp measurementUnit =
      MassFlowRateUnits.KilogramPerSecond "Measurement unit of sensor signal" 
      annotation( Evaluate=true);
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
          string="m_dot", 
          style(color=0)), 
        Line(points=[-70, 0; -100, 0]), 
        Line(points=[70, 0; 100, 0]), 
        Line(points=[0, -70; 0, -90]), 
        Text(extent=[-132, 144; 108, 84], string="%name"), 
        Text(
          extent=[-42, 0; 44, -70], 
          style(color=0), 
          string="m_dot")), 
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
    signal.signal[1] =
      if measurementUnit == MassFlowRateUnits.KilogramPerSecond then m_dot
      else if measurementUnit == MassFlowRateUnits.TonPerHour then m_dot*3.6
      else m_dot; // else should be Undefined
  end MassFlowSensor;
  
  package TemperatureUnits 
    "Type, constants and menu choices for temperature units, as temporary solution until enumerations will be available"
    
    annotation( preferedView="text");
    
    extends Modelica.Icons.Library;
    constant Integer Kelvin=1;
    constant Integer Celsius=2;
    constant Integer Fahrenheit=3;
    constant Integer Rankine=4;
    type Temp 
      "Temporary type of TemperatureUnits with choices for menus (until enumerations will be available)"
      
      extends Integer;
      annotation( Evaluate=true,
                 choices(
                         choice=Modelica_Fluid.Sensors.TemperatureUnits.Kelvin
                         "K",
                         choice=Modelica_Fluid.Sensors.TemperatureUnits.Celsius
                         "°C",
                         choice=Modelica_Fluid.Sensors.TemperatureUnits.Fahrenhit
                         "°F",
                         choice=Modelica_Fluid.Sensors.TemperatureUnits.Rankine
                         "°Rk"));
    end Temp;
  end TemperatureUnits;

  package PressureUnits 
    "Type, constants and menu choices for pressure units, as temporary solution until enumerations will be available"
    
    annotation( preferedView="text");
    
    extends Modelica.Icons.Library;
    constant Integer Pascal=1;
    constant Integer Bar=2;
    constant Integer KiloPascal=3;
    constant Integer MegaPascal=4;
    type Temp 
      "Temporary type of PressureUnits with choices for menus (until enumerations will be available)"
      
      extends Integer;
      annotation( Evaluate=true,
                 choices(
                         choice=Modelica_Fluid.Sensors.PressureUnits.Pascal
                         "Pa",
                         choice=Modelica_Fluid.Sensors.PressureUnits.Bar
                         "bar",
                         choice=Modelica_Fluid.Sensors.PressureUnits.KiloPascal
                         "kPa",
                         choice=Modelica_Fluid.Sensors.PressureUnits.MegaPascal
                         "MPa"));
    end Temp;
  end PressureUnits;

  package MassFlowRateUnits 
    "Type, constants and menu choices for mass flow rate units, as temporary solution until enumerations will be available"
    
    annotation( preferedView="text");
    
    extends Modelica.Icons.Library;
    constant Integer KilogramPerSecond=1;
    constant Integer TonPerHour=2;
    type Temp 
      "Temporary type of MassFlowRateUnits with choices for menus (until enumerations will be available)"
      
      extends Integer;
      annotation( Evaluate=true,
                 choices(
                         choice=Modelica_Fluid.Sensors.MassFlowRateUnits.KilogramPerSecond
                         "kg/s",
                         choice=Modelica_Fluid.Sensors.MassFlowRateUnits.TonPerHour
                         "t/h"));
    end Temp;
  end MassFlowRateUnits;

end Sensors;
