package Sensors 
  "Ideal sensor components (provide the variables in a fluid connector as signals)" 
  extends Modelica.Icons.Library;
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
  
  model Density "Ideal density sensor" 
    extends Interfaces.PartialAbsoluteSensor;
    extends Modelica.Icons.RotationalSensor;
    Medium.BaseProperties medium(p=port.p, h=port.h, X=port.X);
    parameter Modelica_Fluid.Types.DensityUnits signalUnit="kg/m3" 
      "Unit of density at output signal d";
    Modelica.Blocks.Interfaces.RealOutput d(unit = signalUnit) 
      "Density in port medium" annotation (extent=[100,-10; 120,10]);
    
  annotation (
    Diagram(
        Line(points=[70,0; 100,0], style(color=3, rgbcolor={0,0,255}))),
    Icon(
      Text(extent=[-126,160; 138,98],   string="%name"),
        Line(points=[70,0; 100,0], style(color=3, rgbcolor={0,0,255})),
        Text(
          extent=[212,-51; 52,-103],
          style(color=0),
          string="%signalUnit"),
        Line(points=[0,-70; 0,-100])),
    Documentation(info="<HTML>
<p>
This component monitors the density at the medium of the fluid port. The sensor is 
ideal, i.e., it does not influence the fluid and provides the value
of the port in the desired unit.
</p>
</HTML>
"));
  equation 
    if signalUnit == "kg/m3" then
       d = medium.d;
    elseif signalUnit == "g/cm3" then
       d = medium.d*1.e-3;
    else
       assert(false, "parameter signalUnit = \"" + signalUnit + "\" but must be " +
                     "kg/m3, or g/cm3");
    end if;
  end Density;
  
  model Pressure "Ideal pressure sensor" 
    import SI = Modelica.SIunits;
    extends Interfaces.PartialAbsoluteSensor;
    extends Modelica.Icons.RotationalSensor;
    parameter Modelica_Fluid.Types.PressureUnits signalUnit="Pa" 
      "Unit of pressure at output signal p";
    Modelica.Blocks.Interfaces.RealOutput p(unit = signalUnit) 
      "Pressure at port" 
      annotation (extent=[100,-10; 120,10]);
    
    annotation (
      Diagram(Line(points=[0, -70; 0, -100]),
        Line(points=[70,0; 100,0], style(color=3, rgbcolor={0,0,255}))),
      Icon(
        Line(points=[70,0; 100,0], style(color=3, rgbcolor={0,0,255})),
        Line(points=[0, -70; 0, -100]),
        Text(extent=[-159,134; 150,69],   string="%name"),
        Text(
          extent=[206,-60; 46,-112],
          style(color=0),
          string="%signalUnit")),
      Documentation(info="<HTML>
<p>
This component monitors the absolute pressure at its fluid port. The sensor is 
ideal, i.e., it does not influence the fluid and provides the value
of the port in the desired unit.
</p>
</HTML>
"),   Coordsys(grid=[1,1], component=[20,20]));
  equation 
    if signalUnit == "Pa" then
       p = port.p;
    elseif signalUnit == "kPa" then
       p = port.p*1e-3;
    elseif signalUnit == "bar" then
       p = SI.Conversions.to_bar(port.p);
    elseif signalUnit == "MPa" then
       p = port.p*1e-6;
    else
       assert(false, "parameter signalUnit = \"" + signalUnit + "\" but must be " +
                     "Pa, kPa, bar, or MPa");
    end if;
  end Pressure;
  
  model Temperature "Ideal temperature sensor" 
    import SI = Modelica.SIunits;
    extends Interfaces.PartialAbsoluteSensor;
    Medium.BaseProperties medium(p=port.p, h=port.h, X=port.X);
    parameter Modelica_Fluid.Types.TemperatureUnits signalUnit="K" 
      "Unit of temperature at output signal T";
    Modelica.Blocks.Interfaces.RealOutput T(unit = signalUnit) 
      "Temperature in port medium" 
                                 annotation (extent=[100,-10; 120,10]);
    
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
          string="%signalUnit"),
        Text(extent=[-126,160; 138,98],   string="%name")),
      Documentation(info="<HTML>
<p>
This component monitors the temperature at the medium of the fluid port. The sensor is 
ideal, i.e., it does not influence the fluid and provides the value
of the port in the desired unit.
</p>
</HTML>
"));
  equation 
    if signalUnit == "K" then
       T = medium.T;
    elseif signalUnit == "degC" then
       T = SI.Conversions.to_degC(medium.T);
    elseif signalUnit == "degF" then
       T = SI.Conversions.to_degF(medium.T);
    elseif signalUnit == "degR" then
       T = SI.Conversions.to_degRk(medium.T);
    else
       assert(false, "parameter signalUnit = \"" + signalUnit + "\" but must be " +
                     "K, degC, degF, or degR");
    end if;
  end Temperature;
  
  model MassFlowRate "Ideal sensor for mass flow rate" 
    extends Interfaces.PartialFlowRateSensor;
    extends Modelica.Icons.RotationalSensor;
    parameter Modelica_Fluid.Types.MassFlowRateUnits signalUnit="kg/s" 
      "Unit of mass flow rate at output signal m_flow";
      Modelica.Blocks.Interfaces.RealOutput m_flow(unit = signalUnit) 
      "mass flow rate from port_a to port_b" annotation (extent=[-10,-120; 10,
          -100], rotation=-90);
    
  annotation (
    Diagram(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100])),
    Icon(
      Text(extent=[-140,147; 140,82],   string="%name"),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Text(
          extent=[178,-81; 18,-133],
          style(color=0),
          string="%signalUnit"),
        Line(points=[0,-70; 0,-100]),
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255}))),
    Documentation(info="<HTML>
<p>
This component monitors the mass flow rate flowing from port_a to port_b. 
The sensor is ideal, i.e., it does not influence the fluid and provides the value
of the mass flow rate in the desired unit.
</p>
</HTML>
"));
  equation 
    if signalUnit == "kg/s" then
       m_flow = port_a.m_flow;
    elseif signalUnit == "t/h" then
       m_flow = port_a.m_flow*1.e3;
    else
       assert(false, "parameter signalUnit = \"" + signalUnit + "\" but must be " +
                     "kg/s, or t/h");
    end if;
  end MassFlowRate;
  
  model VolumeFlowRate "Ideal sensor for volume flow rate" 
    extends Interfaces.PartialFlowRateSensor;
    extends Modelica.Icons.RotationalSensor;
    Medium.BaseProperties medium(p=port_a.p, h=port_a.h, X=port_a.X);
    parameter Modelica_Fluid.Types.VolumeFlowRateUnits signalUnit="m3/s" 
      "Unit of volume flow rate at output signal V_flow";
    Modelica.Blocks.Interfaces.RealOutput V_flow(unit = signalUnit) 
      "volume flow rate from port_a to port_b" annotation (extent=[-10,-120; 10,
          -100], rotation=-90);
    
  annotation (
    Diagram(
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[0,-70; 0,-100])),
    Icon(
      Text(extent=[-126,160; 138,98],   string="%name"),
        Text(
          extent=[188,-71; 28,-123],
          style(color=0),
          string="%signalUnit"),
        Line(points=[0,-70; 0,-100]),
        Line(points=[-100,0; -70,0], style(color=69, rgbcolor={0,128,255})),
        Line(points=[70,0; 100,0], style(color=69, rgbcolor={0,128,255}))),
    Documentation(info="<HTML>
<p>
This component monitors the volume flow rate flowing from port_a to port_b. 
The sensor is ideal, i.e., it does not influence the fluid and provides the value
of the volume flow rate in the desired unit.
</p>
</HTML>
"));
  equation 
    if signalUnit == "m3/s" then
       V_flow = port_a.m_flow/medium.d;
    elseif signalUnit == "l/s" then
       V_flow = port_a.m_flow*1.e-3/medium.d;
    else
       assert(false, "parameter signalUnit = \"" + signalUnit + "\" but must be " +
                     "m3/s, or l/s");
    end if;
  end VolumeFlowRate;
end Sensors;
