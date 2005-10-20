model TestPortVolumes 
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  Modelica_Fluid.Components.PortVolume PortVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[-40,-6; -20,14]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource1(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-82,-6; -62,14]);
  Modelica_Fluid.Sources.FixedAmbient_phX Sink1(
                                         p=101325, redeclare package Medium = 
               Medium) 
    annotation (extent=[86,-8; 60,16]);
  Modelica_Fluid.Components.PortVolume PortVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[-6,-6; 14,14]);
  Sensors.Temperature Tin(redeclare package Medium = Medium) 
    annotation (extent=[-60,34; -40,54]);
  Sensors.Temperature Tout(redeclare package Medium = Medium) 
    annotation (extent=[28,34; 48,54]);
equation 
  connect(FlowSource1.port, PortVolume1.port) 
    annotation (points=[-61,4; -30,4], style(color=69, rgbcolor={0,127,255}));
  connect(PortVolume1.port, PortVolume2.port) 
    annotation (points=[-30,4; 4,4], style(color=69, rgbcolor={0,127,255}));
  annotation (Diagram);
  connect(PortVolume2.port, Sink1.port) 
    annotation (points=[4,4; 58.7,4],
                                    style(color=69, rgbcolor={0,127,255}));
  connect(Tout.port, PortVolume2.port) annotation (points=[38,33; 38,4; 4,4],
      style(color=69, rgbcolor={0,127,255}));
  connect(Tin.port, FlowSource1.port) annotation (points=[-50,33; -52,33; -52,4;
        -61,4], style(color=69, rgbcolor={0,127,255}));
end TestPortVolumes;
