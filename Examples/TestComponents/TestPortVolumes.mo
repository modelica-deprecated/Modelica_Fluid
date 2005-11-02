model TestPortVolumes 
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  Modelica_Fluid.Components.PortVolume PortVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[-30,-10; -10,10]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource1(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-100,-10; -80,10]);
  Modelica_Fluid.Sources.FixedAmbient_phX Sink1(
                                         p=101325, redeclare package Medium = 
               Medium) 
    annotation (extent=[100,-10; 80,10]);
  Modelica_Fluid.Components.PortVolume PortVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[10,-10; 30,10]);
  Sensors.Temperature Tin(redeclare package Medium = Medium) 
    annotation (extent=[-60,-10; -40,10]);
  Sensors.Temperature Tout(redeclare package Medium = Medium) 
    annotation (extent=[40,-10; 60,10]);
equation 
  connect(PortVolume1.port, PortVolume2.port) 
    annotation (points=[-20,0; 20,0],style(color=69, rgbcolor={0,127,255}));
  annotation (Diagram);
  connect(FlowSource1.port, Tin.port_a) 
    annotation (points=[-79,0; -61,0], style(color=69, rgbcolor={0,127,255}));
  connect(Tin.port_b, PortVolume1.port) 
    annotation (points=[-39,0; -20,0], style(color=69, rgbcolor={0,127,255}));
  connect(PortVolume2.port, Tout.port_a) 
    annotation (points=[20,0; 39,0], style(color=69, rgbcolor={0,127,255}));
  connect(Tout.port_b, Sink1.port) 
    annotation (points=[61,0; 79,0], style(color=69, rgbcolor={0,127,255}));
end TestPortVolumes;
