model TestMixingVolumes 
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  annotation (Diagram, experiment(StopTime=10));
  Modelica_Fluid.Volumes.MixingVolume MixingVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[-30,-40; -10,-20]);
  
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource2(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-100,-40; -80,-20]);
  Modelica_Fluid.Volumes.MixingVolume MixingVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[10,-40; 30,-20]);
  Modelica_Fluid.Sensors.Temperature Tmix_in(
                                         redeclare package Medium = Medium) 
    annotation (extent=[-60,-20; -40,0]);
  Modelica_Fluid.Sensors.Temperature Tmix_out(
                                          redeclare package Medium = Medium) 
    annotation (extent=[40,-20; 60,0]);
  Modelica_Fluid.Pipes.BaseClasses.PortVolume PortVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[-30,10; -10,30]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource1(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-100,10; -80,30]);
  Modelica_Fluid.Pipes.BaseClasses.PortVolume PortVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[10,10; 30,30]);
  Modelica_Fluid.Sensors.Temperature Tport_in(
                                          redeclare package Medium = Medium) 
    annotation (extent=[-60,30; -40,50]);
  Modelica_Fluid.Sensors.Temperature Tport_out(
                                           redeclare package Medium = Medium) 
    annotation (extent=[40,30; 60,50]);
  Modelica_Fluid.Sources.FixedBoundary_phX Sink1(   p=101325, redeclare package
      Medium = Medium,
    h=Medium.h_default) 
    annotation (extent=[100,10; 80,30]);
  Modelica_Fluid.Sources.FixedBoundary_phX Sink2(   p=101325, redeclare package
      Medium = Medium,
    h=Medium.h_default) 
    annotation (extent=[100,-40; 80,-20]);
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[-100,-100; -80,-80]);
equation 
  connect(MixingVolume1.port_b, MixingVolume2.port_a) annotation (points=[-10,-30;
        9.8,-30],                                     style(color=69, rgbcolor=
          {0,127,255}));
  connect(PortVolume1.port,PortVolume2. port) 
    annotation (points=[-20,20; 20,20],
                                     style(color=69, rgbcolor={0,127,255}));
  connect(PortVolume2.port, Sink1.port) 
    annotation (points=[20,20; 80,20], style(color=69, rgbcolor={0,127,255}));
  connect(PortVolume2.port, Tport_out.port) annotation (points=[20,20; 50,20; 50,
        30], style(color=69, rgbcolor={0,127,255}));
  connect(FlowSource1.port, PortVolume1.port) annotation (points=[-80,20; -20,
        20], style(color=69, rgbcolor={0,127,255}));
  connect(FlowSource1.port, Tport_in.port) annotation (points=[-80,20; -50,20;
        -50,30], style(color=69, rgbcolor={0,127,255}));
  connect(FlowSource2.port, MixingVolume1.port_a) annotation (points=[-80,-30;
        -30.2,-30], style(color=69, rgbcolor={0,127,255}));
  connect(FlowSource2.port, Tmix_in.port) annotation (points=[-80,-30; -50,-30;
        -50,-20], style(color=69, rgbcolor={0,127,255}));
  connect(MixingVolume2.port_b, Sink2.port) annotation (points=[30,-30; 80,-30],
      style(color=69, rgbcolor={0,127,255}));
  connect(MixingVolume2.port_b, Tmix_out.port) annotation (points=[30,-30; 50,
        -30; 50,-20], style(color=69, rgbcolor={0,127,255}));
end TestMixingVolumes;
