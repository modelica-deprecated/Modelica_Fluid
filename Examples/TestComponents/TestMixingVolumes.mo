model TestMixingVolumes 
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  annotation (Diagram, experiment(StopTime=10));
  Modelica_Fluid.Components.MixingVolume MixingVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    initOption=Modelica_Fluid.Types.InitTypes.InitialValues) 
                 annotation (extent=[-30,-70; -10,-50]);
  
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource2(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-100,-70; -80,-50]);
  Modelica_Fluid.Components.MixingVolume MixingVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    initOption=Modelica_Fluid.Types.InitTypes.InitialValues) 
                 annotation (extent=[10,-70; 30,-50]);
  Sensors.Temperature Tmix_in(redeclare package Medium = Medium) 
    annotation (extent=[-60,-70; -40,-50]);
  Sensors.Temperature Tmix_out(redeclare package Medium = Medium) 
    annotation (extent=[40,-70; 60,-50]);
  Components.PortVolume PortVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    initOption=Modelica_Fluid.Types.InitTypes.InitialValues) 
                 annotation (extent=[-30,10; -10,30]);
  Sources.PrescribedMassFlowRate_hX FlowSource1(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-100,10; -80,30]);
  Components.PortVolume PortVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    initOption=Modelica_Fluid.Types.InitTypes.InitialValues) 
                 annotation (extent=[10,10; 30,30]);
  Sensors.Temperature Tport_in(redeclare package Medium = Medium) 
    annotation (extent=[-60,10; -40,30]);
  Sensors.Temperature Tport_out(redeclare package Medium = Medium) 
    annotation (extent=[40,10; 60,30]);
  Sources.FixedAmbient_phX Sink1(        p=101325, redeclare package Medium = 
               Medium) 
    annotation (extent=[100,10; 80,30]);
  Sources.FixedAmbient_phX Sink2(        p=101325, redeclare package Medium = 
               Medium) 
    annotation (extent=[100,-70; 80,-50]);
equation 
  connect(MixingVolume1.port_b, MixingVolume2.port_a) annotation (points=[-10,-60; 
        9.8,-60],                                     style(color=69, rgbcolor=
          {0,127,255}));
  connect(PortVolume1.port,PortVolume2. port) 
    annotation (points=[-20,20; 20,20],
                                     style(color=69, rgbcolor={0,127,255}));
  connect(FlowSource1.port, Tport_in.port_a) annotation (points=[-79,20; -61,20], 
      style(color=69, rgbcolor={0,127,255}));
  connect(Tport_in.port_b, PortVolume1.port) annotation (points=[-39,20; -20,20], 
      style(color=69, rgbcolor={0,127,255}));
  connect(PortVolume2.port, Tport_out.port_a)
    annotation (points=[20,20; 39,20], style(color=69, rgbcolor={0,127,255}));
  connect(Tport_out.port_b, Sink1.port)
    annotation (points=[61,20; 79,20], style(color=69, rgbcolor={0,127,255}));
  connect(FlowSource2.port, Tmix_in.port_a) annotation (points=[-79,-60; -61,
        -60], style(color=69, rgbcolor={0,127,255}));
  connect(Tmix_in.port_b, MixingVolume1.port_a) annotation (points=[-39,-60; 
        -30.2,-60], style(color=69, rgbcolor={0,127,255}));
  connect(MixingVolume2.port_b, Tmix_out.port_a) annotation (points=[30,-60; 39,
        -60], style(color=69, rgbcolor={0,127,255}));
  connect(Tmix_out.port_b, Sink2.port) annotation (points=[61,-60; 79,-60], 
      style(color=69, rgbcolor={0,127,255}));
end TestMixingVolumes;
