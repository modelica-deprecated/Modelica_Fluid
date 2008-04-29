within Modelica_Fluid.Test.TestComponents.Sensors;
model TestOnePortSensors1 
  import Modelica_Fluid;
  package Medium = Modelica.Media.Water.StandardWater;
  parameter Real D_a = 0.1;
  parameter Real D_b = 0.2;
  parameter Real A_rel = (D_a/D_b)^2;
  parameter Real zeta =  (1 - A_rel)^2;
  
  annotation (Diagram, experiment(StopTime=25, Algorithm="Dassl"), 
    experimentSetupOutput);
  Modelica_Fluid.Volumes.MixingVolume volume1(
    V=1e-3,
    use_T_start=false,
    redeclare package Medium = Medium,
    h_start=1e5,
    initType=Modelica_Fluid.Types.Init.InitialValues, 
    p_start=101325) 
                 annotation (extent=[-30,30; -10,50]);
  
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource2(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium,
    useFlowRateInput=true) 
                   annotation (extent=[-68,30; -48,50]);
  Modelica_Fluid.Sensors.TemperatureOnePort Tmix1(redeclare package Medium = 
        Medium) 
    annotation (extent=[0,60; 20,80]);
  Modelica_Fluid.Sources.FixedBoundary_phX sink1(             redeclare package
      Medium = Medium, 
    h=5e4, 
    p=101325) 
    annotation (extent=[100,30; 80,50]);
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[-100,-100; -80,-80]);
  Modelica.Blocks.Sources.Ramp ramp(
    height=2,
    duration=20,
    offset=-1) annotation (extent=[-100,30; -80,50]);
  Modelica_Fluid.Volumes.MixingVolume volume2(
    V=1e-3,
    use_T_start=false,
    redeclare package Medium = Medium,
    h_start=1e5, 
    initType=Modelica_Fluid.Types.Init.InitialValues, 
    p_start=101325) 
                 annotation (extent=[-32,-30; -12,-10]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource1(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium,
    useFlowRateInput=true) 
                   annotation (extent=[-68,-30; -48,-10]);
  Modelica_Fluid.Sources.FixedBoundary_phX sink2(             redeclare package
      Medium = Medium, 
    h=5e4, 
    p=101325) 
    annotation (extent=[100,-30; 80,-10]);
  Modelica_Fluid.Sensors.TemperatureTwoPort Tmix2(redeclare package Medium = 
        Medium) annotation (extent=[0,-30; 20,-10]);
  Modelica_Fluid.PressureLosses.SimpleGenericOrifice orifice1(
    redeclare package Medium = Medium, 
    diameter=D_a, 
    zeta=zeta) annotation (extent=[40,30; 60,50]);
  Modelica_Fluid.PressureLosses.SimpleGenericOrifice orifice2(
    redeclare package Medium = Medium, 
    zeta=zeta, 
    diameter=D_a) annotation (extent=[40,-30; 60,-10]);
equation 
  connect(FlowSource2.port, volume1.port_a)       annotation (points=[-48,40; 
        -30.2,40],  style(color=69, rgbcolor={0,127,255}));
  connect(Tmix1.port, volume1.port_b) annotation (points=[10,60; 10,40; -10,40], 
      style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(ramp.y, FlowSource2.m_flow_in) annotation (points=[-79,40; -74,40; 
        -74,46; -67.3,46], style(
      color=74, 
      rgbcolor={0,0,127}, 
      smooth=0));
  connect(FlowSource1.port, volume2.port_a)       annotation (points=[-48,-20; 
        -32.2,-20], style(color=69, rgbcolor={0,127,255}));
  connect(ramp.y, FlowSource1.m_flow_in) annotation (points=[-79,40; -76,40; 
        -76,-14; -67.3,-14], style(
      color=74, 
      rgbcolor={0,0,127}, 
      smooth=0));
  connect(Tmix2.port_a, volume2.port_b) annotation (points=[0,-20; -12,-20], 
      style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(orifice1.port_a, volume1.port_b) annotation (points=[40,40; -10,40], 
      style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(orifice1.port_b,sink1. port) annotation (points=[60,40; 80,40], style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(orifice2.port_a, Tmix2.port_b) annotation (points=[40,-20; 20,-20], 
      style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(orifice2.port_b,sink2. port) annotation (points=[60,-20; 80,-20], 
      style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
end TestOnePortSensors1;
