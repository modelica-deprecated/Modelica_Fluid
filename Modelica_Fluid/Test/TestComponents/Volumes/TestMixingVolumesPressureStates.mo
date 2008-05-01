within Modelica_Fluid.Test.TestComponents.Volumes;
model TestMixingVolumesPressureStates 
  "Test case where in one of the mixing volumes a pressure state appears" 
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  annotation (Diagram, experiment(StopTime=10));
  Modelica_Fluid.Volumes.MixingVolume MixingVolume1(
    V=1e-3,
    redeclare package Medium = Medium,
    p_start=ambient.default_p_ambient,
    use_T_start=true,
    T_start=ambient.default_T_ambient,
    initType=Modelica_Fluid.Types.Init.SteadyState) 
                 annotation (extent=[-30,30; -10,50]);
  
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource2(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-100,30; -80,50]);
  Modelica_Fluid.Volumes.MixingVolume MixingVolume2(
    V=1e-3,
    p_start=ambient.default_p_ambient,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    initType=Modelica_Fluid.Types.Init.NoInit) 
                 annotation (extent=[10,30; 30,50]);
  Modelica_Fluid.Sensors.TemperatureOnePort Tmix_in(
                                         redeclare package Medium = Medium) 
    annotation (extent=[-58,68; -38,88]);
  Modelica_Fluid.Sensors.TemperatureOnePort Tmix_out(
                                          redeclare package Medium = Medium) 
    annotation (extent=[30,68; 50,88]);
  Modelica_Fluid.Sources.FixedBoundary_phX Sink2(   p=101325, redeclare package
      Medium = Medium,
    h=Medium.h_default) 
    annotation (extent=[100,30; 80,50]);
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[-100,-100; -80,-80]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity simpleGenericOrifice2(
    redeclare package Medium = Medium,
    diameter=0.2,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Laminar,
    length=1)                   annotation (extent=[50,30; 70,50]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity simpleGenericOrifice1(
    redeclare package Medium = Medium,
    diameter=0.2,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Laminar,
    length=1)                   annotation (extent=[-70,30; -50,50]);
equation 
  connect(simpleGenericOrifice2.port_b, Sink2.port) annotation (points=[70,40;
        80,40], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(MixingVolume2.port_b, simpleGenericOrifice2.port_a) annotation (
      points=[30,40; 50,40], style(
      color=0,
      rgbcolor={0,0,0},
      smooth=0));
  connect(MixingVolume2.port_b, Tmix_out.port) annotation (points=[30,40; 50,40;
        50,68; 40,68], style(
      color=0,
      rgbcolor={0,0,0},
      smooth=0));
  connect(FlowSource2.port,simpleGenericOrifice1. port_a) annotation (points=[
        -80,40; -70,40], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(simpleGenericOrifice1.port_b, MixingVolume1.port_a) annotation (
      points=[-50,40; -30.2,40], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(Tmix_in.port, MixingVolume1.port_a) annotation (points=[-48,68; -46,
        68; -46,40; -30.2,40], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(MixingVolume1.port_b, MixingVolume2.port_a) annotation (points=[-10,40; 
        9.8,40],      style(
      color=0,
      rgbcolor={0,0,0},
      smooth=0));
end TestMixingVolumesPressureStates;
