model TestPortVolumes2 
  "Test case where in one of the mixing volumes a pressure state appears" 
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  annotation (Diagram, experiment(StopTime=10));
  Modelica_Fluid.Volumes.PortVolume MixingVolume1(
    V=1e-3,
    redeclare package Medium = Medium,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    p_start=ambient.default_p_ambient,
    use_T_start=true,
    T_start=ambient.default_T_ambient) 
                 annotation (extent=[-30,30; -10,50]);
  
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource2(
    m_flow_out=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-100,30; -80,50]);
  Modelica_Fluid.Volumes.PortVolume MixingVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    initType=Modelica_Fluid.Types.Init.NoInit) 
                 annotation (extent=[10,30; 30,50]);
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
    length=1,
    p_a_start=ambient.default_p_ambient,
    p_b_start=ambient.default_p_ambient,
    T_start=ambient.default_T_ambient) 
                                annotation (extent=[50,30; 70,50]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity simpleGenericOrifice1(
    redeclare package Medium = Medium,
    diameter=0.2,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Laminar,
    length=1,
    T_start=ambient.default_T_ambient,
    p_a_start=ambient.default_p_ambient,
    p_b_start=ambient.default_p_ambient) 
                                annotation (extent=[-70,30; -50,50]);
equation 
  connect(simpleGenericOrifice2.port_b, Sink2.port) annotation (points=[70,40;
        80,40], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(FlowSource2.port,simpleGenericOrifice1. port_a) annotation (points=[
        -80,40; -70,40], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(simpleGenericOrifice1.port_b, MixingVolume1.port) annotation (points=
        [-50,40; -20,40], style(color=69, rgbcolor={0,127,255}));
  connect(MixingVolume1.port, MixingVolume2.port)
    annotation (points=[-20,40; 20,40], style(color=69, rgbcolor={0,127,255}));
  connect(MixingVolume2.port, simpleGenericOrifice2.port_a)
    annotation (points=[20,40; 50,40], style(color=69, rgbcolor={0,127,255}));
end TestPortVolumes2;
