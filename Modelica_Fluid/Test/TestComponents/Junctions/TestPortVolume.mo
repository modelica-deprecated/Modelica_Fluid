model TestPortVolume 
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  
  annotation (Diagram);
  Modelica_Fluid.Sources.FixedBoundary_pTX source2(
    T=278.15,
    p=5e5,
    redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
    annotation (extent=[80,-30; 100,-10], rotation=180);
  Modelica_Fluid.Sources.FixedBoundary_pTX source3(
    T=283.15,
    p=2e5,
    redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
    annotation (extent=[20,60; 40,80], rotation=270);
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[-100,80; -80,100]);
  Modelica_Fluid.Sources.PrescribedBoundary_pTX source1(          p=5e5,
      redeclare package Medium = Modelica.Media.Air.DryAirNasa,
    T=ambient.default_T_ambient, 
    usePressureInput=true) 
    annotation (extent=[-48,-30; -28,-10]);
  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    height=-6.5e5,
    offset=7e5) annotation (extent=[-90,-24; -70,-4]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe1(
    redeclare package Medium = Modelica.Media.Air.DryAirNasa, 
    length=1, 
    diameter=0.1)                      annotation (extent=[-12,-30; 8,-10]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe2(redeclare package 
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1)                      annotation (extent=[70,-30; 50,-10]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe3(redeclare package 
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1) 
    annotation (extent=[20,34; 40,14], rotation=90);
  Modelica_Fluid.Volumes.PortVolume portVolume(
    redeclare package Medium = Modelica.Media.Air.DryAirNasa, 
    V=20e-6, 
    initType=Modelica_Fluid.Types.Init.InitialValues, 
    p_start=1e5) annotation (extent=[20,-10; 40,-30]);
equation 
  connect(ramp.y, source1.p_in) annotation (points=[-69,-14; -50,-14], style(
      color=74,
      rgbcolor={0,0,127},
      gradient=2,
      fillColor=69,
      rgbfillColor={0,128,255}));
  connect(source1.port, pipe1.port_a) 
                                     annotation (points=[-28,-20; -12,-20],
      style(color=69, rgbcolor={0,127,255}));
  connect(pipe3.port_a, source3.port) annotation (points=[30,34; 30,47; 30,60; 
        30,60], style(color=69, rgbcolor={0,127,255}));
  connect(pipe2.port_a, source2.port) annotation (points=[70,-20; 80,-20], 
      style(color=69, rgbcolor={0,127,255}));
  connect(pipe1.port_b, portVolume.port)
    annotation (points=[8,-20; 30,-20], style(color=69, rgbcolor={0,127,255}));
  connect(pipe3.port_b, portVolume.port)
    annotation (points=[30,14; 30,-20], style(color=69, rgbcolor={0,127,255}));
  connect(pipe2.port_b, portVolume.port) annotation (points=[50,-20; 30,-20], 
      style(color=69, rgbcolor={0,127,255}));
end TestPortVolume;
