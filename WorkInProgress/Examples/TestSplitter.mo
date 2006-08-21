model TestSplitter 
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  Modelica_Fluid.Junctions.Splitter splitter(
                               redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (extent=[20,-30; 40,-10]);
  Modelica_Fluid.Pipes.LumpedPipe pipe1(
    length=1,
    diameter=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[-10,-30; 10,-10], rotation=0);
  annotation (Diagram);
  Modelica_Fluid.Pipes.LumpedPipe pipe2(
    length=1,
    diameter=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[52,-30; 72,-10], rotation=180);
  Modelica_Fluid.Pipes.LumpedPipe pipe3(
    length=1,
    diameter=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[20,12; 40,32], rotation=270);
  Modelica_Fluid.Sources.FixedAmbient_pTX source2(
    T=278.15,
    p=5e5,
    redeclare package Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[80,-30; 100,-10], rotation=180);
  Modelica_Fluid.Sources.FixedAmbient_pTX source3(
    T=283.15,
    p=2e5,
    redeclare package Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[20,60; 40,80], rotation=270);
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[-100,80; -80,100]);
  Modelica_Fluid.Sources.PrescribedAmbient_pTX source1(           p=5e5,
      redeclare package Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[-40,-30; -20,-10]);
  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    height=-6.5e5,
    offset=7e5) annotation (extent=[-90,-24; -70,-4]);
equation 
  connect(source3.port, pipe3.port_a) annotation (points=[30,60; 30,53; 30,53; 
        30,46; 30,32; 30,32], style(
      color=69,
      rgbcolor={0,127,255},
      gradient=2,
      fillColor=69,
      rgbfillColor={0,128,255}));
  connect(source2.port, pipe2.port_a) annotation (points=[80,-20; 72,-20],
      style(
      color=69,
      rgbcolor={0,127,255},
      gradient=2,
      fillColor=69,
      rgbfillColor={0,128,255}));
  connect(pipe2.port_b, splitter.port_2) annotation (points=[52,-20; 41,-20],
      style(
      color=69,
      rgbcolor={0,127,255},
      gradient=2,
      fillColor=69,
      rgbfillColor={0,128,255}));
  connect(pipe3.port_b, splitter.port_3) annotation (points=[30,12; 30,1.5; 30,
        1.5; 30,-9], style(
      color=69,
      rgbcolor={0,127,255},
      gradient=2,
      fillColor=69,
      rgbfillColor={0,128,255}));
  connect(pipe1.port_b, splitter.port_1) annotation (points=[10,-20; 19,-20],
      style(
      color=69,
      rgbcolor={0,127,255},
      gradient=2,
      fillColor=69,
      rgbfillColor={0,128,255}));
  connect(source1.port, pipe1.port_a) annotation (points=[-20,-20; -10,-20],
      style(
      color=69,
      rgbcolor={0,127,255},
      gradient=2,
      fillColor=69,
      rgbfillColor={0,128,255}));
  connect(ramp.y, source1.p_in) annotation (points=[-69,-14; -42,-14], style(
      color=74,
      rgbcolor={0,0,127},
      gradient=2,
      fillColor=69,
      rgbfillColor={0,128,255}));
end TestSplitter;
