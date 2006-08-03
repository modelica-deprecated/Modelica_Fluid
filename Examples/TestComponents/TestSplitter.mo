model TestSplitter 
  extends Modelica.Icons.Example;
  
  Modelica_Fluid.Junctions.Splitter splitter(redeclare package Medium = 
        Modelica.Media.Air.DryAirNasa)      annotation (extent=[20,-30; 40,-10]);
  annotation (Diagram);
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
  Pipes.LumpedPipe pipe(redeclare package Medium = 
        Modelica.Media.Air.DryAirNasa) annotation (extent=[-12,-30; 8,-10]);
  Pipes.LumpedPipe pipe1(redeclare package Medium = 
        Modelica.Media.Air.DryAirNasa) annotation (extent=[50,-30; 70,-10]);
  Pipes.LumpedPipe pipe2(redeclare package Medium = 
        Modelica.Media.Air.DryAirNasa) 
    annotation (extent=[20,14; 40,34], rotation=90);
equation 
  connect(ramp.y, source1.p_in) annotation (points=[-69,-14; -42,-14], style(
      color=74,
      rgbcolor={0,0,127},
      gradient=2,
      fillColor=69,
      rgbfillColor={0,128,255}));
  connect(source1.port, pipe.port_a) annotation (points=[-20,-20; -12,-20],
      style(color=69, rgbcolor={0,127,255}));
  connect(pipe.port_b, splitter.port_1) 
    annotation (points=[8,-20; 19,-20], style(color=69, rgbcolor={0,127,255}));
  connect(pipe1.port_b, source2.port) annotation (points=[70,-20; 80,-20],
      style(color=69, rgbcolor={0,127,255}));
  connect(splitter.port_2, pipe1.port_a) annotation (points=[41,-20; 50,-20],
      style(color=69, rgbcolor={0,127,255}));
  connect(pipe2.port_b, source3.port) annotation (points=[30,34; 30,47; 30,60;
        30,60], style(color=69, rgbcolor={0,127,255}));
  connect(pipe2.port_a, splitter.port_3) 
    annotation (points=[30,14; 30,-9], style(color=69, rgbcolor={0,127,255}));
end TestSplitter;