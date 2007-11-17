model TestJunctionIdeal2 
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
    annotation (extent=[-40,-30; -20,-10]);
  Modelica.Blocks.Sources.Ramp ramp(
    duration=1,
    height=-6.5e5,
    offset=7e5) annotation (extent=[-90,-24; -70,-4]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe(redeclare package 
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1)                      annotation (extent=[-8,-30; 12,-10]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe1(redeclare package 
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1)                      annotation (extent=[48,-30; 68,-10]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe2(redeclare package 
      Medium = 
        Modelica.Media.Air.DryAirNasa,
    length=1,
    diameter=0.1) 
    annotation (extent=[20,0; 40,20],  rotation=90);
equation 
  connect(ramp.y, source1.p_in) annotation (points=[-69,-14; -42,-14], style(
      color=74,
      rgbcolor={0,0,127},
      gradient=2,
      fillColor=69,
      rgbfillColor={0,128,255}));
  connect(source1.port, pipe.port_a) annotation (points=[-20,-20; -8,-20],
      style(color=69, rgbcolor={0,127,255}));
  connect(pipe1.port_b, source2.port) annotation (points=[68,-20; 80,-20],
      style(color=69, rgbcolor={0,127,255}));
  connect(pipe2.port_b, source3.port) annotation (points=[30,20; 30,40; 30,60; 
        30,60], style(color=69, rgbcolor={0,127,255}));
  connect(pipe.port_b, pipe1.port_a) annotation (points=[12,-20; 48,-20], style(
        color=69, rgbcolor={0,127,255}));
  connect(pipe.port_b, pipe2.port_a) annotation (points=[12,-20; 30,-20; 30,0], 
      style(color=69, rgbcolor={0,127,255}));
end TestJunctionIdeal2;
