within FluidSandbox.Examples.StandardConnectionSemantics;
model IdealJunction 
  extends Icons.Example;
  
  FluidSandbox.Junctions.IdealJunction idealJunction(redeclare package 
      FluidInterface = FluidInterface, redeclare package Medium = Medium) 
  annotation (extent=[-20,-20; 0,0]);
  FluidSandbox.PressureLosses.WallFriction pipeFriction1(
    length=1,
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Medium,
    provide_p_a=false,
    provide_p_b=false,
    provide_T_a=false,
    provide_T_b=false,
    provide_m_flow_ab=false, 
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.Detailed, 
    diameter=0.2, 
    G=1e3) 
  annotation (extent=[-50,-20; -30,0]);
  
  FluidSandbox.PressureLosses.WallFriction pipeFriction2(
    length=1,
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Medium,
    provide_p_a=false,
    provide_p_b=false,
    provide_T_a=false,
    provide_T_b=false,
    provide_m_flow_ab=false, 
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.Detailed, 
    diameter=0.2, 
    G=1e3) 
  annotation (extent=[10,-20; 30,0]);
  
  FluidSandbox.PressureLosses.WallFriction pipeFriction3(
    length=1,
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Medium,
    provide_p_a=false,
    provide_p_b=false,
    provide_T_a=false,
    provide_T_b=false,
    provide_m_flow_ab=false, 
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.Detailed, 
    diameter=0.2, 
    G=1e3) 
  annotation (extent=[-20,10; 0,30], rotation=90);
  
  FluidSandbox.Sources.PrescribedBoundary_pTX_A prescribedBoundary_pTX_A(
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Medium, 
    usePressureInput=true, 
    useTemperatureInput=true, 
    p=100000, 
    T=300) 
  annotation (extent=[-80,-20; -60,0]);
  FluidSandbox.Sources.PrescribedBoundary_pTX_A prescribedBoundary_pTX_A1(
    p=1e5,
    usePressureInput=true,
    T=340,
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Medium) 
  annotation (extent=[-20,40; 0,60], rotation=270);
  FluidSandbox.Sources.PrescribedBoundary_pTX_A prescribedBoundary_pTX_A2(
    p=1.01e5,
    T=320,
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Medium) 
  annotation (extent=[60,-20; 40,0]);
  Modelica.Blocks.Sources.Sine sine(
    phase=0,
    freqHz=2,
    amplitude=0.05e5,
    offset=1e5) 
            annotation (extent=[-60,60; -40,80]);
  Modelica.Blocks.Sources.Sine sine1(
    freqHz=2,
    amplitude=0.05e5,
    offset=1e5, 
    phase=0.01396263401595464) 
            annotation (extent=[-100,40; -80,60], rotation=270);
  Modelica.Blocks.Sources.Sine sine2(
    amplitude=10, 
    freqHz=1, 
    phase=0.00174532925199433, 
    offset=300) 
            annotation (extent=[-100,-60; -80,-40], rotation=90);
equation 
  connect(pipeFriction1.port_b, idealJunction.port_1) 
                                                    annotation (points=[-30,
      -10; -20,-10], style(color=69, rgbcolor={0,127,255}));
  connect(idealJunction.port_3, pipeFriction3.port_a) 
  annotation (points=[-10,0; -10,10], style(color=69, rgbcolor={0,127,255}));
  connect(idealJunction.port_2, pipeFriction2.port_a) 
  annotation (points=[0,-10; 10,-10], style(color=69, rgbcolor={0,127,255}));
  connect(prescribedBoundary_pTX_A1.p_in, sine.y) 
                                                annotation (points=[-4,62; -4,
        70; -39,70],
                   style(color=74, rgbcolor={0,0,127}));
  connect(prescribedBoundary_pTX_A.port, pipeFriction1.port_a) 
                                                             annotation (
    points=[-60,-10; -50,-10], style(color=69, rgbcolor={0,127,255}));
  connect(pipeFriction3.port_b, prescribedBoundary_pTX_A1.port) 
                                                              annotation (
    points=[-10,30; -10,35; -10,40; -10,40], style(color=69, rgbcolor={0,127,
        255}));
  connect(pipeFriction2.port_b, prescribedBoundary_pTX_A2.port) 
                                                              annotation (
    points=[30,-10; 40,-10], style(color=69, rgbcolor={0,127,255}));
  annotation (Diagram, experiment(StopTime=10));
  connect(sine1.y, prescribedBoundary_pTX_A.p_in) annotation (points=[-90,39; 
        -90,-4; -82,-4], style(
      color=74, 
      rgbcolor={0,0,127}, 
      smooth=0));
  connect(sine2.y, prescribedBoundary_pTX_A.T_in) annotation (points=[-90,-39; 
        -90,-10; -82,-10], style(
      color=74, 
      rgbcolor={0,0,127}, 
      smooth=0));
end IdealJunction;
