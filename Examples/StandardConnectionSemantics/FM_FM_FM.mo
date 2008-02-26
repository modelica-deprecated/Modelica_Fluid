within FluidSandbox.Examples.StandardConnectionSemantics;
model FM_FM_FM 
  extends Icons.Example;
  FluidSandbox.Sources.PrescribedBoundary_pTX_A prescribedBoundary_pTX_A(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[70,-20; 90,0], rotation=180);
  FluidSandbox.Sources.PrescribedBoundary_pTX_A prescribedMassFlowRate_TX_A(
    redeclare package Medium = Medium,
    T=320,
    usePressureInput=true,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[-90,-20; -70,0]);
  FluidSandbox.PressureLosses.WallFriction pipeFriction1(
    redeclare package Medium = Medium,
    length=2,
    diameter=0.1,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent,
    redeclare package FluidInterface = FluidInterface, 
    provide_p_a=false, 
    provide_p_b=false, 
    provide_T_a=false, 
    provide_T_b=false, 
    provide_m_flow_ab=false)                  annotation (extent=[-50,-20; -30,
        0]);
  
  Modelica.Blocks.Sources.Sine sine(
    phase=0,
    freqHz=2,
    amplitude=0.05e5,
    offset=1e5) 
              annotation (extent=[-72,20; -92,40], rotation=0);
  PressureLosses.WallFrictionAA pipeFriction2(
    redeclare package Medium = Medium,
    length=2,
    diameter=0.1,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent,
    redeclare package FluidInterface = FluidInterface, 
    provide_p_a=false, 
    provide_p_b=false, 
    provide_T_a=false, 
    provide_T_b=false, 
    provide_m_flow_ab=false)                  annotation (extent=[-10,-20; 10,0]);
  
  FluidSandbox.PressureLosses.WallFriction pipeFriction3(
    redeclare package Medium = Medium,
    length=2,
    diameter=0.1,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent,
    redeclare package FluidInterface = FluidInterface, 
    provide_p_a=false, 
    provide_p_b=false, 
    provide_T_a=false, 
    provide_T_b=false, 
    provide_m_flow_ab=false)                  annotation (extent=[30,-20; 50,0]);
  
equation 
  connect(sine.y, prescribedMassFlowRate_TX_A.p_in) 
    annotation (points=[-93,30; -96,30; -96,-4; -92,-4],
                                       style(color=74, rgbcolor={0,0,127}));
  annotation (Diagram);
  connect(prescribedMassFlowRate_TX_A.port, pipeFriction1.port_a) annotation (
      points=[-70,-10; -50,-10], style(color=69, rgbcolor={0,127,255}));
  connect(pipeFriction1.port_b, pipeFriction2.port_a) annotation (points=[-30,
        -10; -10,-10], style(color=69, rgbcolor={0,127,255}));
  connect(pipeFriction2.port_b, pipeFriction3.port_a) annotation (points=[10,
        -10; 30,-10], style(color=69, rgbcolor={0,127,255}));
  connect(pipeFriction3.port_b, prescribedBoundary_pTX_A.port) annotation (
      points=[50,-10; 60,-10; 60,-10; 70,-10], style(color=69, rgbcolor={0,127,
          255}));
end FM_FM_FM;
