within FluidSandbox.Examples.StandardConnectionSemantics;
model FM 
  extends Icons.Example;
  FluidSandbox.Sources.PrescribedBoundary_pTX_A prescribedBoundary_pTX_A(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[40,-20; 60,0], rotation=180);
  FluidSandbox.Sources.PrescribedBoundary_pTX_A prescribedMassFlowRate_TX_A(
    redeclare package Medium = Medium,
    T=320,
    usePressureInput=true,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[-40,-20; -20,0]);
  FluidSandbox.PressureLosses.WallFriction pipeFriction(
    redeclare package Medium = Medium,
    length=2,
    diameter=0.1,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent,
    redeclare package FluidInterface = FluidInterface) 
                                              annotation (extent=[0,-20; 20,0]);
  
  Modelica.Blocks.Sources.Sine sine(
    phase=0,
    freqHz=2,
    amplitude=0.05e5,
    offset=1e5) 
              annotation (extent=[-80,-14; -60,6]);
equation 
  connect(pipeFriction.port_a, prescribedMassFlowRate_TX_A.port)   annotation (
      points=[0,-10; -20,-10], style(color=69, rgbcolor={0,127,255}));
  connect(pipeFriction.port_b, prescribedBoundary_pTX_A.port) 
    annotation (points=[20,-10; 30,-10; 30,-10; 40,-10],
                                      style(color=69, rgbcolor={0,127,255}));
  connect(sine.y, prescribedMassFlowRate_TX_A.p_in) 
    annotation (points=[-59,-4; -42,-4],
                                       style(color=74, rgbcolor={0,0,127}));
end FM;
