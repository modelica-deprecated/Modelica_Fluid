within FluidSandbox.Examples.StandardConnectionSemantics;
model FM_CV_FM 
  extends Icons.Example;
  FluidSandbox.Sources.PrescribedBoundary_pTX_A source2(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[50,-10; 70,10],rotation=180);
  FluidSandbox.PressureLosses.WallFriction pipeFriction1(
    redeclare package Medium = Medium,
    length=2,
    diameter=0.1,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent,
    redeclare package FluidInterface = FluidInterface) 
                                              annotation (extent=[-40,-10; -20,
        10], rotation=0);
  
  FluidSandbox.Volumes.Volume mixingVolume(
    V=0.1,
    redeclare package Medium = Medium,
    n_a=1,
    n_b=1,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    p_start=1e5,
    use_T_start=true,
    T_start=340,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[-10,-10; 10,10]);
  FluidSandbox.PressureLosses.WallFriction pipeFriction2(
    redeclare package Medium = Medium,
    length=2,
    diameter=0.1,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent,
    redeclare package FluidInterface = FluidInterface) 
                                              annotation (extent=[20,10; 40,-10],
      rotation=180);
  
  FluidSandbox.Sources.PrescribedBoundary_pTX_A prescribedMassFlowRate_TX_A(
    redeclare package Medium = Medium,
    T=320,
    usePressureInput=true,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[-70,-10; -50,10]);
  Modelica.Blocks.Sources.Sine sine1(
    phase=0,
    freqHz=2,
    amplitude=0.05e5,
    offset=1e5) 
              annotation (extent=[-90,30; -70,50], rotation=270);
equation 
  connect(pipeFriction1.port_b, mixingVolume.port_a[1]) annotation (points=[-20,0;
        -10,0],            style(color=69, rgbcolor={0,127,255}));
  connect(pipeFriction2.port_b, mixingVolume.port_b[1]) annotation (points=[20,
        -1.22461e-015; 15,-1.22461e-015; 15,0; 10,0],
                                     style(color=69, rgbcolor={0,127,255}));
  connect(source2.port, pipeFriction2.port_a)                  annotation (
      points=[50,1.22461e-015; 47.5,1.22461e-015; 47.5,0; 45,0; 45,1.22461e-015;
        40,1.22461e-015],                                    style(color=69,
        rgbcolor={0,127,255}));
  connect(sine1.y, prescribedMassFlowRate_TX_A.p_in) annotation (points=[-80,29;
        -80,6; -72,6],       style(color=74, rgbcolor={0,0,127}));
  connect(prescribedMassFlowRate_TX_A.port, pipeFriction1.port_a) annotation (
     points=[-50,0; -40,0],   style(color=69, rgbcolor={0,127,255}));
end FM_CV_FM;
