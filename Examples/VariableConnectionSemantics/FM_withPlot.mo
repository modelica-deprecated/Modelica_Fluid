model FM_withPlot 
  "Boundary condition (pressure) - flow model - boundary condition (pressure), with plot of pressure drop correlation" 
  extends Icons.Example;
  // Visualization
  parameter SI.Pressure plot_dp_max=5000;
  parameter SI.MassFlowRate plot_m_flow_max=1.6;
  parameter Integer plot_n_points=100;
protected 
  parameter SI.MassFlowRate plot_m_flow[plot_n_points]=linspace(
        -plot_m_flow_max,
        plot_m_flow_max,
        plot_n_points);
  SI.Pressure plot_dp[plot_n_points]=pipeFriction.WallFriction.pressureLoss_m_flow(
        plot_m_flow,
        pipeFriction.d_designDirection,
        pipeFriction.d_nonDesignDirection,
        pipeFriction.eta_designDirection,
        pipeFriction.eta_nonDesignDirection,
        pipeFriction.length,
        pipeFriction.diameter,
        roughness=pipeFriction.roughness,
        m_flow_small=0.01);
public 
  UserInteraction.Outputs.SpatialPlot2 spatialPlot(
    minX1=-plot_m_flow_max,
    maxX1=plot_m_flow_max,
    minY1=-plot_dp_max,
    maxY1=plot_dp_max,
    x1=plot_m_flow,
    y1=plot_dp,
    x2={pipeFriction.m_flow,pipeFriction.m_flow},
    y2={-plot_dp_max,plot_dp_max}) 
    annotation (extent=[-30,20; 50,100]);
  
  // Components
  FluidSandbox.Sources.PrescribedBoundary_pTX_A prescribedBoundary_pTX_A(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[60,-10; 80,10],rotation=180);
  FluidSandbox.Sources.PrescribedBoundary_pTX_A prescribedMassFlowRate_TX_A(
    redeclare package Medium = Medium,
    usePressureInput=true,
    redeclare package FluidInterface = FluidInterface,
    T=380) annotation (extent=[-60,-10; -40,10]);
  
  FluidSandbox.PressureLosses.WallFriction pipeFriction(
    redeclare package Medium = Medium,
    length=2,
    diameter=0.1,
    redeclare package FluidInterface = FluidInterface,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent)
                                              annotation (extent=[0,-10; 20,
        10]);
  
  Modelica.Blocks.Sources.Sine sine(
    phase=0,
    freqHz=2,
    amplitude=0.05e5,
    offset=1e5) 
              annotation (extent=[-100,-4; -80,16]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics1(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa, redeclare package 
      FluidInterface = FluidInterface) if FluidInterface.usesNewConnectionSemantics annotation (extent=[-15,-3;
        -25,-1]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics2(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa, redeclare package 
      FluidInterface = FluidInterface) if        FluidInterface.usesNewConnectionSemantics annotation (extent=[35,-3;
        45,-1]);
  annotation (Diagram);
  
equation 
  connect(sine.y, prescribedMassFlowRate_TX_A.p_in) 
    annotation (points=[-79,6; -62,6], style(color=74, rgbcolor={0,0,127}));
  connect(prescribedMassFlowRate_TX_A.port, semantics1.port_b) 
    annotation (points=[-40,0; -32,0; -32,-2; -25,-2],
                                       style(color=69, rgbcolor={0,127,255}));
  connect(semantics1.port_a, pipeFriction.port_a) 
    annotation (points=[-15,-2; -8,-2; -8,0; 0,0],
                                     style(color=69, rgbcolor={0,127,255}));
  connect(semantics2.port_a, pipeFriction.port_b) 
    annotation (points=[35,-2; 28,-2; 28,0; 20,0],
                                     style(color=69, rgbcolor={0,127,255}));
  connect(semantics2.port_b, prescribedBoundary_pTX_A.port) annotation (points=[45,-2;
        52,-2; 52,1.22461e-015; 60,1.22461e-015],
                                  style(color=69,
        rgbcolor={0,127,255}));
  
  // Use plain connections is no new semantics are required for this approach
  if not FluidInterface.usesNewConnectionSemantics then
    connect(prescribedMassFlowRate_TX_A.port, pipeFriction.port_a) 
    annotation (points=[-40,0; -32,0; -32,2; -8,2; -8,0; 0,0],
                                     style(color=69, rgbcolor={0,127,255}));
    connect(pipeFriction.port_b, prescribedBoundary_pTX_A.port) 
                                                              annotation (
      points=[20,0; 28,0; 28,2; 52,2; 52,0; 58,0; 58,1.22461e-015; 60,
          1.22461e-015],                          style(color=69,
        rgbcolor={0,127,255}));
  end if;
end FM_withPlot;
