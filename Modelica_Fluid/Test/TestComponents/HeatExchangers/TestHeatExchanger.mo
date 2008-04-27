within Modelica_Fluid.Test.TestComponents.HeatExchangers;
model TestHeatExchanger 
  
extends Modelica.Icons.Example;
  
//replaceable package Medium = Modelica.Media.Water.StandardWater;
package Medium = Modelica.Media.Incompressible.Examples.Essotherm650;
  Modelica_Fluid.HeatExchangers.BasicHX HEX(
    c_wall=500,
    use_T_start=true,
    T_start_2=300,
    n=20,
    length=2,
    mflow_start_1=0.2,
    static=false,
    mflow_start_2=0.2,
    Twall_start=300,
    redeclare model HeatTransfer_1 = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha (
         alpha0=1000),
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    use_eta_nominal=true,
    k_wall=100,
    T_start_1=304,
    dT=10,
    initType=Modelica_Fluid.Types.Init.SteadyStateHydraulic,
    s_wall=0.005,
    Ac_1=4.5e-4,
    Ac_2=4.5e-4,
    P_1=0.075,
    P_2=0.075,
    Ah_1=0.075*2*20,
    Ah_2=0.075*2*20,
    d_wall=900,
    redeclare package Medium_1 = 
        Medium,
    redeclare package Medium_2 = 
        Medium)                annotation (extent=[-26,-14; 34,46]);
  
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient2(
    p=1e5,
    T=280,
    redeclare package Medium = Medium)                              annotation (extent=[82,-28;
        62,-8]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient1(
    p=1e5,
    T=300,
    redeclare package Medium = Medium)                              annotation (extent=[82,24;
        62,44]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate2(
    m_flow=0.2,
    T=360,
    redeclare package Medium = Medium,
    useFlowRateInput=true,
    useTemperatureInput=false,
    useCompositionInput=false) 
                annotation (extent=[-66,24; -46,44]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate1(
    T=300,
    m_flow=0.5,
    redeclare package Medium = Medium) 
                 annotation (extent=[-66,-10; -46,10]);
  annotation (Diagram, experiment(StopTime=100, Tolerance=1e-005));
  Modelica.Blocks.Sources.Ramp Ramp1(
    startTime=50,
    duration=5,
    height=-1,
    offset=0.5)   annotation (extent=[-100,24; -80,44]);
  inner Modelica_Fluid.Ambient ambient 
                                   annotation (extent=[60,70; 80,90]);
equation 
  connect(massFlowRate2.port, HEX.port_a2)            annotation (points=[-46,34; 
        -40,34; -40,29.8; -29,29.8],     style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=70,
      rgbfillColor={0,63,125},
      fillPattern=1));
  connect(massFlowRate1.port, HEX.port_a1)            annotation (points=[-46,0;
        -40,0; -40,15.4; -29,15.4], style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=70,
      rgbfillColor={0,63,125},
      fillPattern=1));
  connect(HEX.port_b1, ambient1.port)            annotation (points=[37,15.4;
        48.5,15.4; 48.5,34; 62,34], style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=70,
      rgbfillColor={0,63,125},
      fillPattern=1));
  connect(HEX.port_b2, ambient2.port)            annotation (points=[37,2.2; 
        49.5,2.2; 49.5,-18; 62,-18], style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=70,
      rgbfillColor={0,63,125},
      fillPattern=1));
  connect(Ramp1.y, massFlowRate2.m_flow_in) annotation (points=[-79,34; -74,34;
        -74,40; -65.3,40], style(color=74, rgbcolor={0,0,127}));
end TestHeatExchanger;
