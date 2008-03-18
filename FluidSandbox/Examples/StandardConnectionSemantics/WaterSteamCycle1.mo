within FluidSandbox.Examples.StandardConnectionSemantics;
model WaterSteamCycle1 
  "Example adapted from PowerFluid: Basic steam power plant" 
  extends Icons.Example;
  
  import Modelica.SIunits.Conversions.*;
  
  HeatTransfer.EvaporatingVessel evaporator(
    cp_D=500,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_D=100e3,
    V_t=50,
    V_l_start=33,
    redeclare package FluidInterface = FluidInterface, 
    provide_p=false, 
    provide_T=false, 
    initType=Modelica_Fluid.Types.Init.InitialValues, 
    p_start=100000)                                      annotation (extent=[-24,-20; -5,0]);
  annotation (
    Diagram,
    Coordsys(
      extent=[-100,-100; 200,100],
      grid=[1, 1],
      component=[20, 20],
      scale=0.1),
    experiment(StopTime=7200));
  Modelica.Thermal.HeatTransfer.PrescribedHeatFlow furnace 
    annotation (extent=[-24.5,-50; -4.5,-30],rotation=90);
  Modelica.Blocks.Continuous.PI controller(       k=10, T=60) 
    annotation (extent=[-51,33; -65,47]);
  Modelica.Blocks.Math.Feedback feedback 
    annotation (extent=[-26,30; -46,50]);
  Modelica.Blocks.Sources.Constant levelSetPoint(k=33) 
    annotation (extent=[-42,58; -28,72]);
  Modelica.Blocks.Interfaces.RealOutput T_S(redeclare type SignalType = 
        Real (unit="degC")) 
    annotation (extent=[20,66; 28,74]);
  Modelica.Blocks.Interfaces.RealOutput p_S(redeclare type SignalType = 
        Real (unit="bar")) 
    annotation (extent=[40,30; 48,38]);
  Modelica.Blocks.Interfaces.RealOutput qm_S(redeclare type SignalType = 
        Modelica.SIunits.MassFlowRate) 
    annotation (extent=[40,6; 48,14],      rotation=0);
  Modelica.Blocks.Interfaces.RealOutput V_l(redeclare type SignalType = 
        Modelica.SIunits.Volume) 
    annotation (extent=[-8,26; 0,34]);
public 
  Modelica.Blocks.Math.Gain MW2W(k=1e6) 
    annotation (extent=[-33,-65.5; -23,-54.5], rotation=0);
  Modelica.Blocks.Nonlinear.Limiter limiter(uMin=0, uMax=500) 
    annotation (extent=[-85,33; -71,47], rotation=180);
  
  Modelica.Blocks.Sources.TimeTable fuelTable(table=[0,0; 5400,1000; 7210,
        1000]) 
              annotation (extent=[-80,-70; -60,-50]);
  Sources.ControlledPump pump(              redeclare package Medium = 
        Modelica.Media.Water.StandardWater,
    provide_p_b=false,
    provide_T_b=true,
    provide_m_flow_ab=false,
    provide_p_a=true,
    redeclare package FluidInterface = FluidInterface,
    provide_T_a=false) 
    annotation (extent=[-75,-20; -55,0]);
  Modelica.Blocks.Interfaces.RealOutput T_W(redeclare type SignalType = 
        Real (unit="degC")) 
    annotation (extent=[-37,-35; -27,-25]);
  Turbomachinery.TurbineStage turbine(       redeclare package Medium = 
        Modelica.Media.Water.StandardWater, K_t=0.001,
    provide_p_b=false,
    provide_m_flow_ab=true,
    provide_p_a=true,
    provide_T_b=true,
    provide_T_a=true,
    redeclare package FluidInterface = FluidInterface) 
    annotation (extent=[85,26; 105,46]);
  HeatTransfer.EvaporatingVessel condenser(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    V_t=100,
    m_D=10e3,
    cp_D=500,
    redeclare package FluidInterface = FluidInterface, 
    provide_p=false, 
    provide_T=false, 
    p_start=10000) 
              annotation (extent=[141,-70; 161,-90]);
  Modelica.Thermal.HeatTransfer.PrescribedHeatFlow furnace1 
    annotation (extent=[141,-35; 161,-55],   rotation=90);
  Modelica.Blocks.Math.Feedback feedback1 
    annotation (extent=[88,-35; 108,-15], rotation=0);
  Modelica.Blocks.Sources.Constant pressureSetPoint(k=0.05) 
    annotation (extent=[60,-47; 74,-33]);
  Modelica.Blocks.Continuous.PI controller1(T=60, k=1e9) 
    annotation (extent=[119,-32; 133,-18]);
  Modelica.Mechanics.Rotational.ConstantSpeed load(w_fixed=-50) 
    annotation (extent=[163,32.5; 148,47.5]);
  Modelica.Blocks.Interfaces.RealOutput T_LP(
                                            redeclare type SignalType = 
        Real (unit="degC")) 
    annotation (extent=[165,-18; 173,-10]);
public 
  Modelica.Blocks.Math.Gain Pa2bar(k=1e-5) 
    annotation (extent=[93,-65.5; 103,-54.5],  rotation=90);
public 
  Modelica.Blocks.Math.Gain Pa2bar1(k=1e-5) 
    annotation (extent=[60,60; 50,70],         rotation=0);
  Modelica.Blocks.Math.UnitConversions.To_degC partialConversionBlock 
    annotation (extent=[130,0; 150,20]);
  Modelica.Blocks.Math.UnitConversions.To_degC partialConversionBlock1 
    annotation (extent=[60,80; 40,100]);
equation 
  connect(controller.u,feedback.y) 
    annotation (points=[-49.6,40; -45,40], style(rgbcolor={0,0,127}));
  connect(feedback.u2,      evaporator.V) 
    annotation (points=[-36,32; -36,20; -10.7,20; -10.7,1],
                                          style(rgbcolor={0,0,127}));
  connect(evaporator.V, V_l) 
    annotation (points=[-10.7,1; -10.7,30; -4,30], style(rgbcolor={0,0,127}));
  connect(MW2W.y,furnace.Q_flow)       annotation (points=[-22.5,-60; -14.5,-60;
        -14.5,-50],        style(rgbcolor={0,0,127}));
  connect(controller.y, limiter.u) annotation (points=[-65.7,40; -69.6,40],
      style(color=74, rgbcolor={0,0,127}));
  connect(levelSetPoint.y, feedback.u1) annotation (points=[-27.3,65; -20,65; 
        -20,40; -28,40],     style(color=74, rgbcolor={0,0,127}));
  connect(pressureSetPoint.y, feedback1.u1) annotation (points=[74.7,-40; 82,
        -40; 82,-25; 90,-25],
                 style(color=74, rgbcolor={0,0,127}));
  connect(feedback1.y, controller1.u) annotation (points=[107,-25; 117.6,-25],
              style(
      color=74,
      rgbcolor={0,0,127},
      gradient=3,
      fillColor=1,
      rgbfillColor={255,0,0}));
  connect(controller1.y, furnace1.Q_flow) annotation (points=[133.7,-25; 151,
        -25; 151,-35],     style(
      color=74,
      rgbcolor={0,0,127},
      gradient=3,
      fillColor=1,
      rgbfillColor={255,0,0}));
  connect(turbine.flange, load.flange) 
    annotation (points=[105,39; 126.5,39; 126.5,40; 148,40],
                                         style(color=0, rgbcolor={0,0,0}));
  connect(fuelTable.y, MW2W.u) 
                             annotation (points=[-59,-60; -34,-60],
                         style(color=74, rgbcolor={0,0,127}));
  connect(pump.T_b, T_W) annotation (points=[-54,-5; -45,-5; -45,-30; -32,-30],
      style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(limiter.y, pump.massFlowRate) annotation (points=[-85.7,40; -90,40;
        -90,20; -65,20; -65,-1], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(feedback1.u2, Pa2bar.y) annotation (points=[98,-33; 98,-53.95], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(Pa2bar.u, pump.p_a) annotation (points=[98,-66.6; 98,-84; -94,-84;
        -94,-2; -76,-2], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(turbine.m_flow_ab, qm_S) annotation (points=[89,47; 89,52; 60,52; 60,
        4; 35,4; 35,10; 44,10], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(turbine.p_a, Pa2bar1.u) annotation (points=[84,44; 75,44; 75,65; 61,
        65], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(Pa2bar1.y, p_S) annotation (points=[49.5,65; 34,65; 34,34; 44,34],
      style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(turbine.T_b, partialConversionBlock.u) annotation (points=[106,41;
        117,41; 117,10; 128,10], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(partialConversionBlock.y, T_LP) annotation (points=[151,10; 160,10;
        160,-14; 169,-14], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(partialConversionBlock1.y, T_S) annotation (points=[39,90; 10,90; 10,
        70; 24,70], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(turbine.T_a, partialConversionBlock1.u) annotation (points=[84,41; 80,
        41; 80,90; 62,90], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(evaporator.port_b[1], turbine.port_a) annotation (points=[-5,-10; 70,
        -10; 70,36; 85,36], style(
      color=69, 
      rgbcolor={0,127,255}, 
      thickness=2, 
      smooth=0));
  connect(pump.port_b, evaporator.port_a[1]) annotation (points=[-55,-10; -24,
        -10], style(
      color=69, 
      rgbcolor={0,127,255}, 
      thickness=2, 
      smooth=0));
  connect(condenser.port_a[1], pump.port_a) annotation (points=[141,-80; -90,
        -80; -90,-10; -75,-10], style(
      color=69, 
      rgbcolor={0,127,255}, 
      thickness=2, 
      smooth=0));
  connect(condenser.port_b[1], turbine.port_b) annotation (points=[161,-80; 190,
        -80; 190,36; 105,36], style(
      color=69, 
      rgbcolor={0,127,255}, 
      thickness=2, 
      smooth=0));
  connect(evaporator.heatPort, furnace.port) annotation (points=[-14.5,-20; 
        -14.5,-30], style(
      color=42, 
      rgbcolor={191,0,0}, 
      smooth=0));
  connect(furnace1.port, condenser.heatPort) annotation (points=[151,-55; 151,
        -70], style(
      color=42, 
      rgbcolor={191,0,0}, 
      smooth=0));
end WaterSteamCycle1;
