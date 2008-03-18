within FluidSandbox.Examples.StandardConnectionSemantics;
model WaterSteamCycle2 
  "Example adapted from PowerFluid: Simple steam power plant, including superheating and turbine control valves (pressure loss model in superheater is different from PowerFluid)." 
  extends Icons.Example;
  
  import Modelica.SIunits.Conversions.*;
  
  HeatTransfer.EvaporatingVessel evaporator(
    cp_D=500,
    m_D=100e3,
    V_t=50,
    V_l_start=33,
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    provide_p=true,
    provide_T=true,
    p_start=100000)     annotation (extent=[-30,0; -10,20]);
  annotation (
    Diagram,
    Coordsys(
      extent=[-100,-100; 230,100],
      grid=[1, 1],
      component=[20, 20],
      scale=0.1),
    experiment(StopTime=7200));
  Modelica.Thermal.HeatTransfer.PrescribedHeatFlow furnace1 
    annotation (extent=[-30,-35; -10,-15],   rotation=90);
  Modelica.Blocks.Continuous.PI controller(       k=10, T=60) 
    annotation (extent=[-39,43; -53,57]);
  Modelica.Blocks.Math.Feedback feedback 
    annotation (extent=[-14,40; -34,60]);
  Modelica.Blocks.Sources.Constant levelSetPoint(k=33) 
    annotation (extent=[-40,68; -26,82]);
  Modelica.Blocks.Interfaces.RealOutput T_S(redeclare type SignalType = 
        Real (unit="degC")) 
    annotation (extent=[79.5,46; 87.5,54]);
  Modelica.Blocks.Interfaces.RealOutput p_S(redeclare type SignalType = 
        Real (unit="bar")) 
    annotation (extent=[72.5,71; 80.5,79]);
  Modelica.Blocks.Interfaces.RealOutput qm_S(redeclare type SignalType = 
        Modelica.SIunits.MassFlowRate) 
    annotation (extent=[17.5,26; 25.5,34], rotation=0);
  Modelica.Blocks.Interfaces.RealOutput V_l(redeclare type SignalType = 
        Modelica.SIunits.Volume) 
    annotation (extent=[-10,26; -2,34]);
public 
  Modelica.Blocks.Math.Gain MW2W1(k=0.8e6) 
    annotation (extent=[-25,-52.5; -15,-41.5], rotation=90);
  Modelica.Blocks.Nonlinear.Limiter limiter(uMin=0, uMax=500) 
    annotation (extent=[-73,43; -59,57], rotation=180);
  
  Modelica.Blocks.Sources.TimeTable fuelTable(table=[0,0; 5400,1000; 7210,
        1000]) 
              annotation (extent=[-79.5,-70; -59.5,-50]);
  Sources.ControlledPump pump(
    provide_p_a=true,
    provide_p_b=false,
    provide_T_a=false,
    provide_T_b=true,
    provide_m_flow_ab=false,
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[-79.5,0; -59.5,20]);
  Modelica.Blocks.Interfaces.RealOutput T_W(redeclare type SignalType = 
        Real (unit="degC")) 
    annotation (extent=[-42.5,-30; -32.5,-20]);
  Turbomachinery.TurbineStageAA turbineStage1(medium_designDirection(h(start=2.6e6), p(start=10000)),
    K_t=0.001,
    G=1e3,
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    provide_p_a=false,
    provide_p_b=false,
    provide_T_a=false,
    provide_T_b=false,
    provide_m_flow_ab=false) 
    annotation (extent=[140,10; 160,30]);
  HeatTransfer.EvaporatingVessel condenser(
    V_t=100,
    m_D=10e3,
    cp_D=500,
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    provide_p=false,
    provide_T=true, 
    p_start=10000) 
              annotation (extent=[175,-70; 195,-90]);
  Modelica.Thermal.HeatTransfer.PrescribedHeatFlow cooling 
    annotation (extent=[175,-40; 195,-60],   rotation=90);
  Modelica.Blocks.Math.Feedback feedback1 
    annotation (extent=[125,-40; 145,-20],rotation=0);
  Modelica.Blocks.Sources.Constant pressureSetPoint(k=0.05) 
    annotation (extent=[103,-37; 117,-23]);
  Modelica.Blocks.Continuous.PI controller1(T=60, k=1e9) 
    annotation (extent=[153,-37; 167,-23]);
  Turbomachinery.TurbineStage turbineStage2(medium_designDirection(p(start=10000), h(start=2.6e6)),
    K_t=0.001,
    G=1e3,
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    provide_p_a=false,
    provide_p_b=false,
    provide_T_a=false,
    provide_T_b=false,
    provide_m_flow_ab=false) 
    annotation (extent=[169,10; 189,30]);
  HeatTransfer.HeatedPipe superHeater(V_lumped=25, T_start=from_degC(100),
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    provide_m_flow_a=true,
    provide_p_a=false,
    provide_T_a=false,
    provide_p_b=true,
    provide_T_b=true,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent,
    length=5,
    diameter=0.15,
    initType=Modelica_Fluid.Types.Init.InitialValues) 
    annotation (extent=[35.5,0; 55.5,20]);
  
  Modelica.Thermal.HeatTransfer.PrescribedHeatFlow furnace2 
    annotation (extent=[35.5,-35; 55.5,-15], rotation=90);
public 
  Modelica.Blocks.Math.Gain MW2W2(k=0.2e6) 
    annotation (extent=[40.5,-52.5; 50.5,-41.5],
                                             rotation=90);
  Valves.ValveLinear bypass(               Kv=1e-4,
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    provide_p_a=false,
    provide_p_b=false,
    provide_T_a=false,
    provide_T_b=false,
    provide_m_flow_ab=false) 
    annotation (extent=[150,40; 170,60]);
  Modelica.Blocks.Sources.TimeTable bypassTable(table=[0,1; 3600,0; 7210,0]) 
              annotation (extent=[130,70; 150,90]);
  Modelica.Blocks.Interfaces.RealOutput T_E(redeclare type SignalType = 
        Real (unit="degC")) 
    annotation (extent=[32.5,46; 40.5,54]);
  Modelica.Blocks.Interfaces.RealOutput p_E(redeclare type SignalType = 
        Real (unit="bar")) 
    annotation (extent=[32.5,76; 40.5,84]);
  Modelica.Blocks.Interfaces.RealOutput T_C(redeclare type SignalType = 
        Real (unit="degC")) 
    annotation (extent=[223,-54; 231,-46]);
  Modelica.Thermal.HeatTransfer.HeatCapacitor furnace2_mass(C=1e7) 
    annotation (extent=[10,-15; 30,-35]);
  Valves.ValveLinear control(               Kv=1e-4,
    redeclare package FluidInterface = FluidInterface,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    provide_p_a=false,
    provide_p_b=false,
    provide_T_a=false,
    provide_T_b=false,
    provide_m_flow_ab=false) 
    annotation (extent=[110,10; 130,30]);
  Modelica.Blocks.Sources.TimeTable controlTable(table=[0,0; 1800,0; 3600,1;
        7210,1], offset=0) 
              annotation (extent=[90,70; 110,90]);
  Modelica.Mechanics.Rotational.ConstantSpeed load(w_fixed=-50) 
    annotation (extent=[225.5,45.5; 210.5,60.5]);
  Junctions.IdealJunctionAAB idealJunction(redeclare package FluidInterface = 
        FluidInterface, redeclare package Medium = 
        Modelica.Media.Water.StandardWater) 
    annotation (extent=[80,10; 100,30], rotation=270);
  Junctions.IdealJunctionAAB idealJunction1(redeclare package FluidInterface = 
        FluidInterface, redeclare package Medium = 
        Modelica.Media.Water.StandardWater) 
    annotation (extent=[216,10; 196,30], rotation=270);
  Modelica.Blocks.Math.UnitConversions.To_degC partialConversionBlock1 
    annotation (extent=[28,42; 12,58], rotation=180);
  Modelica.Blocks.Math.UnitConversions.To_bar partialConversionBlock2 
    annotation (extent=[28,72; 12,88], rotation=180);
  Modelica.Blocks.Math.UnitConversions.To_degC partialConversionBlock3 
    annotation (extent=[-44,-18; -60,-2], rotation=270);
public 
  Modelica.Blocks.Math.Gain Pa2bar(k=1e-5) 
    annotation (extent=[93,-65.5; 103,-54.5],  rotation=90);
  Modelica.Blocks.Math.UnitConversions.To_degC partialConversionBlock4 
    annotation (extent=[226,-78; 210,-62], rotation=90);
  Modelica.Blocks.Math.UnitConversions.To_degC partialConversionBlock5 
    annotation (extent=[83,24; 67,40], rotation=90);
  Modelica.Blocks.Math.UnitConversions.To_bar partialConversionBlock6 
    annotation (extent=[68,52; 52,68], rotation=90);
equation 
  connect(controller.u,feedback.y) 
    annotation (points=[-37.6,50; -33,50], style(rgbcolor={0,0,127}));
  connect(evaporator.V, V_l) 
    annotation (points=[-16,21; -16,30; -6,30],    style(rgbcolor={0,0,127}));
  connect(MW2W1.y, furnace1.Q_flow)    annotation (points=[-20,-40.95; -20,
        -35],              style(rgbcolor={0,0,127}));
  connect(controller.y, limiter.u) annotation (points=[-53.7,50; -57.6,50],
      style(color=74, rgbcolor={0,0,127}));
  connect(fuelTable.y, MW2W1.u) 
                             annotation (points=[-58.5,-60; -20,-60; -20,
        -53.6],                                                     style(
        color=74, rgbcolor={0,0,127}));
  connect(levelSetPoint.y, feedback.u1) annotation (points=[-25.3,75; -10,75; 
        -10,50; -16,50],     style(color=74, rgbcolor={0,0,127}));
  connect(controller1.y, cooling.Q_flow)  annotation (points=[167.7,-30;
        185,-30; 185,-40], style(
      color=74,
      rgbcolor={0,0,127},
      gradient=3,
      fillColor=1,
      rgbfillColor={255,0,0}));
  connect(superHeater.heatPort, furnace2.port) annotation (points=[45.5,0;
        45.5,-15],       style(
      color=42,
      rgbcolor={191,0,0},
      gradient=3,
      fillColor=1,
      rgbfillColor={255,0,0}));
  connect(MW2W2.y, furnace2.Q_flow) annotation (points=[45.5,-40.95; 45.5,
        -39; 44.5,-39; 44.5,-38; 45.5,-38; 45.5,-35],
                                         style(
      color=74,
      rgbcolor={0,0,127},
      gradient=3,
      fillColor=1,
      rgbfillColor={255,0,0}));
  connect(MW2W1.u, MW2W2.u) 
                           annotation (points=[-20,-53.6; -20,-60; 45.5,-60;
        45.5,-53.6],
                   style(
      color=74,
      rgbcolor={0,0,127},
      gradient=3,
      fillColor=1,
      rgbfillColor={255,0,0}));
  connect(bypassTable.y, bypass.opening) annotation (points=[151,80; 160,80;
        160,59], style(color=74, rgbcolor={0,0,127}));
  connect(controlTable.y, control.opening) annotation (points=[111,80; 120,
        80; 120,29], style(color=74, rgbcolor={0,0,127}));
  connect(turbineStage2.flange, load.flange) 
    annotation (points=[189,23; 198.25,23; 198.25,53; 210.5,53],
                                         style(color=0, rgbcolor={0,0,0}));
  connect(turbineStage1.flange, turbineStage2.flange) 
    annotation (points=[160,23; 189,23], style(color=0, rgbcolor={0,0,0}));
  connect(pressureSetPoint.y, feedback1.u1) annotation (points=[117.7,-30; 127,
        -30], style(color=74, rgbcolor={0,0,127}));
  connect(feedback1.y, controller1.u) annotation (points=[144,-30; 151.6,
        -30], style(color=74, rgbcolor={0,0,127}));
  connect(furnace2_mass.port, furnace2.port) annotation (points=[20,-15;
        45.5,-15], style(color=42, rgbcolor={191,0,0}));
  connect(feedback.u2, evaporator.V) annotation (points=[-24,42; -24,30;
        -16,30; -16,21], style(color=74, rgbcolor={0,0,127}));
  connect(furnace1.port, evaporator.heatPort) annotation (points=[-20,-15; -20,
        0], style(
      color=42,
      rgbcolor={191,0,0},
      smooth=0));
  connect(condenser.heatPort, cooling.port) annotation (points=[185,-70; 185,
        -60], style(
      color=42,
      rgbcolor={191,0,0},
      smooth=0));
  connect(control.port_b, turbineStage1.port_a) annotation (points=[130,20; 140,
        20], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(turbineStage1.port_b, turbineStage2.port_a) annotation (points=[160,
        20; 169,20], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(idealJunction.port_2, superHeater.port_b) annotation (points=[90,10;
        55.5,10], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(idealJunction.port_3, control.port_a) annotation (points=[100,20; 105,
        20; 105,20; 110,20], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(idealJunction.port_1, bypass.port_a) annotation (points=[90,30; 90,50;
        150,50], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(turbineStage2.port_b, idealJunction1.port_3) annotation (points=[189,20; 
        192.5,20; 192.5,20; 196,20],     style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(idealJunction1.port_1, bypass.port_b) annotation (points=[206,30; 206,
        50; 170,50], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(partialConversionBlock1.y, T_E) annotation (points=[28.8,50; 36.5,50],
      style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(partialConversionBlock2.y, p_E) annotation (points=[28.8,80; 36.5,80],
      style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(T_W, partialConversionBlock3.y) annotation (points=[-37.5,-25; -52.25,
        -25; -52.25,-18.8; -52,-18.8], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(pump.T_b, partialConversionBlock3.u) annotation (points=[-58.5,15; 
        -52,15; -52,-0.4], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(Pa2bar.u, pump.p_a) annotation (points=[98,-66.6; 98,-84; -94,-84;
        -94,18; -80.5,18],
                         style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(feedback1.u2,Pa2bar. y) annotation (points=[135,-38; 135,-50.975; 98,
        -50.975; 98,-53.95],                                              style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(partialConversionBlock4.y, T_C) annotation (points=[218,-61.2; 218,
        -50; 227,-50], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(limiter.y, pump.massFlowRate) annotation (points=[-73.7,50; -80,50;
        -80,30; -69.5,30; -69.5,19], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(superHeater.m_flow_a, qm_S) annotation (points=[34.5,18; 13,18; 13,30;
        21.5,30], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(superHeater.T_b, partialConversionBlock5.u) annotation (points=[56.5,
        14; 75,14; 75,22.4], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(partialConversionBlock5.y, T_S) annotation (points=[75,40.8; 75,50;
        83.5,50], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(partialConversionBlock6.u, superHeater.p_b) annotation (points=[60,
        50.4; 60,16; 56.5,16], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(partialConversionBlock6.y, p_S) annotation (points=[60,68.8; 60,75;
        76.5,75], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(evaporator.p_sensor, partialConversionBlock2.u) annotation (points=[
        -31,18; 3,18; 3,80; 10.4,80], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(evaporator.T_sensor, partialConversionBlock1.u) annotation (points=[
        -31,15; 6,15; 6,50; 10.4,50], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(condenser.T_sensor, partialConversionBlock4.u) annotation (points=[
        174,-85; 218,-85; 218,-79.6], style(
      color=74, 
      rgbcolor={0,0,127}, 
      smooth=0));
  connect(pump.port_b, evaporator.port_a[1]) annotation (points=[-59.5,10; -30,
        10], style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(evaporator.port_b[1], superHeater.port_a) annotation (points=[-10,10; 
        35.3,10], style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(pump.port_a, condenser.port_a[1]) annotation (points=[-79.5,10; -90,
        10; -90,-80; 175,-80], style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(idealJunction1.port_2, condenser.port_b[1]) annotation (points=[206,
        10; 206,-80; 195,-80], style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
end WaterSteamCycle2;
