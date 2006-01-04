package Tanks "Examples with Tanks" 
  extends Modelica.Icons.Library;
  
  model TwoTanksSimpleWater 
    import Modelica.SIunits.Conversions.*;
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      Coordsys(grid=[1, 1], component=[20, 20]),
      experiment(StopTime=50),
      Documentation(info="<html>
<p>
Two tanks with different initial levels are connected by a 
horizontal pipe. After about 30 s, the levels of the
two tanks are identical and the system is in steady state.
The two tanks have different initial temperatures. Since the
water is flowing from tank1 to tank2, the temperature of tank2
changes until water from tank1 is flowing into tank2.
</p>
</html>"),
      experimentSetupOutput);
    Components.Tank Tank1(
      area=1,
      T_start=from_degC(50),
      level_start=3,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initOption=Modelica_Fluid.Types.Init.InitialValues,
      V0=0.1,
      pipeArea=0.01) 
      annotation (extent=[-70,20; -50,40]);
    Components.Tank Tank2(
      area=1,
      T_start=from_degC(100),
      level_start=1,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initOption=Modelica_Fluid.Types.Init.InitialValues,
      V0=0.1,
      pipeArea=0.01) 
      annotation (extent=[10,20; 30,40]);
    Components.PressureDropPipe shortPipe1(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      frictionType=Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[-29,0; -9,20]);
    inner Components.FluidOptions fluidOptions 
      annotation (extent=[-100,-100; -80,-80]);
  equation 
    connect(Tank1.port, shortPipe1.port_a) 
      annotation (points=[-60,19; -60,10; -29,10],      style(color=69));
    connect(shortPipe1.port_b, Tank2.port) annotation (points=[-9,10; 20,10; 20,
          19], style(color=69, rgbcolor={0,127,255}));
  end TwoTanksSimpleWater;
  
  model ThreeTanksSimpleWater 
    import Modelica.SIunits.Conversions.*;
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      Coordsys(grid=[1, 1], component=[20, 20]),
      experiment(StopTime=50),
      Documentation(info="<html>
<p>
This example demonstrates the mixing of simple (constant
property) liquid water model
between three tanks with different temperature.
The water model is incompressible.
The same tank system with a compressible water model
is provided in ThreeTanksIF97.
</p>
</html>"),
      experimentSetupOutput);
    Components.Tank Tank1(
      area=1,
      T_start=from_degC(50),
      level_start=3,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initOption=Modelica_Fluid.Types.Init.InitialValues,
      V0=0.1,
      pipeArea=0.01) 
      annotation (extent=[-90,20; -70,40]);
    Components.Tank Tank2(
      area=1,
      T_start=from_degC(100),
      level_start=1,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initOption=Modelica_Fluid.Types.Init.InitialValues,
      H0=0.5,
      V0=0.1,
      pipeArea=0.01) 
      annotation (extent=[-10, 20; 10, 40]);
    Components.PressureDropPipe shortPipe1(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      frictionType=Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[-50, -30; -30, -10]);
    Components.Tank Tank3(
      area=1,
      T_start=from_degC(20),
      level_start=2,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initOption=Modelica_Fluid.Types.Init.InitialValues,
      V0=0.1,
      pipeArea=0.01) 
      annotation (extent=[70,20; 90,40]);
    Components.PressureDropPipe shortPipe3(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      frictionType=Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[30, -30; 50, -10]);
    Components.PressureDropPipe shortPipe2(
      m_flow_nominal=1000,
      dp_nominal=from_bar(0.1),
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      frictionType=Types.FrictionTypes.ConstantLaminar,
      medium_b(T(stateSelect=StateSelect.avoid)),
      medium_a(T(stateSelect=StateSelect.avoid))) 
      annotation (extent=[-10, -10; 10, 10], rotation=-90);
    inner Components.FluidOptions fluidOptions 
      annotation (extent=[-100,-100; -80,-80]);
  equation 
    connect(Tank1.port, shortPipe1.port_a) 
      annotation (points=[-80,19; -80,-20; -50,-20],    style(color=69));
    connect(shortPipe3.port_b, Tank3.port) 
      annotation (points=[50,-20; 80,-20; 80,19],    style(color=69));
    connect(Tank2.port, shortPipe2.port_a) 
      annotation (points=[0,19; 0,10; -6.12303e-016,10],    style(color=69));
    connect(shortPipe1.port_b, shortPipe3.port_a) 
      annotation (points=[-30,-20; 30,-20],   style(color=69));
    connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
          6.12303e-016,-10; 0,-10; 0,-20; 30,-20],     style(color=69));
  end ThreeTanksSimpleWater;
  
  model ThreeTanksIF97 
    import Modelica.SIunits.Conversions.*;
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      Coordsys(grid=[1, 1], component=[20, 20]),
      experiment(StopTime=50),
      Documentation(info="<html>
<p>
This example is the same as the \"ThreeTanksSimpleWater\"
model. The only difference is that the very detailed,
compressible medium model WaterIF97 is used.
</p>
</html>"),
      experimentSetupOutput);
    Components.Tank Tank1(
      area=1,
      T_start=from_degC(50),
      level_start=3,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      initOption=Modelica_Fluid.Types.Init.InitialValues,
      pipeArea=0.01) 
      annotation (extent=[-90, 20; -70, 40]);
    Components.Tank Tank2(
      area=1,
      level_start=1,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      initOption=Modelica_Fluid.Types.Init.InitialValues,
      T_start=from_degC(90),
      pipeArea=0.01) 
      annotation (extent=[-10,20; 10,40]);
    Components.PressureDropPipe shortPipe1(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      frictionType=Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
      annotation (extent=[-50, -30; -30, -10]);
    Components.Tank Tank3(
      area=1,
      T_start=from_degC(20),
      level_start=2,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      initOption=Modelica_Fluid.Types.Init.InitialValues,
      pipeArea=0.01) 
      annotation (extent=[70, 20; 90, 40]);
    Components.PressureDropPipe shortPipe3(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      frictionType=Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
      annotation (extent=[30, -30; 50, -10]);
    Components.PressureDropPipe shortPipe2(
      m_flow_nominal=1000,
      dp_nominal=from_bar(0.1),
      frictionType=Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
      annotation (extent=[-10,-10; 10,10],   rotation=-90);
    inner Components.FluidOptions fluidOptions 
      annotation (extent=[-100,-100; -80,-80]);
  equation 
    connect(Tank1.port, shortPipe1.port_a) 
      annotation (points=[-80,19; -80,-20; -50,-20],    style(color=69));
    connect(shortPipe3.port_b, Tank3.port) 
      annotation (points=[50,-20; 80,-20; 80,19],    style(color=69));
    connect(Tank2.port, shortPipe2.port_a) 
      annotation (points=[0,19; 0,10; -6.12303e-016,10],    style(color=69));
    connect(shortPipe1.port_b, shortPipe3.port_a) 
      annotation (points=[-30,-20; 30,-20],   style(color=69));
    connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
          6.12303e-016,-10; 0,-10; 0,-20; 30,-20],     style(color=69));
  end ThreeTanksIF97;
  
  model ThreeTanksIF97SteadyState 
    import Modelica.SIunits.Conversions.*;
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      Coordsys(grid=[1, 1], component=[20, 20]),
      experiment(StopTime=50),
      Documentation(info="<html>
<p>
This example is the same as ThreeTanks97. The only difference
is that the system starts in steady state, i.e., with constant
(mixing) temperature.
</p>
</html>"),
      experimentSetupOutput);
    Components.Tank Tank1(
      area=1,
      T_start=from_degC(50),
      level_start=3,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      pipeArea=0.01) 
      annotation (extent=[-90, 20; -70, 40]);
    Components.Tank Tank2(
      area=1,
      level_start=1,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      T_start=from_degC(90),
      pipeArea=0.01) 
      annotation (extent=[-10,20; 10,40]);
    Components.PressureDropPipe shortPipe1(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      frictionType=Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
      annotation (extent=[-50, -30; -30, -10]);
    Components.Tank Tank3(
      area=1,
      T_start=from_degC(20),
      level_start=2,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      pipeArea=0.01) 
      annotation (extent=[70, 20; 90, 40]);
    Components.PressureDropPipe shortPipe3(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      frictionType=Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
      annotation (extent=[30, -30; 50, -10]);
    Components.PressureDropPipe shortPipe2(
      m_flow_nominal=1000,
      dp_nominal=from_bar(0.1),
      frictionType=Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
      annotation (extent=[-10, -10; 10, 10], rotation=-90);
    inner Components.FluidOptions fluidOptions(default_initOption=
          Modelica_Fluid.Types.Init.SteadyState) 
      annotation (extent=[-100,-100; -80,-80]);
  equation 
    connect(Tank1.port, shortPipe1.port_a) 
      annotation (points=[-80,19; -80,-20; -50,-20],    style(color=69));
    connect(shortPipe3.port_b, Tank3.port) 
      annotation (points=[50,-20; 80,-20; 80,19],    style(color=69));
    connect(Tank2.port, shortPipe2.port_a) 
      annotation (points=[0,19; 0,10; -6.12303e-016,10],    style(color=69));
    connect(shortPipe1.port_b, shortPipe3.port_a) 
      annotation (points=[-30,-20; 30,-20],   style(color=69));
    connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
          6.12303e-016,-10; 0,-10; 0,-20; 30,-20],     style(color=69));
  end ThreeTanksIF97SteadyState;
end Tanks;
