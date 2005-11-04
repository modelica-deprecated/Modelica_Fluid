package Tanks "Examples with Tanks" 
  extends Modelica.Icons.Library;
  
  model ThreeTanksOneLiquid 
    import Modelica.SIunits.Conversions.*;
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      Coordsys(grid=[1, 1], component=[20, 20]),
      experiment(StopTime=5),
      Documentation(info="<html>
<p>
This example demonstrates the mixing of a single substance flow
between three tanks with different temperature.
</p>
<p>
When comparing with Examples.Tanks.ThreeTanks, it is demonstrated
that it is easy in this case to switch between single and multiple
substance flow. The only differences between the examples ThreeTanks
and ThreeTanksOneLiquid are that the Medium is different and that
for the model ThreeTanksOneLiquid the default values are
used for the initial mass fractions.
</p>
</html>"));
    Components.Tank Tank1(
      area=1,
      T_start=from_degC(50),
      level_start=3,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initOption=Modelica_Fluid.Types.InitTypes.InitialValues,
      V0=0.1) 
      annotation (extent=[-90,21; -70,41]);
    Components.Tank Tank2(
      area=1,
      T_start=from_degC(100),
      level_start=1,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initOption=Modelica_Fluid.Types.InitTypes.InitialValues,
      H0=0.5,
      V0=0.1) 
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
      initOption=Modelica_Fluid.Types.InitTypes.InitialValues,
      V0=0.1) 
      annotation (extent=[71,21; 91,41]);
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
  equation 
    connect(Tank1.port, shortPipe1.port_a) 
      annotation (points=[-80,20; -80,-20; -51,-20],    style(color=69));
    connect(shortPipe3.port_b, Tank3.port) 
      annotation (points=[51,-20; 81,-20; 81,20],    style(color=69));
    connect(Tank2.port, shortPipe2.port_a) 
      annotation (points=[0,19; 0,11; -6.73533e-016,11],    style(color=69));
    connect(shortPipe1.port_b, shortPipe3.port_a) 
      annotation (points=[-29, -20; 29, -20], style(color=69));
    connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
          6.73533e-016,-11; 0,-11; 0,-20; 29,-20],     style(color=69));
  end ThreeTanksOneLiquid;
  
  model ThreeTanksWithPortVolume 
    import Modelica.SIunits.Conversions.*;
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      Coordsys(grid=[1, 1], component=[20, 20]),
      experiment(StopTime=5),
      Documentation(info="<html>
<p>
This example demonstrates the mixing of a single substance flow
between three tanks with different temperatures. The difference to
example \"ThreeTanksOneLiquid\" is that the port where the three
pipes are connected together contains a volume now, in order that
the mixing of the pipe flows is modelled more realistically.
</p>
</html>"));
    Components.Tank Tank1(
      area=1,
      T_start=from_degC(50),
      level_start=3,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initOption=Modelica_Fluid.Types.InitTypes.InitialValues) 
      annotation (extent=[-90, 20; -70, 40]);
    Components.Tank Tank2(
      area=1,
      level_start=1,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initOption=Modelica_Fluid.Types.InitTypes.InitialValues,
      T_start=from_degC(90)) 
      annotation (extent=[-10, 20; 10, 40]);
    Components.PressureDropPipe shortPipe1(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      frictionType=Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[-50, -40; -30, -20]);
    Components.Tank Tank3(
      area=1,
      T_start=from_degC(20),
      level_start=2,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initOption=Modelica_Fluid.Types.InitTypes.InitialValues) 
      "level(fixed=true)" 
      annotation (extent=[70,20; 90,40]);
    Components.PressureDropPipe shortPipe3(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      frictionType=Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[29, -40; 49, -20]);
    Components.PressureDropPipe shortPipe2(
      m_flow_nominal=1000,
      dp_nominal=from_bar(0.1),
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      frictionType=Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[-10, -10; 10, 10], rotation=-90);
    Utilities.PortVolume junctionVolume(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      V=1.e-4,
      T_start=from_degC(50.0),
      initOption=Modelica_Fluid.Types.InitTypes.InitialValues) 
                               annotation (extent=[-10, -40; 10, -20]);
  equation 
    connect(Tank1.port, shortPipe1.port_a) 
      annotation (points=[-80, 19; -80, -30; -51, -30], style(color=69));
    connect(shortPipe3.port_b, Tank3.port) 
      annotation (points=[50,-30; 80,-30; 80,19],    style(color=69));
    connect(Tank2.port, shortPipe2.port_a) 
      annotation (points=[0,19; 0,11; -6.73533e-016,11],    style(color=69));
    connect(shortPipe2.port_b, junctionVolume.port) annotation (points=[
          6.73533e-016,-11; 0,-11; 0,-30],    style(color=69));
    connect(shortPipe1.port_b, junctionVolume.port) 
      annotation (points=[-29, -30; 0, -30], style(color=69));
    connect(shortPipe3.port_a, junctionVolume.port) 
      annotation (points=[28, -30; 0, -30], style(color=69));
  end ThreeTanksWithPortVolume;
  
  model ThreeTanksIF97 
    import Modelica.SIunits.Conversions.*;
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      Coordsys(grid=[1, 1], component=[20, 20]),
      experiment(StopTime=5),
      Documentation(info="<html>
<p>
This example demonstrates the mixing of a single substance flow
between three tanks with different temperature.
</p>
<p>
When comparing with Examples.Tanks.ThreeTanks, it is demonstrated
that it is easy in this case to switch between single and multiple
substance flow. The only differences between the examples ThreeTanks
and ThreeTanksOneLiquid are that the Medium is different and that
for the model ThreeTanksOneLiquid the default values are
used for the initial mass fractions.
</p>
</html>"));
    Components.Tank Tank1(
      area=1,
      T_start=from_degC(50),
      level_start=3,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      initOption=Modelica_Fluid.Types.InitTypes.InitialValues) 
      annotation (extent=[-90, 20; -70, 40]);
    Components.Tank Tank2(
      area=1,
      level_start=1,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      initOption=Modelica_Fluid.Types.InitTypes.InitialValues,
      T_start=from_degC(90)) 
      annotation (extent=[-11,20; 9,40]);
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
      initOption=Modelica_Fluid.Types.InitTypes.InitialValues) 
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
  equation 
    connect(Tank1.port, shortPipe1.port_a) 
      annotation (points=[-80, 19; -80, -20; -51, -20], style(color=69));
    connect(shortPipe3.port_b, Tank3.port) 
      annotation (points=[51, -20; 80, -20; 80, 19], style(color=69));
    connect(Tank2.port, shortPipe2.port_a) 
      annotation (points=[-1,19; -1,11; -6.73533e-016,11],  style(color=69));
    connect(shortPipe1.port_b, shortPipe3.port_a) 
      annotation (points=[-29, -20; 29, -20], style(color=69));
    connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
          6.73533e-016,-11; 0,-11; 0,-20; 29,-20],     style(color=69));
  end ThreeTanksIF97;
  
  model ThreeTanksIF97SteadyHydraulic 
    import Modelica.SIunits.Conversions.*;
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      Coordsys(grid=[1, 1], component=[20, 20]),
      experiment(StopTime=5),
      Documentation(info="<html>
<p>
This example demonstrates the mixing of a single substance flow
between three tanks with different temperature.
</p>
<p>
When comparing with Examples.Tanks.ThreeTanks, it is demonstrated
that it is easy in this case to switch between single and multiple
substance flow. The only differences between the examples ThreeTanks
and ThreeTanksOneLiquid are that the Medium is different and that
for the model ThreeTanksOneLiquid the default values are
used for the initial mass fractions.
</p>
</html>"));
    Components.Tank Tank1(
      area=1,
      T_start=from_degC(50),
      level_start=3,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      initOption=Modelica_Fluid.Types.InitTypes.SteadyState) 
      annotation (extent=[-90, 20; -70, 40]);
    Components.Tank Tank2(
      area=1,
      level_start=1,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      T_start=from_degC(90),
      initOption=Modelica_Fluid.Types.InitTypes.SteadyState) 
      annotation (extent=[-11,20; 9,40]);
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
      initOption=Modelica_Fluid.Types.InitTypes.SteadyState) 
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
  equation 
    connect(Tank1.port, shortPipe1.port_a) 
      annotation (points=[-80, 19; -80, -20; -51, -20], style(color=69));
    connect(shortPipe3.port_b, Tank3.port) 
      annotation (points=[51, -20; 80, -20; 80, 19], style(color=69));
    connect(Tank2.port, shortPipe2.port_a) 
      annotation (points=[-1,19; -1,11; -6.73533e-016,11],  style(color=69));
    connect(shortPipe1.port_b, shortPipe3.port_a) 
      annotation (points=[-29, -20; 29, -20], style(color=69));
    connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
          6.73533e-016,-11; 0,-11; 0,-20; 29,-20],     style(color=69));
  end ThreeTanksIF97SteadyHydraulic;
end Tanks;
