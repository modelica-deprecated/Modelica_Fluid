package TestComponents 
  "Test components (this package will be removed for the final version)" 
  
  extends Icons.Library;
  
  model TestShortPipe "Test ShortPipe component" 
    import Modelica.SIunits.Conversions.*;
    
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      experiment(StopTime=1),
      Coordsys(grid=[1, 1], component=[20, 20]));
    Components.ShortPipe ShortPipe3(
      from_dp=false,
      dp_nominal=from_bar(0.1),
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.DetailedFriction,
      diameter=20,
      roughness=2e-5,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-2, 0; 18, 20]);
    Sources.MassFlowSource MassFlowSource3(T_ambient=from_degC(30),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-41, 0; -21, 20]);
    Modelica_Fluid.Sources.FixedAmbient ambient3(T_ambient=from_degC(15),
        redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[59, 0; 39, 20]);
    Modelica.Blocks.Sources.Ramp ramp1(
      height={6},
      duration={3},
      offset={-3}) annotation (extent=[-80, 0; -60, 20]);
  equation 
    connect(MassFlowSource3.port, ShortPipe3.port_a) 
      annotation (points=[-20, 10; -3, 10], style(color=69));
    connect(ShortPipe3.port_b, ambient3.port) 
      annotation (points=[19, 10; 38, 10], style(color=69));
    connect(ramp1.outPort, MassFlowSource3.m_dot) 
      annotation (points=[-59, 10; -43, 10], style(color=3));
  end TestShortPipe;
  
  model TestLongPipe "Test ShortPipe componet" 
    import Modelica.SIunits.Conversions.*;
    
    extends Modelica.Icons.Example;
    Modelica_Fluid.Sources.FixedAmbient ambient(T_ambient=from_degC(15),
        redeclare package Medium = Modelica_Media.Air.SimpleAir) 
      annotation (extent=[60, 0; 40, 20]);
    annotation (Diagram, experiment(StopTime=4));
    Components.LongPipe LongPipe(
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Air.SimpleAir,
      dp_nominal=from_bar(0.01),
      diameter=0.05,
      length=1) annotation (extent=[0, 0; 20, 20]);
    
    Sources.MassFlowSource MassFlowSource1(T_ambient=from_degC(30),
        redeclare package Medium = Modelica_Media.Air.SimpleAir) 
      annotation (extent=[-40, 0; -20, 20]);
    Modelica.Blocks.Sources.Ramp ramp(
      height={6},
      duration={3},
      offset={-3}) annotation (extent=[-80, 0; -60, 20]);
  equation 
    connect(LongPipe.port_b, ambient.port) 
      annotation (points=[21, 10; 39, 10], style(color=69));
    connect(MassFlowSource1.port, LongPipe.port_a) 
      annotation (points=[-19, 10; -1, 10], style(color=69));
    connect(ramp.outPort, MassFlowSource1.m_dot) 
      annotation (points=[-59, 10; -42, 10], style(color=3));
  end TestLongPipe;
  
  model TestCheckStateSelection 
    "Check whether for the choosen medium the expected states are selected" 
    
    import Modelica.SIunits.Conversions.*;
    
    Modelica_Fluid.Interfaces.JunctionVolume junctionVolume1(
      V=1.e-4,
      T_start=from_degC(50.0),
      redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
      annotation (extent=[-10, 90; 10, 70], rotation=180);
    extends Modelica.Icons.Example;
    Modelica_Fluid.Sources.FixedAmbient fixedAmbient1a(
      p_ambient=from_bar(1.0),
      redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
      T_ambient=from_degC(10)) annotation (extent=[-80, 70; -60, 90]);
    Modelica_Fluid.Examples.Components.ShortPipe shortPipe1a(
      m_dot_nominal=10,
      redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[-40, 70; -20, 90]);
    annotation (
      Diagram(
        Text(
          extent=[-100, 102; -52, 92],
          string="SimpleLiquidWater",
          style(color=58)),
        Text(
          extent=[-113, 55; -54, 47],
          style(color=58),
          string="SimpleAir"),
        Text(
          extent=[-113, 5; -54, -3],
          style(color=58),
          string="DetailedAir")),
      Documentation(info="<html>
<p>
This model is used to select a medium model and then translate.
The log should give information about the selected states.
They should be the ones choosen as StateSelect.prefer in the
corresponding medium model.
</p>
</html>
"),   Coordsys(grid=[1, 1], component=[20, 20]));
    
    Modelica_Fluid.Examples.Components.ShortPipe shortPipe1b(
      m_dot_nominal=10,
      redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[20, 70; 40, 90]);
    Modelica_Fluid.Sources.FixedAmbient fixedAmbient1b(
      p_ambient=from_bar(0.95),
      T_ambient=from_degC(30),
      redeclare package Medium = Modelica_Media.Water.SimpleLiquidWater) 
      annotation (extent=[80, 70; 60, 90]);
    
    Modelica_Fluid.Interfaces.JunctionVolume junctionVolume2(
      redeclare package Medium = Modelica_Media.Air.SimpleAir,
      p_start=from_bar(0.96),
      T_start=from_degC(50.0),
      V=0.1) annotation (extent=[-11, 40; 9, 20], rotation=180);
    
    Modelica_Fluid.Sources.FixedAmbient fixedAmbient2a(
      p_ambient=from_bar(1.0),
      T_ambient=from_degC(10),
      redeclare package Medium = Modelica_Media.Air.SimpleAir) 
      annotation (extent=[-80, 20; -60, 40]);
    Modelica_Fluid.Examples.Components.ShortPipe shortPipe2a(
      redeclare package Medium = Modelica_Media.Air.SimpleAir,
      dp_nominal=from_bar(0.01),
      m_dot_nominal=1,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[-40, 20; -20, 40]);
    
    Modelica_Fluid.Examples.Components.ShortPipe shortPipe2b(
      redeclare package Medium = Modelica_Media.Air.SimpleAir,
      dp_nominal=from_bar(0.01),
      m_dot_nominal=1,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[20, 20; 40, 40]);
    
    Modelica_Fluid.Sources.FixedAmbient fixedAmbient2b(
      p_ambient=from_bar(0.95),
      redeclare package Medium = Modelica_Media.Air.SimpleAir,
      T_ambient=from_degC(10)) annotation (extent=[80, 20; 60, 40]);
    
    Modelica_Fluid.Interfaces.JunctionVolume junctionVolume3(
      p_start=from_bar(0.96),
      T_start=from_degC(50.0),
      V=0.1,
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[-11, -10; 9, -30], rotation=180);
    Modelica_Fluid.Sources.FixedAmbient fixedAmbient3a(
      p_ambient=from_bar(1.0),
      T_ambient=from_degC(10),
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[-80, -30; -60, -10]);
    Modelica_Fluid.Examples.Components.ShortPipe shortPipe3a(
      dp_nominal=from_bar(0.01),
      m_dot_nominal=1,
      redeclare package Medium = Modelica_Media.Air.DetailedAir,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[-40, -30; -20, -10]);
    Modelica_Fluid.Examples.Components.ShortPipe shortPipe3b(
      dp_nominal=from_bar(0.01),
      m_dot_nominal=1,
      redeclare package Medium = Modelica_Media.Air.DetailedAir,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[20, -30; 40, -10]);
    Modelica_Fluid.Sources.FixedAmbient fixedAmbient3b(
      p_ambient=from_bar(0.95),
      T_ambient=from_degC(10),
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[80, -30; 60, -10]);
  equation 
    connect(fixedAmbient1a.port, shortPipe1a.port_a) 
      annotation (points=[-59, 80; -41, 80], style(color=69));
    connect(shortPipe1b.port_b, fixedAmbient1b.port) 
      annotation (points=[41, 80; 59, 80], style(color=69));
    connect(shortPipe1b.port_a, junctionVolume1.port) 
      annotation (points=[19, 80; 0, 80], style(color=69));
    connect(shortPipe1a.port_b, junctionVolume1.port) 
      annotation (points=[-19, 80; 0, 80], style(color=69));
    connect(fixedAmbient2a.port, shortPipe2a.port_a) 
      annotation (points=[-59, 30; -41, 30], style(color=69));
    connect(shortPipe2b.port_b, fixedAmbient2b.port) 
      annotation (points=[41, 30; 59, 30], style(color=69));
    connect(shortPipe2b.port_a, junctionVolume2.port) 
      annotation (points=[19, 30; -1, 30], style(color=69));
    connect(shortPipe2a.port_b, junctionVolume2.port) 
      annotation (points=[-19, 30; -1, 30], style(color=69));
    connect(fixedAmbient3a.port, shortPipe3a.port_a) 
      annotation (points=[-59, -20; -41, -20], style(color=69));
    connect(shortPipe3b.port_b, fixedAmbient3b.port) 
      annotation (points=[41, -20; 59, -20], style(color=69));
    connect(shortPipe3b.port_a, junctionVolume3.port) 
      annotation (points=[19, -20; -1, -20], style(color=69));
    connect(shortPipe3a.port_b, junctionVolume3.port) 
      annotation (points=[-19, -20; -1, -20], style(color=69));
  end TestCheckStateSelection;
  
  model TwoVolumesAir 
    "This example shows the difference of a JunctionVolume and of a FixedComponentVolume" 
    
    import SI = Modelica.SIunits;
    parameter SI.Volume V=0.05 "Size of volume";
    
    Interfaces.JunctionVolume junctionVolume1(
      T_start=from_degC(50.0),
      V=V/2,
      redeclare package Medium = Modelica_Media.Air.DetailedAir,
      initType=Modelica_Media.Interfaces.PartialMedium.Choices.Init.
          InitialStates) annotation (extent=[-30, 40; -10, 20], rotation=180);
    
    import Modelica.SIunits.Conversions.*;
    
    extends Modelica.Icons.Example;
    Sources.FixedAmbient fixedAmbient1(
      p_ambient=from_bar(1.0),
      T_ambient=from_degC(10),
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[-100, 20; -80, 40]);
    Components.ShortPipe shortPipe1(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[-60, 20; -40, 40]);
    annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
    Components.ShortPipe shortPipe3(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[40, 20; 60, 40]);
    Sources.FixedAmbient fixedAmbient3(
      T_ambient=from_degC(20),
      p_ambient=from_bar(0.95),
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[100, 20; 80, 40]);
    Sources.FixedAmbient fixedAmbient2(
      p_ambient=from_bar(1.0),
      T_ambient=from_degC(20),
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[-100, 60; -80, 80]);
    Components.ShortPipe shortPipe2(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[-60, 60; -40, 80]);
    Interfaces.JunctionVolume junctionVolume2(
      T_start=from_degC(50.0),
      V=V/2,
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[0, 40; 20, 20], rotation=180);
    Interfaces.JunctionVolume junctionVolume3(
      T_start=from_degC(50.0),
      V=V,
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[-30, -50; -10, -70], rotation=180);
    Sources.FixedAmbient fixedAmbient4(
      p_ambient=from_bar(1.0),
      T_ambient=from_degC(10),
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[-100, -70; -80, -50]);
    Components.ShortPipe shortPipe4(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[-60, -70; -40, -50]);
    Components.ShortPipe shortPipe5(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[0, -70; 20, -50]);
    Sources.FixedAmbient fixedAmbient5(
      T_ambient=from_degC(20),
      p_ambient=from_bar(0.95),
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[60, -70; 40, -50]);
    Sources.FixedAmbient fixedAmbient6(
      p_ambient=from_bar(1.0),
      T_ambient=from_degC(20),
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[-100, -30; -80, -10]);
    Components.ShortPipe shortPipe6(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Air.DetailedAir) 
      annotation (extent=[-60, -30; -40, -10]);
  equation 
    connect(fixedAmbient1.port, shortPipe1.port_a) 
      annotation (points=[-79, 30; -61, 30], style(color=69));
    connect(shortPipe3.port_b, fixedAmbient3.port) 
      annotation (points=[61, 30; 79, 30], style(color=69));
    connect(fixedAmbient2.port, shortPipe2.port_a) 
      annotation (points=[-79, 70; -61, 70], style(color=69));
    connect(shortPipe2.port_b, junctionVolume1.port) 
      annotation (points=[-39, 70; -20, 70; -20, 30], style(color=69));
    connect(shortPipe1.port_b, junctionVolume1.port) 
      annotation (points=[-39, 30; -20, 30], style(color=69));
    connect(junctionVolume2.port, shortPipe3.port_a) 
      annotation (points=[10, 30; 39, 30], style(color=69));
    connect(junctionVolume1.port, junctionVolume2.port) 
      annotation (points=[-20, 30; 10, 30], style(color=69));
    connect(fixedAmbient4.port, shortPipe4.port_a) 
      annotation (points=[-79, -60; -61, -60], style(color=69));
    connect(shortPipe5.port_b, fixedAmbient5.port) 
      annotation (points=[21, -60; 39, -60], style(color=69));
    connect(shortPipe5.port_a, junctionVolume3.port) 
      annotation (points=[-1, -60; -20, -60], style(color=69));
    connect(fixedAmbient6.port, shortPipe6.port_a) 
      annotation (points=[-79, -20; -61, -20], style(color=69));
    connect(shortPipe6.port_b, junctionVolume3.port) 
      annotation (points=[-39, -20; -20, -20; -20, -60], style(color=69));
    connect(shortPipe4.port_b, junctionVolume3.port) 
      annotation (points=[-39, -60; -20, -60], style(color=69));
  end TwoVolumesAir;
  
  model TwoVolumesDetailedWater 
    "This example shows the difference of a JunctionVolume and of a FixedComponentVolume" 
    
    import SI = Modelica.SIunits;
    parameter SI.Volume V=1.e-4 "Size of volume";
    
    Interfaces.JunctionVolume junctionVolume1(
      T_start=from_degC(50.0),
      V=V/2,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-30, 40; -10, 20], rotation=180);
    
    import Modelica.SIunits.Conversions.*;
    
    extends Modelica.Icons.Example;
    Sources.FixedAmbient fixedAmbient1(
      p_ambient=from_bar(1.0),
      T_ambient=from_degC(10),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-100, 20; -80, 40]);
    Components.ShortPipe shortPipe1(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-60, 20; -40, 40]);
    annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
    Components.ShortPipe shortPipe3(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[40, 20; 60, 40]);
    Sources.FixedAmbient fixedAmbient3(
      T_ambient=from_degC(20),
      p_ambient=from_bar(0.95),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[100, 20; 80, 40]);
    Sources.FixedAmbient fixedAmbient2(
      p_ambient=from_bar(1.0),
      T_ambient=from_degC(20),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-100, 60; -80, 80]);
    Components.ShortPipe shortPipe2(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-60, 60; -40, 80]);
    Interfaces.JunctionVolume junctionVolume2(
      T_start=from_degC(50.0),
      V=V/2,
      redeclare package Medium = Modelica_Media.Water.IF97_ph,
      initType=Modelica_Media.Interfaces.PartialMedium.Choices.Init.
          InitialStates) annotation (extent=[0, 40; 20, 20], rotation=180);
    
    Interfaces.JunctionVolume junctionVolume3(
      T_start=from_degC(50.0),
      V=V,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-30, -50; -10, -70], rotation=180);
    Sources.FixedAmbient fixedAmbient4(
      p_ambient=from_bar(1.0),
      T_ambient=from_degC(10),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-100, -70; -80, -50]);
    Components.ShortPipe shortPipe4(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-60, -70; -40, -50]);
    Components.ShortPipe shortPipe5(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[0, -70; 20, -50]);
    Sources.FixedAmbient fixedAmbient5(
      T_ambient=from_degC(20),
      p_ambient=from_bar(0.95),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[60, -70; 40, -50]);
    Sources.FixedAmbient fixedAmbient6(
      p_ambient=from_bar(1.0),
      T_ambient=from_degC(20),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-100, -30; -80, -10]);
    Components.ShortPipe shortPipe6(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-60, -30; -40, -10]);
  equation 
    connect(fixedAmbient1.port, shortPipe1.port_a) 
      annotation (points=[-79, 30; -61, 30], style(color=69));
    connect(shortPipe3.port_b, fixedAmbient3.port) 
      annotation (points=[61, 30; 79, 30], style(color=69));
    connect(fixedAmbient2.port, shortPipe2.port_a) 
      annotation (points=[-79, 70; -61, 70], style(color=69));
    connect(shortPipe2.port_b, junctionVolume1.port) 
      annotation (points=[-39, 70; -20, 70; -20, 30], style(color=69));
    connect(shortPipe1.port_b, junctionVolume1.port) 
      annotation (points=[-39, 30; -20, 30], style(color=69));
    connect(junctionVolume2.port, shortPipe3.port_a) 
      annotation (points=[10, 30; 39, 30], style(color=69));
    connect(junctionVolume1.port, junctionVolume2.port) 
      annotation (points=[-20, 30; 10, 30], style(color=69));
    connect(fixedAmbient4.port, shortPipe4.port_a) 
      annotation (points=[-79, -60; -61, -60], style(color=69));
    connect(shortPipe5.port_b, fixedAmbient5.port) 
      annotation (points=[21, -60; 39, -60], style(color=69));
    connect(shortPipe5.port_a, junctionVolume3.port) 
      annotation (points=[-1, -60; -20, -60], style(color=69));
    connect(fixedAmbient6.port, shortPipe6.port_a) 
      annotation (points=[-79, -20; -61, -20], style(color=69));
    connect(shortPipe6.port_b, junctionVolume3.port) 
      annotation (points=[-39, -20; -20, -20; -20, -60], style(color=69));
    connect(shortPipe4.port_b, junctionVolume3.port) 
      annotation (points=[-39, -60; -20, -60], style(color=69));
  end TwoVolumesDetailedWater;
  
  model OneVolumeDetailedWater 
    "This example shows the difference of a JunctionVolume and of a FixedComponentVolume" 
    
    import SI = Modelica.SIunits;
    parameter SI.Volume V=1.e-4 "Size of volume";
    
    import Modelica.SIunits.Conversions.*;
    
    extends Modelica.Icons.Example;
    annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
    Interfaces.JunctionVolume junctionVolume3(
      T_start=from_degC(50.0),
      V=V,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-10, 0; 10, -20], rotation=180);
    Sources.FixedAmbient fixedAmbient4(
      p_ambient=from_bar(1.0),
      T_ambient=from_degC(10),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-80, -20; -60, 0]);
    Components.ShortPipe shortPipe4(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-40, -20; -20, 0]);
    Components.ShortPipe shortPipe5(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[20, -20; 40, 0]);
    Sources.FixedAmbient fixedAmbient5(
      T_ambient=from_degC(20),
      p_ambient=from_bar(0.95),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[80, -20; 60, 0]);
    Sources.FixedAmbient fixedAmbient6(
      p_ambient=from_bar(1.0),
      T_ambient=from_degC(20),
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-80, 20; -60, 40]);
    Components.ShortPipe shortPipe6(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica_Media.Water.IF97_ph) 
      annotation (extent=[-40, 20; -20, 40]);
  equation 
    connect(fixedAmbient4.port, shortPipe4.port_a) 
      annotation (points=[-59, -10; -41, -10], style(color=69));
    connect(shortPipe5.port_b, fixedAmbient5.port) 
      annotation (points=[41, -10; 59, -10], style(color=69));
    connect(shortPipe5.port_a, junctionVolume3.port) 
      annotation (points=[19, -10; 0, -10], style(color=69));
    connect(fixedAmbient6.port, shortPipe6.port_a) 
      annotation (points=[-59, 30; -41, 30], style(color=69));
    connect(shortPipe6.port_b, junctionVolume3.port) 
      annotation (points=[-19, 30; 0, 30; 0, -10], style(color=69));
    connect(shortPipe4.port_b, junctionVolume3.port) 
      annotation (points=[-19, -10; 0, -10], style(color=69));
  end OneVolumeDetailedWater;
  
  model TwoVolumesApproximationWater 
    "This example shows the difference of a JunctionVolume and of a FixedComponentVolume" 
    
    import SI = Modelica.SIunits;
    parameter SI.Volume V=1.e-4 "Size of volume";
    
    Interfaces.JunctionVolume junctionVolume1(
      T_start=from_degC(50.0),
      V=V/2,
      p_start=from_bar(100.0),
      initType=Modelica_Media.Interfaces.PartialMedium.Choices.Init.
          InitialStates,
    redeclare package Medium = Modelica_Media.Water.IF97LocalApproximation_ph) 
      annotation (extent=[-12,40; 8,20],     rotation=180);
    
    import Modelica.SIunits.Conversions.*;
    
    extends Modelica.Icons.Example;
    Sources.FixedAmbient fixedAmbient1(
      p_ambient=from_bar(100.0),
      T_ambient=from_degC(110),
    redeclare package Medium = Modelica_Media.Water.IF97LocalApproximation_ph) 
      annotation (extent=[-100, 20; -80, 40]);
    
    Components.ShortPipe shortPipe1(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      T_start=400.0,
    redeclare package Medium = Modelica_Media.Water.IF97LocalApproximation_ph) 
      annotation (extent=[-60, 20; -40, 40]);
    
    annotation (Diagram, Documentation(info="<html>
<p>
This example shows two ways of merging flows. In the upper part
a volume in the connection port is modeled. At the initial time the
fluid in this volume has a higher temperature as the fluid
flowing form the ambient. Therefore, im shortPipe the temperature
will need some time to arrive at the mixing temperature.
</p>
<p>
In the lower part the volume in the connection port is neglected
and ideal mixing takes place. Since there is a steady state mass
flow and no (existing) delay mechanism is modeled, the mixing
is instantaneous. The temperature is the same as in the upper part.
</p>
</html>
"));
    Components.ShortPipe shortPipe3(
      m_dot_nominal=10,
      frictionType=Modelica_Fluid.Examples.Types.FrictionTypes.ConstantLaminar,
      T_start=400,
    redeclare package Medium = Modelica_Media.Water.IF97LocalApproximation_ph) 
      annotation (extent=[40, 20; 60, 40]);
    
    Sources.FixedAmbient fixedAmbient3(
      p_ambient=from_bar(98),
      T_ambient=from_degC(120),
    redeclare package Medium = Modelica_Media.Water.IF97LocalApproximation_ph) 
      annotation (extent=[100, 20; 80, 40]);
    
  equation 
    connect(fixedAmbient1.port, shortPipe1.port_a) 
      annotation (points=[-79, 30; -61, 30], style(color=69));
    connect(shortPipe3.port_b, fixedAmbient3.port) 
      annotation (points=[61, 30; 79, 30], style(color=69));
    connect(shortPipe1.port_b, junctionVolume1.port) 
      annotation (points=[-39,30; -2,30],    style(color=69));
    connect(junctionVolume1.port, shortPipe3.port_a) annotation (points=[-2,
          30; 39,30], style(color=69, rgbcolor={0,127,255}));
  end TwoVolumesApproximationWater;
  
end TestComponents;
